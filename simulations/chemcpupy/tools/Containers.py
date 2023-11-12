#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 19:01:33 2018

@author: jacobrosenstein
"""


from chemcpupy import MixtureList,CompoundList
#from chemcpupy.tools.misc import *
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import csv
import os
import pubchempy as pcp
from pyteomics import mass


# all volumes are in nL



# ======================================================================        



class Container(MixtureList):
    
    def __init__(self,*args,**kwargs):        
        """ A Container is a list of mixtures in a well plate. It 
        inherits from MixtureList. The difference is that a Container also
        stores a location with each mixture.
        """
        self._rows = kwargs.get('rows',1)
        self._cols = kwargs.get('cols',1)
        self.shape = (self._rows,self._cols)
        self._volume_max = kwargs.get('volume_max',1e38)
        self._volume_min = kwargs.get('volume_min',0)
        self._volume_increment = kwargs.get('volume_increment', 1)
        self._autofill = kwargs.get('autofill',True)
        
        if kwargs.get('mixture_json_file') is not None:
            self._autofill = False
        
        MixtureList.__init__(self,**kwargs)
        self._file_subtype='WellPlate-1.0'
    
        if self._autofill:
            for r in range(self._rows):    
                for c in range(self._cols):
                    self._mixture_list.append(
                            {'mixture_name':position_to_lettergrid((r,c)),
                             'compound_list':CompoundList(),
                             'position':(r,c),
                             'volume':0
                             })

    def __getitem__(self,key):
        """ In a Container you can retrieve a mixture either by 
        myplate['mixturename'] or by myplate[row,col].
        """

        if type(key)==tuple or type(key)==list:
            row,col=tuple(key)
            #print('testing key',key,row,col,self._rows,self._cols)
            assert(row<self._rows)
            assert(col<self._cols)
            return [x for x in self._mixture_list if tuple(x['position'])==(row,col)][0]
        else:
            return MixtureList.__getitem__(self,key) 
    
    def _find_empty_location(self):
        filled_positions = [tuple(x['position']) for x in self._mixture_list]        
        for r in range(self._rows):
                for c in range(self._cols):
                    if (r,c) not in filled_positions:
                        return (r,c)
        
        raise Exception('No empty WellPlate positions available.')


    def count_filled(self):
        return len([x for x in self._mixture_list if len(x['compound_list'])>0])

    def list_filled_positions(self):
        return [tuple(x['position']) for x in self._mixture_list if len(x['compound_list'])>0]

    def list_empty_positions(self):
        return [tuple(x['position']) for x in self._mixture_list if len(x['compound_list'])==0]
    
    def add_compounds_to_location(self,position,compound_list,volume=0):
        if type(position) is str:
            position = lettergrid_to_position(position)
            
        mixture_name = [x['mixture_name'] for x in self._mixture_list \
                        if tuple(x['position'])==position][0]
        self.add_compounds_to_mixture(mixture_name,copy.deepcopy(compound_list),\
                                      volume=volume, position=position)
        self.add_volume(position,volume)
        
    def add_new_mixture(self,mixture_name,compound_list,position=None,volume=0):
        """ Add a mixture to a Container at a specified location.
        If no location is specified, it will be put at the first empty location.
        """
        # if a certain volume exceeds the capacity of one well,
        # transfer the mixture to multiple wells
        # we also need to make sure that we won't transfer past the dead volume
        volume_range = self._volume_max - self._volume_min
        num_well_needed = 1 if (volume==0) else int(volume/volume_range) + (volume % volume_range > 0)
        for i in range(num_well_needed):
            if position==None or i > 0:
                position = self._find_empty_location()

            row,col=position
            assert(row<self._rows)
            assert(col<self._cols)

            if position in [tuple(x['position']) for x in self._mixture_list]:
                raise Exception('there is already a mixture at position %s' % (position,))

            mycompounds = copy.deepcopy(compound_list)
            if volume==0:
                transfer_volume=0
            else:
                if volume > volume_range:
                    transfer_volume = self._volume_max
                    volume -= volume_range
                else:
                    transfer_volume = volume + self._volume_min
                    
                formervolume = sum([c['volume'] for c in mycompounds])

                if (formervolume > 0):
                    for c in mycompounds:
                        c['volume'] = transfer_volume*c['volume']/formervolume
                else:
                    for c in mycompounds:
                        c['volume'] = transfer_volume/len(mycompounds)

            self._mixture_list.append({'mixture_name':mixture_name,
                                       'compound_list':mycompounds,
                                       'position':tuple(position),
                                       'volume':transfer_volume
                                       })

    def add_compounds_from_bounding_box(self,compounds,volume=100e-6):
        for c in compounds:
            bounds = c['bounding box'].split(' ')
            includerows = range(ord(bounds[0])-ord('A'),ord(bounds[1])-ord('A')+1)
            includecols = range(int(bounds[2])-1,int(bounds[3]))
            for row in includerows:
                for col in includecols:
                    self.add_compounds_to_location((row,col),
                                                     CompoundList(compound_list=(c,)),
                                                     volume=volume
                                                     )
        
    def export_csv_review(self,filename,typefilter=None):
        if typefilter is None:
                positions = self.list_filled_positions();
        else:
                positions = self.get_positions_by_compound_type(typefilter);                
        compound_positions = [];
        compound_masses = [];
        compound_names = [];
        compound_ids = [];
        compound_type = [];
        compound_volume = [];
        well_volume = [];
        for pos in positions:
            info = self.get_compounds(pos,['mass','name','cid','type','volume']); 
            vol = self.get_volume(pos);
            for i in info:
                current_type = i[3];
                if (typefilter is None) or ( (typefilter is not None) and (current_type==typefilter) ):
                    compound_positions.append(position_to_lettergrid(pos));                
                    compound_masses.append(i[0]);
                    compound_names.append(i[1]);
                    compound_ids.append(i[2]);  
                    compound_type.append(i[3]);  
                    compound_volume.append(i[4]);  
                    well_volume.append(vol);                                                  
        fields = ['position', 'mass', 'id','name', 'type', 'volume', 'well_volume'];
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fields)
            writer.writeheader()
            for n in range(0,len(compound_positions)):
                rowdata = {'position':compound_positions[n],
                           'mass':compound_masses[n],
                           'id':compound_ids[n],
                           'name':compound_names[n],
                           'type':compound_type[n],
                           'volume':compound_volume[n],                           
                           'well_volume':well_volume[n]};
                writer.writerow(rowdata);
        print('Wrote Container CSV review: '+filename)       


    def export_array_review(self,typefilter=None):
        if typefilter is None:
                positions = self.list_filled_positions();
        else:
                positions = self.get_positions_by_compound_type(typefilter);                
        compound_positions = [];
        compound_masses = [];
        compound_names = [];
        compound_ids = [];
        compound_type = [];
        compound_volume = [];
        well_volume = [];
        for pos in positions:
            info = self.get_compounds(pos,['mass','name','cid','type','volume']); 
            vol = self.get_volume(pos);
            for i in info:
                current_type = i[3];
                if (typefilter is None) or ( (typefilter is not None) and (current_type==typefilter) ):
                    compound_positions.append(position_to_lettergrid(pos));   
                    compound_masses.append(i[0]);
                    compound_names.append(i[1]);
                    compound_ids.append(i[2]);  
                    compound_type.append(i[3]);  
                    compound_volume.append(i[4]);  
                    well_volume.append(vol);                                                                      
        return {'position':compound_positions, 'mass':compound_masses, 'id':compound_ids,'name':compound_names, 'type':compound_type, 'volume':compound_volume, 'well_volume':well_volume};
         
        
    def load_platesheet(self,filename,fill_from_pubchem=True,fill_from_csv=False,grid='letter'):
        """ Load the Container from a specified PlateSheet .csv file
        """
        if not os.path.isfile(filename):
            raise Exception('The specified PlateSheet file does not exist.')
            
        df = pd.read_csv(filename, 
                         skip_blank_lines=True,
                         skiprows=5)
        df.Volume *= 1e9    # convert from L to nL
        df = df[df.Volume>0]   # remove entries with zero volume
        
        for index,row in df.iterrows():
                            
            myproperties = {'type':row.Type,
                            'name':row.Name,}
            
            
            if grid=='maldi': 
                # convert to lettergrid for compatibility with below methods
                row.Start = position_to_lettergrid(maldigrid_to_position(row.Start));
                row.End = position_to_lettergrid(maldigrid_to_position(row.End));
            #elif grid=='letter':
            #    # do nothing            
            
            if fill_from_csv:
                fill_from_pubchem = False;
                myproperties['mass'] = row.Mass   
                myproperties['formula'] = row.Formula
                myproperties['concentration'] = row.Concentration
                myproperties['description'] = row.Description
                myproperties['solvent'] = row.Solvent
                myproperties['id'] = row.CID
                        
            if fill_from_pubchem:
                print('Fetch from PubChem: ',int(row.CID))
                m = pcp.Compound.from_cid(int(row.CID))
                myproperties['formula'] = m.molecular_formula.replace('+','')
                myproperties['structure'] = m.isomeric_smiles
                #myproperties['mass'] = m.monoisotopic_mass;                
                myproperties['mass'] = mass.calculate_mass(formula=myproperties['formula'])
                           
            if row.End is np.nan:
                # no End position was specified. Put the compound at the Start.
                pos = lettergrid_to_position(row.Start)
                self.add_compounds_to_location(pos,
                                               CompoundList(cid_list=(row.CID,)),
                                               volume=row.Volume
                                               )
                self.__getitem__(pos)['compound_list'].add_properties_by_cid({row.CID:myproperties})
            else:
                # put the compound in all positions spanning Start to End
                for pos in rectangle_positions(row.Start,row.End):
                    #print('adding to',pos)
                    self.add_compounds_to_location(pos,
                                                   CompoundList(cid_list=(row.CID,)),
                                                   volume=row.Volume
                                                   )
                    self.__getitem__(pos)['compound_list'].add_properties_by_cid({row.CID:myproperties})            
            

    def get_volume(self,key):
        """ Changethe volume of a mixture. The key can either be the mixture_name
        or (row,col).
        """
        found = False
        volume = 0;
        for x in self._mixture_list:
            if (x['mixture_name']==key) or (tuple(x['position'])==tuple(key)):
                found = True
                volume = x['volume'];            
        if (not found):
            raise Exception('Key %s is not found. Task aborted.'
                                    % (str(key)))
        return volume
    
    def add_volume(self,key,add_volume):
        """ Changethe volume of a mixture. The key can either be the mixture_name
        or (row,col).
        """
        found = False
        for x in self._mixture_list:
            if (x['mixture_name']==key) or (tuple(x['position'])==tuple(key)):
                # We already check these in Task.run(). No need to check again.
                '''
                current_position = x['position']
                if (add_volume < 0) and (x['volume'] + add_volume < self._volume_min):
                    raise Exception('Position %s is depleted. Task aborted.'
                                    % (current_position,))
                elif (add_volume > 0) and (x['volume'] + add_volume > self._volume_max):
                    raise Exception('Position %s is full. Task aborted.'
                                    % (current_position,))
                else:
                '''
                found = True
                x['volume'] += add_volume
        if (not found):
            raise Exception('Key %s is not found. Task aborted.'
                                    % (str(key)))

    def shape(self):
        return (self._rows, self._cols)

    def get_param(self, param):
        if param == 'volume_max':
            return self._volume_max
        elif param == 'volume_min':
            return self._volume_min
        elif param == 'volume_increment':
            return self._volume_increment
        #elif param == 'shape':
        #    return (self._rows, self._cols)
        else:
            raise Exception('Unknown parameter: %s'% (param,))
    
    def __str__(self):
        return 'WellPlate: description=%s, wells=%u, shape=%s, filled_wells=%u' % (self._description,
                                                                                   np.prod(self.shape),
                                                                                   self.shape,
                                                                                   len(self))

    def pretty_print(self):
        print(self)
        for x in self._mixture_list:
            if len(x['compound_list'])>0:
                print('------------')
                print('Name:',x['mixture_name'])
                print('Position:',x['position'])
                print('Volume:',x['volume'])
                print(x['compound_list'].describe_properties())
            else:
                print('(empty) ',x['position'])
    
    def pretty_print_volumes(self):
        for r in range(self._rows):    
            myrow=''
            for c in range(self._cols):
                try:
                    #unit is nL
                    myrow += '{0:10d}'.format(int(self[r,c]['volume']))
                except:
                    myrow += '     _____'    
#                myrow += ' %d'%(1e6*self[r,c]['volume'])
            print(myrow)

    
    def graphical_print(self,subtitle='',
                        colorfunction=None,
                        markermap={},
                        edgecolormap={},
                        facecolormap={},
                        textfunction=None):
                            
        f = plt.figure(figsize=(0.5*self._cols,0.5*self._rows))
        legend_dict = {}

        
        for col in range(self._cols):
            for row in range(self._rows):
                if textfunction is not None:
                    mytext = textfunction(self[row,col])
                
                if textfunction is None or mytext is None:
                    mytext = '%s %u' % ('ABCDEFGHIJKLMNOP'[row],col+1)
                    
                if colorfunction is None:
                    colorindex = 0
                else:
                    colorindex = colorfunction(self[row,col])
                
                legend_dict[colorindex]=facecolormap.get(colorindex,'none')
                plt.plot(col,row,
                         markermap.get(colorindex,'o'),
                         markersize=14,
                         markeredgecolor=edgecolormap.get(colorindex,'k'),
                         markerfacecolor=facecolormap.get(colorindex,'none'))
                plt.text(col,row,
                         mytext,
                         fontsize=8,
                         horizontalalignment='center',
                         verticalalignment='center')
                
        plt.gca().xaxis.tick_top()
        plt.tick_params(
            axis='both',
            which='both',
            bottom=False,
            top=False,   
            left=False,
            right=False)
        plt.gca().set_xticks(np.arange(self._cols))
        plt.gca().set_xticklabels(map(str,1+np.arange(self._cols)))
        plt.gca().set_yticks(np.arange(self._rows))
        plt.gca().set_yticklabels(rowletters)
        plt.gca().invert_yaxis()
        plt.text(1,self._rows+0.5,
                 subtitle,
                 fontsize=16)
        # print row numbers on the right side
        ylims = plt.ylim()
        ax2=plt.gca().twinx()
        ax2.set_yticks(np.arange(self._rows))
        ax2.set_yticklabels(1+np.arange(self._rows))
        ax2.invert_yaxis()
        ax2.set_ylim(ylims)
        
        # DRAW A COLOR LEGEND
        patchList = []
        for key in legend_dict:
            data_key = mpatches.Patch(color=legend_dict[key], label=key)
            patchList.append(data_key)        
        plt.legend(handles=patchList,
                   loc='upper left', 
                   bbox_to_anchor=(1.02, 1))

        plt.show()
        return f

        
    def graphical_print_volumes(self,units='nL'):
        if units is 'nL':
            textfunction=lambda w: '%2.0f' % (w['volume'])
        elif units is 'uL':
            textfunction=lambda w: '%2.1f' % (w['volume']/1e3)
        elif units is 'mL':
            textfunction=lambda w: '%2.2f' % (w['volume']/1e6)
        else:
            raise Exception('unknown volume units')
            
        self.graphical_print(subtitle='%s: volume (%s)' % (self._description,units),
                             textfunction=textfunction,
                             colorfunction=lambda w: w['volume']>0,
                             facecolormap={True:'lightblue'}
                             )

    def graphical_print_compounds_per_well(self):        
        max_compounds = max([len(w['compound_list']) for w in self])
        def mytextfunction(w):
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(enumerate(plt.cm.tab20(np.linspace(0, 1, 2*max_compounds))))
        myfacecolors[0]='none'
        self.graphical_print(subtitle=self._description+': number of compounds',
                     textfunction=mytextfunction,
                     colorfunction=lambda w: len(w['compound_list']),
                     facecolormap=myfacecolors
                     )

    def graphical_print_transfer_group(self):        
        all_groups = set([x.get('transfer group','') for x in self._mixture_list])
        def mytextfunction(w):
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(zip(all_groups,plt.cm.tab20(np.linspace(0,1,len(all_groups)))))
        myfacecolors['']='none'
        self.graphical_print(subtitle=self._description+': color by transfer group',
                     textfunction=mytextfunction,
                     colorfunction=lambda w: w.get('transfer group',''),
                     facecolormap=myfacecolors
                     )
        print('colors correspond to the following transfer groups: ',all_groups)
        
    def graphical_print_by_type(self):        
        all_groups = set([str(x['compound_list'].list_all('type')) for x in self._mixture_list])
        def mytextfunction(w):
            return ''
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(zip(all_groups,plt.cm.tab20(np.linspace(0,1,len(all_groups)))))
        myfacecolors['']='none'
        myfacecolors['[]']='none'
        myfacecolors["['acid']"] = [0.9, 0.9, 0.2, 1.0]
        myfacecolors["['base']"] = [0.1, 0.1, 0.7, 1.0]
        # print(myfacecolors)
        # print(self._description+': color by mass')
        # for w in self._mixture_list:
        #     print(w)
        #     print(str(w['compound_list'].list_all('type')))
        def clrfn(w):
            if len(w['compound_list']) > 1:
                if w['mixture_name'] in ['O1', 'P2']:
                    return "['base']"
                elif w['mixture_name'] in ['O2', 'P1']:
                    return "['acid']"
            return str(w['compound_list'].list_all('type'))
        # self.graphical_print(subtitle=self._description+': color by type',
        return self.graphical_print(subtitle='',
                     textfunction=mytextfunction,
                     colorfunction=clrfn,
                     facecolormap=myfacecolors
                     )
        print('colors correspond to the following groups: ',all_groups)
        

    def graphical_print_by_mass(self):        
        all_groups = set([str(x['compound_list'].list_all('mass')) for x in self._mixture_list])
        def mytextfunction(w):
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(zip(all_groups,plt.cm.tab20(np.linspace(0,1,len(all_groups)))))
        myfacecolors['']='none'
        myfacecolors['[]']='none'
        self.graphical_print(subtitle=self._description+': color by mass',
                     textfunction=mytextfunction,
                     colorfunction=lambda w: str(w['compound_list'].list_all('mass')),
                     facecolormap=myfacecolors
                     )
        print('colors correspond to the following groups: ',all_groups)
        
    def graphical_print_by_id(self):        
        all_groups = set([str(x['compound_list'].list_all('cid')) for x in self._mixture_list])
        def mytextfunction(w):
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(zip(all_groups,plt.cm.tab20(np.linspace(0,1,len(all_groups)))))
        myfacecolors['']='none'
        myfacecolors['[]']='none'
        self.graphical_print(subtitle=self._description+': color by PubChem cID',
                     textfunction=mytextfunction,
                     colorfunction=lambda w: str(w['compound_list'].list_all('cid')),
                     facecolormap=myfacecolors
                     )
        print('colors correspond to the following groups: ',all_groups)
        
        
    def graphical_print_by_mass_with_type_filter(self,typefilter=''):   
        all_groups = [];
        def get_label(mixlist):
            values = list(mixlist['compound_list'].list_all('mass'));
            types = list(mixlist['compound_list'].list_all('type'));
            if typefilter in types:
                index = types.index(typefilter) # assume type appears only once
                label = str(values[index]); 
            else:
                label = '';
            return label;
        for x in self._mixture_list:
            all_groups.append(get_label(x));    
        all_groups = set(all_groups);
        def mytextfunction(w):
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(zip(all_groups,plt.cm.tab20(np.linspace(0,1,len(all_groups)))))
        myfacecolors['']='none'
        myfacecolors['[]']='none'
        self.graphical_print(subtitle=self._description+': color by mass (filtered by type: '+str(typefilter)+')',
                     textfunction=mytextfunction,
                     colorfunction=lambda w: get_label(w),
                     facecolormap=myfacecolors
                     )
        print('colors correspond to the following groups: ',all_groups)
        
    def graphical_print_by_id_with_type_filter(self,typefilter=''):   
        all_groups = [];
        def get_label(mixlist):
            values = list(mixlist['compound_list'].list_all('cid'));
            types = list(mixlist['compound_list'].list_all('type'));
            if typefilter in types:
                index = types.index(typefilter) # assume type appears only once
                label = str(values[index]); 
            else:
                label = '';
            return label;
        for x in self._mixture_list:
            all_groups.append(get_label(x));    
        all_groups = set(all_groups);
        def mytextfunction(w):
            if len(w['compound_list'])>0:
                return str(len(w['compound_list']))
            else:
                return ''
        myfacecolors = dict(zip(all_groups,plt.cm.tab20(np.linspace(0,1,len(all_groups)))))
        myfacecolors['']='none'
        myfacecolors['[]']='none'
        self.graphical_print(subtitle=self._description+': color by id (filtered by type: '+str(typefilter)+')',
                     textfunction=mytextfunction,
                     colorfunction=lambda w: get_label(w),
                     facecolormap=myfacecolors
                     )
        print('colors correspond to the following groups: ',all_groups)
                
                
    def get_compounds(self,location,properties):
        if type(properties) is str:
            return self.__getitem__(location)['compound_list'].list_all(properties)
        if type(properties) is list:
            return list(zip(*[self.get_compounds(location,p) for p in properties]))
                
    def get_positions_by_compound_type(self,type_filter,grid=None):        
        positions = [];
        for row in range(0,self._rows):
            for col in range(0,self._cols):
                compound_info = self.get_compounds((row,col),['type']); 
                compound_types = [item for compound_type in compound_info for item in compound_type]; # flatten list of compound types
                if compound_types != []:
                    compound_position = position_to_grid((row,col),grid=grid);  # return position in tuple integer position form
                    #compound_position = position_to_lettergrid((row,col)); # changed on 2018/11/16 C.Arcadia                                    
                    if type_filter in compound_types:                        
                            positions.append(compound_position);
        return positions;
        
    def get_contents(self,locations=[],properties=['mass','name','cid','type'],type_filter=None,grid=None): # ['mass','name','cid','type','volume']
        # formats the contents of a container into a dictionary of lists containing info about the compounds (which can be type filtered) at each of the specified positions in a source plate
        if not isinstance(locations,(list,tuple)):
            locations = [locations];
        if locations == []: # if not given locations assume we want the entire container
            locations = self.list_filled_positions();                   
        contents = dict();
        contents['position'] = [];
        contents['well_volume'] = [];
        for prop in properties:
            contents[prop] = [];
        for location in locations:
            compound_type = self.get_compounds(location,['type']);             
            compound_info = self.get_compounds(location,properties); 
            well_volume = self.get_volume(location);
            for m in range(0,len(compound_info)): 
                if (type_filter is None) or (type_filter == compound_type[m][0]):
                    for n in range(0,len(properties)):
                        value = compound_info[m][n]; 
                        prop = properties[n];
                        contents[prop].append(value); 
                    contents['position'].append(position_to_grid(location,grid=grid));
                    contents['well_volume'].append(well_volume);                    
        return contents
        
 
    

# ======================================================================


class WellPlate1536LDV(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=32
        kwargs['cols']=48
        kwargs['volume_max'] = 5500 #nL
        kwargs['volume_min'] = 1000 #nL
        Container.__init__(self,*args,**kwargs)

class WellPlate384LDV(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=16
        kwargs['cols']=24
        kwargs['volume_max'] = 12000 #nL
        kwargs['volume_min'] = 2500 #nL
        Container.__init__(self,*args,**kwargs)

class WellPlate384PP(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=16
        kwargs['cols']=24
        kwargs['volume_max'] = 65000 #nL
        kwargs['volume_min'] = 15000 #nL
        Container.__init__(self,*args,**kwargs)

class MaldiPlate1536(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=32
        kwargs['cols']=48
        kwargs['volume_max'] = 500 #nL
        kwargs['volume_min'] = 0 #nL
        Container.__init__(self,*args,**kwargs)

class MaldiPlate384(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=16
        kwargs['cols']=24
        kwargs['volume_max'] = 1000 #nL
        kwargs['volume_min'] = 0 #nL
        Container.__init__(self,*args,**kwargs)
        
class WellPlate384(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=16
        kwargs['cols']=24
        kwargs['volume_max'] = 120000 #nL
        kwargs['volume_min'] = 10000 #nL
        Container.__init__(self,*args,**kwargs)
        
class WellPlate96(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=8
        kwargs['cols']=12
        kwargs['volume_max'] = 400000 #nL
        kwargs['volume_min'] = 100000 #nL
        Container.__init__(self,*args,**kwargs)

class WellPlate24(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=4
        kwargs['cols']=6
        kwargs['volume_max'] = 3400000 #nL
        kwargs['volume_min'] = 450000 #nL
        Container.__init__(self,*args,**kwargs)

class microtube15(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=kwargs.get('rows',1)
        kwargs['cols']=kwargs.get('cols',1)
        kwargs['volume_max'] = 1500000 #nL
        kwargs['volume_min'] = 50000 #nL
        Container.__init__(self,*args,**kwargs)

class microtube20(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=kwargs.get('rows',1)
        kwargs['cols']=kwargs.get('cols',1)
        kwargs['volume_max'] = 2000000 #nL
        kwargs['volume_min'] = 70000 #nL
        Container.__init__(self,*args,**kwargs)

class tube15(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=kwargs.get('rows',1)
        kwargs['cols']=kwargs.get('cols',1)
        kwargs['volume_max'] = 15000000 #nL
        kwargs['volume_min'] = 250000 #nL
        Container.__init__(self,*args,**kwargs)

class tube50(Container):
    def __init__(self,*args,**kwargs):
        kwargs['rows']=kwargs.get('rows',1)
        kwargs['cols']=kwargs.get('cols',1)
        kwargs['volume_max'] = 50000000 #nL
        kwargs['volume_min'] = 500000 #nL
        Container.__init__(self,*args,**kwargs)
        

        
 # ======================================================================

# Note: many of the methods below were copied from "misc.py" since they belong in a file that is dedicated to container related operations

# position converters 

rowletters = [chr(x) for x in range(ord('A'),ord('Z')+1)] + ['A' + chr(x) for x in range(ord('A'),ord('Z')+1)];    

def position_to_lettergrid(position):
    """ Returns a string corresponding to the numerical position. Zero indexed.
    Rows are letters, columns are numbers. (0,0) returns 'A1'. (1,3) returns 'B4'.
    """
    if type(position)==tuple:
        return rowletters[position[0]]+str(position[1]+1)
    if type(position)==list:
        return [position_to_lettergrid(x) for x in position]
    # not recognized. assume it does not need to be converted, return as-is
    return position

def lettergrid_to_position(key):
    """ Returns a zero indexed numerical position tuple.
    Rows are letters, columns are numbers. 'A1' returns (0,0). 'B4' returns (1,3).
    """        
    if type(key)==str:
        getletters = ''.join([x for x in key if not x.isdigit()])
        getnumbers = ''.join([x for x in key if x.isdigit()])
        return ( rowletters.index(getletters), int(getnumbers)-1 )
    if type(key)==list:
        return [lettergrid_to_position(x) for x in key]
    # not recognized. assume it does not need to be converted, return as-is
    return key

def position_to_maldigrid(position):
    """ Returns a string corresponding to the MALDI position. Zero indexed.
    Rows are letters, columns are numbers. (0,0) returns 'X01Y01'. (1,3) returns 'X02Y04'.
    """
    if type(position)==tuple:
        return 'X%02dY%02d' % (position[1]+1,position[0]+1)
    if type(position)==list:
        return [position_to_maldigrid(x) for x in position]
    # not recognized. assume it does not need to be converted, return as-is
    return position

def maldigrid_to_position(key):
    """ Returns a zero indexed numerical position tuple.
    Rows are Y numbers, columns are X numbers. 'X01Y01' returns (0,0). 'X02Y04' returns (1,3).
    """        
    if type(key)==str:
        Xind = key.find('X');
        Yind = key.find('Y');
        if Xind<Yind:
            X = int(key[Xind+1:Yind]) - 1;
            Y = int(key[Yind+1:]) - 1;
        else:
            Y = int(key[Yind+1:Xind]) - 1;
            X = int(key[Xind+1:]) - 1;
        return ( Y, X )
    if type(key)==list:
        return [lettergrid_to_position(x) for x in key]
    # not recognized. assume it does not need to be converted, return as-is
    return key

def position_to_grid(position,grid=None):
    if grid=='letter':
        key = position_to_lettergrid(position);                        
    elif grid=='maldi':
        key = position_to_maldigrid(position);
    else:
        key = position; # assume input is in desired form (which can be in position form)
    return key

def grid_to_position(key,grid=None):
    if grid=='letter':
        position = lettergrid_to_position(key);                        
    elif grid=='maldi':
        position = maldigrid_to_position(key);
    else:
        position = key; # assume input is in desired form (which can be in grid/key form)
    return position
      
# -------------------------------------- 
    
def bounds_to_positions(bounds):
    """ Converts a list of bounds into a list of explicit positions.
    The bounds can be a list of single well positions or letter grids,
    or a list of (upper-left,lower-right) tuple pairs.
    
    """
    myposlist = []
    for b in bounds:
        if type(b) is str:   
            # one lettergrid
            mypos = lettergrid_to_position(b)
            includerows = [mypos[0],]
            includecols = [mypos[1],]
        elif type(b) is tuple and type(b[0]) is int:
            # one position
            includerows = [mypos[0],]
            includecols = [mypos[1],]            
        elif type(b) is tuple and type(b[0]) is str: 
            # rectanglar selection of lettergrids
            myposA = lettergrid_to_position(b[0])
            myposB = lettergrid_to_position(b[1])
            includerows = range(myposA[0],myposB[0]+1);
            includecols = range(myposA[1],myposB[1]+1);
        elif type(b) is tuple and type(b[0]) is tuple: 
            # rectanglar selection of positions
            myposA = b[0]
            myposB = b[1]
            includerows = range(myposA[0],myposB[0]+1);
            includecols = range(myposA[1],myposB[1]+1);
        else:                    
            print(type(b),type(b[0]),b)
            raise Exception('unsupported well-plate bounds')
        
        for row in includerows:
            for col in includecols:
                myposlist.append( (row,col) )
    
    return myposlist

# --------------------------------------

def rectangle_positions(topleft,botright):
    return bounds_to_positions( [ (topleft,botright) ] )

def rotate_positions_180(positions,plate):
    rows,cols = plate.shape
    return [(rows-1-x,cols-1-y) for x,y in lettergrid_to_position(positions)]
    
# --------------------------------------

def unique_rows(positions):
    positions = lettergrid_to_position(positions)
    return set([r for r,c in positions])

def unique_cols(positions):
    positions = lettergrid_to_position(positions)
    return set([c for r,c in positions])

def isolate_row(positions,row):
    positions = lettergrid_to_position(positions)
    return [(r,c) for r,c in positions if r==row]

def isolate_col(positions,col):
    positions = lettergrid_to_position(positions)
    return [(r,c) for r,c in positions if c==col]

def separate_by_row(positions):
    positions = lettergrid_to_position(positions)
    return [isolate_row(positions,r) for r in unique_rows(positions)]
    
def separate_by_col(positions):
    positions = lettergrid_to_position(positions)
    return [isolate_col(positions,c) for c in unique_cols(positions)]

def separate_by_count(positions,count):
    positions = lettergrid_to_position(positions)
    return [positions[i:i+count] for i in range(0,len(positions),count)]