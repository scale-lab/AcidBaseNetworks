
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 19:01:33 2018

@author: jacobrosenstein
"""


from chemcpupy import MixtureList,CompoundList
from chemcpupy.tools.misc import *
import numpy as np
import copy
import matplotlib.pyplot as plt



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
                
                plt.plot(col,row,
                         markermap.get(colorindex,'o'),
                         markersize=24,
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
        plt.show()
        return plt

        
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
        
    def get_compounds(self,location,properties):
        if type(properties) is str:
            return self.__getitem__(location)['compound_list'].list_all(properties)
        if type(properties) is list:
            return list(zip(*[self.get_compounds(location,p) for p in properties]))

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
        kwargs['rows']=16
        kwargs['cols']=24
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

    