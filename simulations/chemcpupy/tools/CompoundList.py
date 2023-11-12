# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 10:52:39 2018

@author: jacobrosenstein

This is a basic class for storing a list of compounds 
with their IDs (PubChem CID, InChI Key) and 
properties (mass, formula, InChI, and custom properties).

Requires pubchempy and pyteomics.
"""

import json
import csv
from datetime import datetime
import os

import pubchempy as pcp
from pyteomics import mass

from pprint import pprint


class CompoundList:
    
    # ----------------------------------------------
    def __init__(self, *args, **kwargs):
        """A CompoundList stores a list of compounds with their properties
        and unique identifiers like PubChem CID and InChI Key.
    
        Args:
           None.
    
        Kwargs:
           cid_list (list): PubChem CIDs to add to CompoundList.
           inchikey_list (list): InChI Keys to add to the CompoundList.
           compound_json_file (string): An existing CompoundList .JSON file to add.
           cid_csv_file (string): A CSV file with CIDs in the first column.
           name_csv_file (string): (Not recommended.) A CSV file with compound names or CAS numbers in the first column.
    
        Returns:
            None.
            
        Raises:
           None.
    
        Example uses:
    
        mylist1 = CompoundList(cid_list=(11160305,97233,68250,))
        mylist2 = CompoundList(compound_json_file='isocyanides.json')
    
        """        
        self._compound_list = []
        self._description = kwargs.pop('description','(no description)')
        self._file_type='CompoundList-1.0'
        self._file_subtype=''
        
        # add compounds from keyword-specified sources
        self.add_compounds(**kwargs)
    
    
    # ----------------------------------------------
    def set_description(self,description):
        """Sets a text description of the CompoundList which will be stored
        in the JSON file. Use this to describe in English what the list is for.
    
        Args:
           description (string): Plain text description of the compound list.
    
        Kwargs:
            None.
            
        Returns:
            None.
            
        Raises:
           None.
    
        Example use:
    
        mylist1.set_description('This is a list of BOC-protected amino acids
                                 to be used as carboxylic acids in Ugi synthesis.')
    
        """         
        self._description=description        

    # ----------------------------------------------
    def add_compounds(self, *args, **kwargs):
        """Adds compounds to an existing CompoundList.
    
        Args:
           None.
    
        Kwargs:
           compound_list (list): A direct list of compounds to include.
           cid_list (list): PubChem CIDs to add to CompoundList.
           inchikey_list (list): InChI Keys to add to the CompoundList.
           compound_json_file (string): An existing CompoundList .JSON file to add.
           cid_csv_file (string): A CSV file with CIDs in the first column.
           name_csv_file (string): (Not recommended.) A CSV file with compound names or CAS numbers in the first column.
           csv_col (int):  Optional. The CSV file column number to use. Default = 0 (first column)
           csv_properties (list):  Optional. A list of pairs of (property_name,csv_column)
        Returns:
            None.
            
        Raises:
           None.
    
        Example uses:
    
        mylist1.add_compounds(cid_list=(11160305,97233,68250,))
        mylist2.add_compounds(compound_json_file='isocyanides.json')
    
        """          
        _compound_list = kwargs.get('compound_list',None)
        _cid_list = kwargs.get('cid_list',None)    
        _inchikey_list = kwargs.get('inchikey_list',None)    
        _compoundjsonfile = kwargs.get('compound_json_file',None)
        _cid_csv_file = kwargs.get('cid_csv_file',None)    
        _name_csv_file = kwargs.get('name_csv_file',None)        
        _csv_col = kwargs.get('csv_col',0)
        _csv_properties = kwargs.get('csv_properties',[])
        _total_volume = kwargs.get('total_volume',0)
        _other_properties = kwargs.get('other_properties',[])

        if _compound_list is not None:
            for c in _compound_list:
                found = self._find(c)
                if found:
                    self._compound_list[found[0]]['volume'] = self._compound_list[found[0]].get('volume',0)+\
                                                                _total_volume/len(_compound_list)
                    for p,v in _other_properties:
                        self._compound_list[found[0]][p]=v
                else:
                    self._compound_list.append(c)
                    self._compound_list[-1]['volume'] = _total_volume/len(_compound_list)
                    for p,v in _other_properties:
                        self._compound_list[-1][p]=v
            
        if _cid_list is not None:
            # add provided CIDs
            for cid in _cid_list:
                self._compound_list.append( { 'cid' : cid } )
                self._compound_list[-1]['volume'] = _total_volume/len(_cid_list)
                
        if _inchikey_list is not None:
            # add provided InChIKeys
            for inchikey in _inchikey_list:
                self._compound_list.append( { 'inchikey' : inchikey } )
                self._compound_list[-1]['volume'] = _total_volume/len(_inchi_list)
                                
        if _compoundjsonfile is not None:
            print('Loading compounds from ',_compoundjsonfile)
            with open(_compoundjsonfile,'r') as f:
                filecontents = json.load(f)
                self._compound_list.extend(filecontents['compound_list'])
                if self._description == '(no description)':
                    self._description = ''
                self._description = self._description + filecontents['description'] + ' / '

        if _cid_csv_file is not None:
            # add a list of CIDs from a CSV file. 
            my_reader = csv.reader(open(_cid_csv_file, newline=''), 
                                   delimiter=',')
            for row in my_reader:
                self._compound_list.append( { 'cid' : row[_csv_col] } )
                self._compound_list[-1]['volume'] = _total_volume
                for p in _csv_properties:
                    self._compound_list[-1].update( { p[0] : row[p[1]] } )
        
        if _name_csv_file is not None:
            # add a list of named compounds from a CSV file. 
            my_reader = csv.reader(open(_name_csv_file, newline=''), 
                                   delimiter=',')
            for row in my_reader:
                self._compound_list.append( { 'name' : row[_csv_col] } )
                self._compound_list[-1]['volume'] = _total_volume
                for p in _csv_properties:
                    self._compound_list[-1].update( { p[0] : row[p[1]] } )
        
    # ----------------------------------------------        
    def add_properties_by_cid(self, cids_and_properties):
        """Adds properties to compounds in a CompoundList. These can be any
        properties you want to store. You can choose new property names.
    
        Args:
           cids_and_properties (dict): This is a dictionary indexed by CID.
                                       The contents are themselves dictionaries 
                                       associating property names with values.
                                       See example.
    
        Example uses:
    
            mylist1.add_properties_by_cid(
                {11101:{'mass':1234.5678,'custom_label':'something special'}},
                 83589:{'supplier':'Aldrich'}})   
    
        """          
        for cid in cids_and_properties.keys():
            for compound in self._compound_list:
                if compound.get('cid')==cid:
                    compound.update(cids_and_properties[cid])


    # ----------------------------------------------        
    def save_json(self,tofile):
        """Saves the CompoundList to a .JSON file.
    
        Args:
           tofile (string): The desired .JSON file name to write.
    
        Example use:
    
            mylist1.save_json('testjson1.json')    
    
        """          
        print('Saving compounds to JSON file: ',tofile)
        myinfo = {'compound_list':self._compound_list,
                  'created':str(datetime.now()),
                  'file type':self._file_type,
                  'file subtype':self._file_subtype,
                  'description':str(self._description),                  
                  }        
        with open(tofile, 'w') as f:
            json.dump(myinfo, f, indent=2)
    
    # ----------------------------------------------    
    # use pubchem & pyteomics to fill out whole properties list
    def auto_fill_properties_from_pubchem(self,overwrite=False):
        """Retrieves properties from PubChem to fill out all compounds 
        where possible.
    
        Kwargs:
           overwrite (bool): If overwrite=True, existing properties will be overwritten.
    
        Example use:
    
            mylist1=CompoundList(cid_list=(11160305,97233,68250,))
            mylist1.auto_fill_properties_from_pubchem(overwrite=False)
    
        """         

        for i in range(len(self._compound_list)):
            
            m = None
            
            if 'cid' in self._compound_list[i].keys():
                print('Fetch from PubChem (CID)',self._compound_list[i]['cid'])
                m = pcp.Compound.from_cid(self._compound_list[i]['cid'])
                                            
            elif 'inchikey' in self._compound_list[i].keys():
                print('Fetch from PubChem (InChIKey)',self._compound_list[i]['inchikey'])
                c = pcp.get_compounds(self._compound_list[i]['inchikey'], 
                                      'inchikey')
                if len(c)>0:
                    m = c[0]
                else:
                    print('WARNING: no PubChem match found.')

            elif 'name' in self._compound_list[i].keys():
                print('Fetch from PubChem (name)',self._compound_list[i]['name'])
                c = pcp.get_compounds(self._compound_list[i]['name'], 
                                      'name')
                if len(c)>0:
                    m = c[0]
                else:
                    print('WARNING: no PubChem match found.')


            if m is not None:                    
                if overwrite or ('formula' not in self._compound_list[i].keys()):
                    self._compound_list[i]['formula'] = m.molecular_formula
                    
                if overwrite or ('mass' not in self._compound_list[i].keys()):
                    self._compound_list[i]['mass'] = \
                            mass.calculate_mass(formula=m.molecular_formula)
    
                if overwrite or ('inchikey' not in self._compound_list[i].keys()):
                    self._compound_list[i]['inchikey'] = m.inchikey

                if overwrite or ('inchi' not in self._compound_list[i].keys()):
                    self._compound_list[i]['inchi'] = m.inchi

                if overwrite or ('cid' not in self._compound_list[i].keys()):
                    self._compound_list[i]['cid'] = m.cid
                
                if overwrite or ('smiles' not in self._compound_list[i].keys()):
                    self._compound_list[i]['canonical_smiles'] = m.canonical_smiles


    # ----------------------------------------------
    def calculate_compositions(self,overwrite=False):
        """Calculates atomic compositions. Requires 'formula' property.
    
        Kwargs:
           overwrite (bool): Overwrite existing compositions.
    
        Example use:
    
            mylist1=CompoundList(cid_list=(11160305,97233,68250,))
            mylist1.auto_fill_properties_from_pubchem()
            mylist1.calculate_compositions()
    
        """    
        for k in range(len(self)):
            if overwrite or ('composition' not in self._compound_list[k].keys()):
                if ('formula' in self._compound_list[k].keys()):
                    self._compound_list[k]['composition']=mass.Composition(formula=self._compound_list[k]['formula'])                            
                else:
                    print('WARNING: formula not included. Skipping composition calculation.')
            
    # ----------------------------------------------
    def calculate_isotopes(self,isotope_threshold=0.005,overall_threshold=0.01,verbose=False):
        """Calculates expected isotopic masses. Requires 'composition' property.
    
        Kwargs:
           isotope_threshold (float): Threshold for isotopic abundance of single atoms.
           overall_threshold (float): Threshold for isotopic abundance of overall composition.
    
        Example use:
    
            mylist1=CompoundList(cid_list=(11160305,97233,68250,))
            mylist1.auto_fill_properties_from_pubchem()
            mylist1.calculate_compositions()
            mylist1.calculate_isotopes()
    
        """    
        for k in range(len(self)):
            if verbose:
                print('%d/%d: calculate isotopes' % (k,len(self)))
            if ('composition' in self._compound_list[k].keys()):
                composition = self._compound_list[k]['composition']
                isotope_masses = []
                isotope_abundance = []
                for iso_comp,iso_abund in mass.isotopologues(composition=composition,
                                                report_abundance=True,
                                                isotope_threshold=isotope_threshold,
                                                overall_threshold=overall_threshold):
                    if verbose:
                        print('  ',mass.calculate_mass(composition=iso_comp),iso_abund)
                    isotope_masses.append(mass.calculate_mass(composition=iso_comp))
                    isotope_abundance.append(iso_abund)
                self._compound_list[k]['isotope_masses']=isotope_masses
                self._compound_list[k]['isotope_abundance']=isotope_abundance
            else:
                print('WARNING: composition not included. Skipping isotope calculation.')
            

    # ----------------------------------------------
    def calculate_ion_mz(self):
        """Calculates m/z of expected ions M+, MH+, MNa+, MK+. Requires 'composition' property.
    
        Kwargs:
           None.
    
        Example use:
    
    
        """    
        for k in range(len(self)):
            if ('composition' in self._compound_list[k].keys()):
                composition = self._compound_list[k]['composition']
                
                self._compound_list[k]['mz'] = self._compound_list[k]['mass']
                
                self._compound_list[k]['mz_H+'] = mass.calculate_mass(
                        composition=(composition+mass.Composition(formula='H')))
                
                self._compound_list[k]['mz_Na+'] = mass.calculate_mass(
                        composition=(composition+mass.Composition(formula='Na')))
                
                self._compound_list[k]['mz_K+'] = mass.calculate_mass(
                        composition=(composition+mass.Composition(formula='K')))
                
            else:
                print('WARNING: composition not included. Skipping ion m/z calculation.')
            

    # ----------------------------------------------
    def __str__(self):
        """The string of a CompoundList includes its length and description.
    
        Example use:
    
            >> mylist1=CompoundList(cid_list=(11160305,97233,68250,))
            >> mylist1.set_description('This is an example list.')
            >> print(mylist1)
    
        """        
        mystr = 'CompoundList (length=%d) ' % (len(self._compound_list))
        mystr = mystr + str(self._description)
        return mystr
    
    # ----------------------------------------------
    def __len__(self):
        """The length of a CompoundList is the number of compounds it contains.        
        """        
        return len(self._compound_list)

    def __getitem__(self, key):
        """You can retrieve compounds from a CompoundList by their CID or InChIKey.
        
        Args:
            key (int or str): The CID or InChIKey of the compound.
            
        Example use:
            >> x = my_compound_list[68250]
            >> y = my_compound_list['CMWKITSNTDAEDT-UHFFFAOYSA-N']            
            
        """
        if type(key) is int:
            return [x for x in self._compound_list if x['cid']==key][0]
        elif type(key) is str:
            return [x for x in self._compound_list if x['inchikey']==key][0]            
        else:
            raise TypeError

    def __iter__(self):
        self._iter_index = 0
        return self

    def __next__(self):
        if self._iter_index >= len(self):
            raise StopIteration
        else:
            self._iter_index += 1
            return self._compound_list[self._iter_index-1]
            
        
    def __contains__(self,item):
        """ A CompoundList contains another CompoundList if 
        there are matches for all of its elements.
        """
        for x in item:
            findmatch = [y for y in self._compound_list if self.same_compound(x,y)]
            if len(findmatch)==0:
                return False            
        return True

    def _find(self,item):
        """ Returns the _compound_list indices which match another compound.
        """            
        return [i for i,c in enumerate(self._compound_list) if self.same_compound(c,item)]

        
        
    def same_compound(self,x,y):
        """ Returns true of the two compounds are the same.
        A match requires that: 
            - all cid_xxxx parameters match
            - all inchikey_xxxxx parameters match

        Inherited classes of CompoundList may extend this function to other criteria, if
        CIDs or InChIKeys are not available.
        """
        
        for k in y.keys():            
            if 'cid' in k and x.get(k,'n/a')!=y.get(k):
                return False                    
            if 'inchikey' in k and x.get(k,'n/a')!=y.get(k):
                return False
        
        return True
        
        

    def set_concentration(self,cid,concentration):
        """ Set the concentration of one chemical in the list.
        """
        self.add_properties_by_cid( {cid:{'concentration':concentration}} )

        
    # ----------------------------------------------
    def describe_properties(self,properties=None):
        """Produces a detailed enumeration of the compounds, and a
        given list of properties. 
    
        Kwargs:
            properties (list):  A list of properties to describe. If no list is
                                provided, then all properties will be included.
        
        Returns:
            A string listing all of the compounds and their properties.
                
        Example use:
    
            >> mylist1 = CompoundList(compound_json_file='isocyanides.json')
            >> print(mylist1.describe_properties(('cid','mass')))
        """        
        
        mystr = ''
        for c in self._compound_list:
            mystr = mystr + 'Compound: \n'
            
            if properties is not None:
                for p in properties:
                    mystr = mystr + ' - ' + p + ' = ' 
                    mystr = mystr + str(c.get(p,'[None]')) + '\n'
                    
            else:
                for p in c.keys():
                    mystr = mystr + ' - ' + p + ' = ' + str(c[p]) + '\n'                    
                    
                    
        return mystr

    # ----------------------------------------------
    def list_all(self,prop='cid'):
        """Returns a list of the specified property value.
    
        Args:
            prop (string):  A compound property. e.g. 'cid','mass','inchi',etc.
        
        Returns:
            A list of the values of the given property for all compounds.
                
        Example use:
    
            >> mymasses = mylist2.list_all('mass')
        """        
        
        return [x.get(prop,None) for x in self]
    
     # ----------------------------------------------
    def get_subset(self,prop,val):
        """Returns a subset of the CompoundList which matches the given property and value.
    
        Args:
            prop (string):  property name to match.
            val (anything):  property value to match
        
        Returns:
            A new CompoundList with only the matching compounds
                
        Example use:
    
            >> aldehydes = ugi_components.get_subset('type','aldehyde')
        """        
        
        new_compound_list =  CompoundList()
        new_compound_list._compound_list.extend([x for x in self._compound_list if x[prop]==val])
        return new_compound_list


    def get_absent_compounds(self,present_compounds):
        """Returns the absent subset of the CompoundList, compared to the provided CompoundList
        """
        new_compound_list =  CompoundList()
        new_compound_list._compound_list.extend([x for x in self._compound_list if (x,) not in present_compounds])
        return new_compound_list


    def remove_if_zero_volume(self):
        self._compound_list = [x for x in self._compound_list if x['volume']>0]



# EXAMPLES OF HOW TO USE THIS CLASS
if __name__ == '__main__':
    
    print('Creating a compound list.')
    mylist1=CompoundList(cid_list=(11160305,97233,68250,))
    
    mylist1.set_description('This is a test compound list to demonstrate the class.')

    print('Adding other compounds to the list.')
    mylist1.add_compounds(cid_list=(11101,83589))
    mylist1.add_compounds(inchikey_list=('ZGEAYXBIJAYBKA-UHFFFAOYSA-N',))    
    
    print('Adding custom properties to the list.')
    mylist1.add_properties_by_cid({11101:{'mass':1234.5678,'custom_label':'something special'}})   

    print('Auto-fill properties from PubChem.')
    mylist1.auto_fill_properties_from_pubchem(overwrite=False)

    print('Save the list to a .json file.')    
    mylist1.save_json('testjson1.json')    
    
    print('Load the .json file into a new compound list, and add another compound')
    mylist2 = CompoundList(compound_json_file='testjson1.json')
    mylist2.add_compounds(cid_list=(16211045,7504,))    
    mylist2.auto_fill_properties_from_pubchem(overwrite=False)

    print('Summary')
    print(mylist2)
    
    print('List specific properties')
    print(mylist2.describe_properties(properties=('cid','mass')))
    
    print('List all properties')
    print(mylist2.describe_properties())
    
    print('List of the masses')
    print(mylist2.list_all('mass'))
    
def load_compounds_from_file(title, csv_file, json_file):

    if os.path.isfile(json_file):
        compound_list = CompoundList(compound_json_file=json_file)
    else:
        # Create a list of compounds from a CSV file with CAS Numbers and other labels
        compound_list = CompoundList(
                            description=title,
                            name_csv_file=csv_file,
                            csv_col=2,
                            csv_properties=( ('type',0),('common name',1),('location',3) ) 
                            )
        # assumes CSV is formated with the following columns: type, common name, location
        compound_list.auto_fill_properties_from_pubchem()
    return compound_list
