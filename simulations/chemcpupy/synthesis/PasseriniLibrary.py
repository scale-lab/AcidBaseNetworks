# -*- coding: utf-8 -*-
"""
Created on Tues Sep 11 13:32:12 2018

@author: jacobrosenstein & chrisarcadia

A PasseriniLibrary stores a list of Passerini molecules. It inherits from the CompoundList class. 
"""


import json
import csv
from datetime import datetime

import pubchempy as pcp
from pyteomics import mass

from chemcpupy import CompoundList, Containers
#from chemcpupy.tools.misc import *


class PasseriniLibrary(CompoundList):
    
    # ----------------------------------------------
    def __init__(self, *args, **kwargs):
        """An PasseriniLibrary stores a list of Passerini molecules. It inherits from the
        CompoundList class. 
    
        Args:
           None.
    
        Kwargs:
            description (string): Optional. Text description of the library.
    
        """        
        
        
        self._description = kwargs.pop('description','(no description)')
        self._file_subtype='PasseriniLibrary-1.0'
        CompoundList.__init__(self,*args,**kwargs)
        
    def add_reagents(self, *args, **kwargs):
        """Add reagents for Passerini reactions.
    
        Args:
           None.
    
        Kwargs:
            all_components (CompoundList):  A list of all components, identified by 'type'.
            isocyanides (CompoundList):   A list of isocyanides.
            aldehydes(CompoundList):   A list of aldehydes.
            carboxylic_acids (CompoundList):   A list of carboxylic acids.
            
        Returns:
            None.
            
        Raises:
           None.
    
        Example uses:
    
            TBD    
        """        
        if 'all_components' in kwargs:
            c = kwargs.get('all_components')
            self.isocyanides = c.get_subset('type','isocyanide')
            self.aldehydes = c.get_subset('type','aldehyde')
            self.carboxylic_acids = c.get_subset('type','carboxylic acid')
            
        else:
            self.isocyanides = kwargs.get('isocyanides',CompoundList())
            self.aldehydes = kwargs.get('aldehydes',CompoundList())
            self.carboxylic_acids = kwargs.get('carboxylic_acids',CompoundList())
        self._build_params = kwargs.get('build_params',{})


    def synthesize_library(self,verbose=False):
        """Create all possible Passerini molecules from the reagents.
       
        Kwargs:
           Verbose (bool): Default=False. Include detailed annotations of the 4 components for each Passerini compound. 
       
        """          
        
        for iso in self.isocyanides:
            for ald in self.aldehydes:
                for car in self.carboxylic_acids:

                    passerini_composition = mass.Composition(formula=''.join(
                                                        (iso['formula'],
                                                         ald['formula'],
                                                         car['formula']))
                                                        )

                    passerini_mass = mass.calculate_mass(composition=passerini_composition)
                    
                    self._compound_list.append( {'type':'passerini',
                                                 'composition':passerini_composition,
                                                 'mass':passerini_mass,
                                                 'cid':0, # indicate not in PubChem database
                                                 'name':'Passerini Product (Ald:'+str(ald['cid'])+'-'+'Car:'+str(car['cid'])+'-'+'Iso:'+str(iso['cid'])+')',                                                 
                                                 'cid_isocyanide':iso['cid'],
                                                 'cid_aldehyde':ald['cid'],
                                                 'cid_carboxylic_acid':car['cid'],
                                                 })

                    
                    if verbose:
                        print(len(self),
                              iso['formula'],
                              ald['formula'],
                              car['formula'],
                              passerini_mass)

                    for newparam in self._build_params:
                        self._compound_list[-1].update({newparam:''})
                        for subparam in self._build_params[newparam]:
                            if subparam[0]=='isocyanide':
                                self._compound_list[-1][newparam] += iso[subparam[1]]
                            if subparam[0]=='aldehyde':
                                self._compound_list[-1][newparam] += ald[subparam[1]]
                            if subparam[0]=='carboxylic acid':
                                self._compound_list[-1][newparam] += car[subparam[1]]
                                
                        
                        
                        

    def same_compound(self,x,y):
        """ Returns true of the two compounds are the same.
        A match can be from: 
           - all 4 Passerini component CIDs matching
        """
                        
        if ((x.get('cid_aldehyde',-1)==y.get('cid_aldehyde')) 
            and (x.get('cid_carboxylic_acid',-1)==y.get('cid_carboxylic_acid'))
            and (x.get('cid_isocyanide',-1)==y.get('cid_isocyanide'))):
            return True

        return CompoundList.same_compound(self,x,y)




    def generate_passerini_products(my_well_plate,reaction_yield=0.95):
        for w in my_well_plate:
            c = w['compound_list']
            for iso in c.get_subset('type','isocyanide'):
                for ald in c.get_subset('type','aldehyde'):
                    for carb in c.get_subset('type','carboxylic acid'):
                        my_passerini_lib = PasseriniLibrary(description=w['mixture_name'])
                        my_passerini_lib.add_reagents(isocyanides=c.get_subset('type','isocyanide'),
                                                aldehydes=c.get_subset('type','aldehyde'),
                                                carboxylic_acids=c.get_subset('type','carboxylic acid')
                                                )
                        my_passerini_lib.synthesize_library()
                        my_passerini_lib.calculate_isotopes()
                        
                        myvol=w['volume']
                        my_well_plate.add_compounds_to_location(w['position'],
                                                                c,
                                                                volume = -reaction_yield*myvol)
                        
                        my_well_plate.add_compounds_to_location(w['position'],
                                                                my_passerini_lib,
                                                                volume=reaction_yield*myvol)
            c.remove_if_zero_volume()
                            

    def list_passerini_products(my_well_plate,bounds=None):
        my_passerinis = CompoundList()
        
        if bounds is None:
            for w in my_well_plate:
                more_passerinis = w['compound_list'].get_subset('type','passerini')
                my_passerinis.add_compounds(compound_list=more_passerinis,
                                      other_properties=[ ('position',w['position']),
                                                         ('location',Containers.position_to_lettergrid(w['position'])) 
                                                         ])
        else:
            myposlist = Containers.bounds_to_positions(bounds)
            
            for row,col in myposlist:
                more_passerinis = my_well_plate[row,col]['compound_list'].get_subset('type','passerini')
                my_passerinis.add_compounds(compound_list=more_passerinis,
                                      other_properties=[ ('position',(row,col)),
                                                         ('location',Containers.position_to_lettergrid((row,col))) 
                                                         ])
    
        return my_passerinis



# EXAMPLES OF HOW TO USE THIS CLASS
if __name__ == '__main__':
    
    print('Creating an Passerini library.')
    
    passerini_components = CompoundList(
                            description='Reagents for 102 Passerini Demo, Feb 2018',
                            name_csv_file='examples/demo1_ugi_reagents1.csv',
                            csv_col=2,
                            csv_properties=( ('type',0),('common name',1) ) 
                            )    
    passerini_components.auto_fill_properties_from_pubchem()
    
    my_passerini_lib = PasseriniLibrary(description='my first Passerini library')
    my_passerini_lib.add_reagents(isocyanides=passerini_components.get_subset('type','isocyanide'),
                            aldehydes=passerini_components.get_subset('type','aldehyde'),
                            carboxylic_acids=passerini_components.get_subset('type','carboxylic acid')
                            )
    my_passerini_lib.synthesize_library()
    print(str(my_passerini_lib))
    
    