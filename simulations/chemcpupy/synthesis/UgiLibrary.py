# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 20:26:05 2018

@author: jacobrosenstein

An UgiLibrary stores a list of Ugi molecules. It inherits from the CompoundList class. 
"""


import json
import csv
from datetime import datetime

import pubchempy as pcp
from pyteomics import mass

from chemcpupy import CompoundList, Containers
#from chemcpupy.tools.misc import *


class UgiLibrary(CompoundList):
    
    # ----------------------------------------------
    def __init__(self, *args, **kwargs):
        """An UgiLibrary stores a list of Ugi molecules. It inherits from the
        CompoundList class. 
    
        Args:
           None.
    
        Kwargs:
            description (string): Optional. Text description of the library.
    
        """        
        
        
        self._description = kwargs.pop('description','(no description)')
        self._file_subtype='UgiLibrary-1.0'
        CompoundList.__init__(self,*args,**kwargs)
        
    def add_reagents(self, *args, **kwargs):
        """Add reagents for Ugi reactions.
    
        Args:
           None.
    
        Kwargs:
            all_components (CompoundList):  A list of all components, identified by 'type'.
            isocyanides (CompoundList):   A list of isocyanides.
            amines (CompoundList):   A list of amines.
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
            self.amines = c.get_subset('type','amine')
            self.aldehydes = c.get_subset('type','aldehyde')
            self.carboxylic_acids = c.get_subset('type','carboxylic acid')
            
        else:
            self.isocyanides = kwargs.get('isocyanides',CompoundList())
            self.amines = kwargs.get('amines',CompoundList())
            self.aldehydes = kwargs.get('aldehydes',CompoundList())
            self.carboxylic_acids = kwargs.get('carboxylic_acids',CompoundList())
        self._build_params = kwargs.get('build_params',{})


    def synthesize_library(self,verbose=False):
        """Create all possible Ugi molecules from the reagents.
       
        Kwargs:
           Verbose (bool): Default=False. Include detailed annotations of the 4 components for each Ugi compound. 
       
        """          
        
        for iso in self.isocyanides:
            for ami in self.amines:
                for ald in self.aldehydes:
                    for car in self.carboxylic_acids:

                        ugi_composition = mass.Composition(formula=''.join(
                                                            (iso['formula'],
                                                             ami['formula'],
                                                             ald['formula'],
                                                             car['formula']))
                                                            ) - mass.Composition(formula='H2O')

                        ugi_mass = mass.calculate_mass(composition=ugi_composition)
                        
                        self._compound_list.append( {'type':'ugi',
                                                     'composition':ugi_composition,
                                                     'mass':ugi_mass,
                                                     'cid':0, # indicate not in PubChem database
                                                     'name':'Ugi Product (Ami:'+str(ami['cid'])+'-'+'Ald:'+str(ald['cid'])+'-'+'Car:'+str(car['cid'])+'-'+'Iso:'+str(iso['cid'])+')',
                                                     'cid_isocyanide':iso['cid'],
                                                     'cid_aldehyde':ald['cid'],
                                                     'cid_carboxylic_acid':car['cid'],
                                                     'cid_amine':ami['cid']
                                                     })

                        
                        if verbose:
                            print(len(self),
                                  iso['formula'],
                                  ami['formula'],
                                  ald['formula'],
                                  car['formula'],
                                  ugi_mass)
    
                        for newparam in self._build_params:
                            self._compound_list[-1].update({newparam:''})
                            for subparam in self._build_params[newparam]:
                                if subparam[0]=='isocyanide':
                                    self._compound_list[-1][newparam] += iso[subparam[1]]
                                if subparam[0]=='aldehyde':
                                    self._compound_list[-1][newparam] += ald[subparam[1]]
                                if subparam[0]=='carboxylic acid':
                                    self._compound_list[-1][newparam] += car[subparam[1]]
                                if subparam[0]=='amine':
                                    self._compound_list[-1][newparam] += ami[subparam[1]]
                                    
                            
                            
                        

    def same_compound(self,x,y):
        """ Returns true of the two compounds are the same.
        A match can be from: 
           - all 4 Ugi component CIDs matching
        """
                        
        if ((x.get('cid_aldehyde',-1)==y.get('cid_aldehyde')) 
            and (x.get('cid_amine',-1)==y.get('cid_amine'))
            and (x.get('cid_carboxylic_acid',-1)==y.get('cid_carboxylic_acid'))
            and (x.get('cid_isocyanide',-1)==y.get('cid_isocyanide'))):
            return True

        return CompoundList.same_compound(self,x,y)




    def generate_ugi_products(my_well_plate,reaction_yield=0.95):
        for w in my_well_plate:
            c = w['compound_list']
            for iso in c.get_subset('type','isocyanide'):
                for am in c.get_subset('type','amine'):
                    for ald in c.get_subset('type','aldehyde'):
                        for carb in c.get_subset('type','carboxylic acid'):
                            my_ugi_lib = UgiLibrary(description=w['mixture_name'])
                            my_ugi_lib.add_reagents(isocyanides=c.get_subset('type','isocyanide'),
                                                    amines=c.get_subset('type','amine'),
                                                    aldehydes=c.get_subset('type','aldehyde'),
                                                    carboxylic_acids=c.get_subset('type','carboxylic acid')
                                                    )
                            my_ugi_lib.synthesize_library()
                            my_ugi_lib.calculate_isotopes()
                            
                            myvol=w['volume']
                            my_well_plate.add_compounds_to_location(w['position'],
                                                                    c,
                                                                    volume = -reaction_yield*myvol)
                            
                            my_well_plate.add_compounds_to_location(w['position'],
                                                                    my_ugi_lib,
                                                                    volume=reaction_yield*myvol)
            c.remove_if_zero_volume()
                            

    def list_ugi_products(my_well_plate,bounds=None):
        my_ugis = CompoundList()
        
        if bounds is None:
            for w in my_well_plate:
                more_ugis = w['compound_list'].get_subset('type','ugi')
                my_ugis.add_compounds(compound_list=more_ugis,
                                      other_properties=[ ('position',w['position']),
                                                         ('location',Containers.position_to_lettergrid(w['position'])) 
                                                         ])
        else:
            myposlist = Containers.bounds_to_positions(bounds)
            
            for row,col in myposlist:
                more_ugis = my_well_plate[row,col]['compound_list'].get_subset('type','ugi')
                my_ugis.add_compounds(compound_list=more_ugis,
                                      other_properties=[ ('position',(row,col)),
                                                         ('location',Containers.position_to_lettergrid((row,col))) 
                                                         ])
    
        return my_ugis



# EXAMPLES OF HOW TO USE THIS CLASS
if __name__ == '__main__':
    
    print('Creating an Ugi library.')
    
    ugi_components = CompoundList(
                            description='Reagents for 102 Ugi Demo, Feb 2018',
                            name_csv_file='examples/demo1_ugi_reagents1.csv',
                            csv_col=2,
                            csv_properties=( ('type',0),('common name',1) ) 
                            )    
    ugi_components.auto_fill_properties_from_pubchem()
    
    my_ugi_lib = UgiLibrary(description='my first Ugi library')
    my_ugi_lib.add_reagents(isocyanides=ugi_components.get_subset('type','isocyanide'),
                            amines=ugi_components.get_subset('type','amine'),
                            aldehydes=ugi_components.get_subset('type','aldehyde'),
                            carboxylic_acids=ugi_components.get_subset('type','carboxylic acid')
                            )
    my_ugi_lib.synthesize_library()
    print(str(my_ugi_lib))
    
    