#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 13:26:41 2018

@author: jacobrosenstein
"""

import json
from datetime import datetime
import numpy as np
import copy
import chemcpupy as ccpu

class MixtureList:
    
    def __init__(self,*args,**kwargs):
        self._mixture_list = []
        self._description = kwargs.pop('description','')   
        self._file_type='MixtureList-1.0'
        self._file_subtype=''
        
        # add compounds from keyword-specified sources
        self.add_mixtures(**kwargs)
                
    def set_description(self,description):
        self._description=description
        
    def add_mixtures(self,*args,**kwargs):
        _mixture_list = kwargs.get('mixture_list',None) 
        _mixturejsonfile = kwargs.get('mixture_json_file',None)
        
        if _mixture_list is not None:
            self._mixture_list.extend(_mixture_list)
                                
        if _mixturejsonfile is not None:
            print('Loading mixtures from ',_mixturejsonfile)
            with open(_mixturejsonfile,'r') as f:
                filecontents = json.load(f)
                self._mixture_list.extend(filecontents['mixture_list'])
                self._description = self._description + filecontents['description'] + ' / '
            for m in self._mixture_list:
                m['compound_list']=ccpu.CompoundList(compound_list=m['compound_list']['_compound_list'])

        
    def add_mixture(self,mixture_name,compound_list):
        self._mixture_list.append({'mixture_name':mixture_name,
                                   'compound_list':compound_list
                                   })
        
    def add_compounds_to_mixture(self,mixture_name,compound_list,volume=0, position=None):
        for i in range(len(self)):
            if self._mixture_list[i]['mixture_name']==mixture_name:
                if (position is None) or (self._mixture_list[i]['position']==position):
                    self._mixture_list[i]['compound_list'].add_compounds(compound_list=compound_list,
                                                                     total_volume=volume)
        
    def __getitem__(self,key):
        """ You can retrieve a mixture by name.
        """
        if type(key)==str:
            return [x for x in self._mixture_list if x['mixture_name']==key][0]
        else:
            raise Exception('MixtureList only supports indexing by name')

    def __iter__(self):
        self._iter_index = 0
        return self

    def __next__(self):
        if self._iter_index >= len(self):
            raise StopIteration
        else:
            self._iter_index += 1
            return self._mixture_list[self._iter_index-1]
        
    def __len__(self):
        """ The length of a MixtureList is the number of mixtures it contains.
        """
        return len(self._mixture_list)        
        
    def __getitem__(self,key):
        """ You can retrieve a mixture by name.
        """
        return [x for x in self._mixture_list if x['mixture_name']==key][0]            

    def __iter__(self):
        self._iter_index = 0
        return self

    def __next__(self):
        if self._iter_index >= len(self):
            raise StopIteration
        else:
            self._iter_index += 1
            return self._mixture_list[self._iter_index-1]
        
    def __len__(self):
        """ The length of a MixtureList is the number of mixtures it contains.
        """
        return len(self._mixture_list)        
        
    def save_json(self,tofile):        
        
        # save compound lists as dictionaries to be JSON compatible
        m_temp = copy.deepcopy(self._mixture_list)
        for m in m_temp:
            m.update( {'mixture_size':len(m['compound_list'])} )
            m.update( {'compound_list':m['compound_list'].__dict__})
        
        print('Saving compounds to JSON file: ',tofile)
        myinfo = {'mixture_list':m_temp,
                  'created':str(datetime.now()),
                  'file type':self._file_type,
                  'file subtype':self._file_subtype,
                  'description':str(self._description),                  
                  }        
        with open(tofile, 'w') as f:
            json.dump(myinfo, f, indent=2)
            
        
    def get_mixture_names(self):
        return [m['mixture_name'] for m in self._mixture_list]
    
    def get_compound_list(self,mixture_name):
        for m in self._mixture_list:
            if m.get('mixture_name')==mixture_name:
                return m['compound_list']
        return []
            
    def get_mixture_masses(self,mixture_name):
        return np.array([c['mass'] for c in self.get_compound_list(mixture_name)])
        
    def calculate_isotopes(self):
        for m in self._mixture_list:
            cl = CompoundList(compound_list=m['compound_list'])
            cl.calculate_isotopes()
            m['compound_list'] = cl._compound_list
    
    def get_mixture_isotope_masses(self,mixture_name):
        return np.array([c['isotope_masses'] for c in self.get_compound_list(mixture_name)])
    
if __name__=="__main__":
    
    mixturefile = 'demos_and_data/ugi_demo102_Feb2018/demo102_mixture_list.json'
    mymixtures = MixtureList(mixture_json_file=mixturefile)
    
    
    
    
    