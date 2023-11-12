#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 10:00:03 2018

@author: jacobrosenstein
"""

import numpy as np
import chemcpupy as ccpu

# This file is for miscellaenous methods.
# Note:
#   To avoid creating a messy code base please refrain from adding final methods to this file.
#   It is much better to have methods kept in a class/file with other related ones, than here.

# ======================================================================        
# HELPER FUNCTIONS

def binary_image_to_positions(image,offset=(0,0)):
    """ Converts a binary image (numpy bool array) into a list of positions
    where the value is 1.
    """
    return list(map(tuple,(np.argwhere(image)+np.array(offset)).tolist()))

# alias
image2pos=binary_image_to_positions


# =========================================================
# COMMON COMPOUNDS

WATER = ccpu.CompoundList(cid_list=[962,])
DMSO = ccpu.CompoundList(cid_list=[679,])
HCCA = ccpu.CompoundList(cid_list=[5328791,])




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LEGACY METHODS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



# ======================================================================        
# DATA FUNCTIONS 

# IMPORTANT NOTE:
# the below functions are depreciated and kept here for legacy reasons
# please instead use Plating.py to access these and related methods  

def make_permutation(mydata,seed=12345):
    np.random.seed(seed)
    perm = np.random.permutation(mydata.size)
    antiperm = perm.argsort()
    permdata = mydata[perm]
    return perm,antiperm,permdata


# ======================================================================        
# POSITION FUNCTIONS 

# IMPORTANT NOTE:
# the below functions are depreciated and kept here for legacy reasons
# please instead use Containers.py to access these and related methods  

#rowletters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 
#              'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 
#              'Y', 'Z', 'AA', 'AB', 'AC', 'AD', 'AE', 'AF']
rowletters = [chr(x) for x in range(ord('A'),ord('Z')+1)] + ['A' + chr(x) for x in range(ord('A'),ord('Z')+1)];    

# -----------------------------------
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

# alias
p2l=position_to_lettergrid

# -----------------------------------
def lettergrid_to_position(key):
    """ Returns a zero indexed numerical position tuple.
    Rows are letters, columns are numbers. 'A1' returns (0,0). 'B4 returns (1,3).
    """        
    if type(key)==str:
        getletters = ''.join([x for x in key if not x.isdigit()])
        getnumbers = ''.join([x for x in key if x.isdigit()])
        return ( rowletters.index(getletters), int(getnumbers)-1 )
    
    if type(key)==list:
        return [lettergrid_to_position(x) for x in key]

    # not recognized. assume it does not need to be converted, return as-is
    return key

# alias
l2p=lettergrid_to_position

# -----------------------------------
def position_to_maldigrid(position):
    """ Returns a string corresponding to the MALDI position. Zero indexed.
    Rows are letters, columns are numbers. (0,0) returns 'A1'. (1,3) returns 'B4'.
    """
    if type(position)==tuple:
        return 'X%02dY%02d' % (position[1]+1,position[0]+1)

    if type(position)==list:
        return [position_to_maldigrid(x) for x in position]

    # not recognized. assume it does not need to be converted, return as-is
    return position

# alias
p2m=position_to_maldigrid

# -----------------------------------
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

# --------------------------------------

def rotate_positions_180(positions,plate):
    rows,cols = plate.shape
    return [(rows-1-x,cols-1-y) for x,y in l2p(positions)]
    
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
