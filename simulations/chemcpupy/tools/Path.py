#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: Path.py
@author: chrisarcadia
@description: file and directory path related methods
@created: Thu May 31 10:00:03 2018
"""

import os

# file handling

def get_dir(file_name): 
    # get the directory that contains the specified file 
    path = os.path.dirname(sanitize(file_name));
    return path

def make_temp_dir(directory_name='__temporary__'): 
    # make temporary directory
    path = sanitize(directory_name);
    if not os.path.isdir(path):
        os.mkdir(path)
    return path

def sanitize(path):    
    # sanitize path slashes (setting them for the current os and removing duplicates)
    path = os.path.normpath(path);
    return path;

def join(path_parts):    
    # join parts of a path (using os specific slashes)
    path = '';
    if isinstance(path_parts, (list, tuple,)):
        path = path_parts[0];            
        for n in range(1,len(path_parts)):
            path = path + os.sep + path_parts[n];
    elif isinstance(path_parts, (str,)):
        path = path_parts;   
    else: # provided invalid input
        path = '';
    path = sanitize(path);        
    return path;

