# -*- coding: utf-8 -*-
"""
Created on Thu May 10 12:17:52 2018

@author: Eamonn
"""


import chemcpupy as ccpu
from pprint import pprint

import os
from sys import platform
import numpy as np
import pylab as pl

from pyteomics import mass


if not os.path.isdir('__temporary__'):
    os.mkdir('__temporary__')


#mydir = r"M:\molinfo\mass_spectra\Bruker_20180511\Serial Test"  #/Serial.d
mydir = r"M:\molinfo\mass_spectra\Bruker_20180511\Plate1"

# Get all serial folder names
all_analyses = ccpu.BrukerMS.list_bruker_analyses(mydir)

for this_analysis in all_analyses: # For each serial file

    # Find the index associated with spot number for the serial file
    spot_index  = ccpu.BrukerMS.find_spot_number_index(mydir,this_analysis)

    # Summarize the data and print the summary details
    summary = ccpu.BrukerMS.summarize_serial_file(mydir,this_analysis,spot_index)
    print('"' + this_analysis + '" contains ' +str(max(summary.keys()))+ ' unique spectra starting at ' + summary[1])
    
    # Generate a dictionary of the mz and intensity arrays for each well
    dataset = ccpu.BrukerMS.read_serial_file(mydir,this_analysis,summary)
    
    # Write the serial file contents as mat files in a serial-unique folder
    #serial_dir = mydir + '\\' + this_analysis.partition('.')[0]
    #ccpu.BrukerMS.mat_write(dataset,summary,serial_dir)
    #ccpu.BrukerMS.txt_write(dataset,summary,'.')
    