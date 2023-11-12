# -*- coding: utf-8 -*-
"""
@title: test_convert_raw.py
@description: Example Python script to test the raw data converter and to test reloading the converted raw data
@author: chrisarcadia 
@created: 2018/10/26
"""

import Bruker
import matplotlib.pyplot as pyplot
import h5py
import numpy

# Convert a data file
try:
    input_filename = r'C:\Users\ChrisTow\Desktop\Examples\single\Mix_1_100_1.d';
    output_filename = r'C:\Users\ChrisTow\Desktop\Converted\single\Mix_1_100_1.hdf5';        
    settings = Bruker.get_measurement_settings(input_filename);
    positions = settings['positions'];    
    Bruker.convert_settings(output_filename, settings);
    Bruker.convert_data_raw(input_filename,output_filename,settings);
    converted = 1;     
except:
    converted = 0; # do nothing if file already exists
    
# Load the converted data
hf = h5py.File(output_filename, 'r');
raw = hf.get('raw/');
index = 0;
time = numpy.array(raw['time']);
signal = numpy.array(raw['signal'][index,:]);
info = {}
for k in raw.attrs.keys():
    info.update({k:raw.attrs[k]});
hf.close();
print('Data Info:')
print(info);

# Plot one of the time series
pyplot.plot(time,signal, 'bo', markersize=1)
pyplot.ylabel('Signal')
pyplot.xlabel('Time [s]')
pyplot.grid(True)
pyplot.show()
