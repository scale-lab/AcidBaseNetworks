# -*- coding: utf-8 -*-
"""
@title: test_convert_spectra.py
@description: Example Python script to test the spectra data converter and to test reloading the converted spectra data
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
    Bruker.convert_data_spectra(input_filename,output_filename,settings);
    converted = 1;     
except:
    converted = 0; # do nothing if file already exists
    
# Load the converted data
hf = h5py.File(output_filename, 'r');
spec = hf.get('spectra/');
index = 0;
mz = numpy.array(spec['mass_to_charge'][index,:]);
signal = numpy.array(spec['signal'][index,:]);
info = {}
for k in spec.attrs.keys():
    info.update({k:spec.attrs[k]});
hf.close();
print('Data Info:')
print(info);

# Plot one of the spectra
pyplot.plot(mz,signal, 'bo', markersize=1)
pyplot.ylabel('Signal')
pyplot.xlabel('M/Z')
pyplot.grid(True)
pyplot.show()




    
       
        
    