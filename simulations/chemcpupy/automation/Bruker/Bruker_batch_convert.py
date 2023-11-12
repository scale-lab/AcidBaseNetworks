# -*- coding: utf-8 -*-
"""
@title: Bruker_batch_convert.py
@description: script that converts each Bruker '.d' file in a directory into an '.hdf5' file containing data series and measurement metadata
@author: chrisarcadia 
@created: 2018/10/29
"""
import Bruker
import timeit
import os, sys

# Provide a source folder or data file 
#path_source = os.path.join(os.getcwd(),'resource','single_position_MALDI.d');
# Examples:
#   path_source = r'C:\Users\ChrisTow\Desktop\Examples\single\Mix_1_100_1.d'; # a single (single spot) data file 
#   path_source = r'C:\Users\ChrisTow\Desktop\Examples\multiple'; # a folder with multiple (single spot) data files
#path_source = r'C:\Users\ChrisTow\Desktop\Examples\maldi'; # a folder with a multi-spot MALDI data file
path_source = r'Z:\group\mass_spectra\0165 Library 13\Averaging and Resolutions Tests'

# Set a destination folder
path_destin = r'C:\Users\ChrisTow\Downloads\Converted'; # place to write converted data to
#path_destin = os.path.join(os.getcwd(),'__temporary__'); # place to write converted data to

# Prepare destination folder
path_destin = Bruker.prep_output_dir(path_source,path_destin);

# Get list of available '.d' files
file_list = Bruker.list_bruker_files(path_source);
#path_source = ;  # location of raw data
#file_list = [r'Z:\group\mass_spectra\0165 Library 13\e0165p01t03.d']; # manually specify multiple files

# Conversion options
# time and mz
export_compression = None; # compression to use for export {None, 'gzip', ...}
export_mz_domain = False; # export mass/charge domain data
export_mz_domain_resampled = True; # export mass/charge domain data after resampling all spectra to the same mz grid
export_mz_domain_resampled_sum = True; # during export, accumulate the aggregate spectra by summing all intensity vectors
export_time_domain = False; # export time domain data
export_all_settings_to_MAT = True; # save an exhaustive copy of file settings to a seperate MAT (MATLAB file)
export_every_mz_vector = False; # export all mass/charge vectors, as opposed to a single one
export_statistics = True; # during export, compute statistical measures of spectra intensity 
export_reduced = True; # export lossy reduced/compressed version of aMALDIms for portability
# time only
#export_compression = None; # compression to use for export {None, 'gzip', ...}
#export_mz_domain = False; # export mass/charge domain data
#export_mz_domain_resampled = False; # export mass/charge domain data after resampling all spectra to the same mz grid
#export_mz_domain_resampled_sum = False; # during export, accumulate the aggregate spectra by summing all intensity vectors
#export_time_domain = True; # export time domain data
#export_all_settings_to_MAT = True; # save an exhaustive copy of file settings to a seperate MAT (MATLAB file)
#export_every_mz_vector = False; # export all mass/charge vectors, as opposed to a single one
#export_statistics = False; # during export, compute statistical measures of spectra intensity 

# Convert each file        
conversion_times = [];
print('starting conversion of "' + path_source + '"')
for file in file_list:
    
    # Start timer
    toc = timeit.default_timer(); 
    
    # Convert files to HDF5 (a self contained single universally readable data file, with essential settings)
    name = Bruker.get_file_parts(file)['name'];      
    fileout = os.path.join(path_destin,name);
    Bruker.convert_data(file, fileout+'.hdf5', 
                       include_spectra = export_mz_domain, 
                       include_respectra = export_mz_domain_resampled, 
                       include_raw = export_time_domain, 
                       copy_settings_to_mat = export_all_settings_to_MAT,
                       compression = export_compression,
                       save_each_mz_vector = export_every_mz_vector,
                       get_statistics = export_statistics,
                       collect_sum = export_mz_domain_resampled_sum); # this will write the data to '.hdf5' file                            
        
    # End timer
    tic = timeit.default_timer(); 
    conversion_times = conversion_times + [(tic-toc)/60]; # [min]
    print('converted "' + name + '" (in ' + '{0:3.1f}'.format(conversion_times[-1])  + ' minutes)')     
    
    # Convert HDF5 file to mat (reduced/compressed MATLAB version of aMALDIms) 
    if export_reduced:
        try:
            import matlab.engine 
            # for the above package, you need a licensed version of MATLAB installed and must install the MATLAB API for Python: https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
            # to do so run the following in the terminal (replace the matlab root folder with your install directory)
            # > cd "/Applications/MATLAB_R2019a.app/extern/engines/python" # or for windows: # > cd "C:\Program Files\MATLAB\R2019b\extern\engines\python"
            # > python setup.py install          
        except:
            print('The library "matlab.engine" could not be imported. Please check that it is installed.')     
        if   export_mz_domain_resampled and ('matlab.engine' in sys.modules):
            MATLAB = matlab.engine.start_matlab()
            amaldims = MATLAB.aMALDIms('');
            [compressfiles] = MATLAB.compress_multiple_respectra_by_threshold(amaldims, [fileout+'.hdf5'],
                'Algorithm','PeakThreshold', # compression method
                'Threshold',30,              # signal cutoff value
                'FindPeaks',True,            # also save all spectra peaks
                'NormalizeIntensity',True,   # intensity (false) or snr (true)
                'IncludeRawData',False);     # save raw time-series data)
            print('reduced all files')     
        
print('finished conversion, see "' + path_destin + '"')
