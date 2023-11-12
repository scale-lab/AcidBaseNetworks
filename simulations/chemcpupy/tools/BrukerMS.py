# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 11:38:51 2018

@author: jrosenstein


Functions to read Bruker .d data analysis directories. This uses Windows 
Component Object Model (COM) interfaces and will not run on other operating systems.

This requires installing Bruker CompassXtract. 
Tested with CompassXtract 3.2.201 (64 bit) on Windows 7
and python 3.6 in Anaconda
"""

from sys import platform

if platform != "win32":
    raise Exception('Bruker CompassXtract is Windows-only')



import os
import comtypes.client
import numpy as np
import pylab as pl
import scipy.io as sio

from scipy.ndimage.filters import maximum_filter1d,median_filter

# ==========================================================================

SPECTRUMTYPE_LINE = 0
SPECTRUMTYPE_PROFILE = 1

SPECTRUMTYPE_UNCALIBRATED = False
SPECTRUMTYPE_RECALIBRATED = True



def list_bruker_analyses(mydir):
    # returns all .d analyses in the specified directory
    myanalyses = [ name for name in os.listdir(mydir) if (os.path.splitext(name)[1]=='.d') ]
    return myanalyses


###########################################################
# Function to ensure the existence of a directory

def ensure_dir(saveDir):
    os.path.exists(saveDir)
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
        
        
############################################################        

def read_bruker_analysis(mydir,myanalysis,get_all_parameters=True):
    
    MSAnalysis = comtypes.client.CreateObject("EDAL.MSAnalysis")
    #MSAnalysis = comtypes.client.CreateObject("EDAL.IMSAnalysis2")
    MSAnalysis.Open(os.path.join(mydir,myanalysis))

    count = MSAnalysis.GetTypeInfoCount()
    
    
    other_parameters={}
    for i in range(count):
        #print('i=',i)
        try:
            sp = MSAnalysis.MSSpectrumCollection(i+1)
            mz,intensity,points = sp.GetMassIntensityValues(SPECTRUMTYPE_PROFILE)
            #mz,intensity,points = sp.GetRecalibratedMassIntensityValues(SPECTRUMTYPE_PROFILE,
            #                                                            SPECTRUMTYPE_UNCALIBRATED)
            print('Length: ',points)
        
            if get_all_parameters:
                j=1
                while j is not -1:
                    try:
                        p=sp.MSSpectrumParameterCollection(j)
                        other_parameters[p.ParameterName]=p.ParameterValue
                    except:
                        j=-1
                    else:
                        j+=1
        except:
            #print('exception for i=',i)
            print('ERROR: problem reading Bruker log file with CompassXtract')
            
    return { 'm/z array':np.array(mz,dtype='float32'),
             'intensity array':np.array(intensity,dtype='float32'),
             'points':points,
             'Analysis':myanalysis,
             'AnalysisName':MSAnalysis.AnalysisName,
             'AnalysisDateTimeIsoString':MSAnalysis.AnalysisDateTimeIsoString,
             'SampleName':MSAnalysis.SampleName,
             'MethodName':MSAnalysis.MethodName,
             'AnalysisDescription':MSAnalysis.AnalysisDescription,
             'other_parameters':other_parameters
             }





###########################################################
# A function to read and output a dictionary of all spectra in a serial file

def read_serial_file(mydir,myanalysis,summary):
    
    MSAnalysis = comtypes.client.CreateObject("EDAL.MSAnalysis")
    MSAnalysis.Open(os.path.join(mydir,myanalysis))
    
    dataset={}
    for i in range(len(summary)):
        wellName = summary[i + 1]
        print(i + 1, wellName)
        sp = MSAnalysis.MSSpectrumCollection(i + 1)
        mz,intensity,points = sp.GetMassIntensityValues(1)
        dataset[wellName+ ' m/z array'] = np.array(mz,dtype='float32')
        dataset[wellName+ ' intensity array'] = np.array(intensity,dtype='float32')
            
    return dataset

###########################################################
# A function to read and output a dictionary of all spectra in a serial file

def read_serial_file_to_array(mydir,
                              myanalysis,
                              summary,
                              mz_interp=None,
                              rollingmedian=None, # optional:subtract a rolling median
                              normalizeSNR=False,
                              maxspots=None,
                              maxwindow=None,
                              maldi_positions=None):

    if maldi_positions is not None:
        # follow specified positions instead of the raw summary
        summary = dict([(i+1,maldi_positions[i]) for i in range(len(maldi_positions))])
        print('following specified MALDI positions ',len(maldi_positions))
    
    if maxspots is not None:
        summary = dict([(i,summary[i]) for i in list(range(1,maxspots+1))])
        print('limiting maximum spots to ',maxspots)


    MSAnalysis = comtypes.client.CreateObject("EDAL.MSAnalysis")
    MSAnalysis.Open(os.path.join(mydir,myanalysis))

    dataset = {}
    for i in range(len(summary)):
        wellName = summary[i + 1]
        print(i + 1, wellName)
        sp = MSAnalysis.MSSpectrumCollection(i + 1)
        mz,intensity,points = sp.GetMassIntensityValues(1)

        if mz_interp is not None:
            intensity = np.interp(mz_interp,mz,intensity)
            mz = mz_interp
            
        if i==0:
            numspots = len(summary)
            mzlen = len(mz)
            dataset['m/z'] = np.array(mz,dtype=np.float32)
            dataset['intensity array'] = np.zeros([numspots,mzlen],dtype=np.float32)
        
        dataset['intensity array'][i,:] = np.array(intensity,dtype='float32')
        
        if rollingmedian is not None:
            dataset['intensity array'][i,:] = dataset['intensity array'][i,:] - median_filter(dataset['intensity array'][i,:], rollingmedian)
            
        if normalizeSNR:
            filterthresh = np.percentile(dataset['intensity array'][i,:],95)    # find the 95th percentile amplitude
            stdev = np.std(dataset['intensity array'][i,dataset['intensity array'][i,:]<filterthresh])   # exclude outliers from STDEV
            dataset['intensity array'][i,:] = dataset['intensity array'][i,:] / stdev
            
        if maxwindow is not None:
            # APPLY SLIDING MAXIMUM FILTER
            dataset['intensity array'][i,:] = maximum_filter1d(dataset['intensity array'][i,:], maxwindow)
                
    return dataset


###########################################################
# Function to summarize and error check well numbers / spectra in a serial file

def summarize_serial_file(mydir,myanalysis,spot_index):
    
    MSAnalysis = comtypes.client.CreateObject("EDAL.MSAnalysis")
    MSAnalysis.Open(os.path.join(mydir,myanalysis))
    count = MSAnalysis.GetTypeInfoCount()
    
    analysis = {}
    analysis_summary = {}
    
    i = 1
    while True:
        try:
            sp = MSAnalysis.MSSpectrumCollection(i)     
            p=sp.MSSpectrumParameterCollection(spot_index)
            stringID = p.ParameterValue.partition(':')[0]
            analysis[i] = stringID # generate a dictionary of well IDs                                 
            i+=1
        except:
            
            # Check for double recordings of the same well in the serial
            for key,value in analysis.items():
                if value not in analysis_summary.values():
                    analysis_summary[key] = value
            numReps = len(analysis) - len(analysis_summary);
            if numReps>=1:
                print('WARNING: ' +str(myanalysis) + ' contains ' +str(numReps)+ ' duplicate spectra:')         
                print(str(set(analysis.items()).symmetric_difference(set(analysis_summary.items()))))
            return analysis_summary
            break
      
###########################################################
# Function to find the parameter index of 'Spot Number'

def find_spot_number_index(mydir,myanalysis):  
    
    MSAnalysis = comtypes.client.CreateObject("EDAL.MSAnalysis")
    MSAnalysis.Open(os.path.join(mydir,myanalysis))
    count = MSAnalysis.GetTypeInfoCount()
    sp = MSAnalysis.MSSpectrumCollection(1)
    
    other_parameters={}
    j=1
    while j is not -1:
        try:
            p=sp.MSSpectrumParameterCollection(j)
            if p.ParameterName == 'Spot Number':
                spot_index = j
        except:
                j=-1
        else:
                j+=1
    return spot_index
        
############################################################## 
# Function to write a serial dataset to single .mat files

def mat_write(dataset,summary,serial_dir):
    saveDir = serial_dir + '_matOnly'
    ccpu.BrukerMS.ensure_dir(saveDir)
    
    for wellName in summary.values(): 
        os.chdir(saveDir)
        saveName = 'rawData_'  + wellName
        rawData = {'mz': dataset[wellName+ ' m/z array'], 'intensity': dataset[wellName+ ' intensity array']}
        sio.savemat(saveName+'.mat', {saveName : rawData})
        
        
############################################################## 
# Function to write a serial dataset to series of .txt files

def txt_write(dataset,summary,serial_dir):
    #saveDir = serial_dir + '_txtOnly'
    #ccpu.BrukerMS.ensure_dir(saveDir)
    saveDir = serial_dir
    
    for wellName in summary.values(): 
        #os.chdir(saveDir)
        saveName = 'rawData_'  + wellName
        rawData = {'mz': dataset[wellName+ ' m/z array'], 'intensity': dataset[wellName+ ' intensity array']}
        np.savetxt(saveName+'_mz.txt', rawData['mz'])
        np.savetxt(saveName+'_intensity.txt', rawData['intensity'])
        
        
  
        
        
        
        
        

# EXAMPLE

if __name__=="__main__":
    
    import numpy as np
    import pylab as pl

    #mydir = r"M:/molinfo/mass_spectra/Bruker_20180218_102Ugi"
    mydir = r"M:/molinfo/mass_spectra/Bruker_20180504/SerialAcquisitionTest/serial_file"  #/Serial.d
    all_analyses = list_bruker_analyses(mydir)
    
    
    for this_analysis in all_analyses:
        summary = summarize_serial_bruker_analysis(mydir,this_analysis)
        print([x['Spot Number'] for x in summary])
        
#        
#    for this_analysis in all_analyses:
#        results = read_bruker_analysis(mydir,this_analysis)
#        
#        pl.plot(results['m/z array'],results['intensity array'])
#        pl.xlim(200,1000)
#        pl.title(results['AnalysisName'])
#        pl.show()

