# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 15:06:37 2018

@author: labuser
"""


#import chemcpupy as ccpu
#from pprint import pprint

import os
#import numpy as np
#import pylab as pl

#from pyteomics import mass


if not os.path.isdir('__temporary__'):
    os.mkdir('__temporary__')
    
    
    
    

#mydir = r"M:\molinfo\mass_spectra\Bruker_20180227_102Ugi_individual"
#mydir = r"M:\molinfo\mass_spectra\Bruker_20180218_102Ugi"
#mydir = r"M:\molinfo\mass_spectra\Bruker_20180216_102Ugi"
mydir = r"M:\molinfo\mass_spectra\Bruker_20171115_25Ugi"

#all_analyses = sorted(ccpu.BrukerMS.list_bruker_analyses(mydir))

bruker_analyses = [ name for name in os.listdir(mydir) if (os.path.splitext(name)[1]=='.d') ]
        
for this_analysis in bruker_analyses:
    
    myanalysis = this_analysis
    
    # remove spaces from the filename
    if ' ' in this_analysis:
        myanalysis = this_analysis.replace(' ','_')
        print('Changing analysis name to underscores instead of spaces: ',myanalysis)
        
        os.rename(os.path.join(mydir,this_analysis),
                  os.path.join(mydir,myanalysis))
    
    
    # requires ProteoWizard MSConvert tool
    # http://proteowizard.sourceforge.net/
    proteowizard_path = '\"' + r"C:\Program Files\ProteoWizard\ProteoWizard 3.0.11806\msconvert.exe" + '\"'
    
    run_str = proteowizard_path+' '+os.path.join(mydir,myanalysis)+' --outdir '+mydir+' --mzML >> __temporary__/proteowizard_log.txt'
    print(run_str)
    
    if os.system(run_str)==0:
        print('Success')
    else:
        print('MSCONVERT ERROR')
    
    
    
