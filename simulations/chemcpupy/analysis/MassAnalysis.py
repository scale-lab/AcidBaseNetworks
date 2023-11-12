# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:45:16 2018

@author: jacobrosenstein

A collection of functions to analyze the masses of lists of compounds.
"""

import os
import numpy as np
import chemcpupy as ccpu
from pyteomics import mass, mzml, auxiliary
import pylab as pl


# ----------------------------------------------
def MassStats(cl):
    """Generate summary statistics describing the masses in a CompoundList
   
    Args:
       cl (CompoundList): The CompoundList to analyze.

    Raises:
        Exception
        
    """
    all_masses = cl.list_all('mass')
    
    unique_masses = set(all_masses)
    sorted_unique_masses = sorted(unique_masses)
    mass_spacing = [sorted_unique_masses[x+1]-sorted_unique_masses[x] for x in range(len(sorted_unique_masses)-1)]
    
    if all_masses.count(None)>0:
        raise Exception('Tried to generate MassStats, but mass not defined.')
    
    mystats = { 'count_all':len(all_masses),
                'mean':np.mean(all_masses),
                'min':np.min(all_masses),
                'max':np.max(all_masses),
                'count_unique':len(unique_masses),
                'min_spacing':min(mass_spacing),
                'max_spacing':max(mass_spacing),
                }

    return mystats



def rolling_window(a, window):
    # function by Erik Rigtorp, www.rigtorp.se
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def rolling_maximum(intensity,window_halfwidth=13):
    """ Helper function, calculates a rolling maximum intensity.
    """
    max_local_intensity=np.max(rolling_window(intensity,2*window_halfwidth+1),-1)
    return max_local_intensity



def get_nearest_maxima(mz,intensity,mass_list,ion_list=('H',),window_halfwidth=13):
    """ Finds the maximum intensities near a given list of masses.
    """
    max_local_intensity = rolling_maximum(intensity,window_halfwidth)
    window_mz = mz[(window_halfwidth):-(window_halfwidth)]
        
    maxima={}
    maxima['total'] = np.zeros_like(mass_list)
    for ion in ion_list:
        maxima[ion] = (np.interp(mass_list + mass.calculate_mass(formula=ion),
                                 window_mz,
                                 max_local_intensity))
    
        maxima['total'] += maxima[ion]
    
    return maxima


def plot_expected_mz(mass_list,ion_list):
    """ Plotting helper function.
    """
    for ion in ion_list:
        pl.plot(mass_list + mass.calculate_mass(formula=ion), -5e5*np.ones((len(mass_list),1)), '.', label=('expected %s'%(ion,)))
    
    
def plot_meas_intensities(mass_list,ion_list,measured):
    """ Plotting helper function.
    """
    for ion in ion_list:
        pl.plot(mass_list + mass.calculate_mass(formula=ion), 
                measured[ion], 
                'o', 
                markerfacecolor='None',
                markersize=12,
                label=('measured %s'%(ion,)))

def plot_meas_intensity_histograms(measured,ion_list,labelprefix,color='r'):
    """ Plotting helper function
    """
    pl.xscale('log')
    mybins = np.logspace(3,10,num=12)
    bincenters = 0.5*(mybins[1:]+mybins[:-1])
    
    for ion in ion_list:
        n,bins=np.histogram(measured[ion],mybins)    
        pl.plot(bincenters,n/len(measured[ion]),color,label=(labelprefix+' M'+ion))
        #pl.plot(bincenters,n,color,label=(labelprefix+' M'+ion))
    

def calc_absent_masses(present_masses,all_masses,digits=5):
    """ Given a list of present masses, and a list of all masses, 
    calculates the absent masses.
    """
    absent_masses = np.array(list(set(np.round(all_masses,digits)).difference(np.round(present_masses,digits))))
    return absent_masses


def list_mzml(mydir):
    """ Lists all .mzML analyses in the specified directory
    """
    mymzml = [ name for name in os.listdir(mydir) if (os.path.splitext(name)[1]=='.mzML') ]
    return mymzml

def read_mzml_analysis(mydir,myanalysis):
    """ Temporarily only supports 1 scan per .mzML file
    """
    with mzml.read(os.path.join(mydir,myanalysis),huge_tree=True) as reader:
        return reader["scan=1"]
