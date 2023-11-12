#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chrisarcadia 
@created: 2018/03/21
"""

#-------------- Readme

# run from chemicalcpu root

# for simplicity this simple example does not require knowledge of nor does it 
# currently utilize the proposed Tool classes (nor any others from 'chemcpupy')

#-------------- Setup

# packages
import os
import numpy
import pubchempy 
import matplotlib
#import chemcpupy 
#import mpldatacursor

# temporary folder
pathtemp = '__temporary__';
if not os.path.isdir(pathtemp):
    os.mkdir(pathtemp)
    
# make inline plots vector graphics instead of raster graphics
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('pdf', 'svg')

# specify plot container for IPython
# %matplotlib auto # {auto, inline}

    
#################################### User Input
    

# compound source samples (identified by PubChem CID)
    
# unique/distinguishable molecules for data storage
moleculeCID           = [767,284,176]; 
moleculeSolventCID    = [962,962,962];
moleculeVolumeInitial = [1e-3,1e-3,1e-3];
moleculeContainerInitial = ['1.5mL test tube','1.5mL test tube','1.5mL test tube'];
moleculeConcentration = [1,1,1]; 
moleculeCount = len(moleculeCID);

# universal solvent for which to dilute molecules
solventCID            = [962];
solventVolumeInitial  = [45e-3];
moleculeContainerInitial = ['50mL falcon tube'];   

# tasks
tranfer0Volume = [25e-6]; # transfer 0: from solvent to wells
tranfer1Volume = [10e-6]; # transfer 1: from source to image
tranfer2VolumeMax = [10e-6]; # transfer 2: from image to pool (multiply this by weights to get transfer2Volume)

# binary image data (currently generates dummy weights & images)
    
# to write
imageDimension = [3,3]; #[16,16];
image = [];
for n in range(moleculeCount) :
    image.append((numpy.round(numpy.random.rand(imageDimension[0],imageDimension[1]))).astype(int)); 
imageDimension = image[0].shape;

    
# query
weights = 2*(numpy.random.rand(imageDimension[0],imageDimension[1]))-1;
weights = numpy.round(weights*100)/100;

# preview images & weights
previewImages = 1;
previewWeights = 1;

#################################### Compound Information

# download PubChem compound information
moleculeName = [];
moleculeCompound = [];
moleculeMolarMass = [];
moleculeSMILES = [];
moleculeFormula = [];
for n in range(moleculeCount) :
    moleculeCompound.append(pubchempy.Compound.from_cid(moleculeCID[n]));
    moleculeName.append(moleculeCompound[n].synonyms[0]);
    moleculeMolarMass.append(moleculeCompound[n].molecular_weight);
    moleculeSMILES.append(moleculeCompound[n].isomeric_smiles);
    moleculeFormula.append(moleculeCompound[n].molecular_formula);
    
solventCompound = pubchempy.Compound.from_cid(solventCID);
solventName = solventCompound.synonyms[0];
solventMolarMass = (solventCompound.molecular_weight);
solventSMILES = (solventCompound.isomeric_smiles);
solventFormula = (solventCompound.molecular_formula);


#################################### Visualization
        
# print 2D numpy array as image
#def formattingForHover(**kwargs):
#    i = kwargs['i'];
#    j = kwargs['j'];
#    z = kwargs['z'];
#    return '(i, j) = ({:}, {:})\nz = {:.02g}'.format(i,j,z)
def plotImage( title, celldata ):
    
    data = celldata;
    title = title;
    fig, ax = matplotlib.pyplot.subplots()
    ax.set_title(title,loc='center')    
    ax.set_aspect('equal', adjustable='box')
    cmap ='binary_r'; # coolwarm cool binary binary_r ('_r' gives reversed colormap)
    
    
    ax.xaxis.set_ticks(numpy.arange(0, imageDimension[0], 1))
    ax.yaxis.set_ticks(numpy.arange(0, imageDimension[1], 1))
    ax.imshow(data,cmap=matplotlib.pyplot.get_cmap(cmap)); #,extent=[0,imageDimension[0], 0,imageDimension[1]])
#    mpldatacursor.datacursor(hover=True, bbox=dict(alpha=1, fc='w'), arrowprops=dict(arrowstyle='simple', fc='white', alpha=0.5), formatter=formattingForHover); # 'i, j = {i}, {j}\nz = {z:.02g}'.format)
    
#    rows = range(imageDimension[0]);
#    columns = range(imageDimension[1]);
#    headercolor = "#ffffff";
#    ax.axis('tight')
#    ax.axis('off')    
#    tableplot = ax.table(cellText=data, cellColours=None, cellLoc='center',
#                              rowLabels=rows, rowColours=[headercolor]*imageDimension[0], rowLoc='center',
#                              colLabels=columns, colColours=[headercolor]*imageDimension[1], colLoc='center', colWidths=None, 
#                              loc='center', bbox=None)
#    #tableplot.auto_set_font_size(False);
#    #tableplot.set_fontsize(18);
#    #tableplot.scale(3,3);    
#    datamin = numpy.min(data);
#    datamax = numpy.max(data);
#    colorMap = matplotlib.pyplot.get_cmap(cmap) 
#    colorNorm  = matplotlib.colors.Normalize(vmin=datamin, vmax=datamax)
#    scalarMap = matplotlib.cm.ScalarMappable(norm=colorNorm, cmap=colorMap)
#    for x in range(0, imageDimension[0]):
#        for y in range(0, imageDimension[1]):
#            color = scalarMap.to_rgba(data[x,y]); 
#            tableplot._cells[(x+1, y)].set_facecolor(color);   
       
    
    matplotlib.pyplot.show()
    
    return;
    
# preview images
if previewImages:
    for n in range(moleculeCount) :
        plotImage('Image ' + str(n) + ' : ' + moleculeName[n], image[n])
    
# preview weights
if previewWeights:
    plotImage('Weights', weights)

#################################### 

# write image from source material

# weight (by volume) and pool image pixels (positive & negative weights seperate)



# incomplete 

