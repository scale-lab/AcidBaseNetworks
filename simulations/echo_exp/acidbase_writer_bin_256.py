#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chris3arcadia 
@description: template script to generate transfer lists for acid/base computing
@created: Thr Feb 28 2019
"""

# Libraries
import numpy as np
import chemcpupy
from chemcpupy import Plating, Containers
from acidbase_network import AcidBaseNetwork
import matplotlib.pyplot as plt
import os

# Set paths
path_script = chemcpupy.Path.get_dir(__file__) # get current script's directory
path_temp = chemcpupy.Path.make_temp_dir() # make temporary directory

# Script options
# show_plates = True # plot images of plate contents
show_plates = True # plot images of plate contents
simulate_seperately = True # simulate tasklists seperately

# Demo information
demo_name = 'acidbase_bin_16x16_demo'
demo_file = 'demo_export.mat'
demo_source_platesheet = 'demo_source_platesheet.csv'
data_dir = os.path.join(path_script, 'data')
kernel_file = os.path.join(data_dir, 'kernel_mem_bin_256.txt')
img_name = 'img_zero_71'
# img_name = 'img_one_91'
sample_image1 = os.path.join(data_dir, 'bin_images_256/' + img_name + '.txt')
# sample_image1 = os.path.join(data_dir, 'bin_images_256/img_zero_71.txt')

# zero_1.txt -> [[4, -4], [-20, 20]] -> (AB, BA) -> (1, -1)
# one_1.txt -> [[-34, 34], [18, -18]] -> (BA, AB) -> (-1, 1)

def read_file(filepath):
    with open(filepath, "r") as f:
        return f.read()


# Global Variables
min_vol_tranfer = 2.5 # minimum transferrable volume [nL]
max_vol_transfer = 10000 # max transferrable volume[nL]
vol1536_data = 4500 # volume to transfer from source to data plate [nL]
vol_pool = 150 # volume to transfer to pool from a data well
source_max_transfers = 5 # maximum transfers from each source well

indicator_transfer_volume = 100 # Add 100nl to each well in the data plate


grid384_height = 16
grid384_width = 24

grid1536_height = 32
grid1536_width = 48

img_width = 16
img_height = 16

img_np = np.reshape([int(float(val)) for val in read_file(sample_image1).splitlines()], (img_width, img_height))
# p = plt.imshow(img_np, cmap='gray')

# import sys
# sys.exit(0)


# Transfer volumes:
# Min: 2.5nL (multiples)
# Max: 400 nL
# Recommended: 20-60nL

## PP:
# 20-50uL (inside the well)

# LDV
# 4-12uL (inside the well)

# Stock solutions: 60uL

# Each pooling well (PP) will have 400nL x 64 = 25.6uL
# Each source well (PP) will feed upto 5 (LDV) data wells (5uL each)
# 256 data wells = 26 well acid & 26 well base


# Helper functions
#   simulate/run transfer
simulate = lambda x : x.run(verbose=False,
           volume_increment=2.5,
           enforce_volume_limits=False,
           robot='echo') 
#   convert number to string
num2str = lambda x, places=3: ("{:0."+str(places)+"f}").format(x)
# convert volume to an Echo transferrable volume
sanitizeVolume = lambda x: np.round(np.array(x) / min_vol_tranfer)*min_vol_tranfer
# convert weight to a Echo transferrable volume
weight2volume = lambda x: sanitizeVolume(max_vol_transfer * np.array(x))


###############################################################################    
print('\n\n########## SOURCE PLATE ########################################\n')
###############################################################################    
                
# Load source plate
source_plate = chemcpupy.WellPlate384PP(description='source')    
source_plate.load_platesheet(os.path.join(path_script, demo_source_platesheet),fill_from_csv=True,grid='letter')

# Visualize the plate 
if show_plates:    
    print('\n'+'masses')        
    x = source_plate.graphical_print_by_type().savefig('src.png')
    print('\n'+'volumes')            
    source_plate.graphical_print_volumes(units='uL')
    print('')

# Get info about the source compounds
acid_positions = source_plate.get_positions_by_compound_type('acid',grid='letter')
base_positions = source_plate.get_positions_by_compound_type('base',grid='letter')
indicator_positions = source_plate.get_positions_by_compound_type('indicator',grid='letter')
acid_info = source_plate.get_contents(locations=acid_positions,properties=['mass','name','cid','concentration','type','volume'])
acid_conc = acid_info['concentration'][0]
acid_volume = acid_info['volume'][0]
acid_cids = acid_info['cid'][0]
base_info = source_plate.get_contents(locations=base_positions,properties=['mass','name','cid','concentration','type','volume'])
base_conc = base_info['concentration'][0]
base_volume = base_info['volume'][0]
base_cids = base_info['cid'][0]
print('\n'+'Source acid concentration: ' + num2str(acid_conc) + 'M, well count: ' + str(len(acid_positions)) + ', volume per well: ' + num2str(acid_volume*1e-3) + 'uL')     
print('\n'+'Source base concentration: ' + num2str(base_conc) + 'M, well count: ' + str(len(base_positions)) + ', volume per well: ' + num2str(acid_volume*1e-3) + 'uL')     

###############################################################################    
print('\n\n########## DATA PLATE ##########################################\n')      
###############################################################################  

     
# Initialize plate positions
data_plate = chemcpupy.WellPlate1536LDV(description='data')
available_positions = data_plate.list_empty_positions()
used_positions = []


# Initialize AcidBaseNetwork
network = AcidBaseNetwork(num_neurons=2, weights_per_neuron=img_width*img_height, img_width=img_width, img_height=img_height)
network.load_weights(kernel_file)
network.load_image(sample_image1)

# Write the data wells
tasklist_data = chemcpupy.TaskList(description='Write Data')
# generate source and destination positions
pos_list_data_source = []
pos_list_data_destin = []
data_positions = []      
data_values = []
# write data to wells as an image
acid_source_index = 0
base_source_index = 0
indicator_source_index = 0

# print(acid_positions)
# print(base_positions)
# print("xxxxxxxxxxx")

pos_list_data_source, pos_list_data_destin = network.generate_weighted_data_plate(
                                    acid_positions=acid_positions,
                                    base_positions=base_positions,
                                    water_positions=[],
                                    source_max_transfers=source_max_transfers,
                                    starting_letter='A',
                                    starting_index=1,
                                    max_row=grid1536_height,
                                    max_col=grid1536_width,
                                    transfer_unit_vol=(vol1536_data * 1e-9),
                                    acid_conc=(acid_conc * 1e-3),
                                    base_conc=(base_conc * 1e-3))

# Add indidcator
pos_list_indicator_source, pos_list_indicator_destin = [], []
# indicator_pos = 0
# for d in pos_list_data_destin:
#     pos_list_indicator_source.append(indicator_positions[indicator_pos])
#     pos_list_indicator_destin.append(d)
#     indicator_pos = (indicator_pos + 1) % len(indicator_positions)


used_positions = used_positions + pos_list_data_destin + pos_list_indicator_destin
Plating.update_free_positions(available_positions,used_positions)                                                        
#    create the transfer task
task = chemcpupy.TransferTask(from_plate=source_plate,
                        from_positions=pos_list_data_source + pos_list_indicator_source,
                        to_plate=data_plate,
                        to_positions=pos_list_data_destin + pos_list_indicator_destin,
                        transfer_volumes=[vol1536_data] * len(pos_list_data_destin) + [indicator_transfer_volume] * len(pos_list_indicator_destin),
                        transfer_group_label='write')
tasklist_data.add(task)
# run the transfer task
print('\n'+'Simulating transfers')
simulate(tasklist_data)
tasklist_data.summarize()
if show_plates:
    print('\n\t'+'source plate:')
    source_plate.graphical_print_by_type()
    print('\n\t'+'data plate:')    
    data_plate.graphical_print_by_type().savefig('data.png')
    data_plate.graphical_print_volumes(units='uL')
chemcpupy.Echo.write_Echo_csv_picklist(tasklist_data, path_temp + '/' + demo_name + '-1_write_data-Echo_picklist' + img_name + '.csv')

###############################################################################    
print('\n\n########## pooling PLATE ########################################\n')
###############################################################################            

# Initialize plate positions
available_positions = source_plate.list_empty_positions()
used_positions = []


# Apply weights

pos_list_data_source, pos_list_data_destin = network.generate_summation_plate(source_image=pos_list_data_destin, starting_letter='O', starting_index=1, max_row=grid384_height, max_col=grid384_width)
expected_output = network.get_expected_outputs()

used_positions = used_positions + pos_list_data_destin     
Plating.update_free_positions(available_positions,used_positions)                                                        
#    create the transfer task
task = chemcpupy.TransferTask(from_plate=data_plate,
                        from_positions=pos_list_data_source,
                        to_plate=source_plate,
                        to_positions=pos_list_data_destin,
                        transfer_volumes=[vol_pool]*len(pos_list_data_destin),
                        transfer_group_label='pooling')

tasklist_data = chemcpupy.TaskList(description='Apply pooling')
tasklist_data.add(task)
# run the transfer task
print('\n'+'Simulating transfers')
print("Expected outputs", expected_output)
simulate(tasklist_data)
tasklist_data.summarize()     
if show_plates:    
    print('\n\t'+'pooling plate:')
    source_plate.graphical_print_by_type().savefig('pool.png')
    source_plate.graphical_print_volumes(units='uL')
chemcpupy.Echo.write_Echo_csv_picklist(tasklist_data, path_temp + '/' + demo_name + '-2_write_pooling-Echo_picklist' + img_name + '.csv')

