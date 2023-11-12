#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chris3arcadia 
@description: template script to generate transfer lists for acid/base computing
@created: Thr Feb 28 2019
"""

# Libraries
import numpy
import chemcpupy
from chemcpupy import Plating, Containers
from acidbase_network import AcidBaseNetwork
import os
import sys

dirname = os.path.dirname(__file__)


datatype = sys.argv[1]

assert datatype in ['bin', '3bit'] 

imgsize = int(sys.argv[2])

assert(imgsize in [8, 12, 16, 28])

num_classes = 2
if len(sys.argv) >= 3:
	num_classes = int(sys.argv[3])


# Set paths
path_script = chemcpupy.Path.get_dir(__file__) # get current script's directory
path_temp = chemcpupy.Path.make_temp_dir() # make temporary directory

# Script options
# show_plates = True # plot images of plate contents
show_plates = False # plot images of plate contents
simulate_seperately = True # simulate tasklists seperately
silent = True

data_path = os.path.join(dirname, './datasets/dataset_{}'.format(datatype))
# Demo information
demo_name = 'acidbase_demo2022'
demo_file = 'demo_export.mat'
demo_source_platesheet = 'platesheet.csv'
dataset_path = data_path + '/img{}'.format(imgsize)
# kernel_file = os.path.join(data_path, 'kernel_mem_{}_{}.txt'.format(num_classes, imgsize))
kernel_file = os.path.join(data_path, 'kernel_mem_{}_{}.txt'.format(num_classes, imgsize))
neurons_file = os.path.join(data_path, 'neurons_{}_{}.txt'.format(num_classes, imgsize))
labels_file = os.path.join(data_path, 'labels_{}_{}.txt'.format(num_classes, imgsize))

dataset_imgs = os.listdir(dataset_path)



num_map = {
	0: 'zero',
	1: 'one',
	2: 'two',
	3: 'three',
	4: 'four',
	5: 'five',
	6: 'six',
	7: 'seven',
	8: 'eight',
	9: 'nine'
}
invalid_cnt = 0
counters = {}
correct = {}
tf_match = {}
print_every = 50
img_prefix = {}
for digit in range(num_classes):
	counters[digit] = 0
	correct[digit] = 0
	tf_match[digit] = 0
	img_prefix[digit] = 'img_' + num_map[digit]

dataset_labels = {}
for img in dataset_imgs:
    lbl = None
    for dig, name in img_prefix.items():
        if img.startswith(name):
            lbl = dig
            break
    if lbl is not None:
        dataset_labels[img] = lbl

tf_outputs = {}
tf_labels = {}

with open(neurons_file, mode='r') as neurons_file_f:
    with open(labels_file, mode='r') as labels_file_f:
        neurons_file_lines = neurons_file_f.read().splitlines()
        labels_file_lines = labels_file_f.read().splitlines()
        for i, n_line in enumerate(neurons_file_lines):
            l_line = labels_file_lines[i]
            n_parts = n_line.split(',')
            lbl = l_line.split(',')[1]
            img_file = n_parts[0]
            n_outputs = [float(outp) for outp in n_parts[1:]]
            tf_labels[img_file] = int(lbl)
            tf_outputs[img_file] = n_outputs


# sample_image1 = os.path.join(data_dir, 'bin_images_256/img_zero_71.txt')

# zero_1.txt -> [[4, -4], [-20, 20]] -> (AB, BA) -> (1, -1)
# one_1.txt -> [[-34, 34], [18, -18]] -> (BA, AB) -> (-1, 1)


# Global Variables
min_vol_tranfer = 2.5 # minimum transferrable volume [nL]
max_vol_transfer = 10000 # max transferrable volume[nL]
vol1536_data = 2000 # volume to transfer from source to data plate [nL]
vol_pool = 200 # volume to transfer to pool from a data well
source_max_transfers = 5 # maximum transfers from each source well

indicator_transfer_volume = 100 # Add 100nl to each well in the data plate


grid384_height = float('inf')
grid384_width = float('inf')

grid1536_height = float('inf')
grid1536_width = float('inf')

img_width = imgsize
img_height = imgsize
imglevel = [-1, 1, 3, 5]


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
sanitizeVolume = lambda x: numpy.round(numpy.array(x) / min_vol_tranfer)*min_vol_tranfer
# convert weight to a Echo transferrable volume
weight2volume = lambda x: sanitizeVolume(max_vol_transfer * numpy.array(x))

num_ones = 0
num_zeroes = 0

def print_summary():
    sum_correct = 0
    sum_data = 0
    for dig, corr in correct.items():
        if counters[dig] > 0:
            print("Correct", num_map[dig], ':', corr, 'Accuracy:', round(100 * (corr/counters[dig]), 2))
        else:
            print("Correct", num_map[dig], ':', corr, 'Accuracy:', 'N/A')
        sum_correct = sum_correct + corr
        sum_data = sum_data + counters[dig]
    print("Invalid Output", invalid_cnt)
    print("Accuracy", round(100 *((sum_correct) / (sum_data)), 2))


data_cnt = 1


wrong_cnt = 0

use_tf_label = True

for img, lbl in dataset_labels.items():
    # print('evaluating', img, '(', (data_cnt + 1), '/', len(dataset_labels), ')')


    data_cnt = data_cnt + 1
    if use_tf_label:
        lbl = tf_labels[img]

    
    
    
    counters[lbl] = counters[lbl] + 1

    sample_image1 = os.path.join(dataset_path, img)

    if not silent:
        ###############################################################################    
        print('\n\n########## SOURCE PLATE ########################################\n')
        ###############################################################################    
                    
    # Load source plate
    source_plate = chemcpupy.WellPlate384PP(description='source')    
    source_plate.load_platesheet(os.path.join(path_script,demo_source_platesheet),fill_from_csv=True,grid='letter')

    # Visualize the plate 
    if show_plates:    
        print('\n'+'masses')        
        source_plate.graphical_print_by_mass()
        print('\n'+'volumes')            
        source_plate.graphical_print_volumes(units='uL')
        print('')


    # Get info about the source compounds
    acid_positions = source_plate.get_positions_by_compound_type('acid',grid='letter')
    base_positions = source_plate.get_positions_by_compound_type('base',grid='letter')
    water_positions = source_plate.get_positions_by_compound_type('water',grid='letter')
    indicator_positions = source_plate.get_positions_by_compound_type('indicator',grid='letter')
    acid_info = source_plate.get_contents(locations=acid_positions,properties=['mass','name','cid','concentration','type','volume'])
    acid_conc = acid_info['concentration'][0]
    acid_volume = acid_info['volume'][0]
    acid_cids = acid_info['cid'][0]
    base_info = source_plate.get_contents(locations=base_positions,properties=['mass','name','cid','concentration','type','volume'])
    base_conc = base_info['concentration'][0]
    base_volume = base_info['volume'][0]
    base_cids = base_info['cid'][0]
    if not silent:
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
    network = AcidBaseNetwork(num_neurons=num_classes, weights_per_neuron=img_width*img_height, img_width=img_width, img_height=img_height)
    network.load_weights(kernel_file)
    if datatype == '3bit':
        network.load_image(sample_image1, img_levels=[-4, -3, -2, -1, 1, 2, 3, 4])
    else:
        network.load_image(sample_image1)

    # Write the data wells
    if not silent:
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

    pos_list_data_source, pos_list_data_destin = network.generate_weighted_data_plate(
                                        acid_positions=acid_positions,
                                        base_positions=base_positions,
                                        water_positions=water_positions,
                                        source_max_transfers=source_max_transfers,
                                        starting_letter='A',
                                        starting_index=1,
                                        max_row=grid1536_height,
                                        max_col=grid1536_width,
                                        transfer_unit_vol=(vol1536_data * 1e-9),
                                        acid_conc=(acid_conc * 1e-3),
                                        base_conc=(base_conc * 1e-3))

    pos_list_indicator_source, pos_list_indicator_destin = [], []


    used_positions = used_positions + pos_list_data_destin + pos_list_indicator_destin
    Plating.update_free_positions(available_positions,used_positions)                 
    if not silent:
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
        chemcpupy.Echo.write_Echo_csv_picklist(tasklist_data, path_temp + '/' + demo_name + '-1_write_data-Echo_picklist.csv')
    if show_plates:
        print('\n\t'+'source plate:')
        source_plate.graphical_print_by_mass()
        print('\n\t'+'data plate:')    
        data_plate.graphical_print_by_mass()
        data_plate.graphical_print_volumes(units='uL')

    if not silent:
        ###############################################################################    
        print('\n\n########## pooling PLATE ########################################\n')
        ###############################################################################            

    # Initialize plate positions
    available_positions = source_plate.list_empty_positions()
    used_positions = []


    # Apply weights

    pos_list_data_source, pos_list_data_destin = network.generate_summation_plate(source_image=pos_list_data_destin, starting_letter='O', starting_index=1, max_row=grid384_height, max_col=grid384_width)
    expected_output = network.get_expected_outputs()



    invalid = False
    already_found = False

    neurons = []
    for e in expected_output:
        if e[0] < 7.0 and e[1] > 7.0:
            neurons.append(1) 
            if already_found:
                invalid = True
            else:
                already_found = True
        elif (e[0] > 7.0 and e[1] < 7.0) or (e[0] == 7.0 and e[1] == 7.0):
            neurons.append(0)
        else:
            neurons.append(-1)
            invalid = True

    is_correct = True
    for n, neur in enumerate(neurons):
        if (neur == 1 and tf_outputs[img][n] <= 0) or (neur == 0 and tf_outputs[img][n] > 0):
            is_correct = False
            break
    if not invalid and is_correct:
        correct[lbl] = correct[lbl] + 1
    else:
        if invalid:
            invalid_cnt = invalid_cnt + 1
        wrong_cnt = wrong_cnt + 1


    used_positions = used_positions + pos_list_data_destin     
    Plating.update_free_positions(available_positions,used_positions)                                                        
    #    create the transfer task
    if not silent:
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
        chemcpupy.Echo.write_Echo_csv_picklist(tasklist_data, path_temp + '/' + demo_name + '-2_write_pooling-Echo_picklist.csv')
    if show_plates:    
        print('\n\t'+'pooling plate:')
        source_plate.graphical_print_by_mass()
        source_plate.graphical_print_volumes(units='uL')
    if not silent and data_cnt > 0 and ((data_cnt % print_every) == 0):
        print('evaluating', img, '(', (data_cnt + 1), '/', len(dataset_labels), ')')
        print_summary()

print("===========")

print_summary()