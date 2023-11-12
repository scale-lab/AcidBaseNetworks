#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 09:49:38 2018

@author: jacobrosenstein
"""

import csv

echo_csv_fieldnames = ['Source Plate Name',
                      'Source Plate Type',
                      'Source Well',
                      'Destination Plate Name',
                      'Destination Plate Type',
                      'Destination Well',
                      'Destination Well X Offset',
                      'Destination Well Y Offset',
                      'Delay',
                      'Transfer Volume']


def write_Echo_csv_picklist(mytasklist,CSVfilename):
        
    with open(CSVfilename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=echo_csv_fieldnames)
    
        writer.writeheader()
        for t in mytasklist._task_list:
            for r in _TransferTask_to_fields(t):
                writer.writerow(r)

    print('Wrote Echo CSV picklist:',CSVfilename)


def _TransferTask_to_fields(mytask):
    myfields = []
    
    num_transfers,from_plate,from_positions,to_plate,to_positions,transfer_volumes,pre_transfer_delays = mytask.unpack_transfers()
    
    for from_pos,to_pos,vol,delay in zip(from_positions,to_positions,transfer_volumes,pre_transfer_delays):
        myfields.append( {
                          'Source Plate Name' : from_plate._description,
                          'Source Plate Type' : '384PP_DMSO2',
                          'Source Well' : position_to_lettergrid(from_pos),
                          'Destination Plate Name' : to_plate._description,
                          'Destination Plate Type' : '384PP_Dest',
                          'Destination Well' : position_to_lettergrid(to_pos),
                          'Destination Well X Offset' : str(0),
                          'Destination Well Y Offset' : str(0),
                          'Delay' : round(delay), # round to nearest millisecond to be read properly                                                   
                          'Transfer Volume' : vol
                          } )
        
    return myfields


def position_to_lettergrid(position):
    # Returns a string corresponding to the given (zero indexed) numerical position. Rows are letters and columns are numbers, such that (1,3) returns 'B4'.
    rowletters = [chr(x) for x in range(ord('A'),ord('Z')+1)] + ['A' + chr(x) for x in range(ord('A'),ord('Z')+1)];    
    if type(position)==tuple:
        return rowletters[position[0]]+str(position[1]+1)
    elif type(position)==list:
        return [position_to_lettergrid(x) for x in position]
    else:
        position = position; # not recognized. assume it does not need to be converted
    return position


# todo: add option to sort by source and then destination (or the reverse) - to minimize stage movements




