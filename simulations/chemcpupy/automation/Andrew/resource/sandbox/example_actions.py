#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 03:51:43 2018

@author: chrisarcadia
"""

import Andrew
from lxml import etree

# alert action
alert = Andrew.Action(  position = 1,
                        type = 'alert',
                        level = 0,
                        message = 'alert! alert!' );
# user action
user = Andrew.Action(   position = 1,
                        type = 'user',
                        level = 0,
                        message = 'waiting for user',
                        expected_time = 60 ); # seconds

# incubation action                     
incubate = Andrew.Action( position = 1,
                          type = 'incubate',
                          level = 0,
                          message = 'incubating samples',
                          duration = 300, # seconds
                          position_offset = 1); 
                         
# pipetting action                                                                      
water = Andrew.StockSolution( id = '35651B84-7CF4-CAC7-59A9-43E3E7399E96',
                              name = 'water',
                              description = 'deionized water',
                              color = 16766720,
                              concentration = 1,
                              concentration_unit = 'a.u.');
                             
dye = Andrew.StockSolution(   id = '8ED34505-69B1-5BC7-A13D-40B01B09E97D',
                              name = 'red dye',
                              description = 'red food coloring diluted in water',                              
                              color = 16738740,
                              concentration = 1,
                              concentration_unit = 'a.u.');                                           
plate384well = Andrew.LibraryConsumable(     shortname='MTP384',
                                             type = 'microplate384',
                                             preparation = 'group',
                                             version = 2,
                                             id = 'microplate384',
                                             name = '384-well plate',
                                             display_name = 'MTP384',
                                             description = 'Generic microplate with 384 wells. This consumable describes typical square-shaped wells of capacity below 120 uL - with exclusion of 384-wells PCR microplates',
                                             category = 'Microplates',
                                             volume_dead = 10000, # [nL]
                                             volume_max = 120000,  # [nL]
                                             size_columns = 24,
                                             size_rows = 16);
 
tube50mL = Andrew.LibraryConsumable(     shortname='TB50',
                                         type = 'tube50',
                                         preparation = 'individual',
                                         version = 2,
                                         id = 'tube50',
                                         name = '50 mL conical centrifuge tube',
                                         display_name = 'TUBE 50',
                                         description = 'Plastic 50 mL tube, conical bottom (without skirt) and screwable cap. Also called Falcon tube or centrifuge tube',
                                         category = 'Tubes',
                                         volume_dead =  500000, # [nL]
                                         volume_max = 50000000,  # [nL]
                                         size_columns = 1,
                                         size_rows = 1 );   
plateTarget = Andrew.Consumable(    id = 'F9DB029D-9778-1EC4-7EA4-3F1BC3F0D649',
                                    gui_position_x = 326,
                                    gui_position_y = 110,
                                    library_consumable = plate384well,                           
                                    label = '384-well plate',
                                    comment = 'a canvas on which to draw');                                
tubeDye = Andrew.Consumable(        id = '8CF4F670-48DD-E62D-EE18-3F1BDA813E43',
                                    gui_position_x = 143,
                                    gui_position_y = 90,
                                    library_consumable = tube50mL,                           
                                    label = '50mL Falcon Tube (red dye)',
                                    comment = '');                                 
tubeDye.setSolution(    positions=[[1,1]],
                        volumes = [-1], # volume in μL
                        stock_solutions = [dye] );                    
tubeWater = Andrew.Consumable(      id = '516A731A-2BC7-0321-742B-3F2C542953FF',
                                    gui_position_x = 140,
                                    gui_position_y = 228,
                                    library_consumable = tube50mL,                           
                                    label = '50mL Falcon Tube (water)',
                                    comment = '');                                          
tubeWater.setSolution(          positions=[[1,1]],
                                volumes = [-1], # volume in μL
                                stock_solutions = [water] );                                                                                                                                                      
transfer = dict();
transfer['source']     = {'consumable' : tubeDye,
                           'positions'  : [ [1,1] ]
                          };
transfer['destination'] = {'consumable' : plateTarget, 
                            'positions'  : [ [2,1], [2,2], [2,3], [2,4], [2,6], [2,7], [2,8], [2,9], [2,11], [2,12], [2,13], [2,15], [2,19], [2,21], [2,24], [3,1], [3,4], [3,6], [3,9], [3,11], [3,13], [3,15], [3,19], [3,21], [3,24], [4,1], [4,2], [4,3], [4,6], [4,9], [4,11], [4,13], [4,15], [4,17], [4,19], [4,21], [4,22], [4,24], [5,1], [5,4], [5,6], [5,7], [5,8], [5,11], [5,13], [5,15], [5,17], [5,19], [5,21], [5,23], [5,24], [6,1], [6,4], [6,6], [6,9], [6,11], [6,13], [6,15], [6,17], [6,19], [6,21], [6,24], [7,1], [7,2], [7,3], [7,4], [7,6], [7,9], [7,11], [7,12], [7,13], [7,16], [7,18], [7,21], [7,24], [10,1], [10,2], [10,3], [10,4], [10,6], [10,7], [10,8], [10,9], [10,11], [10,12], [10,13], [10,14], [10,16], [10,17], [10,18], [10,19], [10,21], [10,22], [10,23], [10,24], [11,1], [11,4], [11,6], [11,9], [11,11], [11,14], [11,16], [11,19], [11,21], [11,24], [12,1], [12,4], [12,6], [12,9], [12,11], [12,14], [12,16], [12,19], [12,21], [12,24], [13,1], [13,4], [13,6], [13,7], [13,8], [13,9], [13,11], [13,12], [13,13], [13,16], [13,17], [13,18], [13,19], [13,21], [13,22], [13,23], [13,24], [14,1], [14,4], [14,6], [14,9], [14,11], [14,14], [14,16], [14,21], [14,24], [15,1], [15,2], [15,3], [15,4], [15,6], [15,9], [15,11], [15,14], [15,16], [15,21], [15,24] ]
                           };
transfer['volume'] = [25000];
transfer['params'] = {
                         'tip'         : {'no_blow_out': False, 'change_before':True, 'change_between':False, 'verify':False, 'touch_off':False,  'filter':False},
                         'precision'   : {'required':False}, 
                         'operation'   : {'speed':'normal', 'alicots':'-1', 'type':'forward'},
                         'viscosity'   : {'type':'normal'},
                         'cushion'     : {'top':False, 'bottom':False},
                         'source'      : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}},
                         'destination' : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}}
                      };
pipette = Andrew.Action(  
                          position    = 1,
                          type        = 'pipette',
                          level       = 0,
                          message     = 'pipetting',
                          source      = transfer['source'],
                          destination = transfer['destination'],
                          volume      = transfer['volume'],
                          params      = transfer['params'],
                      );                                                 
                  
# all actions               
actions = [user, alert, incubate, pipette];
#print('A[2]=' + str((actions[2]).position) + '-' + str((actions[2]).type) + '  A[3]=' + str((actions[3]).position) + '-' + str((actions[3]).type) + '  N=' + str(len(actions)))
Andrew.Action.reposition(actions);  # sort actions according to how the actions are listed (taking account of incubation steps)
#print('A[2]=' + str((actions[2]).position) + '-' + str((actions[2]).type) + '  A[3]=' + str((actions[3]).position) + '-' + str((actions[3]).type) + '  A[4]=' + str((actions[4]).position) + '-' + str((actions[4]).type) + '  N=' + str(len(actions)))
x = Andrew.Action.toXML(actions[4])
print(etree.tostring(x, encoding = "unicode", pretty_print=True))
            