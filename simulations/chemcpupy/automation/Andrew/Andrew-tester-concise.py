# -*- coding: utf-8 -*-
"""
@title: a more concise tester script for the Andrew class
@description: writes "BROWN DARPA" in a 384 well plate (based on 'example_protocol_2-unencrypted.anp')
@author: chrisarcadia 
@created: 2018/04/03
"""

import Andrew

# create a new protocol
recipe = Andrew.Protocol();

# set protocol information
recipe.name = 'Brown University DARPA';
recipe.author = 'Chris Arcadia';
recipe.author_email = 'christopher_arcadia@brown.edu';


# define the stock solutions

water = Andrew.StockSolution( name = 'water',
                              description = 'deionized water',
                              concentration = 1,
                              concentration_unit = 'a.u.');
                             
dye = Andrew.StockSolution(   name = 'red dye',
                              description = 'red food coloring diluted in water',                              
                              concentration = 1,
                              concentration_unit = 'a.u.');  
                                         
recipe.stock_solutions = [water,dye];

# define the library of consumables (objects help in Andrew's dominos)

plate384well = Andrew.LibraryConsumable.fromStandard('microplate384');                                                                                         

tube50mL = Andrew.LibraryConsumable.fromStandard('tube50');                                                                                                              
                                    
recipe.library_of_consumables = [plate384well, tube50mL];       

# set the contents of the consumables

plateTarget = Andrew.Consumable(    library_consumable = plate384well,                           
                                    label = '384-well plate',
                                    gui_position_x = 326,
                                    gui_position_y = 110,
                                    comment = 'a canvas on which to draw');
                                

tubeDye = Andrew.Consumable(        library_consumable = tube50mL,                           
                                    label = '50mL Falcon Tube (red dye)',
                                    gui_position_x = 143,
                                    gui_position_y = 90);     
                            
tubeDye.setSolution(    positions=[[1,1]],
                        volumes = [-1], # volume in nL
                        stock_solutions = [dye] );

                    
tubeWater = Andrew.Consumable(      library_consumable = tube50mL,                           
                                    label = '50mL Falcon Tube (water)',
                                    gui_position_x = 140,
                                    gui_position_y = 228);            
                              
tubeWater.setSolution(          positions=[[1,1]],
                                volumes = [-1], # volume in nL
                                stock_solutions = [water] );
                    
recipe.consumables = [plateTarget, tubeDye, tubeWater];

# create a list of actions to perform on the consumables' contents

### action 1 : pipette from falcon tube of dye to many wells on a plate (text foreground)
transfer1 = dict();
transfer1['source']     = {'consumable' : tubeDye,
                           'positions'  : [ [1,1] ]
                          };
transfer1['destination'] = {'consumable' : plateTarget, 
                            'positions'  : [ [2,1], [2,2], [2,3], [2,4], [2,6], [2,7], [2,8], [2,9], [2,11], [2,12], [2,13], [2,15], [2,19], [2,21], [2,24], [3,1], [3,4], [3,6], [3,9], [3,11], [3,13], [3,15], [3,19], [3,21], [3,24], [4,1], [4,2], [4,3], [4,6], [4,9], [4,11], [4,13], [4,15], [4,17], [4,19], [4,21], [4,22], [4,24], [5,1], [5,4], [5,6], [5,7], [5,8], [5,11], [5,13], [5,15], [5,17], [5,19], [5,21], [5,23], [5,24], [6,1], [6,4], [6,6], [6,9], [6,11], [6,13], [6,15], [6,17], [6,19], [6,21], [6,24], [7,1], [7,2], [7,3], [7,4], [7,6], [7,9], [7,11], [7,12], [7,13], [7,16], [7,18], [7,21], [7,24], [10,1], [10,2], [10,3], [10,4], [10,6], [10,7], [10,8], [10,9], [10,11], [10,12], [10,13], [10,14], [10,16], [10,17], [10,18], [10,19], [10,21], [10,22], [10,23], [10,24], [11,1], [11,4], [11,6], [11,9], [11,11], [11,14], [11,16], [11,19], [11,21], [11,24], [12,1], [12,4], [12,6], [12,9], [12,11], [12,14], [12,16], [12,19], [12,21], [12,24], [13,1], [13,4], [13,6], [13,7], [13,8], [13,9], [13,11], [13,12], [13,13], [13,16], [13,17], [13,18], [13,19], [13,21], [13,22], [13,23], [13,24], [14,1], [14,4], [14,6], [14,9], [14,11], [14,14], [14,16], [14,21], [14,24], [15,1], [15,2], [15,3], [15,4], [15,6], [15,9], [15,11], [15,14], [15,16], [15,21], [15,24] ]
                           };
transfer1['volume'] = [25000]; # in nL
transfer1['params'] = {
                         'tip'         : {'no_blow_out': False, 'change_before':True, 'change_between':False, 'verify':False, 'touch_off':False,  'filter':False},
                         'precision'   : {'required':False}, 
                         'operation'   : {'speed':'normal', 'alicots':'1', 'type':'forward'},
                         'viscosity'   : {'type':'normal'},
                         'cushion'     : {'top':False, 'bottom':False},
                         'source'      : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}},
                         'destination' : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}}
                      };

action1 = Andrew.Action(  
                          type        = 'pipette',
                          message     = 'filling foreground',
                          source      = transfer1['source'],
                          destination = transfer1['destination'],
                          volume      = transfer1['volume'],
                          params      = transfer1['params'],
                      );               
        
### action 2 : pipette from falcon tube of water to many wells on a plate (text background)
transfer2 = dict();
transfer2['source']     = {'consumable' : tubeWater,
                           'positions'  : [ [1,1] ]
                          };
transfer2['destination'] = {'consumable' : plateTarget, 
                            'positions'  : [ [1,1], [1,2], [1,3], [1,4], [1,5], [1,6], [1,7], [1,8], [1,9], [1,10], [1,11], [1,12], [1,13], [1,14], [1,15], [1,16], [1,17], [1,18], [1,19], [1,20], [1,21], [1,22], [1,23], [1,24], [2,5], [2,10], [2,14], [2,16], [2,17], [2,18], [2,20], [2,22], [2,23], [3,2], [3,3], [3,5], [3,7], [3,8], [3,10], [3,12], [3,14], [3,16], [3,17], [3,18], [3,20], [3,22], [3,23], [4,4], [4,5], [4,7], [4,8], [4,10], [4,12], [4,14], [4,16], [4,18], [4,20], [4,23], [5,2], [5,3], [5,5], [5,9], [5,10], [5,12], [5,14], [5,16], [5,18], [5,20], [5,22], [6,2], [6,3], [6,5], [6,7], [6,8], [6,10], [6,12], [6,14], [6,16], [6,18], [6,20], [6,22], [6,23], [7,5], [7,7], [7,8], [7,10], [7,14], [7,15], [7,17], [7,19], [7,20], [7,22], [7,23], [8,1], [8,2], [8,3], [8,4], [8,5], [8,6], [8,7], [8,8], [8,9], [8,10], [8,11], [8,12], [8,13], [8,14], [8,15], [8,16], [8,17], [8,18], [8,19], [8,20], [8,21], [8,22], [8,23], [8,24], [9,1], [9,2], [9,3], [9,4], [9,5], [9,6], [9,7], [9,8], [9,9], [9,10], [9,11], [9,12], [9,13], [9,14], [9,15], [9,16], [9,17], [9,18], [9,19], [9,20], [9,21], [9,22], [9,23], [9,24], [10,5], [10,10], [10,15], [10,20], [11,2], [11,3], [11,5], [11,7], [11,8], [11,10], [11,12], [11,13], [11,15], [11,17], [11,18], [11,20], [11,22], [11,23], [12,2], [12,3], [12,5], [12,7], [12,8], [12,10], [12,12], [12,13], [12,15], [12,17], [12,18], [12,20], [12,22], [12,23], [13,2], [13,3], [13,5], [13,10], [13,14], [13,15], [13,20], [14,2], [14,3], [14,5], [14,7], [14,8], [14,10], [14,12], [14,13], [14,15], [14,17], [14,18], [14,19], [14,20], [14,22], [14,23], [15,5], [15,7], [15,8], [15,10], [15,12], [15,13], [15,15], [15,17], [15,18], [15,19], [15,20], [15,22], [15,23], [16,1], [16,2], [16,3], [16,4], [16,5], [16,6], [16,7], [16,8], [16,9], [16,10], [16,11], [16,12], [16,13], [16,14], [16,15], [16,16], [16,17], [16,18], [16,19], [16,20], [16,21], [16,22], [16,23], [16,24] ]
                           };
transfer2['volume'] = [20000];
transfer2['params'] = {
                         'tip'         : {'no_blow_out': False, 'change_before':True, 'change_between':False, 'verify':False, 'touch_off':False,  'filter':False},
                         'precision'   : {'required':False}, 
                         'operation'   : {'speed':'normal', 'alicots':'1', 'type':'forward'},
                         'viscosity'   : {'type':'normal'},
                         'cushion'     : {'top':False, 'bottom':False},
                         'source'      : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}},
                         'destination' : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}}
                      };

action2 = Andrew.Action(  
                          type        = 'pipette',
                          message     = 'filling background',
                          source      = transfer2['source'],
                          destination = transfer2['destination'],
                          volume      = transfer2['volume'],
                          params      = transfer2['params'],
                      );      
        
                        
recipe.actions = [action1, action2];
Andrew.Action.reposition(recipe.actions);  # sort actions according to how the actions are listed (taking account of incubation steps by expanding each into two actions)                                  

# preview the recipe XML
print(recipe)

# save the recipe to file
recipe.write('recipe.anp'); # this file can be directely loaded into Andrew Lab


