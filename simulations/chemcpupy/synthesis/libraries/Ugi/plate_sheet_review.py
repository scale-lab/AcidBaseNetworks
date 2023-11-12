#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon August 6 2018

@author: chrisarcadia
@description: script to review a library platesheet and to optionally 
              generate a mass list and Andrew recipe for that library
              
"""

# Import python libraries
import chemcpupy 
import csv
import copy

#%%############################################################################    
print('\n\n########## REVIEW LIBRARY PLATE ################################\n')
###############################################################################    

# Select platesheet
plate_sheet = 'library_4';

# Script options
export_masses = True;
generate_recipe_for_Andrew = True;
generate_recipe_for_Echo = False; # not yet implemented 
visualize_plate = True;

# Make temperorary directory
temp_dir = chemcpupy.Path.make_temp_dir(); 

# Load the platesheet
path = chemcpupy.Path.get_dir(__file__);
filename = path+'/'+plate_sheet+'_platesheet.csv';
plate = chemcpupy.WellPlate384PP(description='library plate')
plate.load_platesheet(filename)
print('\n'+'Loaded "' + filename + '"')
plate_initial = copy.deepcopy(plate); # before reaction

# Generate Ugi Compounds
chemcpupy.UgiLibrary.generate_ugi_products(plate,reaction_yield=1.0)    

# Order to add compounds (by type)
ordered_of_addition = ['amine','aldehyde','carboxylic acid','isocyanide'];      

# Visualize the plate 
if visualize_plate:
        
    print('\n'+'Plotting plate before Ugi reactions')
    print('\n'+' - '+'amines           - monoisotopic masses [g/mol]:')            
    plate_initial.graphical_print_by_id_with_type_filter('amine')
    print('\n'+' - '+'aldehydes        - monoisotopic masses [g/mol]:')            
    plate_initial.graphical_print_by_id_with_type_filter('aldehyde')
    print('\n'+' - '+'carboxylic acids - monoisotopic masses [g/mol]:')            
    plate_initial.graphical_print_by_id_with_type_filter('carboxylic acid')
    print('\n'+' - '+'isocyanides      - monoisotopic masses [g/mol]:')            
    plate_initial.graphical_print_by_id_with_type_filter('isocyanide')
        
    print('\n'+'Plotting plate after Ugi reactions')
    print('\n'+' - '+'compounds [#]:')            
    plate.graphical_print_compounds_per_well()
    print('\n'+' - '+'types:')                
    plate.graphical_print_by_type()    
    print('\n'+' - '+'monoisotopic masses [g/mol]:')        
    plate.graphical_print_by_mass()
    print('\n'+' - '+'volumes [uL]:')            
    plate.graphical_print_volumes(units='uL')
    print('\n'+' - '+'IDs [PubChem]:')            
    plate.graphical_print_by_id()
    
#%%############################################################################    
print('\n\n########## TESTING COMMON FUNCTIONS ############################\n')
###############################################################################        

# Get empty and filled positions
pos_filled = plate.list_filled_positions();
pos_empty = plate.list_empty_positions();

# View contents of a well
well_name = 'A1'; # pos_filled[0];
well_pos = chemcpupy.Containers.lettergrid_to_position(well_name); 
well_contents = plate.get_compounds(well_pos,['type','name','cid','mass']);
print('Contents of well' + well_name + ':')
print(well_contents)
    
#%%############################################################################    
print('\n\n########## EXPORT DATA #########################################\n')
###############################################################################  
    
# Export compound info to CSV
if export_masses:
    
#    # Find wells with a particular compound type and get compound info
#    compound_info = plate.get_contents(type_filter='ugi',grid='letter');
#    compound_positions = compound_info['position'];
#    compound_masses = compound_info['mass'];
#    compound_names = compound_info['name'];
#    compound_ids = compound_info['cid'];
#    
#    # Write to CSV file
#    CSVfilename = temp_dir + '/'+plate_sheet+'_review.csv';     
#    csv_fieldnames = ['position',
#                      'mass'];
#    with open(CSVfilename, 'w', newline='') as csvfile:
#        writer = csv.DictWriter(csvfile, fieldnames=csv_fieldnames)
#        writer.writeheader()
#        for n in range(0,len(compound_positions)):
#            rowdata = {'mass':compound_masses[n],
#                       'position':compound_positions[n],};
#            writer.writerow(rowdata)  
#    print('Wrote compound info CSV to file: ' + CSVfilename)
        
    print('\n'+'Exporting Compound Info')     
    plate_initial.export_csv_review(temp_dir + '/'+plate_sheet+'_reagent'+'_review.csv');
    plate.export_csv_review(temp_dir + '/'+plate_sheet+'_review.csv');    
    plate.export_csv_review(temp_dir + '/'+plate_sheet+'_ugi_review.csv',typefilter='ugi');
            
# or equivalently: 
            
#   plate.export_csv_review(CSVfilename,typefilter='ugi');
    
# Generate recipe to make the library with Andrew
if generate_recipe_for_Andrew:
        
    # Get filled positions
    pos_filled = plate_initial.list_filled_positions();
        
    # Get unique compounds
    reag_ids = [];
    reag_infos= {};
    pos_offset = (1,1); # Andrew's row/column numbering starts at 1, not 0    
    for pos in pos_filled:
       infolist = plate_initial.get_compounds(pos,['type','name','cid']);
       shiftedpos = (pos[0]+pos_offset[0], pos[1]+pos_offset[1]);
       for info in infolist:
           cid = info[2];
           if cid not in reag_infos:
               reag_infos.update({cid:{'name':info[1],'type':info[0],'positions':[shiftedpos]}});
               reag_ids.append(cid);
           else:
               reag_infos[cid]['positions'].append(shiftedpos);
    num_reag = len(reag_ids);
          
    # Initialize the protocol
    recipe = chemcpupy.Andrew.Protocol();
    recipe.name = plate_sheet + ' (from platesheet)';
    recipe.author = 'Chris Arcadia';
    recipe.author_email = 'christopher_arcadia@brown.edu';

    # Define the stock solutions
    reagents = [None] * num_reag;   
    for n in range(0,len(reag_ids)):    
        cid = reag_ids[n];
        cname = reag_infos[cid]['name'];
        ctype = reag_infos[cid]['type'];
        reagents[n] = chemcpupy.Andrew.StockSolution(  
                name = cname,
                description = '(CID: ' + str(cid) + ') ' + str(cname),                              
                concentration = 1,
                concentration_unit = 'a.u.'
            );  
    recipe.stock_solutions = reagents;

    # Define the library of consumables
    plate384well = chemcpupy.Andrew.LibraryConsumable.fromStandard('microplate384');                                                                                         
    tube2mL = chemcpupy.Andrew.LibraryConsumable.fromStandard('microtube20');     
    recipe.library_of_consumables = [plate384well, tube2mL];       

    # Set the contents of the consumables
    counter = 0;
    type_counter = {}; # used to sort consumbale position in the GUI by compound type
    reagentVials = [None] * num_reag;   
    for n in range(0,num_reag):   
        cid = reag_ids[n];
        cname = reag_infos[cid]['name'];
        ctype = reag_infos[cid]['type'];   
        if ctype in type_counter:            
                type_counter[ctype]['xoff'] = type_counter[ctype]['xoff'] + 1;
        else:
                type_counter.update({ctype:{'yoff':counter,'xoff':0}});
                counter = counter + 1;
        reagentVials[n] = chemcpupy.Andrew.Consumable(   
            gui_position_x = 50 + type_counter[ctype]['xoff']*125,
            gui_position_y = 175 + type_counter[ctype]['yoff']*100,
            library_consumable = tube2mL,                       
            label =  cname,
            comment = '2mL Test Tube (#' + str(n+1) + ') ' + cname + ' (CID: ' + str(cid) + ')', 
        );                                                                                      
        reagentVials[n].setSolution(      
            positions= [[1,1]],
            volumes = [-1], 
            stock_solutions = [reagents[n]],
        );                             
    reactionPlate = chemcpupy.Andrew.Consumable(   
            gui_position_x = 50,
            gui_position_y = 25,
            library_consumable = plate384well,                           
            label = 'reaction plate',
            comment = '384-well plate for Ugi reactions'
    );                            
    recipe.consumables = reagentVials + [reactionPlate];


    # Pipetting parameters    
    pipette_single_use_viscous =  {
         'tip'         : {'no_blow_out': False, 'change_before':True, 'change_between':True, 'verify':False, 'touch_off':True,  'filter':False},
         'precision'   : {'required':False}, 
         'operation'   : {'speed':'normal', 'alicots':'1', 'type':'forward'},
         'viscosity'   : {'type':'high'},
         'cushion'     : {'top':True, 'bottom':False},
         'source'      : {'pipette_position':{'height':0, 'type':'liquid', 'policy':'default'}, 'mixing':{'speed':'normal','number':0}},
         'destination' : {'pipette_position':{'height':0, 'type':'bottom', 'policy':'default'}, 'mixing':{'speed':'normal','number':2}}
    };

    # Simulate transfers  
    transfer_actions = [];                     
    for otype in ordered_of_addition:
        type_write = [];        
        for n in range(0,num_reag):   
            cid = reag_ids[n];
            cname = reag_infos[cid]['name'];
            ctype = reag_infos[cid]['type'];  
            if ctype == otype:                                                                            
                transfer = dict();
                transfer['source']     = {'consumable' : reagentVials[n],
                                           'positions'  : [ [1,1] ]
                                          };
                transfer['destination'] = {'consumable' : reactionPlate, 
                                            'positions'  : reag_infos[reag_ids[n]]['positions'],
                                           };
                transfer['volume'] = [10e3]; # [nL]
                transfer['params'] = pipette_single_use_viscous;
                action = chemcpupy.Andrew.Action(  
                      type        = 'pipette',
                      message     = 'writing '+cname,
                      source      = transfer['source'],
                      destination = transfer['destination'],
                      volume      = transfer['volume'],
                      params      = transfer['params'],
                );     
                type_write.append(action);                
                      
        type_write_alert = chemcpupy.Andrew.Action(  
              type        = 'alert',
              message     = 'writing '+otype,
        );
        transfer_actions = transfer_actions + [type_write_alert] + type_write;                

    recipe.actions = transfer_actions;
    chemcpupy.Andrew.Action.reposition(recipe.actions);  # sort actions according to how the actions are listed (taking account of incubation steps by expanding each into two actions)                                  
    
    # Save the recipe to file
    filename = temp_dir + '/'+plate_sheet+'_recipe.anp';
    recipe.write(filename); # this file can be directely loaded into Andrew Lab
    #print(recipe)    
    print('Wrote Andrew recipe to file: ' + filename)

#%%############################################################################    
print('\n\n################################################################\n')
###############################################################################    
