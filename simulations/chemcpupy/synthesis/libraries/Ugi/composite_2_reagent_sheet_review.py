#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur September 05 2018

@author: chrisarcadia
@description: script to review a composite library sourceplatesheet and to optionally 
              generate a mass list and Echo recipe for that library 
              
              Composite Library Style #2
              - each reaction has slots for 1xAmine, 1xAldehyde, 1xCarboxylicAcids, 5xIsocyanides 
              - all combinations (with replacement) are generated              
              - the order of inclusion is Amine, Aldehyde, CarboxylicAcid, Isocyanide_1, Isocyanide_2, Isocyanide_3, Isocyanide_4
              - the volume included (following the above order) is 200, 200, 200, 200, 200, 200, 200 nL
              - second repetition on the plate:
              - the order of inclusion is Isocyanide_1, Isocyanide_2, Isocyanide_3, Isocyanide_4, Amine, Aldehyde, CarboxylicAcid
              - the volume included (following the above order) is 200, 200, 200, 200, 200, 200, 200 nL
                  
"""

# Import python libraries
import chemcpupy
import copy
import math
import scipy.io
from datetime import datetime

#%%############################################################################    
print('\n\n########## REVIEW SOURCE/REAGENT PLATE #########################\n')
###############################################################################    

# Select sourceplatesheet
source_plate_sheet = 'library_9';

# Script options
export_masses = True;
generate_recipe_for_Andrew = False; # not yet implemented 
generate_recipe_for_Echo = True; 
generate_picklist_for_Solarix = True;
export_as_platesheet = False; # not yet implemented
visualize_plate = True;
save_variables = True;

# Make temperorary directory
temp_dir = chemcpupy.Path.make_temp_dir(); 
      
# Load the sourceplatesheet
path = chemcpupy.Path.get_dir(__file__);
filename = path+'/'+source_plate_sheet+'_sourceplatesheet.csv';
plateSource = chemcpupy.WellPlate384PP(description='source plate')
plateSource.load_platesheet(filename)
source_working_volume = 30e3; # [nL] (for 384PP)
print('\n'+'Loaded "' + filename + '"')

# Visualize the source plate 
if visualize_plate:
        
    print('\n'+'Plotting source plate')
    print('\n'+' - '+'compounds [#]:')            
    plateSource.graphical_print_compounds_per_well()
    print('\n'+' - '+'types:')                
    plateSource.graphical_print_by_type()
    print('\n'+' - '+'monoisotopic masses [g/mol]:')        
    plateSource.graphical_print_by_mass()
    print('\n'+' - '+'volumes [uL]:')            
    plateSource.graphical_print_volumes(units='uL')
    print('\n'+' - '+'amines           - monoisotopic masses [g/mol]:')            
    plateSource.graphical_print_by_id_with_type_filter('amine')
    print('\n'+' - '+'aldehydes        - monoisotopic masses [g/mol]:')            
    plateSource.graphical_print_by_id_with_type_filter('aldehyde')
    print('\n'+' - '+'carboxylic acids - monoisotopic masses [g/mol]:')            
    plateSource.graphical_print_by_id_with_type_filter('carboxylic acid')
    print('\n'+' - '+'isocyanides      - monoisotopic masses [g/mol]:')            
    plateSource.graphical_print_by_id_with_type_filter('isocyanide')
    print('\n'+' - '+'solvent      - monoisotopic masses [g/mol]:')            
    plateSource.graphical_print_by_id_with_type_filter('solvent')
    
# Find wells with a particular compound type 
pos_sol = plateSource.get_positions_by_compound_type('solvent');  
pos_mat = plateSource.get_positions_by_compound_type('matrix');  
pos_ali = plateSource.get_positions_by_compound_type('aligner');  
num_sol = len(pos_sol);    
num_mat = len(pos_mat);    
num_ali = len(pos_ali);    
#pos_ami = plateSource.get_positions_by_compound_type('amine');
#pos_ald = plateSource.get_positions_by_compound_type('aldehyde');
#pos_car = plateSource.get_positions_by_compound_type('carboxylic acid');
#pos_iso = plateSource.get_positions_by_compound_type('isocyanide');  
props = ['mass','name','cid','type','formula','structure'];
con_ami = plateSource.get_contents(type_filter='amine',properties=props);
con_ald = plateSource.get_contents(type_filter='aldehyde',properties=props);
con_car = plateSource.get_contents(type_filter='carboxylic acid',properties=props);
con_iso = plateSource.get_contents(type_filter='isocyanide',properties=props);
#cid_ami = con_ami['cid'];
#cid_ald = con_ald['cid'];
#cid_car = con_car['cid'];
#cid_iso = con_iso['cid'];
#uid_ami = list(set(con_ami['cid']));
#uid_ald = list(set(con_ald['cid']));
#uid_car = list(set(con_car['cid']));
#uid_iso = list(set(con_iso['cid']));
#num_ami = len(uid_ami);
#num_ald = len(uid_ald);
#num_car = len(uid_car);
#num_iso = len(uid_iso);
#num_ugi = num_ami*num_ald*num_car*num_iso;
#pos_ami = con_ami['position'];
#pos_ald = con_ald['position'];
#pos_car = con_car['position'];
#pos_iso = con_iso['position'];
def group_content_by_compound(contents):    
    list_of_ids = contents['cid'];
    list_of_positions = contents['position'];     
    ids = list(set(list_of_ids));
    count = len(ids);
    well_counts = []; wells = []; masses = []; names = []; formula = []; structure = [];
    for n in range(0,count):
        indices = [ind for ind in range(len(list_of_ids)) if list_of_ids[ind]==ids[n]];
        well_counts.append(len(indices));
        wells.append([list_of_positions[ind] for ind in indices]);
        masses.append(contents['mass'][indices[0]]);
        names.append(contents['name'][indices[0]]);
        formula.append(contents['formula'][indices[0]]);
        structure.append(contents['structure'][indices[0]]);        
    return ids,names,masses,wells,formula,structure,count,well_counts;
[uid_ami, nam_ami, mas_ami, loc_ami, for_ami, str_ami, num_ami, well_count_ami] = group_content_by_compound(con_ami);
[uid_ald, nam_ald, mas_ald, loc_ald, for_ald, str_ald, num_ald, well_count_ald] = group_content_by_compound(con_ald);
[uid_car, nam_car, mas_car, loc_car, for_car, str_car, num_car, well_count_car] = group_content_by_compound(con_car);
[uid_iso, nam_iso, mas_iso, loc_iso, for_iso, str_iso, num_iso, well_count_iso] = group_content_by_compound(con_iso);
num_ugi = num_ami*num_ald*num_car*num_iso;

# Check if we will be adding alignment positions
add_alignment_spots = (pos_ali!=[]); # reserve several positions on the plate, on the outer perimeter of the plate, for spotting aligner (matrix whose crystals that can be readily seen under the Solarix microscope)

# Collect reagent information
reag = {};
reag['class'] = ['amine','aldehyde','carboxylic acid','isocyanide'];
reag['id'] = [uid_ami,uid_ald,uid_car,uid_iso];
reag['name'] = [nam_ami,nam_ald,nam_car,nam_iso];        
reag['count'] = [num_ami,num_ald,num_car,num_iso];        
reag['position'] = [loc_ami,loc_ald,loc_car,loc_iso];
reag['mass'] = [mas_ami,mas_ald,mas_car,mas_iso];
reag['formula'] = [for_ami,for_ald,for_car,for_iso];
reag['structure'] = [str_ami,str_ald,str_car,str_iso];

# Determine the wells associated with each compound (to handle multiple occurences of a compounds)
#def collect_associated_wells(pos,num,uid,cid):
#    counter = [];    
#    loc = [];
#    for n in range(0,num):
#        ccid = uid[n];
#        indices = [ind for ind in range(len(cid)) if cid[ind]==ccid];
#        counter.append(len(indices));
#        loc.append([pos[ind] for ind in indices]);
#    return loc,counter;
#loc_ami,well_count_ami = collect_associated_wells(pos_ami,num_ami,uid_ami,cid_ami)
#loc_ald,well_count_ald = collect_associated_wells(pos_ald,num_ald,uid_ald,cid_ald)
#loc_car,well_count_car = collect_associated_wells(pos_car,num_car,uid_car,cid_car)
#loc_iso,well_count_iso = collect_associated_wells(pos_iso,num_iso,uid_iso,cid_iso)
#reag['position'] = [loc_ami,loc_ald,loc_car,loc_iso];


#%%############################################################################    
print('\n\n########## MAKE DESTINATION/LIBRARY PLATE ######################\n')
###############################################################################    

# Initialize destination plate
plateDestin = chemcpupy.WellPlate1536LDV(description='destination plate')
pos_empty = plateDestin.list_empty_positions();
if len(pos_empty)<num_ugi:
    raise ValueError('Not enough positions in well plate to make this library.')

# Set aside positions for alignment (done here so that the positions on the well plate match those on the final output, the MALDI plate)
if add_alignment_spots:    
    plate_rows = plateDestin._rows;
    plate_cols = plateDestin._cols;
#    alignment_positions = [(0,0), (plate_rows-1,0), (0,plate_cols-1), (plate_rows-1,plate_cols-1)]; # at the four corners of the plate 
    alignment_positions = [(0, 0), (0, 47), (31, 23)]; # three point alignment : ['X01Y01', 'X48Y01', 'X24Y32']     
    for p_ali in alignment_positions:
        pos_empty.remove(p_ali);
            
# Generate the combinations/positions for the 4-component reagent combinations
# (make all possible combinations of {1xAmine, 1xAldehyde, 1xCarboxylicAcid, 5xIsocyanide} and do so by adding isocyanides in one order and then the reverse)        
num_copies = 1; # number of copies
pos_use_ami = []; # source
pos_use_ald = []; # source
pos_use_car = []; # source
pos_use_iso = []; # source
pos_use_iso2 = []; # source
pos_use_iso3 = []; # source
pos_use_iso4 = []; # source
pos_used = []; # destination
ind_ami = reag['class'].index('amine');
ind_ald = reag['class'].index('aldehyde');    
ind_car = reag['class'].index('carboxylic acid');
ind_iso = reag['class'].index('isocyanide');    
counter_ami = 0;
counter_ald = 0;
counter_car = 0;
counter_iso = 0;
counter_iso2 = 0;
counter_iso3 = 0;
counter_iso4 = 0;
count_ami = min([len(x) for x in loc_ami]);   
count_ald = min([len(x) for x in loc_ald]);                 
count_car = min([len(x) for x in loc_car]);                 
count_iso = min([len(x) for x in loc_iso]);  
count_iso2 = min([len(x) for x in loc_iso]);  
count_iso3 = min([len(x) for x in loc_iso]);  
count_iso4 = min([len(x) for x in loc_iso]);  
#                      
# repitition 1
for n_rep in range(0,num_copies): 
    for n_ami in range(0,num_ami):
        for n_ald in range(0,num_ald):
            for n_car in range(0,num_car):
                for n_iso in range(0,num_iso):
                    for n_iso2 in range(0,num_iso):
                        for n_iso3 in range(0,num_iso):
                            for n_iso4 in range(0,num_iso):                            
                                pos_use_ami.append(reag['position'][ind_ami][n_ami][counter_ami]); 
                                counter_ami=counter_ami+1; 
                                if counter_ami==count_ami: 
                                    counter_ami=0;
                                pos_use_ald.append(reag['position'][ind_ald][n_ald][counter_ald]); 
                                counter_ald=counter_ald+1; 
                                if counter_ald==count_ald: 
                                    counter_ald=0;
                                pos_use_car.append(reag['position'][ind_car][n_car][counter_car]); 
                                counter_car=counter_car+1; 
                                if counter_car==count_car: 
                                    counter_car=0;
                                pos_use_iso.append(reag['position'][ind_iso][n_iso][counter_iso]); 
                                counter_iso=counter_iso+1; 
                                if counter_iso==count_iso: 
                                    counter_iso=0;
                                pos_use_iso2.append(reag['position'][ind_iso][n_iso2][counter_iso2]); 
                                counter_iso2=counter_iso2+1; 
                                if counter_iso2==count_iso2: 
                                    counter_iso2=0;
                                pos_use_iso3.append(reag['position'][ind_iso][n_iso3][counter_iso3]); 
                                counter_iso3=counter_iso3+1; 
                                if counter_iso3==count_iso3: 
                                    counter_iso3=0;
                                pos_use_iso4.append(reag['position'][ind_iso][n_iso4][counter_iso4]); 
                                counter_iso4=counter_iso4+1; 
                                if counter_iso4==count_iso4: 
                                    counter_iso4=0;                                
                                pos_used.append(pos_empty.pop(0));    
# repitition 2
for n_rep in range(0,num_copies): 
    for n_iso in range(0,num_iso):
        for n_iso2 in range(0,num_iso):
            for n_iso3 in range(0,num_iso):
                for n_iso4 in range(0,num_iso):                                            
                    for n_ami in range(0,num_ami):
                        for n_ald in range(0,num_ald):
                            for n_car in range(0,num_car):
                                pos_use_ami.append(reag['position'][ind_ami][n_ami][counter_ami]); 
                                counter_ami=counter_ami+1; 
                                if counter_ami==count_ami: 
                                    counter_ami=0;
                                pos_use_ald.append(reag['position'][ind_ald][n_ald][counter_ald]); 
                                counter_ald=counter_ald+1; 
                                if counter_ald==count_ald: 
                                    counter_ald=0;
                                pos_use_car.append(reag['position'][ind_car][n_car][counter_car]); 
                                counter_car=counter_car+1; 
                                if counter_car==count_car: 
                                    counter_car=0;
                                pos_use_iso.append(reag['position'][ind_iso][n_iso][counter_iso]); 
                                counter_iso=counter_iso+1; 
                                if counter_iso==count_iso: 
                                    counter_iso=0;
                                pos_use_iso2.append(reag['position'][ind_iso][n_iso2][counter_iso2]); 
                                counter_iso2=counter_iso2+1; 
                                if counter_iso2==count_iso2: 
                                    counter_iso2=0;
                                pos_use_iso3.append(reag['position'][ind_iso][n_iso3][counter_iso3]); 
                                counter_iso3=counter_iso3+1; 
                                if counter_iso3==count_iso3: 
                                    counter_iso3=0;
                                pos_use_iso4.append(reag['position'][ind_iso][n_iso4][counter_iso4]); 
                                counter_iso4=counter_iso4+1; 
                                if counter_iso4==count_iso4: 
                                    counter_iso4=0;                                
                                pos_used.append(pos_empty.pop(0));    

# Transfer reagents from source to destination
sort_transfers = False;                            
reagent_count = len(reag['class'])+3; # updated on 2020/02/26 # 2;
reagent_concentration = [0.5]*reagent_count; # [M]
reagent_transfer_volume = [200,200,200,200,200,200,200]; #[200,200,200,160]; #[100,100,100,80]; #[100]*reagent_count; # [nL]                      
reagent_pos = [ pos_use_ami, pos_use_ald, pos_use_car,  pos_use_iso,  pos_use_iso2,  pos_use_iso3,  pos_use_iso4 ];
reagent_class = [ 'amine',  'aldehyde', 'carboxylic acid', 'isocyanide', 'isocyanide', 'isocyanide' , 'isocyanide'];  
reagent_type = [ 'amine 1',  'aldehyde 1', 'carboxylic acid 1', 'isocyanide 1', 'isocyanide 2', 'isocyanide 3', 'isocyanide 4' ];  
reagent_symbol = [ 'Ami1',  'Ald1', 'Car', 'Iso1', 'Iso2', 'Iso3', 'Iso4' ];  
#reagent_num = [num_ami, num_ald, num_car, num_iso, num_car, num_iso];
#reagent_equivalence = numpy.multiply(reagent_concentration,reagent_transfer_volume);       
#reagent_equivalence = reagent_equivalence/max(reagent_equivalence);               
tasklist_combo = chemcpupy.TaskList(description='Create 4 Component Combinations');    
trans = {};
trans['source']=[];    
trans['destination']=[];    
trans['volume']=[];    
#trans['label']=[];  
def sort_related_lists(list_to_sort,list_to_follow):
    if len(list_to_sort) != len(list_to_follow):
        raise ValueError('Related lists are not of equal length.')
    list_sorted, list_followed = zip(*[(sortee, follower) for sortee, follower in sorted(zip(list_to_sort, list_to_follow))]) # sort by source
    list_sorted = list(list_sorted);
    list_followed = list(list_followed);        
    return list_sorted,list_followed;
for n in range(0,reagent_count):
    pos_source = reagent_pos[n];
    pos_destin = pos_used;
    if sort_transfers:
        pos_source_sorted, pos_destin_sorted_by_source = sort_related_lists(pos_source,pos_destin);
        pos_source_for_transfer = pos_source_sorted;
        pos_destin_for_transfer = pos_destin_sorted_by_source;
    else:
        pos_source_for_transfer = pos_source;
        pos_destin_for_transfer = pos_destin;        
    task = chemcpupy.TransferTask(from_plate=plateSource,
                            from_positions=pos_source_for_transfer,
                            to_plate=plateDestin,
                            to_positions=pos_destin_for_transfer,
                            transfer_volumes=[reagent_transfer_volume[n]],
                            transfer_group_label=reagent_type[n]);
    tasklist_combo.add(task)
    trans['source'] = trans['source'] + pos_source_for_transfer;
    trans['destination'] = trans['destination'] + pos_destin_for_transfer;
    trans['volume'] = trans['volume'] + [reagent_transfer_volume[n]]*len(pos_destin_for_transfer);
    #trans['label'] = trans['label'] + [reagent_type[n]];
#product_volume = sum(reagent_transfer_volume);    
used_positions = sorted(list(set(pos_used)));
    
# Simulate the transfers    
print('\n'+'Simulating reagent transfers')   
tasklist_combo.summarize();             
tasklist_combo.run(verbose=False,
       volume_increment=2.5,
       enforce_volume_limits=False,
       robot='echo')

# Generate Ugi Compounds
incubation_time = 24; # [hr]
print('\n'+'Simulating reactions')   
plateDestinInit = copy.deepcopy(plateDestin); # copy of destination plate before reaction
chemcpupy.UgiLibrary.generate_ugi_products(plateDestin,reaction_yield=1.0)  
plateDestinRxn = copy.deepcopy(plateDestin); # copy of destination plate after reaction
print('\n'+'Incubation time: % 0.3f'% (incubation_time) + ' hr.')
#print('\n'+'Product concentration %0.2f mM' % (product_concentration*1e3))

# Generate Other Possible Compounds
plateDestinRxnPasserini = copy.deepcopy(plateDestinInit); 
chemcpupy.PasseriniLibrary.generate_passerini_products(plateDestinRxnPasserini,reaction_yield=1.0) # simulate Passerini reactions

# Visualize the source plate after transfers
if visualize_plate:
    print('\n'+'Plotting source plate after transfers')
    print('\n'+' - '+'volumes [uL]:')            
    plateSource.graphical_print_volumes(units='uL')     

# Visualize the destination plate 
if visualize_plate:        
    print('\n'+'Plotting destination plate before Ugi reactions')
    print('\n'+' - '+'amines           - monoisotopic masses [g/mol]:')            
    plateDestinInit.graphical_print_by_id_with_type_filter('amine')
    print('\n'+' - '+'aldehydes        - monoisotopic masses [g/mol]:')            
    plateDestinInit.graphical_print_by_id_with_type_filter('aldehyde')
    print('\n'+' - '+'carboxylic acids - monoisotopic masses [g/mol]:')            
    plateDestinInit.graphical_print_by_id_with_type_filter('carboxylic acid')
    print('\n'+' - '+'isocyanides      - monoisotopic masses [g/mol]:')            
    plateDestinInit.graphical_print_by_id_with_type_filter('isocyanide')
        
    print('\n'+'Plotting destination plate after Ugi reactions')
    print('\n'+' - '+'compounds [#]:')            
    plateDestin.graphical_print_compounds_per_well()
    print('\n'+' - '+'types:')                
    plateDestin.graphical_print_by_type()    
    print('\n'+' - '+'monoisotopic masses [g/mol]:')        
    plateDestin.graphical_print_by_mass()
    print('\n'+' - '+'volumes [uL]:')            
    plateDestin.graphical_print_volumes(units='uL')

#%%############################################################################    
print('\n\n########## DILUTE DESTINATION/LIBRARY PLATE ####################\n')
###############################################################################    

# Dilute products after incubation period
solvent_transfer_volume = 3200; #1600; # [nL]
num_dest_per_source = math.ceil(len(used_positions)/num_sol);
tasklist_dilute = chemcpupy.TaskList(description='Dilute Products');  
counter = 0;
if num_dest_per_source*solvent_transfer_volume > source_working_volume:
    raise ValueError('Volume taken from solvent well exceeds the working range.')
for n in range(0,len(used_positions),num_dest_per_source):
    task = chemcpupy.TransferTask(from_plate=plateSource,
                            from_positions=pos_sol[counter],
                            to_plate=plateDestin,
                            to_positions=used_positions[n:(n+num_dest_per_source)],
                            transfer_volumes=[solvent_transfer_volume],
                            transfer_group_label='dilution');
    tasklist_dilute.add(task);
    counter = counter + 1; 
#diluted_product_volume = product_volume+solvent_transfer_volume;
#diluted_product_concentration = product_concentration*product_volume/diluted_product_volume;    
     
# Simulate the transfers    
print('\n'+'Simulating dilutions')   
tasklist_dilute.summarize();             
tasklist_dilute.run(verbose=False,
       volume_increment=2.5,
       enforce_volume_limits=False,
       robot='echo')   
#print('\n'+'Product concentration after dilution: %0.2f mM' % (diluted_product_concentration*1e3))

# Visualize the destination plate 
if visualize_plate:
    print('\n'+'Plotting destination plate after dilution')
    print('\n'+' - '+'volumes [uL]:')            
    plateDestin.graphical_print_volumes(units='uL')
    
# plateDestinInit.get_compounds([0,7],['type','name','cid','mass'])
    
#%%############################################################################    
print('\n\n########## SPOT LIBRARY TO MALDI PLATE #########################\n')
###############################################################################    

# Prepare to spot products to MALDI plate 
library_spotting_volume = 20; # [nL]
matrix_spotting_volume = library_spotting_volume; # [nL]
matrix_concentration = 176.21; # [M] (176.21mM, 52.86mM, 26.43mM, or 15.87mM for 10mg HCCA in 0.3mL, 1mL, 2mL, or 3.33mL DMSO) (10e-3/189.17/Xe-3*1e3)    
aligner_spotting_volume = matrix_spotting_volume + library_spotting_volume; # [nL]
plateLibrary = copy.deepcopy(plateDestin); # copy of destination plate 
plateMALDI = chemcpupy.MaldiPlate1536(description='MALDI plate');
tasklist_spot_matrix = chemcpupy.TaskList(description='MALDI Spot Matrix');    
tasklist_spot_sample = chemcpupy.TaskList(description='MALDI Spot Products');
#spotted_product_volume = matrix_spotting_volume + library_spotting_volume;
#spotted_product_concentration = diluted_product_concentration*library_spotting_volume/spotted_product_volume;    
#spotted_matrix_concentration = matrix_concentration*matrix_spotting_volume/spotted_product_volume;    

# Spot aligner to MALDI plate 
if add_alignment_spots:    
    for pos in alignment_positions:
        task = chemcpupy.TransferTask(from_plate=plateSource,
                                from_positions=pos_ali[0],
                                to_plate=plateMALDI,
                                to_positions=pos,
                                transfer_volumes=[aligner_spotting_volume],
                                transfer_group_label='MALDI');
        tasklist_spot_matrix.add(task);
            
# Spot matrix to MALDI plate  
counter = 0;
well_index= 0;
draws_per_matrix_well = math.ceil(len(used_positions)/len(pos_mat));
for pos in used_positions:
    task = chemcpupy.TransferTask(from_plate=plateSource,
                            from_positions=pos_mat[well_index],
                            to_plate=plateMALDI,
                            to_positions=pos,
                            transfer_volumes=[matrix_spotting_volume], 
                            transfer_group_label='matrix');
    tasklist_spot_matrix.add(task)      
    if counter < draws_per_matrix_well:
        counter = counter + 1;    
    else:
        counter = 0;      
        well_index = well_index + 1;
        
# Spot library to MALDI plate 
for pos in used_positions:
    task = chemcpupy.TransferTask(from_plate=plateLibrary,
                            from_positions=pos,
                            to_plate=plateMALDI,
                            to_positions=pos,
                            transfer_volumes=[library_spotting_volume],
                            transfer_group_label='MALDI');
    tasklist_spot_sample.add(task);             
         
# Simulate the transfers    
print('\n'+'Simulating MALDI spotting')   
tasklist_spot = chemcpupy.TaskList(description='MALDI Spot');    
tasklist_spot.merge(tasklist_spot_matrix);
tasklist_spot.merge(tasklist_spot_sample);    
tasklist_spot.summarize();             
tasklist_spot.run(verbose=False,
       volume_increment=2.5,
       enforce_volume_limits=False,
       robot='echo')

#print('\n'+'Matrix concentration after spotting: %0.2f mM' % (spotted_product_concentration*1e3))    
#print('\n'+'Product concentration after spotting: %0.2f mM' % (spotted_product_concentration*1e3))  

# Visualize the source plate after transfers
if visualize_plate:
    print('\n'+'Plotting source plate after spotting')
    print('\n'+' - '+'volumes [uL]:')            
    plateSource.graphical_print_volumes(units='uL')     

# Visualize the library plate after transfers
if visualize_plate:
    print('\n'+'Plotting library plate after spotting')
    print('\n'+' - '+'volumes [uL]:')            
    plateLibrary.graphical_print_volumes(units='uL')     

# Visualize the MALDI plate 
if visualize_plate:
    print('\n'+'Plotting MALDI plate after spotting')
    print('\n'+' - '+'volumes [uL]:')            
    plateMALDI.graphical_print_volumes(units='nL')
    
#%%############################################################################    
print('\n\n########## EXPORT DATA #########################################\n')
###############################################################################    
    
# Export all transfer tasks to a picklist (for Echo)  
if generate_recipe_for_Echo:    
    print('\n'+'Exporting Transfers')  
    chemcpupy.Echo.write_Echo_csv_picklist(tasklist_combo, temp_dir+'/'+source_plate_sheet+'_step_1_combo'+'_Echo_picklist'+'.csv')     
    chemcpupy.Echo.write_Echo_csv_picklist(tasklist_dilute, temp_dir+'/'+source_plate_sheet+'_step_2_dilute'+'_Echo_picklist'+'.csv') 
    chemcpupy.Echo.write_Echo_csv_picklist(tasklist_spot_matrix, temp_dir+'/'+source_plate_sheet+'_step_3a_spot_matrix'+'_Echo_picklist'+'.csv') 
    chemcpupy.Echo.write_Echo_csv_picklist(tasklist_spot_sample, temp_dir+'/'+source_plate_sheet+'_step_3b_spot_sample'+'_Echo_picklist'+'.csv') 

# Export used plate positions to picklist (for mass spec) 
if generate_picklist_for_Solarix:  
    skip_positions = [];
    if add_alignment_spots:    
        skip_positions = alignment_positions;
    print('\n'+'Exporting Positions')    
    chemcpupy.Bruker.write_Bruker_xlsx_picklist(plateMALDI,temp_dir+'/'+source_plate_sheet+'_Bruker_picklist.xlsx', exclude_positions = skip_positions)       

# Export compound info to CSV
if export_masses:
    print('\n'+'Exporting Compound Info')     
    plateSource.export_csv_review(temp_dir + '/'+source_plate_sheet+'_reagent_source'+'_review.csv');
    plateDestinInit.export_csv_review(temp_dir + '/'+source_plate_sheet+'_reagent'+'_review.csv');
    plateDestinRxn.export_csv_review(temp_dir + '/'+source_plate_sheet+'_review.csv');    
    plateDestinRxn.export_csv_review(temp_dir + '/'+source_plate_sheet+'_ugi_review.csv',typefilter='ugi');
    plateDestinRxnPasserini.export_csv_review(temp_dir + '/'+source_plate_sheet+'_passerini'+'_review.csv',typefilter='passerini');

# Convert the sourceplatesheet to platesheet
if export_as_platesheet:
    not_yet_implemented = True;
    
# Save details to MAT
if save_variables:
    print('\n'+'Exporting Variables')         
    export_variables = dict();    
    time = datetime.now();    
    about = {'name':source_plate_sheet, 'filename':filename, 'timestamp':time.strftime("%Y/%m/%d %H:%M:%S")}
    used_options = {'num_copies':num_copies, 'sort_transfers':sort_transfers, 'has_alignment_spots':add_alignment_spots};
    reagent_slots = {'count':reagent_count, 'class':reagent_class, 'symbol':reagent_symbol, 'name':reagent_type, 'concentration':reagent_concentration, 'transfer_volume':reagent_transfer_volume, 'position_map':reagent_pos, 'positions':pos_used};
    #plate_contents = {'source':plateSource.get_contents(), 'destination':plateDestinInit.get_contents()}; #, 'destinationPreReact':plateDestinInit.get_contents(),'destinationPostReact':plateDestinRxn.get_contents(),'MALDI':plateMALDI.get_contents()};    
    review_plates = {'source':plateSource.export_array_review(), 'destinationPreReaction':plateDestinInit.export_array_review(), 'destinationPostReaction':plateDestinRxn.export_array_review(), 'destinationUgiProducts':plateDestinRxn.export_array_review(typefilter='ugi'), 'MALDI':plateMALDI.export_array_review()};    
    export_variables.update({'reagents':reag, 'transfers':trans, 'about':about, 'options':used_options, 'slots': reagent_slots, 'plates':review_plates});    
#    export_variables.update({'reag':reag, 'trans':trans, 'num_copies':num_copies, 'pos_used':pos_used,
#                             'sort_transfers':sort_transfers,
#                             'reagent_count':reagent_count,'reagent_type':reagent_type,'reagent_concentration':reagent_concentration,
#                             'reagent_transfer_volume':reagent_transfer_volume,'reagent_pos':reagent_pos,
#                             })
    scipy.io.savemat(temp_dir + '/'+source_plate_sheet+'_variables.mat', mdict=export_variables);
    print('\n'+'Wrote Variables: ' + temp_dir + '/'+source_plate_sheet+'_variables.mat') 

#%%##############################################################################    
print('\n\n################################################################\n')
###############################################################################    

