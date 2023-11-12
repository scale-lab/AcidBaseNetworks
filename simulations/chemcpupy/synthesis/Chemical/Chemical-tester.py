#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chrisarcadia 
@created: 2018/06/04
"""

# libraries
import Chemical

###############################################################################
# Periodic Table Of Elements
 
## ways of instantiating the periodic table of elements
elements = Chemical.PeriodicTableOfElements();

## reasons to use this class

### to look up element information 
print(elements.name('H')) # get name of the element
print(elements.mass('H')) # get mass of the element
print(elements.protons('H')) # get number of protons (atomic number) in the element
print(elements.neutrons('H')) # get number of neutrons in the element
print(elements.nucleons('H')) # get number of nucleons (mass number) in the element
print()
# given values are for the most abundant isotope of the element
# since this table is intended to be used for mass spec analyses

### to look up compound information from simplified molecular formula (no parentheses allowed, nor charges)
print(elements.formula_mass('H2O')) # get principle monoisotopic mass of the compound
print(elements.formula_weight('H2O')) # get weighted average mass (molecular weight) of the compound
print(elements.formula_mass_number('H2O')) # get mass number of the compound
print()

### to look up particle masses [Da]
print(elements.mass('e')) # electron
print(elements.mass('p')) # proton
print(elements.mass('n')) # neutron
print()

### to get elemental composition formula and mass of a chemical reaction output
stoichiometric_list2 = [{'moles':1,'formula':'CH4'},
                       {'moles':2,'formula':'O2'},
                       {'moles':-2,'formula':'H2O'},]; 
                       # combustion of methane
                       # CH4(g) + 2O2(g) -> CO2(g) + 2H2O(g)   
stoichiometric_list = [{'moles':3,'formula':'H2'},
                       {'moles':1,'formula':'N2'}]; 
                       # sythesis of ammonia
                       # 3H2(g) + N2(g) -> 2NH3(g)    
print(elements.formula_sum(stoichiometric_list2)) # carbon dioxide                       
product_formula = elements.formula_sum(stoichiometric_list); # ammonia
print(product_formula + ' (' + str(elements.formula_mass(product_formula)) +' Da)')  # ammonia
# if you know how many moles there will be of a product, you can get its the molecular formula:
moles_expected = 2; # ammonia
product_molecular_formula = elements.formula_sum(stoichiometric_list,moles_expected); # ammonia
print(product_molecular_formula + ' (' + str(elements.formula_mass(product_molecular_formula)) +' Da)')  # ammonia



###############################################################################
# Compound 

## ways of instantiating a compound

### providing a correctly formatted dictionary with all required properties (will also save additional properties as Compound.info)
properties = {'Type':'salt', 
              'Name':'sodium chloride', 
              'ID [PubChem]':5234, # if the ID is unknown then set this field as: -1 (or any other negative number)
              'Mass [Da]':57.959, # if the mass is unknown then set this field as: -1
              'Molecular Formula':'NaCl'}; # if the molecular formula is unknown then set this field as: ''
              # for all the required fields to be populated, at a minimum either compound ID, formula, or mass must be specified
a = Chemical.Compound(properties);
a.fetch_PubChem_info(); # not required, compound has been fully specfiied

### providing a PubChem CID 
c = Chemical.Compound(5234);
c.fetch_PubChem_info(); # required 

### directly specifying all properties (no need for PubChem)
b = Chemical.Compound();
b.id = 5234;
b.mass = 57.959;
b.name = 'sodium chloride';
b.type = 'salt';
b.formula = 'NaCl';
b.fetch_PubChem_info(); # not required, compound has been fully specfiied

### if you have a chemical whose is not in the PubChem database then set the ID to a nonpositive integer (0,-1, etc.) so that a unique ID (for in house use) can be generated

## View all compound data as dictionary
a.to_dictionary() # as dict
str(a) # as a str from dict

## PubChem data

### looking at PubChem data
a.info['Synonyms'];

### if you want to retreive more properties from PubChem than are supported by the Compound class simply store the PubChem data
pubchemdata = a.fetch_PubChem_info();

###############################################################################
# Compound Data Base

## ways of instantiating a compound data base

solvent = {'Type':'solvent', 'Name':'dimethyl sulfoxide', 'ID [PubChem]':679, 'Mass [Da]':-1, 'Molecular Formula':''};
analyte1 = {'Type':'phenol', 'Name':'2,4,6 tri-tert-butyl phenol', 'ID [PubChem]':12902, 'Mass [Da]':-1, 'Molecular Formula':''};
analyte2 = {'Type':'phenol', 'Name':'2,6 dimethyl phenol', 'ID [PubChem]':11335, 'Mass [Da]':-1, 'Molecular Formula':''};
analyte3 = {'Type':'phenol', 'Name':'4-nitro phenol', 'ID [PubChem]':980, 'Mass [Da]':-1, 'Molecular Formula':''};
dictionary_list = [solvent, analyte1, analyte2, analyte3];
compound_list = [Chemical.Compound(solvent), Chemical.Compound(analyte1), Chemical.Compound(analyte2), Chemical.Compound(analyte3)];
id_list = [679,12902,11335,980];

### from list of Compounds
adb = Chemical.CompoundDataBase('DataBase A',compound_list);
adb.fetch_PubChem_info();

### from list of formatted dictionaries
bdb = Chemical.CompoundDataBase('DataBase B',dictionary_list);
bdb.fetch_PubChem_info();

### from list of PubChem IDs
cdb = Chemical.CompoundDataBase('DataBase C',id_list);
cdb.fetch_PubChem_info();

## merge databases
db1 = Chemical.CompoundDataBase('DataBase 1',[compound_list[1]]);
db2 = Chemical.CompoundDataBase('DataBase 2',[compound_list[2]]);
db3 = Chemical.CompoundDataBase('DataBase 3',[compound_list[3]]);
mdb = Chemical.CompoundDataBase('Composite DataBase');
mdb.merge([db1,db2,db3]);
print(mdb.get_names());

## view all values of a particular field
print(adb.get_names())
print(adb.get_types())
print(adb.get_masses())
print(adb.get_ids())
print(adb.get_values('ID [PubChem]'))
print(adb.get_values('Molecular Formula'))

## lookup compounds by their IDs in the database
compound_found = adb.lookup(11335);
print(compound_found.name)
compounds_found = adb.lookup_multiple([11335,1,980]);
print(compounds_found[1])
print(compounds_found[2].name)

## save database to file
adb.save('example_cdb.json')

## load database from file
a2db = Chemical.CompoundDataBase();
a2db.load('example_cdb.json')
#a2db.fetch_PubChem_info(); # unless you manually fetch from PubChem, loaded data will not be updated from PubChem


###############################################################################
# Solution

## ways of instantiating a chemical solution

# Convention: arbitrary concentrations (such as that for a solvent) should be defined as negative

### adding to empty solution
compound_ids = [679,12902,11335,980];
compound_concentrations = [-1,1,1,1]; # denote solvent by -1
volumes_to_add = [10e-6,1e-6,1e-6,1e-6];
compounds = {'Volume [L]':volumes_to_add, 'ID [PubChem]':compound_ids, 'Concentration [M]':compound_concentrations};
sol = Chemical.Solution();
sol.add(compounds)   
   
### with initial contents of the solution 
initial = {'Volume [L]':10e-6, 'ID [PubChem]':[679], 'Concentration [M]':[-1]}
soli = Chemical.Solution(initial);
# Note: When you provide solution intial properties the volume field is a single float and not a list, as this value is the initial total volume of the solution and the concentrations are given for the solution,
#       Contrast this to when you add solutions to a solution, where the volume field is a list of the desired transfer volume for each compound whose source solution concentration is specified

### adding solutions to empty solution
initial1 = {'Volume [L]':200e-6, 'ID [PubChem]':[679], 'Concentration [M]':[-1]}
initial2 = {'Volume [L]':200e-6, 'ID [PubChem]':[887], 'Concentration [M]':[-1]}
sole = Chemical.Solution();
sol1 = Chemical.Solution(initial1);
sol2 = Chemical.Solution(initial2);
solutions = {'Volume [L]':[1e-6,99e-6], 'Solution':[sol1,sol2]};
sole.add(solutions);
# show final volumes and concentrations 
print(sol1.compounds) # source 1
print(sol1.concentrations) 
print(sol1.volume*1e6) 
print()
print(sol2.compounds) # source 2
print(sol2.concentrations) 
print(sol2.volume*1e6) 
print()
print(sole.compounds) # destination
print(sole.concentrations) 
print(sole.volume*1e6) 
print()
        

###############################################################################
# SimpleReaction

## performing a one-pot sythesis reaction
print('------ before reaction ------')
print()    

### make chemical database
amine = {'Type':'amine', 'Name':'paramethoxybenzlamine', 'ID [PubChem]':75452, 'Mass [Da]':-1, 'Molecular Formula':''};
aldehyde = {'Type':'aldehyde', 'Name':'2-nitrobenzaldehyde', 'ID [PubChem]':11101, 'Mass [Da]':-1, 'Molecular Formula':''};
carboxylicacid = {'Type':'carboxylic acid', 'Name':'Boc-O-benzyl L-beta-homotyrosine', 'ID [PubChem]':2761555, 'Mass [Da]':-1, 'Molecular Formula':''};
isocyanide = {'Type':'isocyanide', 'Name':'methyl isocyanoacetate', 'ID [PubChem]':547815, 'Mass [Da]':-1, 'Molecular Formula':''};
dictionary_list = [amine, aldehyde, carboxylicacid, isocyanide];
database = Chemical.CompoundDataBase('DataBase Demo',dictionary_list);
database.fetch_PubChem_info();
print(database)
print()

### make solution
compound_ids = [75452,11101,2761555,547815];
compound_concentrations = [4,4,4,4]; 
volumes_to_add = [10e-6,10e-6,10e-6,10e-6];
compounds = {'Volume [L]':volumes_to_add, 'ID [PubChem]':compound_ids, 'Concentration [M]':compound_concentrations};
solution = Chemical.Solution();
solution.add(compounds)  
print(solution)  
print()
database_subset_in_solution = solution.get_compounds(database);
database_subset_in_solution.get_names()
print(database_subset_in_solution.get_values('Concentration [M]')) # only databases that are processed from the Solution class have this field
print()

### run simple sythesis reaction
rxn = Chemical.SimpleReaction();
reaction = 'Ugi';
#efficiency = 0.90; # 90% conversion to product
efficiency = 1.00; # 100% conversion to product
product_ID = rxn.react(reaction,efficiency,solution,database); # success status is also the product ID
if not product_ID:
    print('------ reaction failed ------');    
else:    
    print('------ reaction completed ------');   
    print()
    print(rxn.get_equation(reaction));
    print()    
    print('------ after reaction ------')
    print()
    print(database)
    print()
    print(solution)
    print()
    database_subset_in_solution = solution.get_compounds(database);
    print(database_subset_in_solution.get_names())
    ### remove compounds with zero concentration from the solution
    solution.clean_up();
    print()
    print('------ after cleanup ------')
    print()
    print(database)
    print()
    print(solution)
    print()
    database_subset_in_solution = solution.get_compounds(database);
    print(database_subset_in_solution.get_names())
    print()
    print('------ product info ------')
    print()    
    product_concentration = solution.get_concentration(product_ID); 
    product_compound = database_subset_in_solution.lookup(product_ID);
    product_name = product_compound.name;
    product_precursors = product_compound.info['Precursor'];
    print('Name: '+product_name)        
    print('Concentration [M]: '+str(product_concentration))    
    print('ID:' + str(product_ID))    
    print('Precursor IDs [PubChem]: ' + str(product_precursors))  
    
        

###############################################################################
# Container
# #not implemented yet
    
## ways of instantiating a chemical container
    
### start with a natively supported container type
container = Chemical.Container('WellPlate384'); # choose a type and provide it as input
print(container.get_native_types()) # print list of supported types


### provide container information
initial = {'Volume [L]':10e-6, 'ID [PubChem]':[679], 'Concentration [M]':[-1]}
solution = Chemical.Solution(initial);
solutions = [solution];
positions = [(1,1)];
properties = {'Type':'WellPlate384',
              'Solutions':solutions,
              'Positions':positions};
container1 = Chemical.Container(properties); 

### provide container information with custom container type: (volumes must be in liters)
customtype = {'Type':'Tube25mLCustom', 
              'Rows':1, 
              'Columns':1, 
              'Dead Volume [L]':500e-6, 
              'Max Volume [L]':25e-3, 
              'Category':'tube'};                                  
properties = {'Type':customtype,
              'Solutions':solutions,
              'Positions':positions};
container2 = Chemical.Container(properties); 

### editing an empty container
con1 = Chemical.Container();


###### this class is incomplete (will return to at later point)    


###############################################################################
# Work Bench

# a workbench object will contain a list of container and a compound database
# it should be a complete specification of the materials for an experiment 
# and should be a save/load-able file

# import a plate sheet as a workbench

###############################################################################
# todo's

# - finish Container class
# - write WorkBench class (including import from import a PlateSheet file)
# - convert/export transfer tasks in the virtual workspace as physical transfer tasks for Andrew or Echo
# - consider using pyteomics to calculate mass from chemical formula

