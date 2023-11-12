#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: Chemical classes
@description: classes for chemical objects, calculations, and simulations
@author: chrisarcadia 
@created: 2018/06/03
"""

# units
#   monoisotopic mass : Dalton [Da]
#   concentration : moles per liter [mol/L] = [M] 
#   volume : liters [L]
#   IDs : PubChem CID (or unique local)

import os, copy, datetime, random, json, re, pubchempy

###############################################################################

class PeriodicTableOfElements(object):
    
    # Periodic Table of Elements
    #   provides information about each element  
    
    def __init__(self, *args): 
                                        
        # class constants
        self.__table = {   
                        'H':{'protons':1, 'nucleons':1, 'name':'Hydrogen', 'mass':1.007825032, 'weight':1.00794075382579},
                        'He':{'protons':2, 'nucleons':4, 'name':'Helium', 'mass':4.002603254, 'weight':4.00260193199093},
                        'Li':{'protons':3, 'nucleons':7, 'name':'Lithium', 'mass':7.016003437, 'weight':6.940036603255},
                        'Be':{'protons':4, 'nucleons':9, 'name':'Beryllium', 'mass':9.012183065, 'weight':9.012183065},
                        'B':{'protons':5, 'nucleons':11, 'name':'Boron', 'mass':11.00930536, 'weight':10.81102804641},
                        'C':{'protons':6, 'nucleons':12, 'name':'Carbon', 'mass':12, 'weight':12.010735896788},
                        'N':{'protons':7, 'nucleons':14, 'name':'Nitrogen', 'mass':14.003074, 'weight':14.006703207036},
                        'O':{'protons':8, 'nucleons':16, 'name':'Oxygen', 'mass':15.99491462, 'weight':15.9994049247427},
                        'F':{'protons':9, 'nucleons':19, 'name':'Fluorine', 'mass':18.99840316, 'weight':18.99840316},
                        'Ne':{'protons':10, 'nucleons':20, 'name':'Neon', 'mass':19.99244018, 'weight':20.180046383602},
                        'Na':{'protons':11, 'nucleons':23, 'name':'Sodium', 'mass':22.98976928, 'weight':22.98976928},
                        'Mg':{'protons':12, 'nucleons':24, 'name':'Magnesium', 'mass':23.9850417, 'weight':24.305051622827},
                        'Al':{'protons':13, 'nucleons':27, 'name':'Aluminum', 'mass':26.98153853, 'weight':26.98153853},
                        'Si':{'protons':14, 'nucleons':28, 'name':'Silicon', 'mass':27.97692653, 'weight':28.0854987013117},
                        'P':{'protons':15, 'nucleons':31, 'name':'Phosphorus', 'mass':30.973762, 'weight':30.973762},
                        'S':{'protons':16, 'nucleons':32, 'name':'Sulfur', 'mass':31.97207117, 'weight':32.064787401779},
                        'Cl':{'protons':17, 'nucleons':35, 'name':'Chlorine', 'mass':34.96885268, 'weight':35.452937580608},
                        'Ar':{'protons':18, 'nucleons':40, 'name':'Argon', 'mass':39.96238312, 'weight':39.9477985599133},
                        'K':{'protons':19, 'nucleons':39, 'name':'Potassium', 'mass':38.96370649, 'weight':39.0983009135851},
                        'Ca':{'protons':20, 'nucleons':40, 'name':'Calcium', 'mass':39.96259086, 'weight':40.0780225081095},
                        'Sc':{'protons':21, 'nucleons':45, 'name':'Scandium', 'mass':44.95590828, 'weight':44.95590828},
                        'Ti':{'protons':22, 'nucleons':48, 'name':'Titanium', 'mass':47.94794198, 'weight':47.866744962722},
                        'V':{'protons':23, 'nucleons':51, 'name':'Vanadium', 'mass':50.94395704, 'weight':50.941465037425},
                        'Cr':{'protons':24, 'nucleons':52, 'name':'Chromium', 'mass':51.94050623, 'weight':51.9961317554337},
                        'Mn':{'protons':25, 'nucleons':55, 'name':'Manganese', 'mass':54.93804391, 'weight':54.93804391},
                        'Fe':{'protons':26, 'nucleons':56, 'name':'Iron', 'mass':55.93493633, 'weight':55.8451444338659},
                        'Co':{'protons':27, 'nucleons':59, 'name':'Cobalt', 'mass':58.93319429, 'weight':58.93319429},
                        'Ni':{'protons':28, 'nucleons':58, 'name':'Nickel', 'mass':57.93534241, 'weight':58.6933471099477},
                        'Cu':{'protons':29, 'nucleons':63, 'name':'Copper', 'mass':62.92959772, 'weight':63.54603994583},
                        'Zn':{'protons':30, 'nucleons':64, 'name':'Zinc', 'mass':63.92914201, 'weight':65.377782529525},
                        'Ga':{'protons':31, 'nucleons':69, 'name':'Gallium', 'mass':68.9255735, 'weight':69.7230660725936},
                        'Ge':{'protons':32, 'nucleons':74, 'name':'Germanium', 'mass':73.92117776, 'weight':72.627550166039},
                        'As':{'protons':33, 'nucleons':75, 'name':'Arsenic', 'mass':74.92159457, 'weight':74.92159457},
                        'Se':{'protons':34, 'nucleons':80, 'name':'Selenium', 'mass':79.9165218, 'weight':78.959388556298},
                        'Br':{'protons':35, 'nucleons':79, 'name':'Bromine', 'mass':78.9183376, 'weight':79.90352778051},
                        'Kr':{'protons':36, 'nucleons':84, 'name':'Krypton', 'mass':83.91149773, 'weight':83.7979999968875},
                        'Rb':{'protons':37, 'nucleons':85, 'name':'Rubidium', 'mass':84.91178974, 'weight':85.467663596857},
                        'Sr':{'protons':38, 'nucleons':88, 'name':'Strontium', 'mass':87.9056125, 'weight':87.61664446962},
                        'Y':{'protons':39, 'nucleons':89, 'name':'Yttrium', 'mass':88.9058403, 'weight':88.9058403},
                        'Zr':{'protons':40, 'nucleons':90, 'name':'Zirconium', 'mass':89.9046977, 'weight':91.22364159706},
                        'Nb':{'protons':41, 'nucleons':93, 'name':'Niobium', 'mass':92.906373, 'weight':92.906373},
                        'Mo':{'protons':42, 'nucleons':98, 'name':'Molybdenum', 'mass':97.90540482, 'weight':95.959788541188},
                        'Tc':{'protons':43, 'nucleons':97, 'name':'Technetium', 'mass':96.9063667, 'weight':97.9066098687601},
                        'Ru':{'protons':44, 'nucleons':102, 'name':'Ruthenium', 'mass':101.9043441, 'weight':101.06494013916},
                        'Rh':{'protons':45, 'nucleons':103, 'name':'Rhodium', 'mass':102.905498, 'weight':102.905498},
                        'Pd':{'protons':46, 'nucleons':106, 'name':'Palladium', 'mass':105.9034804, 'weight':106.41532750734},
                        'Ag':{'protons':47, 'nucleons':107, 'name':'Silver', 'mass':106.9050916, 'weight':107.868149634557},
                        'Cd':{'protons':48, 'nucleons':114, 'name':'Cadmium', 'mass':113.9033651, 'weight':112.41155783105},
                        'In':{'protons':49, 'nucleons':115, 'name':'Indium', 'mass':114.9038788, 'weight':114.8180866507},
                        'Sn':{'protons':50, 'nucleons':120, 'name':'Tin', 'mass':119.9022016, 'weight':118.71011259491},
                        'Sb':{'protons':51, 'nucleons':121, 'name':'Antimony', 'mass':120.903812, 'weight':121.75978367348},
                        'Te':{'protons':52, 'nucleons':130, 'name':'Tellurium', 'mass':129.9062227, 'weight':127.60312647465},
                        'I':{'protons':53, 'nucleons':127, 'name':'Iodine', 'mass':126.9044719, 'weight':126.9044719},
                        'Xe':{'protons':54, 'nucleons':132, 'name':'Xenon', 'mass':131.9041551, 'weight':131.292761474025},
                        'Cs':{'protons':55, 'nucleons':133, 'name':'Cesium', 'mass':132.905452, 'weight':132.905452},
                        'Ba':{'protons':56, 'nucleons':138, 'name':'Barium', 'mass':137.905247, 'weight':137.326891623585},
                        'La':{'protons':57, 'nucleons':139, 'name':'Lanthanum', 'mass':138.9063563, 'weight':138.905468873713},
                        'Ce':{'protons':58, 'nucleons':140, 'name':'Cerium', 'mass':139.9054431, 'weight':140.115730737836},
                        'Pr':{'protons':59, 'nucleons':141, 'name':'Praseodymium', 'mass':140.9076576, 'weight':140.9076576},
                        'Nd':{'protons':60, 'nucleons':142, 'name':'Neodymium', 'mass':141.907729, 'weight':144.241596031827},
                        'Pm':{'protons':61, 'nucleons':145, 'name':'Promethium', 'mass':144.9127559, 'weight':145.91395045},
                        'Sm':{'protons':62, 'nucleons':152, 'name':'Samarium', 'mass':151.9197397, 'weight':150.36635571193},
                        'Eu':{'protons':63, 'nucleons':153, 'name':'Europium', 'mass':152.921238, 'weight':151.96437812638},
                        'Gd':{'protons':64, 'nucleons':158, 'name':'Gadolinium', 'mass':157.9241123, 'weight':157.25213064688},
                        'Tb':{'protons':65, 'nucleons':159, 'name':'Terbium', 'mass':158.9253547, 'weight':158.9253547},
                        'Dy':{'protons':66, 'nucleons':164, 'name':'Dysprosium', 'mass':163.9291819, 'weight':162.499472819424},
                        'Ho':{'protons':67, 'nucleons':165, 'name':'Holmium', 'mass':164.9303288, 'weight':164.9303288},
                        'Er':{'protons':68, 'nucleons':166, 'name':'Erbium', 'mass':165.9302995, 'weight':167.259082649669},
                        'Tm':{'protons':69, 'nucleons':169, 'name':'Thulium', 'mass':168.9342179, 'weight':168.9342179},
                        'Yb':{'protons':70, 'nucleons':174, 'name':'Ytterbium', 'mass':173.9388664, 'weight':173.054150166317},
                        'Lu':{'protons':71, 'nucleons':175, 'name':'Lutetium', 'mass':174.9407752, 'weight':174.966814957855},
                        'Hf':{'protons':72, 'nucleons':180, 'name':'Hafnium', 'mass':179.946557, 'weight':178.4849787234},
                        'Ta':{'protons':73, 'nucleons':181, 'name':'Tantalum', 'mass':180.9479958, 'weight':180.947875636227},
                        'W':{'protons':74, 'nucleons':184, 'name':'Tungsten', 'mass':183.9509309, 'weight':183.84177754094},
                        'Re':{'protons':75, 'nucleons':187, 'name':'Rhenium', 'mass':186.9557501, 'weight':186.2067045456},
                        'Os':{'protons':76, 'nucleons':192, 'name':'Osmium', 'mass':191.961477, 'weight':190.22485962824},
                        'Ir':{'protons':77, 'nucleons':193, 'name':'Iridium', 'mass':192.9629216, 'weight':192.2160516521},
                        'Pt':{'protons':78, 'nucleons':195, 'name':'Platinum', 'mass':194.9647917, 'weight':195.084456867452},
                        'Au':{'protons':79, 'nucleons':197, 'name':'Gold', 'mass':196.9665688, 'weight':196.9665688},
                        'Hg':{'protons':80, 'nucleons':202, 'name':'Mercury', 'mass':201.9706434, 'weight':200.59916702622},
                        'Tl':{'protons':81, 'nucleons':205, 'name':'Thallium', 'mass':204.9744278, 'weight':204.38341283936},
                        'Pb':{'protons':82, 'nucleons':208, 'name':'Lead', 'mass':207.9766525, 'weight':207.216908063},
                        'Bi':{'protons':83, 'nucleons':209, 'name':'Bismuth', 'mass':208.9803991, 'weight':208.9803991},
                        'Po':{'protons':84, 'nucleons':209, 'name':'Polonium', 'mass':208.9824308, 'weight':209.48265245},
                        'At':{'protons':85, 'nucleons':210, 'name':'Astatine', 'mass':209.9871479, 'weight':210.48732225},
                        'Rn':{'protons':86, 'nucleons':211, 'name':'Radon', 'mass':210.9906011, 'weight':217.67319091566},
                        'Fr':{'protons':87, 'nucleons':223, 'name':'Francium', 'mass':223.019736, 'weight':223.019736},
                        'Ra':{'protons':88, 'nucleons':223, 'name':'Radium', 'mass':223.0185023, 'weight':225.273798825},
                        'Ac':{'protons':89, 'nucleons':227, 'name':'Actinium', 'mass':227.0277523, 'weight':227.0277523},
                        'Th':{'protons':90, 'nucleons':232, 'name':'Thorium', 'mass':232.0380558, 'weight':232.0380558},
                        'Pa':{'protons':91, 'nucleons':231, 'name':'Protactinium', 'mass':231.0358842, 'weight':231.0358842},
                        'U':{'protons':92, 'nucleons':238, 'name':'Uranium', 'mass':238.0507884, 'weight':238.028910461657},
                        'Np':{'protons':93, 'nucleons':236, 'name':'Neptunium', 'mass':236.04657, 'weight':236.5473718},
                        'Pu':{'protons':94, 'nucleons':238, 'name':'Plutonium', 'mass':238.0495601, 'weight':240.722556698112},
                        'Am':{'protons':95, 'nucleons':241, 'name':'Americium', 'mass':241.0568293, 'weight':242.0591053},
                        'Cm':{'protons':96, 'nucleons':243, 'name':'Curium', 'mass':243.0613893, 'weight':245.5665940578},
                        'Bk':{'protons':97, 'nucleons':247, 'name':'Berkelium', 'mass':247.0703073, 'weight':248.0726475},
                        'Cf':{'protons':98, 'nucleons':249, 'name':'Californium', 'mass':249.0748539, 'weight':250.578118975},
                        'Es':{'protons':99, 'nucleons':252, 'name':'Einsteinium', 'mass':252.08298, 'weight':252.08298},
                        'Fm':{'protons':100, 'nucleons':257, 'name':'Fermium', 'mass':257.0951061, 'weight':257.0951061},
                        'Md':{'protons':101, 'nucleons':258, 'name':'Mendelevium', 'mass':258.0984315, 'weight':259.10104075},
                        'No':{'protons':102, 'nucleons':259, 'name':'Nobelium', 'mass':259.10103, 'weight':259.10103},
                        'Lr':{'protons':103, 'nucleons':262, 'name':'Lawrencium', 'mass':262.10961, 'weight':262.10961},
                        'Rf':{'protons':104, 'nucleons':267, 'name':'Rutherfordium', 'mass':267.12179, 'weight':267.12179},
                        'Db':{'protons':105, 'nucleons':268, 'name':'Dubnium', 'mass':268.12567, 'weight':268.12567},
                        'Sg':{'protons':106, 'nucleons':271, 'name':'Seaborgium', 'mass':271.13393, 'weight':271.13393},
                        'Bh':{'protons':107, 'nucleons':272, 'name':'Bohrium', 'mass':272.13826, 'weight':272.13826},
                        'Hs':{'protons':108, 'nucleons':270, 'name':'Hassium', 'mass':270.13429, 'weight':270.13429},
                        'Mt':{'protons':109, 'nucleons':276, 'name':'Meitnerium', 'mass':276.15159, 'weight':276.15159},
                        'Ds':{'protons':110, 'nucleons':281, 'name':'Darmstadtium', 'mass':281.16451, 'weight':281.16451},
                        'Rg':{'protons':111, 'nucleons':280, 'name':'Roentgenium', 'mass':280.16514, 'weight':280.16514},
                        'Cn':{'protons':112, 'nucleons':285, 'name':'Copernicium', 'mass':285.17712, 'weight':285.17712},
                        'Nh':{'protons':113, 'nucleons':284, 'name':'Nihonium', 'mass':284.17873, 'weight':284.17873},
                        'Fl':{'protons':114, 'nucleons':289, 'name':'Flerovium', 'mass':289.19042, 'weight':289.19042},
                        'Mc':{'protons':115, 'nucleons':288, 'name':'Moscovium', 'mass':288.19274, 'weight':288.19274},
                        'Lv':{'protons':116, 'nucleons':293, 'name':'Livermorium', 'mass':293.20449, 'weight':293.20449},
                        'Ts':{'protons':117, 'nucleons':292, 'name':'Tennessine', 'mass':292.20746, 'weight':292.20746},
                        'Og':{'protons':118, 'nucleons':294, 'name':'Oganesson', 'mass':294.21392, 'weight':294.21392},
                       };   
                       # data based on values from NIST database: Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015), Atomic Weights and Isotopic Compositions (version 4.1). [Online] Available: <http://physics.nist.gov/Comp> [2018/06/06]. National Institute of Standards and Technology, Gaithersburg, MD. (see: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl)                    
        # here the <mass = monoisotopic mass> and <weight = molecular weight> 
        # - mass of an element [Da] is the (relative atomic) mass of the principal (most abundant) isotope (monoisotopic value as opposed to the isotopic average). if abundance information was unavailable the mass of the lightest isotope was used.
        # - weight of an element [Da] is the abundance-weighted isotopic average (relative atomic) mass (as in the molecular weight). if abundance information was unavailable the unweighted average mass of all isotopes was used.
        # this dictionary was automatically generated with the script & spreadsheet in the './Chemical/periodic table/' subdirectory
         
        particles = {
                        'e':{'protons':0, 'nucleons':0, 'name':'electron', 'mass':5.48579909070e-4},
                        'p':{'protons':1, 'nucleons':1, 'name':'proton', 'mass':1.007276466879},
                        'n':{'protons':0, 'nucleons':1, 'name':'neutron', 'mass':1.00866491588},                
                    }; 
        self.__table.update(particles);   
        # masses [Da] of each particle was taken from NIST database (see: https://physics.nist.gov/cuu/Constants/index.html)           

    def get_property(self,element,field):
        value = None;
        if element in self.__table:
            value = self.__table[element][field];
        else:
            print('Failed to get property: Not given a valid element symbol.')
        return value;

    def table(self):
        return copy.deepcopy(self.__table);
    
    def mass(self,element):
        return self.get_property(element,'mass');
    
    def weight(self,element):
        return self.get_property(element,'weight');    
    
    def name(self,element):
        return self.get_property(element,'name');
    
    def atomic_number(self,element): # atomic number = number of protons
        return self.get_property(element,'protons');    

    def protons(self,element):
        return self.get_property(element,'protons');
    
    def mass_number(self,element): # mass number = number of nucleons = number of protons + number of neutrons
        return self.get_property(element,'nucleons');    
    
    def nucleons(self,element):
        return self.get_property(element,'nucleons');
    
    def neutrons(self,element):
        return self.get_property(element,'nucleons') - self.get_property(element,'protons');
        
    def formula_parse(self, formula): # assumes simplified molecular formula of uncharged compound (does not support parenthesis "()", reduced or oxidized forms "+"/"-", or spacing)
        atoms = [];
        regex = re.findall(r'([A-Z][a-z]*)(\d*)', formula);
        for expression in regex:
            symbol = expression[0];
            count = expression[1];
            if count=='':
                count = 1;
            else:
                count = int(count);
            atoms.append({'symbol':symbol,'count':count});
        return atoms
        
    def get_formula_property(self,formula,field):
        atoms = self.formula_parse(formula);
        value = 0;
        for atom in atoms:
            value = value + atom['count']*self.__table[atom['symbol']][field];            
        return value;   
        
    def formula_mass(self,formula): # get mass of a compound from its molecular formula
        return self.get_formula_property(formula,'mass');   

    def formula_weight(self,formula): # get weight of a compound from its molecular formula
        return self.get_formula_property(formula,'weight');  
    
    def formula_mass_number(self,formula): # get mass number of a compound from its molecular formula
        return self.get_formula_property(formula,'nucleons');  
            
    def formula_sum(self,stoichiometric_list,*args):  
        # make stoichiometrically specified combination of formulas to get the empirical formula of a product ()
        #   a stoichiometric_list is a list of dictaries each with a 'formula' and 'moles' field 
        #       ex: stoichiometric_list = [{'moles':2,'formula':'H2'},{'moles':1,'formula':'O2'}] 
        #           stoichiometric_list = [{'moles':2,'formula':'KMnO4'},{'moles':3,'formula':'H2O2'}] 
        #   the number of moles for a reagent should be set as positive and for a side product should be set as negative
        #   moles must be integers 
        # if you know how many moles there will be of a product, you can get its correct molecular formula by simply also providing the number of moles as an argument
        formula_sum = ''; 
        elemental_composition = dict();
        for stoichiometric in stoichiometric_list:
            formula = stoichiometric['formula'];
            moles = stoichiometric['moles'];
            atoms = self.formula_parse(formula);
            for atom in atoms:
               symbol = atom['symbol'];
               count = atom['count']*moles;
               if symbol in elemental_composition:
                   elemental_composition[symbol] = elemental_composition[symbol] + count;                   
               else:
                   elemental_composition[symbol] = count;
        #print(elemental_composition)
        for element, count in elemental_composition.items():
            if len(args)>0:
                moles_of_product = args[0];
                count = count//moles_of_product;
            if isinstance(count,int):
                if count>0:
                    if count==1:
                        formula_sum = formula_sum + element;                        
                    else:
                        formula_sum = formula_sum + element + str(count);
            else:
                raise Exception('Failed to sum formulas: Element count was not an integer.'+'('+element+':'+str(count)+')')
        return formula_sum; 
        
###############################################################################
        
class Compound(object):
    
    # Compound
    #   a chemical compound    
    
    def __init__(self, *args): 
                                        
        # class constants
        self.__required_fields = [{'name':'Name','attribute':'name','type':str},
                                  {'name':'Type','attribute':'type','type':str},
                                  {'name':'ID [PubChem]','attribute':'id','type':int},
                                  {'name':'Molecular Formula','attribute':'formula','type':str},
                                  {'name':'Mass [Da]','attribute':'mass','type':(int, float)}];
        self.__auto_fetch_pubchem_info = False; # should available properties be auto fetched from PubChem info                                            
        self.__overwrite_with_pubchem_info = True; # should available properties be updated based on PubChem info          
        self.__calculate_mass_from_formula = True; # should we calculate monoisotopic mass based on the molecular formula of the compound        
        self.__number_of_synonyms_to_keep = 5; # number of compound synonyms to keep
        
        # class properties
        self.id = 0; # unique chemical identifier (when specified this is equivalent to the PubChem ID, for compounds without a provided ID or that are not in the PubChem database this value will be a negative-valued unique local serial number) 
        self.name = ''; 
        self.type = '';    
        self.mass = -1; # monoisotopic mass [g/mol] = [Da]
        self.formula = ''; # simplified linear molecular formula
        self.info = dict(); # additional info about the compound        

        # update properties from given 
        if len(args)>0: 
            given = args[0];
            if isinstance(given, dict): # given dictionary
                self.from_dictionary(given); 
            elif isinstance(given, int): # given PubChem ID
                self.id = given;
            else:
                print('Incorrectly formatted input. Given data not used.')
            if not self.is_valid():
                print('Compound not valid. Please specify compound ID, mass, and/or formula.')
                
        # generate unique ID if ID is nonvalid 
        if not self.has_valid_id():
            self.get_unique_id();  

        # update data from formula
        if self.__calculate_mass_from_formula:
            self.get_mass_from_formula();             
                            
        # update data from PubChem
        if self.__auto_fetch_pubchem_info:
            self.fetch_PubChem_info();
                                    
    def has_valid_id(self):
        # check for valid compound ID (for PubChem)
        return self.id>0;   

    def is_valid(self):
        # check if a compound is sufficiently specified
        return (self.id != 0) or (self.mass != -1) or (self.formula != '');            

    def get_mass_from_formula(self):
        # get mass from the molecular formula
        if self.formula != '':
            elements = PeriodicTableOfElements();
            self.mass = elements.formula_mass(self.formula);
                            
    def fetch_PubChem_info(self):
        # fetch compound data from PubChem (http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest)
        if self.has_valid_id():
            try:
                # get PubChem data
                pubchemdata = pubchempy.Compound.from_cid(self.id);
                                
                # append PubChem data
                self.info['PubChem CID'] = pubchemdata.cid;
                self.info['SMILES'] = pubchemdata.isomeric_smiles; # options: {isomeric_smiles,canonical_smiles}               
                #self.info['Molecular Formula'] = pubchemdata.molecular_formula; # used below
                self.info['Molar Mass [g/mol]'] = pubchemdata.molecular_weight; # The molecular weight is the sum of all atomic weights of the constituent atoms in a compound, measured in g/mol. In the absence of explicit isotope labelling, averaged natural abundance is assumed. If an atom bears an explicit isotope label, 100% isotopic purity is assumed at this location.
                self.info['Monoisotopic Mass [Da]'] = pubchemdata.monoisotopic_mass; # The mass of a molecule, calculated using the mass of the most abundant isotope of each element.             
                self.info['Exact Mass [Da]'] = pubchemdata.exact_mass; # The mass of the most likely isotopic composition for a single molecule, corresponding to the most intense ion/molecule peak in a mass spectrum.            
                self.info['Synonyms'] = [pubchemdata.synonyms[n] for n in range(0,min([self.__number_of_synonyms_to_keep,len(pubchemdata.synonyms)]))];                                
                self.info['InChI'] = {'name':pubchemdata.iupac_name, 'string':pubchemdata.inchi, 'key':pubchemdata.inchikey}; #  IUPAC International Chemical Identifier                                                                           
                self.info['Net Charge'] =  pubchemdata.charge;
                self.info['Complexity'] =  pubchemdata.complexity; # molecular complexity rating obtained with the Bertz/Hendrickson/Ihlenfeldt formula.
                self.info['XLogP'] =  pubchemdata.xlogp; # computationally generated octanol-water partition coefficient or distribution coefficient. XLogP is used as a measure of hydrophilicity or hydrophobicity of a molecule.
                self.info['Polar Surface Area'] =  pubchemdata.tpsa; # Topological polar surface area, computed by the algorithm described in the paper by Ertl et al.                
                # PubChem data fields not currently supported
                #self.info['Coordinate Type'] = pubchemdata.coordinate_type;                
                #self.info['Counts'] =  {'':pubchemdata.atoms,
                #                        '':pubchemdata.bonds,    
                #                        '':pubchemdata.atom_stereo_count,
                #                        '':pubchemdata.bond_stereo_count,
                #                        '':pubchemdata.rotatable_bond_count,
                #                        '':pubchemdata.isotope_atom_count,
                #                        '':pubchemdata.h_bond_acceptor_count,
                #                        '':pubchemdata.h_bond_donor_count,
                #                        '':pubchemdata.heavy_atom_count,
                #                        '':pubchemdata.defined_atom_stereo_count,
                #                        '':pubchemdata.defined_bond_stereo_count,
                #                        '':pubchemdata.covalent_unit_count,
                #                        '':pubchemdata.undefined_atom_stereo_count,
                #                        '':pubchemdata.undefined_bond_stereo_count};   
                #self.info['3D'] =  {'':pubchemdata.effective_rotor_count_3d,
                #                    '':pubchemdata.feature_selfoverlap_3d,
                #                    '':pubchemdata.mmff94_energy_3d,
                #                    '':pubchemdata.mmff94_partial_charges_3d,
                #                    '':pubchemdata.pharmacophore_features_3d,
                #                    '':pubchemdata.multipoles_3d,
                #                    '':pubchemdata.shape_fingerprint_3d,
                #                    '':pubchemdata.shape_selfoverlap_3d,
                #                    '':pubchemdata.conformer_id_3d,
                #                    '':pubchemdata.conformer_rmsd_3d,
                #                    '':pubchemdata.volume_3d};
                
                # overwrite available properties with PubChem data
                if self.__overwrite_with_pubchem_info:
                    if self.name.strip() == '':
                        if len(pubchemdata.synonyms)>0:
                            self.name = pubchemdata.synonyms[0];
                        else:
                            self.name = pubchemdata.iupac_name;
                    self.mass = pubchemdata.monoisotopic_mass;
                    self.formula = pubchemdata.molecular_formula;                            
                    if self.__calculate_mass_from_formula:
                        self.get_mass_from_formula();
                    
                # give all PubChem data
                return pubchemdata;
            except:
                print('Failed to load info from PubChem. Review compound ID ('+str(self.id)+') and ensure internet is accessible.')
        else:
            print('Failed to load info from PubChem: Invalid chemical ID ('+str(self.id)+').')             
    
    def from_dictionary(self,dictionary):
        # convert dictionary to compound
        dictionary = copy.deepcopy(dictionary);
        for field in self.__required_fields:
            if field['name'] in dictionary:
                value = dictionary[field['name']];
                if isinstance(value,field['type']):
                    setattr(self, field['attribute'], value);
                    del dictionary[field['name']];
                else:
                    print('Failed to load compound: Field with incorrect type given.');    
            else:
                print('Failed to load compound: Missing field "'+ field['name'] +'".');    
        self.info = dictionary;
        
    def to_dictionary(self):
        # convert compound to dictionary
        dictionary = self.info;
        for field in self.__required_fields:
            dictionary.update({field['name']:getattr(self, field['attribute'])})
        return dictionary

    def __str__(self):
        return str(self.to_dictionary())     
    
    def get_unique_id(self):
        # provides a unique ID 
        
        # get unique timestamp
        date = (datetime.datetime.now()).strftime("%Y%m%d%H%M%S%f");
        
        # get random number
        lengthRand = 16;
        rand = str(round(random.uniform(0, 1)*pow(10,lengthRand)));
        
        # merge timestamp and random number 
        ID = int(rand + date);
        
        # assign the ID
        self.id = -ID;
        
#    def plot(self):
#        key = 'SMILES';
#        if key in self.info:
#            smiles = self.info[key];
#            molecule = rdkit.Chem.MolFromSmiles(smiles);
#            rdkit.Chem.Draw.MolToMPL(molecule);
        
###############################################################################

class CompoundDataBase(object):
    
    # Compound DataBase
    #   a database of chemical compounds (can save and load '_cdb.json' files)  
    

    def __init__(self, *args): 
                                        
        # class constants
        self.__format_name = 'Chemical CompoundDataBase';
        self.__format_extension = '_cdb.json'; # class file format 
        self.__format_version = '1.0'; # version of the file format supported by this class               
        
        # class properties
        self.version = '1.0'; # version of the file format used in the file
        self.filename = 'my'+self.__format_extension; # name of the file to import data from                   
        self.name = ''; 
        self.data = []; # list of compounds 
        # could have also defined database data as a dictionary where the compound IDs were keys 
        
        # update properties from given 
        for given in args:
            if isinstance(given,list):
                compounds = [];
                if all(isinstance(each, Compound) for each in given): # given list of chemical compounds
                    compounds = copy.deepcopy(given);
                elif all(isinstance(each, dict) for each in given): # given list of formatted dictionaries
                    compounds = self.get_compounds_from_dictionaries(given);                        
                elif all(isinstance(each, int) for each in given): # given list of PubChem IDs
                    compounds = self.get_compounds_from_ids(given);
                else:
                    print('Failed to load given data: Incorrectly formatted list')
                self.append_compounds(compounds); 
            else: 
                if isinstance(given, str): # given a string 
                    self.name = given;
                else:
                    print('Failed to load given data: Not provided a list or name.')

    def __len__(self): 
        return len(self.data);
    
    def __str__(self): 
        return str(self.get_names());
    
    def get_mass_from_formula(self):
        for compound in self.data:         
            compound.get_mass_from_formula();   
        
    def fetch_PubChem_info(self):
        for compound in self.data:         
            compound.fetch_PubChem_info();

    def append_compounds(self,compound_list):
        ids = self.get_ids();
        for compound in compound_list:
            if compound.id not in ids:
                self.data.append(compound);
            else:
                print('Adding compound failed: Compound already exists in the database (ID: '+compound.id+').' );            
                
    def merge(self,database_list): # merge multiple chemical data bases
        for database in database_list:
            self.append_compounds(database.data);                
            
    def get_compounds_from_dictionaries(self,dictionary_list):
        compound_list = [];
        for dictionary in dictionary_list:  
            compound = Compound(dictionary);
            compound_list.append(compound);
        return compound_list; 
            
    def get_compounds_from_ids(self,id_list):
        compound_list = [];        
        for ID in id_list:  
            compound = Compound(ID);
            compound_list.append(compound);
        return compound_list;  
    
    def get_dictionaries(self):
        dictionary_list = [];
        for compound in self.data:  
            dictionary_list.append(compound.to_dictionary());
        return dictionary_list;                
    
    def get_ids(self):
        id_list = [];
        for compound in self.data:  
            id_list.append(compound.id);
        return id_list;    
    
    def get_formulas(self):
        formula_list = [];
        for compound in self.data:  
            formula_list.append(compound.formula);
        return formula_list;       

    def get_masses(self):
        mass_list = [];
        for compound in self.data:  
            mass_list.append(compound.mass);
        return mass_list;   

    def get_names(self):
        name_list = [];
        for compound in self.data:  
            name_list.append(compound.name);
        return name_list;

    def get_types(self):
        type_list = [];
        for compound in self.data:  
            type_list.append(compound.type);
        return type_list;            
            
    def get_values(self,field):
        value_list = [];        
        for dictionary in self.get_dictionaries(): 
            if field in dictionary:
                value_list.append(dictionary[field]);
            else:
                print('Error loading values: Field missing.');
        return value_list;    

    def unique(self): # check if compound ids are unique
        id_list = self.get_ids();
        return len(set(id_list)) == len(id_list)      
    
    def lookup(self,ID): # lookup compound from id
        compound = None;
        id_list = self.get_ids();
        if ID in id_list:
            index = id_list.index(ID);
            compound = self.data[index]; 
        return compound    

    def lookup_multiple(self,ids): # lookup compounds from ids
        compound_list = [];
        id_list = self.get_ids();
        if not isinstance(ids,list): 
            ids = [ids]; # make list from single integer input
        for ID in ids:
            if ID in id_list:
                index = id_list.index(ID);
                compound_list.append(self.data[index]);
            else:
                compound_list.append(None);
        return compound_list;                
    
    def save(self,filename):
        if not os.path.isfile(filename):        
            with open(filename, 'w') as fileID:
                data = {'name':self.name,'data':self.get_dictionaries(),'format':{'name':self.__format_name,'version':self.__format_version},'date':(datetime.datetime.now()).strftime("%Y/%m/%d %H:%M")};
                json.dump(data, fileID, sort_keys=True, indent=4);
        else:
            raise Exception('Saving CompoundDataBase Failed: File already exists.');        

    def load(self,filename): 
        if os.path.isfile(filename):        
            with open(filename, 'r') as fileID:
                data = json.load(fileID);
                self.name = data['name'];
                self.version = data['format']['version'];
                dictionary_list = data['data'];                
                compound_list = [];
                for dictionary in dictionary_list:  
                    compound = Compound(dictionary);
                    compound_list.append(compound);
                self.data = compound_list;
        else:
            raise Exception('Reading CompoundDataBase Failed: File was not found.');

###############################################################################
            
class Solution(object):
    
    # Solution
    #   a chemical solution comprised of compounds (specified by ID) of known concentrations in a known volume
       
    def __init__(self, *args): 
        
        # class constants
        self.__required_fields = [{'name':'Volume [L]','attribute':'volume','type':(int, float)},                                  
                                  {'name':'Concentration [M]','attribute':'concentrations','type':list},
                                  {'name':'ID [PubChem]','attribute':'compounds','type':list}];    
        self.__transfer_fields ={'compounds':['Volume [L]','Concentration [M]','ID [PubChem]'],
                                 'solutions':['Volume [L]','Solution']};                                  
                       
        # class properties
        self.volume = 0; # total solution volume in liters [L]
        self.compounds = []; # list of ids of the compounds in the solution (intended to be used in parallel with a compound database)
        self.concentrations = []; # list of the molar concentrations of the compounds in the solution [M]
        # could have also defined the content of a solution as a list of dictionaries (of compounds, each with fields: concentration & ID)

        # update properties from given 
        if len(args)>0: 
            given = args[0];
            if isinstance(given, dict): # given dictionary
                self.from_dictionary(given); 
            else:
                print('Incorrectly formatted input. Given data not used.')

    def from_dictionary(self,dictionary):
        # convert dictionary to solution
        dictionary = copy.deepcopy(dictionary);
        for field in self.__required_fields:
            if field['name'] in dictionary:
                value = dictionary[field['name']];
                if isinstance(value,field['type']):
                    setattr(self, field['attribute'], value);
                    del dictionary[field['name']];
                else:
                    print('Failed to load solution: Field with incorrect type given.');    
            else:
                print('Failed to load solution: Missing field "'+ field['name'] +'".');    
        
    def to_dictionary(self):
        # convert solution to dictionary
        dictionary = dict();
        for field in self.__required_fields:
            dictionary.update({field['name']:getattr(self, field['attribute'])})
        return dictionary
    
    def __str__(self):
        return str(self.to_dictionary())     
    
    def add(self,content):
        if isinstance(content,dict):                                   
            if all(field in content for field in self.__transfer_fields['solutions']):  # adding solutions
                self.add_solutions(content);
            elif all(field in content for field in self.__transfer_fields['compounds']): # adding compounds
                self.add_compounds(content);
            else:
                print('Failed to add content: Input is missing a field.')        
        else:
            print('Failed to add content: Given input type is unsupported.') 
        
    def add_solutions(self,content):
        # add solution contents 
        
        # format inputs
        volumes = content['Volume [L]'];
        solutions = content['Solution'];
        if isinstance(volumes,list) and isinstance(solutions,list):
                
            # initialize variables
            ids = [];
            masses = []; # for the default units this is in moles [mol] and as such is not technically a mass [g]
            concentrations = [];      
                                                        
            # get final volume (including the destination solution)
            solutions = solutions + [self];                
            volumes = volumes + [self.volume];
            total_volume = sum(volumes);        
                    
            # convert all concentrations to masses and flatten into list
            ids_flat = [];
            masses_flat = [];
            N = len(solutions);
            for n in range(0,N):
                
                # get values
                CID_list = solutions[n].compounds;
                C_list = solutions[n].concentrations;
                V = volumes[n];
                
                # compute mass to be added
                m_list = [V*C for C in C_list];
                
                # collect compound id and mass
                ids_flat = ids_flat + CID_list;
                masses_flat = masses_flat + m_list;
                
                # remove the volume from source solutions
                if n < (N-1):
                    solutions[n].volume = solutions[n].volume - V;   
                                
            # combine masses for identical ids            
            for ID in set(ids_flat):
               indices = [index for index, value in enumerate(ids_flat) if value == ID]; 
               total_mass = 0;
               for index in indices:
                   total_mass = total_mass + masses_flat[index];
               concentration = total_mass/total_volume;
               masses.append(total_mass);   
               concentrations.append(concentration);
               ids.append(ID);
                       
            # update properties
            self.compounds = ids;   
            self.concentrations = concentrations;               
            self.volume = total_volume;  
            
        else:
            print('Failed to add content: Input not a list.')
            
    def add_compounds(self,content):
        # add solution contents 

        # format inputs        
        volumes = content['Volume [L]'];
        compounds = content['ID [PubChem]'];
        concentrations = content['Concentration [M]'];
        if isinstance(volumes,list) and isinstance(compounds,list) and isinstance(concentrations,list):
            
            # initialize variables
            recontent = dict(); # reformat content as volumes and solution
            recontent['Volume [L]'] = volumes;
            
            # construct solutions from inputs
            solutions = [];
            for n in range(0,len(compounds)):                
                initial = {'Volume [L]':volumes[n], 'ID [PubChem]':[compounds[n]], 'Concentration [M]':[concentrations[n]]};
                solutions.append(Solution(initial));  
            recontent['Solution'] = solutions;
            
            # add newly made solutions to current
            self.add_solutions(recontent)
            
            # provide the intermediate solutions and volumes to add of each 
            return recontent;     
        else:
            print('Failed to add content: Input not a list.')
            
    def get_concentration(self,compound):
        concentration = 0;
        if compound in self.compounds:
            concentration = self.concentrations[self.compounds.index(compound)];
        return concentration;            
            
    def get_compounds(self,database):
        compound_list = [];
        for n in range(0,len(self.compounds)):
            ID = self.compounds[n];
            C = self.concentrations[n];
            compound = copy.deepcopy(database.lookup(ID));
            compound.info['Concentration [M]'] = C;
            compound_list.append(compound);
        database_subset = CompoundDataBase(compound_list);
        return database_subset; # compounds in solution given as a compound database
    
    def clean_up(self):
        # remove compounds with a concentration of zero
        compounds = []
        concentrations = [];
        for n in range(0,len(self.compounds)):
            if not self.concentrations[n]==0:
                compounds.append(self.compounds[n]);
                concentrations.append(self.concentrations[n]);
        self.compounds = compounds;   
        self.concentrations = concentrations;
  
############################################################################### 
            
class SimpleReaction(object):
    
    # SimpleReaction
    #   a chemical reaction between compounds in a solution based on simple rules (the known reaction equation) and compound types. This class performs highly simplisitic, rule based reactions and as such uses no form of numerical chemical simulation to determine likelihood or validity or true structure of product.
    
    def __init__(self): 
        
        # class constants
        self.__auto_cleanup = False; # automatically remove solution entries that end up with a concentration of zero                                                   
        # self.__track_in_formula = True; # if set to True reactions will be conducted based on formulas, if False they will be based on mass
        self.__supported_reactions = {
            'Ugi':{
                   'inputs': {'amine':1,'aldehyde':1,'carboxylic acid':1,'isocyanide':1}, # type:moles (for each "input" type there must be a matching single compound in the solution for the reaction to be successful)
                   'side products': {'H2O':1}, # formula:moles
                   'side products IDs': [962], # ID [PubChem]
                   'product': {'ugi':1} # (amide) # type:moles # must list only on product whose fomula is unknown and to be determined (all other products should be listed as side products)
                   # equation: 1 (amine) + 1 (aldehyde) + 1 (carboxylic acid) + 1 (isocyanide) -> 1 (ugi) + 1 H2O
                  } # multi-component reaction in organic chemistry involving an aldehyde, an amine, an isocyanide and a carboxylic acid to form a bis-amide.' <https://en.wikipedia.org/wiki/Ugi_reaction>
        };        
        
    def react(self,reaction,efficiency,solution,database):
        # perform a specific kind of reaction in a solution with a set conversion efficieny (of reagents to products)
        success = False;
        if isinstance(solution,Solution) and isinstance(database,CompoundDataBase) and 0<=efficiency<=1:
            sol = copy.deepcopy(solution); # copy so to not modify these until end
            db = copy.deepcopy(database); # copy so to not modify these until end       
            if reaction in self.__supported_reactions:
                elements = PeriodicTableOfElements();                
                compounds_in_solution = sol.get_compounds(db);
                types = compounds_in_solution.get_types();
                IDs = compounds_in_solution.get_ids();            
                formulas = compounds_in_solution.get_formulas();  
                concentrations = compounds_in_solution.get_values('Concentration [M]');            
                required_types = list(self.__supported_reactions[reaction]['inputs'].keys());
                required_types_moles = list(self.__supported_reactions[reaction]['inputs'].values());
                side_product_formulas = list(self.__supported_reactions[reaction]['side products'].keys());
                side_product_moles = list(self.__supported_reactions[reaction]['side products'].values());
                side_product_IDs = self.__supported_reactions[reaction]['side products IDs'];
                side_product_concentrations = [];
                product_types = list(self.__supported_reactions[reaction]['product'].keys())[0];
                product_moles = list(self.__supported_reactions[reaction]['product'].values())[0];
                product_concentrations = 0;    
                new_concentrations = []
                new_ids = [];
                # check if solution contains expected reagents
                if set(required_types).issubset(types):
                    products_compound_list = [];                
                    reagent_dict = dict();
                    stoichiometric_list = [];
                    counter = 0;
                    for m in range(0,len(required_types)):
                        reqtyp = required_types[m];  
                        moles = required_types_moles[m];   
                        for n in range(0,len(types)):  
                            typ = types[n];
                            ID = IDs[n];
                            formula = formulas[n];                        
                            if typ == reqtyp: # at this point currently assumes the solution contains only a single compound per each type 
                                counter = counter + 1;
                                per_volume_multiplier = concentrations[m]/moles; # this is a constant a constant                               
                                concentrations[m] = (1.00-efficiency)*concentrations[m];                                
                                reagent_dict.update({reqtyp:ID});
                                stoichiometric_list.append({'moles':moles,'formula':formula})
                                new_concentrations.append(concentrations[m]);
                                new_ids.append(ID); 
                    for o in range(0,len(side_product_formulas)):
                        moles =  side_product_moles[o];    
                        ID = side_product_IDs[o];
                        formula = side_product_formulas[o];  
                        mass = elements.formula_mass(formula);       
                        concentration = per_volume_multiplier*moles*efficiency;
                        side_product_concentrations.append(concentration);                                                        
                        stoichiometric_list.append({'moles':-moles,'formula':formula})
                        side_product = Compound({'Type':'side product', 'Name':'', 'ID [PubChem]':ID, 'Mass [Da]':mass, 'Molecular Formula':formula});
                        side_product.fetch_PubChem_info();  
                        products_compound_list.append(side_product);
                        new_concentrations.append(concentration);
                        new_ids.append(ID);                                                        
                    product_formula = elements.formula_sum(stoichiometric_list,product_moles);  
                    product_concentrations = per_volume_multiplier*product_moles*efficiency;
                    product_mass = elements.formula_mass(product_formula);                                
                    product = Compound({'Type':product_types, 'Name':product_types + ' product X', 'ID [PubChem]':-1, 'Mass [Da]':product_mass, 'Molecular Formula':product_formula});
                    product.info['Precursor'] = reagent_dict;
                    products_compound_list.append(product);
                    new_concentrations.append(product_concentrations);
                    new_ids.append(product.id); 
                    
                    # update status
                    success = product.id; # provide status as the product ID (since nonzero integers are boolean Trues)
                    if counter > len(required_types):
                        print('Reaction failed: Multiple compounds in the solution express the same required type.')                    
                        success = False;
                    
                    # update the input objects
                    if success:              
                        
                        # update solutions
                        solution.compounds = new_ids;
                        solution.concentrations = new_concentrations;
                        if self.__auto_cleanup:
                            solution.clean_up();
                        
                        # update database
                        database.append_compounds(products_compound_list);
                    
                else:
                    print('Reaction failed: Could not find all required reagents.')                    
            else:
                print('Reaction failed: Given unsupported reaction type.')
        else:
            print('Reaction failed: Given incorrect input type.')
        return success; 
    
    def get_equation(self,reaction):
        equation = '';
        if reaction in self.__supported_reactions:
            info = self.__supported_reactions[reaction];
            reagents = info['inputs'];
            product = info['product'];
            sideproducts = info['side products'];            
            for key, value in reagents.items():                
                equation = equation + str(value) + '(' + key + ') + ';   
            equation = equation[:-2]
            equation = equation + '-> '; 
            for key, value in product.items():                
                equation = equation + str(value) + '(' + key + ') + ';  
            if len(sideproducts.keys())>0:
                for key, value in sideproducts.items():                
                    equation = equation + str(value) + '(' + key + ') + ';    
                equation = equation[:-2]   
            else:
                equation = equation[:-2]                
        else:
            print('Given unsupported reaction type.')
        return equation;       
        
#     def add_new_reaction():
#          allow user to define new reaction type

        
###############################################################################
    
class Container(object):
    # Container
    # a container of solutions 
        
    def __init__(self, *args): 
        
        # class constants
        self.__enforce_volume_limits = True; # choose whether or not volume limits are enforced/checked where applicable   
        self.__autofill = False; # should the positions of new solutions be automatically chosen to as empty/available positions in the container 
        self.__required_fields = [{'name':'Type','attribute':'type','type':(str,dict)},                                  
                                  {'name':'Solutions','attribute':'solutions','type':list},
                                  {'name':'Positions','attribute':'positions','type':list},                                  
                                  ];    
        self.__type_fields = [{'name':'Type','attribute':'type','type':str},
                              {'name':'Rows','attribute':'rows','type':int},
                              {'name':'Columns','attribute':'columns','type':int},
                              {'name':'Dead Volume [L]','attribute':'volume_dead','type':(int,float)},
                              {'name':'Max Volume [L]','attribute':'volume_max','type':(int,float)},
                              {'name':'Category','attribute':'category','type':str},
                              ];                                  
        self.__native_types = {
            'WellPlate384':{
                   'rows': 16,
                   'columns': 24,
                   'dead volume': 10e-6,
                   'max volume': 120e-6,                    
                   'category': 'microplate',   
                  },
            'WellPlate96':{
                   'rows': 8,
                   'columns': 12,
                   'dead volume': 100e-6,
                   'max volume': 400e-6,                    
                   'category': 'microplate',   
                  },
            'WellPlate24':{
                   'rows': 4,
                   'columns': 6,
                   'dead volume': 450e-6,
                   'max volume': 3.4e-3,                    
                   'category': 'microplate',   
                  },
            'Tube50mL':{
                   'rows': 1,
                   'columns': 1,
                   'dead volume': 500e-6,
                   'max volume': 50e-3,                    
                   'category': 'tube',   
                  },
            'Tube15mL':{
                   'rows': 1,
                   'columns': 1,
                   'dead volume': 250e-6,
                   'max volume': 15e-3,                    
                   'category': 'tube',   
                  },
            'Tube1.5mL':{
                   'rows': 1,
                   'columns': 1,
                   'dead volume': 50e-6,
                   'max volume': 1.5e-3,                    
                   'category': 'microtube',   
                  },
            'Tube2mL':{
                   'rows': 1,
                   'columns': 1,
                   'dead volume': 70e-6,
                   'max volume': 2e-3,                    
                   'category': 'microtube',   
                  },                                   
        };    
                                                                                
        # class properties
        self.type = ''; # type of container
        self.rows = 1; # number of rows in the container
        self.columns = 1; # number of columns in the container        
        self.volume_dead = 0; # dead volume for a well, in liters [L]
        self.volume_max = 1; # max volume allowed in a well, in liters [L]        
        self.category = ''; # category of container {microplate, tube, etc.}
        self.solutions = []; # list of solutions in the container
        self.positions = []; # list of positions of the solutions in the container               
        # could have also defined the solutions of container as a list of dictionaries (of solutions, each with fields: solution & position)        

        # update properties from given 
        if len(args)>0: 
            given = args[0];
            if isinstance(given, dict): # given dictionary
                self.container_from_dictionary(given); 
            elif isinstance(given, str): # given native container type 
                self.container_from_native_type(given);                 
            else:
                print('Incorrectly formatted input. Given data not used.')

    def get_native_types(self):
        return list(self.__native_types.keys());
        
    def container_from_native_type(self,native_type):
        if native_type in self.__native_types:
            presets = self.__native_types[native_type];
            self.type = native_type; 
            self.rows = presets['rows']; 
            self.columns = presets['columns']; 
            self.volume_max = presets['max volume']; 
            self.volume_dead = presets['dead volume']; 
            self.category = presets['category']; 
        else:
            print('Loading container presets failed: Unsupported native container type specified.')                    

    def container_from_custom_type(self,custom_type):  
        for field in self.__type_fields:
            if field['name'] in custom_type:
                value = custom_type[field['name']];
                if isinstance(value,field['type']):
                    setattr(self, field['attribute'], value);
                else:
                    print('Loading container failed: Unexpected value given in field (' + field['name'] + ').')
            else:
                print('Loading container failed: Missing container field (' + field['name'] + ').')
                            
    def container_from_dictionary(self,dictionary):
        # convert dictionary to container
        dictionary = copy.deepcopy(dictionary);
        for field in self.__required_fields:
            if field['name'] in dictionary:
                value = dictionary[field['name']];
                if isinstance(value,field['type']):
                    setattr(self, field['attribute'], value);
                    del dictionary[field['name']];
                else:
                    print('Failed to load solution: Field with incorrect type given.');    
            else:
                print('Failed to load solution: Missing field "'+ field['name'] +'".');  
        self.update_container_type();
            
    def update_container_type(self):        
        # update parameters from current type
        if not self.type == '':
            current_type = copy.deepcopy(self.type);
            if isinstance(current_type,str): # given a native type  
                self.container_from_native_type(current_type);
            if isinstance(current_type,dict): # given a custom type
                self.container_from_custom_type(current_type);
         
         
#    def from_dictionary(self,dictionary):
#        # convert dictionary to container
#        dictionary = copy.deepcopy(dictionary);
#        for field in self.__required_fields:
#            if field['name'] in dictionary:
#                value = dictionary[field['name']];
#                if isinstance(value,field['type']):
#                    setattr(self, field['attribute'], value);
#                    del dictionary[field['name']];
#                else:
#                    print('Failed to load solution: Field with incorrect type given.');    
#            else:
#                print('Failed to load solution: Missing field "'+ field['name'] +'".');    
#        
#    def to_dictionary(self):
#        # convert container to dictionary
#        dictionary = dict();
#        for field in self.__required_fields:
#            dictionary.update({field['name']:getattr(self, field['attribute'])})
#        return dictionary
    
#    list of solutions
 # self.shape = (self._rows,self._cols)
 
#    # def react(self,type):
##    type -> ugi
##    reagent - compound list
##    volume list
##    concentration list
##    conversion efficiency
##    class
##    product - compound list

#    chemical database


###############################################################################
#            
#class WorkBench(object):
    
    # Work Bench
    #   a virtual chemical work bench of containers and a database of compounds (can save and load '_wb.json' files)  


# #not implemented yet
# opens cwb.json files
# comprised of a list of containers and a compound database
