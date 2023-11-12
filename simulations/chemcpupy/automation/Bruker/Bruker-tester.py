# -*- coding: utf-8 -*-
"""
@title: Bruker-tester.py
@description: script to test the database reading methods of the Bruker class
@author: chrisarcadia 
@created: 2018/10/12
"""
import os
from lxml import etree
import sqlite3
import os
import matplotlib.pyplot as pyplot
import Bruker

# Specify data file path
path_script = os.path.dirname(__file__)
path_data = path_script + '/resource/example.d/'; # data file
dir_contents = os.listdir(path_data);
dir_method = [path for path in dir_contents if path[-2:]=='.m'];
path_method = path_data + dir_method[-1] + '/';

# Manually load settings that are in XML files
file1 = path_method + 'apexAcquisition.method';
file2 = path_method + 'submethods.xml';
methodXML = Bruker.import_XML(file1);
submethodXML = Bruker.import_XML(file2);
#print(etree.tostring(methodXML))
method_name = methodXML.find('methodmetadata').find('primarykey').find('methodname').text;
method_report_node = methodXML.find('reportinfo');
def retrieve(head,title,attribute):
    return method_report_node.find('.//section[@title="'+ head +'"]//param[@title="' + title + '"]').attrib[attribute];    
report = {
    'name': method_name,
    'polarity_mode': retrieve('Polarity','Polarity','value'),   
    'laser_power_unit': retrieve('MALDI','Laser Power','unit'),
    'laser_power': float(retrieve('MALDI','Laser Power','value')), 
    'laser_shot_count': float(retrieve('MALDI','Number of Laser Shots','value')),
    'laser_frequency': float(retrieve('MALDI','Laser Frequency','value')),
    'laser_frequency_unit': retrieve('MALDI','Laser Frequency','unit'),  
    'laser_focus': retrieve('MALDI','LaserFocus','value'),       
    'smart_walk_enabled': retrieve('MALDI','MTP_SmartWalkEnable','value'),    
    'smart_walk_pattern': retrieve('MALDI','Smart walk pattern table','value'),        
    'smart_walk_grid_width': float(retrieve('MALDI','Smart walk grid width','value')),
    'smart_walk_grid_width_unit': retrieve('MALDI','Smart walk grid width','unit'),  
    'broadband_mass_unit': retrieve('Acquisition Mass Control','Broadband Low Mass','unit'),                                    
    'broadband_mass_low': retrieve('Acquisition Mass Control','Broadband Low Mass','value'),                                    
    'broadband_mass_high': retrieve('Acquisition Mass Control','Broadband High Mass','value'),                                                                      
    'detection_mode': retrieve('Acquisition Mass Control','Detection Mode','value'), 
}
print(report)

# Manually load settings that are in SQLite files
file3 = path_data + 'parameterChanges.sqlite';
parameter_changes = Bruker.import_SQLite(file3);
#parameter_changes.keys()
scan_Number = parameter_changes['ScanInfo']['ScanNumber'];
scan_TIC = parameter_changes['ScanInfo']['Tic'];
scan_MaxPeak = parameter_changes['ScanInfo']['MaxPeak'];
scan_TimePointMinute = parameter_changes['ScanInfo']['TimeMinutes'];
pyplot.plot(scan_Number, scan_TIC, 'bo', markersize=1)
pyplot.xlabel('Scan Number')
pyplot.ylabel('Total Ion Count')
pyplot.yscale('log')
pyplot.grid(True)
pyplot.show()
pyplot.plot(scan_Number, scan_MaxPeak, 'bo', markersize=1)
pyplot.xlabel('Scan Number')
pyplot.ylabel('Max Peak in Spectrum')
pyplot.yscale('log')
pyplot.grid(True)
pyplot.show()
pyplot.plot(scan_TimePointMinute,scan_Number, 'bo', markersize=1)
pyplot.xlabel('Time [min]')
pyplot.ylabel('Scan Number')
pyplot.grid(True)
pyplot.show()

# Manually load BASIC files
file4 = path_method + 'PULPROG__basic.compiled';
file5 = path_method + 'PULPROG__basic.event_seq';
BASIC_compiled = Bruker.import_NLDF(file4);
BASIC_event_sequence = Bruker.import_NLDF(file5);

# Manually load excite sweep
file6 = path_method + 'ExciteSweep';
excite_sweep = Bruker.import_NLDF_as_num(file6);
pyplot.plot(excite_sweep, 'bo', markersize=1)
pyplot.ylabel('Excite Sweep')
pyplot.grid(True)
pyplot.show()


