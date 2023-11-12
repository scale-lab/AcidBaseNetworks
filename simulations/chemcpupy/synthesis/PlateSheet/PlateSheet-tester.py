# -*- coding: utf-8 -*-

"""
@author: chrisarcadia 
@created: 2018/05/31
"""

# libraries
import PlateSheet

# load PlateSheet file 
filename = 'example_platesheet.csv';
sheet = PlateSheet.PlateSheet();    
sheet.load(filename);

# print properties 
print(sheet.format)
print() # newline
print(sheet.properties)
print() # newline
print(sheet.content)

# print specific content (available data fields: ['Class', 'Name', 'ID [PubChem]', 'Concentration [M]', 'Volume [L]', 'Positions'])
print() # newline
print(sheet.content[0]['Name']) # specific row in a column
print() # newline
print(sheet.list_values('Name')) # all rows in a column
print() # newline
print(sheet.restructure()) # restructure content for more concise representation

# modify the platesheet 
field = 'Name';
for n in range(0,len(sheet)):
    sheet.content[n][field] = sheet.content[n][field].replace(' ','-'); # replace all spaces in names with hyphens

#field = 'Volume [L]';
#for n in range(0,len(sheet)):
#    sheet.content[n][field] = sheet.content[n][field].split(';'); # replace all spaces in names with hyphens

# empty sheet (for a blank sheet)
# sheet.empty();
    
# save platesheet to file
sheet.save(filename.replace('_platesheet.csv','_saved_platesheet.csv'))

##########
