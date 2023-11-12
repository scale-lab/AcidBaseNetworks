#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: PlateSheet class
@description: class for loading "_platesheet.csv" files
@author: chrisarcadia 
@created: 2018/05/31
"""

import os, csv, copy, datetime

# class for handling Plate Sheet (_platesheet.csv) files 
class PlateSheet(object):
    
    def __init__(self):   
        # intialize class instance 
        
        # class constants        
        self.__format_extension = '_platesheet.csv'; # class file format    
        self.__format_version = '1.0'; # version of the file format supported by this class                      
        self.__fields = {'format':['Name','Version'],
                       'properties':['Plate', 'Barcode', 'Title', 'Author', 'Date',	'Description'],                       
                       'content':['Type', 'Name', 'ID [PubChem]', 'Concentration [M]', 'Volume [L]', 'Positions']}; # required fields per each section        
        self.__parsetype = {'single integers':['Barcode','ID [PubChem]'],
                            'delimited floats':['Volume [L]','Concentration [M]'],
                            'delimited strings':['Positions']}; # defines how certain fields are parsed  
        self.__delimiter = ';'; # value delimiter for a multi-valued entry                            
        self.__python_to_excel_row_offset = 2;        
        
        # class properties
        self.version = '1.0'; # version of the file format used in the file
        self.filename = 'my'+self.__format_extension; # name of the file to import data from

        # class sections           
        self.format = dict();
        self.content = []; # [dict()];
        self.properties = dict();                                     
                                                                                         
    def load(self,filename):         
        # load sheet data from csv
        self.filename = filename;        
        if os.path.isfile(self.filename):
            with open(self.filename, 'r') as fileID:
                
                # import the csv data
                filedata = list(csv.reader(fileID));
                
                # print the csv data 
                #print(filedata);  
                
                # read the format
                form = filedata[0][0];
                form = form.split(" v",1);
                del(filedata[0]); # remove it  
                del(filedata[0]); # remove blank row               
                self.format['Name'] = form[0].strip();
                self.format['Version'] = form[1].strip();                
                
                # read the properties
                labels = filedata[0]; 
                del(filedata[0]); # remove it  
                data = [filedata[0]];
                del(filedata[0]); # remove it  
                del(filedata[0]); # remove blank row                
                self.properties = self.read_rows_as_dictionaries(labels,data)[0];
    
                # read the content
                labels = filedata[0]; 
                del(filedata[0]); # remove it  
                data = filedata;
                del(filedata); # remove it             
                self.content = self.read_rows_as_dictionaries(labels,data);
                                                                         
                # check for required fields 
                fields_of_properties = list(self.properties.keys());
                if not set(self.__fields['properties']).issubset(fields_of_properties):
                    raise Exception('Reading PlateSheet Failed: File missing required property.');  
                fields_of_content = list(self.content[0].keys());
                if not set(self.__fields['content']).issubset(fields_of_content):
                    raise Exception('Reading PlateSheet Failed: File missing required content.'); 
                    
                # expand multivalued entries
                self.expand_multivalued_cells();                    
                    
                # check for duplicate IDs
                uniqueIDs = self.check_values_are_unique('ID [PubChem]')
                if not uniqueIDs:
                    raise Exception('Reading PlateSheet Failed: File contains repeated chemical IDs. Please consolidate repeated ID entries into a single entry.');                                        
                                                         
        else:
            raise Exception('Reading PlateSheet Failed: File was not found.');
            
    def empty(self):  
        # empty the platesheet (gives blank sheet)
        # self.format = {field: '' for field in self.__fields['format']}; # keep file format information                                    
        self.content = [{field: '' for field in self.__fields['content']}];
        self.properties = {field: '' for field in self.__fields['properties']};                                     
                        
    def __str__(self):
        return str(self.content);
    
    def __len__(self):
        return len(self.content);
        
    def list_values(self,field):
        # list all values in a given field
        listOfValues = [];
        for row in self.content:
            listOfValues.append(row[field]);
        return listOfValues;
    
    def check_values_are_unique(self,field):
        # check if the values in a given field are unique
        listOfValues = self.list_values(field);
        is_unique = len(set(listOfValues)) == len(listOfValues);
        return is_unique;
        
    def read_rows_as_dictionaries(self,labels,data):
        # read each csv row as a dictionary and append to a list
        dictionaries = [];
        for row in range(0,len(data)):
            rowdata = data[row];
            rowdictionary = dict();
            rowisnotempty = any(entry.split() for entry in rowdata);
            if rowisnotempty:
                for column in range(0,len(labels)):
                    label = labels[column];
                    labelisnotempty = label.strip() != '';
                    if labelisnotempty:
                        celldata = rowdata[column].strip();
                        cellisnotempty = celldata != '';
                        if cellisnotempty:
                            celldata_parsed = self.parse_cell(label, celldata);
                            rowdictionary.update({label:celldata_parsed})
                        else:
                            raise Exception('Reading PlateSheet Failed: File missing a cell entry (column: '+label+', row: '+str(row+self.__python_to_excel_row_offset)+')');
                dictionaries.append(rowdictionary);                                            
        return dictionaries
            
    def restructure(self):
        # restructure content data as a dictionary of lists instead of as a list of dictionaries (for compactness)
        content = dict();
        fields = self.__fields['content'];
        for field in fields:
            content[field] = self.list_values(field);
        return content;
    
    def parse_cell(self, label, celldata):   
        # parse the cell data according to entry type
        celldata_delimited = celldata.split(self.__delimiter);             
        if label in  self.__parsetype['single integers']:
            data = int(celldata);  
        elif label in self.__parsetype['delimited floats']:
            data = [];
            for n in range(0,len(celldata_delimited)):
                data.append(float(celldata_delimited[n]))
        elif label in self.__parsetype['delimited strings']:
            data = [];
            for n in range(0,len(celldata_delimited)):
                data.append(celldata_delimited[n].strip())
        else:
            data = celldata;                          
        return data;
                         
    def unparse_cells(self,dictionaries):  
        # undo parsing of cells for file export
        for row in range(0,len(dictionaries)):
            keys = dictionaries[row].keys()
            for col in keys:
                fielddata = dictionaries[row][col];
                if isinstance(fielddata, list):
                    data = str(fielddata[0]);
                    N = len(fielddata);
                    if N>1:
                        for n in range(1,N):
                            data = data + self.__delimiter + str(' ') + str(fielddata[n]);
                    dictionaries[row][col] = data;   
                else:
                    dictionaries[row][col] = fielddata;   
        return dictionaries; 

    def expand_multivalued_cells(self):
        # expand all single-valued cells to multi-valued cells if any field in that entry/row has multiple values
        for row in range(0,len(self)):   
            
            # get list of fields that can have multi-valued cells
            delimited_fields = self.__parsetype['delimited floats'] + self.__parsetype['delimited strings'];
            
            # get number of values in the largest cell of an entry/row
            length_of_cells = [];
            value_of_cells = [];
            for col in delimited_fields: # self.content[row].keys():
                cell = self.content[row][col];
                if isinstance(cell, list):
                    length = len(cell);
                    value = cell;
                else:
                    length = 1;  
                    value = [cell];
                length_of_cells.append(length);
                value_of_cells.append(value);
            num_values = max(length_of_cells)
            
            # repeat the values of single-valued cells to match the length of the multivalued cells
            for n in range(0,len(delimited_fields)):
                if length_of_cells[n] != num_values:
                    if length_of_cells[n] == 1:
                        col = delimited_fields[n];
                        self.content[row][col] = value_of_cells[n]*num_values;
                    else:
                        raise Exception('Reading PlateSheet Failed: Number of values in a multi-valued entry are in disagreement, please review the sheet (column: '+col+', row: '+str(row+self.__python_to_excel_row_offset)+')');

    def collapse_multivalued_cells(self):
        # convert all multi-valued cells which consist of only a single repeated value into a single-valued cell
        for row in range(0,len(self)):   
            
            # get list of fields that can have multi-valued cells
            delimited_fields = self.__parsetype['delimited floats'] + self.__parsetype['delimited strings'];
                        
            # collapse cells with only a single repeated value
            for col in delimited_fields:
                cell = self.content[row][col];
                if isinstance(cell, list):
                    length = len(cell);
                    value = cell;
                else:
                    length = 1;  
                    value = [cell];
                if length>1 and len(set(value)) == 1:
                    self.content[row][col] = [value[0]];                
        
    def save(self,filename):
        # save class data to csv
        if not os.path.isfile(filename):
            with open(filename, 'w') as fileID:    
                
                # copy current PlateSheet instance to modifiable variable
                sheet = copy.deepcopy(self); 
                
                # update date
                sheet.properties['Date'] = (datetime.datetime.now()).strftime("%Y/%m/%d %H:%M");                                                         
                
                # convert multi-valued cells made of only repeated values into single-valued cells (for conciseness) 
                sheet.collapse_multivalued_cells();
                
                # write format
                form = ['PlateSheet v' + sheet.format['Version'],''];
                writer1 = csv.writer(fileID);              
                writer1.writerow(form)                                              
                
                # write blank
                blank = ['',''];
                writer0 = csv.writer(fileID);              
                writer0.writerow(blank)     
                
                # write properties
                writer2 = csv.DictWriter(fileID, fieldnames=list(sheet.properties.keys()))
                writer2.writeheader()
                for data in self.unparse_cells([sheet.properties]):
                    writer2.writerow(data)                              
                
                # write blank
                writer0.writerow(blank)                              
                
                # write content
                writer3 = csv.DictWriter(fileID, fieldnames=list(sheet.content[0].keys()))
                writer3.writeheader()
                for data in self.unparse_cells(sheet.content):
                    writer3.writerow(data)    
                    
        else:
            raise Exception('Saving PlateSheet Failed: File already exists.');        
                    

    
    