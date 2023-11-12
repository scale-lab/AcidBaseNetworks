#The Plate Sheet Format

*Chris Arcadia* 
2018/05/31

## What is a Plate Sheet?

A **plate sheet** file (suffixed "_platesheet.csv") is a comma-seperated file that unambiguously defines chemical reagents and the desired positions & quantities of those reagents in a well plate. 

A plate sheet is divided into three parts:

1. statement of the used plate **format**
2. plate **properties** and related information
3. definition of plate **content** 

It is important to note that a plate sheet will only properly be read if all of these sections are present and seperated by a blank row. In addition, the header strings/labels for both the properties and content sections must be written *exactly* as specified. For consistency, the fields (headers and their in-column values) should be kept in the shown order, but this is not required for proper read in.

### Format

This single cell entry states the format used in the document, including the version of the format. It is the first line of every plate sheet file. This section will look like:

| PlateSheet v1.0 |
| ---------------------- |
|                             |


### Properties

This section is dedicated to storing global properties for the plate. This two row block can also be used to store any additional plate parameter (its label and value) you would like to include in the file (these will not be used by the program). Currently, there are several required parameters (these are shown in the below table). This section will look like:

| Plate | Barcode | Title   | Author        | Date   | Description           |
| ----- | ------- | ------- | ------------- | ------ | --------------------- |
| 384   | 0       | Example | Chris Arcadia | 6/1/15 | an example platesheet |



### Content

The data specifying the contents of a plate will be stored here. This section will look like this:

| Type           | Name                             | ID [PubChem] | Concentration [M] | Volume [L] | Positions      |
| --------------- | -------------------------------- | ------------ | ----------------- | ---------- | -------------- |
| amine           | paramethoxybenzlamine            | 75452        | 4                 | 1.00E-05   | A1:J17         |
| aldehyde        | 2-nitrobenzaldehyde              | 11101        | 4                 | 1.00E-05   | A1:E17         |
| carboxylic acid | Boc-O-benzyl L-beta-homotyrosine | 2761555      | 4                 | 1.00E-05   | A1:J1          |
| isocyanide      | methyl isocyanoacetate           | 547815       | 4                 | 1.00E-05   | A1:A17; F1:F17 |



Each column of the plate sheet is labelled at its top with parameter name to indicate its contents. This list of parameters (the header) is required to be present for proper loading of a plate sheet. 

#### Type

The type of a chemical. For instance, the type of an Ugi product precursor will either be "aldehyde", "amine", "isocyanide", or "carboxylic acid". For a solvents, such as water or dimethyl sulfide, the type will be "solvent".

#### Name

The name of the chemical. This can be the common name.

#### PubChem CID

The unique identifier of the chemical in the [PubChem](https://pubchem.ncbi.nlm.nih.gov/) database. This is used to automatically load information about the chemical such as molecular mass, CAS number, etc. Each row should specify a distinct chemical with and its unique ID (row entries containing the same chemical are not allowed).  

#### Positions

A desired position may be specified as a single well (such as "B3") or as a rectangular region (bounding box) of wells (such as "B3:B3", "C1:C17", or "B1:O5"). Multiple single wells or region of wells may specified by using a semicolon-space ("; ") delimiter (such as "B3:B8; C1; H1:I2"). This can be used when noncontiguous regions of reagent are needed.

To be clear, here are some examples of how to specify a bounding box region:

|           | Start Row | Start Column | Symbol | End Row | End Column | Wells                              |
| --------- | --------- | ------------ | ------ | ------- | ---------- | ---------------------------------- |
| Example 1 | A         | 1            | :      | A       | 7          | A1; A2; A3; A4; A5; A6; A7;        |
| Example 2 | A         | 1            | :      | C       | 1          | A1; B1; C1;                        |
| Example 3 | A         | 1            | :      | C       | 3          | A1; A2; A3; B1; B2; B3; C1; C2; C3 |



#### Volume

The volume to move per transfer into the specified well positions. This must be specified in liters [L]. Volumes are assigned for transferred on a per region basis (one concentration per region deliminiter ";"), so if you wanted to transfer 10μL to "A1:C1; A2:A9" then you should specify the volume as "10E-6; 10E-6". If you only write a single value (such as "10E-6") as the volume for an entry with multiple regions, it will be assumed that you want to use that volume for all  of the regions (as in "10E-6; 10E-6" will used).

#### Concentration

The concentration of the chemical in the source container. This must be specified as a molar concentration [M = mol/L]. Just as for the transfer volumes, the source concentration is specified on a per region basis.