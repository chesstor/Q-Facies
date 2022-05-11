# Qfacies - A tool for the quantitative interpretation of spatial and temporal distribution of hydrochemical facies

## Overview
Q-Facies is a method for quantitative spatiotemporal analysis of hydrochemical facies based on their spatial distribution in the panels of the Piper diagram

## How does it work?
To make a D-Piper diagram just follow these steps:
### 1. Prepare the file structure:
Keep the file structure as shown in the figure:

![image](https://user-images.githubusercontent.com/12763571/167836289-92ec2385-87f6-425b-8d7d-dbb9c7f3ebf9.png)


 * The input analytical data files must be saved to the *Data* folder. The sample files used in the article to make Figures 2, 3, 4 and 5 can be found in this folder.
 * The Q-Facies plots will be saved to the *Graphics* folder. It currently contains some example images from the paper and the How_to_Q_Facies.pdf reference guide file.
 * The file *main.py* is the Python script to be run to generate the Q-Facies tables and graphics.
 * The file *calculation.py* contains two classes that allow to apply all the Euclidean transformations required and to calculate all of the indices for each panel, as well as to identify outliers.
 * The file *diagram.py* contains all the classes related to the elements of the diagrams. A diagram is always conformed by three panels: Cation, Anion and Diamond, all of them subclassed from the Panel class and bound together with the Diagram class.
 * The file *plot.py* contains all the graphic methods of Q-Facies. 
 * The file *Options.txt*, necessary to run the *main.py* script, contains the contains the different options for executing the program and commented help text with the description of variables and parameters.
 
Do not delete or change any of the above files or folders. If necessary, you can add other folders or files. Note that automatically a *pycache* folder is created, requiered for an optimal execution of the program.
 
The script always reads the options from the *Options.txt* file. If you want to save a set of options, you must rename the file or save it to another location. You will have to name it as the original one (*Options.txt*) or update the new name in line 51 of the *main.py* file for further use.

### 2. Preparing the data file:

The data file must be in ASCII format and saved to the *Data* folder as .txt or .csv.
The data file can have any name.
The structure of the data file should be as follows:

* TAB-separated ASCII file
* It must contain 9 columns:
      
![image](https://user-images.githubusercontent.com/12763571/167825027-c10d289b-b607-4b97-a126-b45ec71b560f.png)

*Group/Date*: Each sample can belong to a different group (e.g., geographical area), or a time series of a single sampling point (e.g., well) can be provided for temporal analysis.<br>
*Eight columns with analytical results*. The ion content should be expressed in mg/L (milligrams per liter) or ppm (parts per million), or directly in percentage milliequivalents (%epm).<br>

### 3. Preparing the _Options_ file:

The *Options.txt* file must be in ASCII format.
It is a self-explanatory file that shows all the options of the program.
It is essential to comply with the file format to avoid execution errors.

The file How_to_Q_Facies.pdf describes, step by step, how to use the Q-Facies program and the meaning of each of the variables and parameters contained in the Options.txt file.

## Authors
All authors belong to IGME Geological Survey of Spain. C/Ríos Rosas 23, 28003 Madrid, Spain

* M. González-Jiménez         miguigonn@gmail.com

* L. Moreno Merino            l.moreno@igme.es

* H. Aguilera Alonso          h.aguilera@igme.es

* A. Romero                   a.romero@igme.es


## Copyright
License: GPLv3 
