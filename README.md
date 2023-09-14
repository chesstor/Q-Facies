# Qfacies - A tool for the quantitative interpretation of spatial and temporal distribution of hydrochemical facies

## Overview
Q-Facies is a method for quantitative spatiotemporal analysis of hydrochemical facies based on their spatial distribution in the panels of the Piper diagram

## How to install

The program could be easily **installed** in a **Windows** system by:

* **Download** and unzip the repositoy [*zip file*](https://github.com/chesstor/Q-Facies/archive/refs/heads/main.zip) in your working directory.

or
* **Cloning** the repository using git*:
    `git clone https://github.com/chesstor/Q-Facies.git`


Q-Facies have been programmed using **Python `3.7`**, and depends mainly on common Python packages.

For installing, running and managing Python packages, we recommend using [Anaconda](https://docs.anaconda.com/free/anaconda/install/windows/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html#windows-installers) distribution. The latter is a lighter release of Anaconda. 

For avoiding conflicts among the dependencies and Python versions while running Q-Facies, it's recommended to create a **new conda environment** with the same Python and packages versions in which it was programmed. You can do this by:

* **Option 1**. [*Manually*](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) from the Anaconda Prompt by running:

    `conda create -n qfacies python=3.7 numpy=1.19.2 pandas=1.3.0 matplotlib=3.3.2 openpyxl git`

	`pip install scikit-learn==0.23.2`     # *This line after activating the env*

    The library `openpyxl` is also included since it's needed for Excel creation and not always included in `pandas` installation. The use of *pip* is highly recommended for sckikit-learn installaion for avoiding a bad installation of the `scipy` dependency (see regarding error [here](https://stackoverflow.com/questions/39020361/python-scipy-module-import-error-due-to-missing-ufuncs-dll)).

or

* **Option 2**. From an [*yml file*](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) using the *qfacies.yml* file included in the repository. This option will create an environment with all the packages and its versions included in the file. The first line of the *yml* file sets the new environment's name. Be free to change it for other env name.  Simply run the following line on the Anaconda Prompt from within the folder that holds the *yml* file:

    `conda env create -f qfacies.yml`

## Dependecies:
The program relies on Python standard libraries.

_ `Numpy 1.19.2`

_ `Pandas 1.3.0`

_ `Matplotlib 3.3.2`

_ `Scikit-learn 0.23.2`

## How does it work?
To make a Q-Facies diagram just follow these steps:
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
All authors belong to the Geological Survey of Spain: Instituto Geológico y Minero de España (IGME), CSIC, Ríos Rosas, 23, 28003 Madrid, Spain

* M. González-Jiménez         miguel.gj@ipna.csic.es

* L. Moreno Merino            l.moreno@igme.es

* H. Aguilera Alonso          h.aguilera@igme.es

* A. Romero                   a.romero@igme.es

## Paper
If you use the program in one of your studies, please cite this paper:
* González-Jiménez, M. G., Aguilera, H., Merino, L. M., & Prados, A. R. (2023). Q-Facies: A tool for the quantitative interpretation of groundwater hydrochemical facies. SoftwareX, 101450.

    https://www.softxjournal.com/article/S2352-7110(23)00146-2/fulltext

## Copyright
License: GNU General Public License v3.0 
