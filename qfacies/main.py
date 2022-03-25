# -*- coding: utf-8 -*-
'''
                               IGME 
                     Geological Survey of Spain

                             Q-Facies
                            Version 1.0
                            01-10-2021
  Authors:
	M. González Jiménez   miguel.gonzalez@igme.es
	L. Moreno             l.moreno@igme.es	
	H. Aguilera           h.aguilera@igme.es
	A. De la Losa         a.delalosa@igme.es
	A. Romero             a.romero@igme.es

If you have some problem with the code or want to propose some modification,
please contact with M. González Jiménez.

===============================================================================
===============================================================================


 ------------------------------------------------------------------------------
 --------------------- DATA ENTRY FILE FORMAT ---------------------------------
 ------------------------------------------------------------------------------
Position:           1         2  3   4   5    6     7    8   9
Variable:  ID_group or Date  Ca  Mg  Na  K   HCO3  CO3  SO4  Cl

ASCII file separated by tabs   (\t)
The decimal separator is point (.)
===============================================================================
===============================================================================
'''

__version__ = "1.0"

import numpy as np
import pandas as pd
from matplotlib import cm

import time, warnings

from plot import Overlap, Plot_all, Evolution_time, Evolution_groups, Skeleton
import plot as qplot
from diagram import Diagram
import calculous as cal

###############################################################################
###################---------- Variables defining ---------#####################
###############################################################################

with open('Options.txt', 'r') as file:
    data = file.readlines()
    text = [i.strip() for i in data if i[0] != ';']
    text = [i for i in text if i != '']
    text = [i for i in text if '=' in i]
    plain_text = ','.join([i for i in text])
    
    try:
        exec(f'kw = dict({plain_text})')
        
    except SyntaxError:
        print("Something is wrong in the Options.txt file \n \
               Check out that non-variables lines  start by ';' \n \
               and the following on the variable lines: \n \
                _Variables should be define by '=' \n \
                _The name of the variable cannot contain spaces. \n \
                _When assigning a string to a variable, this must be \n \
                    enclosed in quotation marks")
        raise
        
    fname = kw['fname']
    way = kw['way']
    [kw.pop(i) for i in ('fname','way')]
              
# =============================================================================
# ================================ CORE =======================================
# =============================================================================

cal.import_skl() if kw['lof'] else None
qplot._extension_graph = kw['extension_graph']

class Dataset:
    ''' Main class of the program'''
    
    def __init__(self, fname, way='by_groups', **kw):
        '''Read the dataset from 'Data' folder
        way:str {'by_groups', 'by_time'} default 'by_groups
            Way of analyze hydrochemical facies evolution.
            'by_groups' -> study the facies of defined groups
            'by_time'   -> study the facies of the temporal series throughout
                           a rolling window. In this case, window params need
                           to be set.'''
        file = ReadDataFile(fname, way, **kw)
        self.fname = fname.split('.')[0]
        self.data = file.get()
        self.way = way
        self.main(self.way, **kw)
    
    def main(self, way, **kw):
        '''Program execution'''
        _g = self.by_groups() if way == 'by_groups' else self.by_time(
                       window_size=kw['window_size'], freq=kw['freq'])
        # Clean the groups of NaN rows.
        self.groups = self._group_filter(_g)
        
        polygon_colors = self._color_diagram(folder_color=kw['folder_color'])
        kw.pop('pol_color') # Delete so as not to overwrite kwargs.
        self.diagrams = [Diagram(i, self.way, self.fname, pol_color=j, **kw)
                         for i,j in zip(self.groups, polygon_colors)]
        
        self.excel() if kw['excel'] else None
        self.plotting(**kw)
    
    def _color_diagram(self, folder_color='Spectral', reverse=False):
        ''' Return a discretized list of a determined colormap with as many 
        elements as groups has to analyze'''
        try:
            folder_color = cm.get_cmap(folder_color)
            folder_color = folder_color.reversed() if not reversed else folder_color
            lsp = np.linspace(0,folder_color.N, len(self.groups)).astype(int)
            return [folder_color(i) for i in lsp]
        
        except:
            return [folder_color for i in range(len(self.groups))]
            
    def _group_filter(self, g):
        '''Clean Nan rows and empty groups (all nan values) and set a threshold:
            groups of less than three points are not valid.'''
        # Drop NaN rows:
        non_NaN_groups = [i.dropna(axis=0, thresh=8) for i in g]
        
        # If there are empty groups: delete them and raise a Warning
        idx = [i.empty for i in non_NaN_groups] # Bool type list of df to remove
        groups = [i for i,j in zip(non_NaN_groups,idx) if j==False]
        del_groups = len(non_NaN_groups) - len(groups)
        
        if del_groups>0:
            warnings.warn('{} empty group(s) has been deleted'.format(del_groups))
        
        # Raise an exception when groups are formed by less than three points:
        fails = len([i for i in groups if i.shape[0]<3])
        assert fails==0, "{a} group(s) have less than three points (minimun \
            requiered). \nTry to define a wider window size (if creating groups\
            by time) or delete those groups conformed by less than 3 analyses \
            (if creating groups by ID_Group).".format(a=fails)
        
        return groups

    def by_groups(self):
        ''' Extract all groups contained in the dataset and calculate the
        facies parameters for one of them. The groups are made from ID-group
        fields differences.'''
        return [self.data.loc[self.data.Group_ID == i] for i in 
                self.data.Group_ID.unique()]
        
    def by_time(self, window_size=10, freq=5):
        ''' Extract all groups contained in the dataset and calculate the
        facies parameters for one of them. Minimum: every three points
        '''
        index = self.data.index.to_numpy()
        assert window_size >3, "A group must be conformed by three or more points."
        assert freq<window_size, 'Window size is wider than slicing frequency.\
             There will be analyses not taken into account'
             
        def rolling_window(array, window_size, freq):
            shape   = (array.shape[0] - window_size + 1, window_size)
            strides = (array.strides[0],) + array.strides
            rolled  = np.lib.stride_tricks.as_strided(array, shape=shape, strides=strides)
            step_rolled = rolled[np.arange(0,shape[0],freq)]
            start, end = step_rolled[:,0], step_rolled[:,-1]
            # Mask array for adding 1 to all end-array elements
            mask = np.ones_like(end).astype(bool)
            mask[-1] = False
            return start, np.where(mask, end+1, end)

        start, end = rolling_window(index, window_size, freq)
        return [self.data.iloc[i:j,:] for i,j in zip(start, end)]
   
    def _info_table(self):
        '''Insert all diagrams' parameters into a DataFrame.'''
        return pd.DataFrame([j for i in self.diagrams for j in i.get_params()])
    
    def _mapa(self):
        ''' Dictionary to map columns names '''
        return {'Area':'Area (%)','Shape':'Shape index (%)',
                'Angle':'Orientation ({})'.format('$^{o}$'),
                'panel':'Panel','Blau':'Blau index (%)',
                'points':'Points','A':'Blau (A)',
                'B':'Blau (B)','C':'Blau (C)', 'D':'Blau (D)',
                'Time':'Time', 'Dispersion':'Standard Distance'}
         
    def excel(self):
        ''' Create an Excel file with all the groups' parameters info and 
        save it at 'Data' folder'''
        df = self._info_table()
        df.rename(mapper=self._mapa(), axis=1, inplace=True)
        df.rename(mapper={'Orientation ($^{o}$)':'Orientation ({})'.format(u"\u00b0")},
                  axis=1, inplace=True)
        df.Panel = df.Panel.map(lambda x: x.capitalize()) # Format the fields
        
        # Create two dataframes
        df1 = df.set_index(['Group','Panel'])
        df1.drop('Dominant',axis=1, inplace=True)
        
        a = df.columns.to_list()
        [a.remove(i) for i in ('Group','Panel')]
        df2 = df.set_index('Group').pivot(columns='Panel',values=a).swaplevel(0,1,
                                    axis=1).sort_index(axis=1)
        
        # Saving the dataframes as two different sheets of the same Excel file
        book = pd.ExcelWriter('Data\{}.xlsx'.format(self.fname))
        df1.to_excel(book, '{}_A'.format(self.fname), float_format="%.2f")
        df2.to_excel(book, '{}_B'.format(self.fname), float_format="%.2f")
        book.save()
         
    def plotting(self, **kw):
        '''Different representation function via 'Plot' module'''
        # Create a Figure for each group if 'Figs' is set to True
        if kw['figs']:
            try:
                import tqdm
                [i.plot() for i in tqdm.tqdm(self.diagrams)]
            except:
                [i.plot() for i in self.diagrams]
        
        # Create an Evolution figure:
        if kw['evolution_fig']:
            if self.way == 'by_groups':
                Evolution_groups(self.diagrams,self._info_table(),
                            self.fname,self.way,self._mapa())
            elif self.way == 'by_time':
                Evolution_time(self.diagrams,self._info_table(),
                            self.fname, self.way, self._mapa())
        # Create an Overlap figure:
        if kw['overlap_fig']:
            Overlap(self.diagrams, self.way, self.fname, **kw) #['overlap_var']
            
        # Create a unique figure with ALL the points:
        if kw['plot_all']:
            if self.way == 'by_time':
                df_all = Diagram(self.data, 'by_time',
                                     self.fname).get_all_points()
                Plot_all('by_time', self.fname, df=df_all, **kw)

            elif self.way == 'by_groups':
                Plot_all('by_groups', self.fname, diagrams=self.diagrams, **kw)

class ReadDataFile:
    '''Read the input data file, check for some decimal errors and convert
    messures from ppm o mg/L values to percentage (epm) ones.'''
    
    def __init__(self, fname, way, space=False, freq=None,
                 transform=True, **kw):
        '''fname: str
                  Filepath of the input data.
           way: str {'by_groups','by_time'}
           space: bool, default False. Optional
                 Only takes place when way set to 'by_time'
           freq: str {following pandas frequency aliases} default None
                 Fore more info: https://pandas.pydata.org/pandas-docs/stable/
                 user_guide/timeseries.html#offset-aliases
           transform: bool, default True
                     Transform ionic data from mg/L (or ppm) to % meq/L.'''

        self.ions=['Ca','Mg','Na','K','HCO3','CO3','SO4','Cl']
        self.cols = dict(by_groups=['Group_ID',*self.ions], by_time=['Time',*self.ions])

        self.way = way
        self.df = pd.read_csv('Data/{}'.format(fname), sep='\t')
        self.df.columns = self.cols[way]
        self.formatt()
        self._converse() if transform == True else self.df
        self.time(space, freq, agregation=kw['agregation'],
                          format=kw['time_format']) if way=='by_time' else None
        
    def get(self):
        '''Return formated input file as a pandas DataFrame.'''
        return self.df
    
    def formatt(self):
        '''Format columns depending on 'way' parameter '''
        head = self.cols[self.way][0]
        
        # Replace '<' charachters for 0, since are considered bellow detection
        # limit.
        self.df[self.ions] = self.df[self.ions].applymap(lambda x:\
             (np.nan if '<' in x else x) if (type(x) == str) else x)
        
        # Points are decimal separator and float values
        self.df[self.ions] = self.df[self.ions].applymap(lambda x:\
             x.replace(',','.') if type(x) == str else x).astype(float)
        

        self.df.fillna(0, inplace=True)
            
        # Treat fields as strings or dates depending on 'way' parameter
        # Time format will be added in via 'time' method.
        if self.way=='by_groups':
            self.df[head] = self.df[head].map(lambda x: str(x))
              
    def time(self, space, freq, agregation='mean', format='%d/%m/%Y'):
        '''Turn time column to pandas datetime format and resample if desired'''
        try:         
            self.df['Time'] = pd.to_datetime(self.df['Time'],
                              format=format, errors='raise')
            self.df.sort_values('Time',inplace=True)
            
        except TypeError:
            raise("There are troubles with date format. Try to format them as \
                   follows: 'dd-MM-YYYY' or 'dd-MM-YYYY 'hours:minutes:seconds' \
                   for higher resolution.")
                   
        # Resample to specifed frequency if desired to be evenly spaced.
        self.resample(freq, agregation) if space else None
        
    def resample(self, freq, agregation):
        '''Resample via pd.resample() '''
        _trans = {'yearly':'Y','monthly':'M','weekly':'W','daily':'D'}
        freq = _trans[freq] if freq in _trans.values() else freq
        self.df = self.df.set_index('Time').resample(freq, origin='start') \
                          .agg(agregation).reset_index()
    
    ###########################################################################
    # ----------- Percentage equivalents (epm) calculation --------------------
    ###########################################################################
    def _converse(self):
        '''Convert mg/L or ppm analyses to epm values (ie, % values)
        '''
        # equivalent number (eq) = substance weight () * (valence (V) / molecular mass(Mm))
        # Conversion Factor (CF) = V/Mm ----> One per molecule
        
        Mm = {'Ca': 40.078, 'Mg': 24.305, 'Na':22.989769, 'K': 39.0983,
              'HCO3': 61.01684, 'CO3': 60.0089, 'SO4': 96.0626, 'Cl':35.453}
        
        V = {'Ca': 2, 'Mg': 2, 'Na':1, 'K': 1, 'HCO3': 1, 'CO3': 2, 'SO4': 2, 'Cl':1}
        
        cf = {i: V[i]/Mm[i] for i in V}
        
        # --------------------------------------------------------------------
        # --------- Convert mass units into equivalent units (eq) ------------
        # --------------------------------------------------------------------
        # CATIONS
        self.df['Ca_eq'] = (self.df.Ca * cf['Ca']).round(3)
        self.df['Mg_eq'] = (self.df.Mg * cf['Mg']).round(3)
        self.df['Na_eq'] = (self.df.Na * cf['Na']).round(3)
        self.df['K_eq']  = (self.df.K  * cf['K'] ).round(3)
        
        # ANIONS:
        self.df['HCO3_eq'] = (self.df.HCO3 * cf['HCO3']).round(3)
        self.df['CO3_eq']  = (self.df.CO3  * cf['CO3'] ).round(3)
        self.df['SO4_eq']  = (self.df.SO4  * cf['SO4'] ).round(3)
        self.df['Cl_eq']   = (self.df.Cl   * cf['Cl']  ).round(3)
        
        self.df.drop(self.ions, axis=1, inplace=True) # Deleting useless fields
        elements = ['{}_eq'.format(i) for i in self.ions]

        # --------------------------------------------------------------------
        # --------- Convert equivalents into percentage units (epm) ----------
        # --------------------------------------------------------------------
        
        # Summatory:     
        self.df['CAT_sum'] = self.df[['Ca_eq','Mg_eq','Na_eq','K_eq']].sum(axis=1)
        self.df['ANI_sum'] = self.df[['HCO3_eq','CO3_eq','SO4_eq','Cl_eq']].sum(axis=1)
        
        # CATIONS
        self.df['Ca_epm']  = ((self.df.Ca_eq/self.df.CAT_sum)*100).round(3)
        self.df['Mg_epm']  = ((self.df.Mg_eq/self.df.CAT_sum)*100).round(3)
        self.df['NaK_epm'] = (((self.df.Na_eq + self.df.K_eq)/self.df.CAT_sum)*100).round(3)
        
        # ANIONS:
        self.df['HCO3CO3_epm'] = (((self.df.HCO3_eq + self.df.CO3_eq)/self.df.ANI_sum)*100).round(3)
        self.df['SO4_epm']      = ((self.df.SO4_eq/self.df.ANI_sum)*100).round(3)
        self.df['Cl_epm']       = ((self.df.Cl_eq/self.df.ANI_sum)*100).round(3)
        
        self.df.drop(elements, axis=1, inplace=True) # Deleting useless fields
        
# =============================================================================
# ========================= CUSTOM WARNING =====================================
# =============================================================================

def format_warning(msg, *args, **kwargs):
    return f'\nWarning: {msg}\n' # ignore everything except the message

warnings.formatwarning = format_warning

warnings.simplefilter('ignore') if kw['ignore_warnings'] else None

# =============================================================================
# ============================== MAIN =========================================
# =============================================================================

def main():
    start = time.time()
    D = Dataset(fname, way=way, **kw)
    end = time.time()
    d = end - start
    print(
    f"Execution time: {round(d//60)} minutes and {round(d%60,2)} seconds.")
    return D
    
if __name__ == '__main__':
    D = main()
    
# =============================================================================
# ======================== END OF THE PROGRAM =================================
# =============================================================================