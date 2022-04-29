# -*- coding: utf-8 -*-
"""
===============================================================================
======================== Module of Q-Facies package ============================
===============================================================================
Contains all the classes related to the elements of the diagrams.
A diagram is always conformed by three panels: Cation, Anion and Diamond (each 
one subclassed from Panel class and associated at Diagram class)

All euclidean transformations requiered for each panel are imported from the
'calculation.py' module.

All plotting methods for visual representation are imported from the 'plot.py' module.

===============================================================================
===============================================================================
"""

import pandas as pd
from plot import Plot_Diagram
from calculation import Indixes, Transform

class Diagram:
    ''' Contains all the information of a classical Piper Diagram'''
    
    def __init__(self, df, way, df_name, pol_color='orange', **kw):
        ''' Create an object for the Cation, Anion and Diamond panels, and store
        them in a tuple.'''
        self.df, self.way, self.kw  = df, way, kw
        self.df_name, self.pol_color = df_name, pol_color
        self.df.dropna(inplace=True)
        #self.df.fillna(0, inplace=True) # same results as dropna
        # Classes agregation:
        self.cation = Cation(df, **kw)
        self.anion = Anion(df, **kw)
        self.diamond = Diamond(df, **kw)
        self.panels  = (self.cation, self.anion, self.diamond)

    def __str__(self):
        ''' Set Group name'''
        if self.way == 'by_groups':
            return self.df['Group_ID'].iloc[0]
        elif self.way == 'by_time':
            return (self.df['Time'].min(), self.df['Time'].max())
        
    def plot(self):
        ''' Plot the Piper diagram with all the information.
            Executed via the 'plot.py' module'''
        Plot_Diagram(self.panels, name=self.__str__(),
                 way=self.way, df_name=self.df_name,
                 pol_color=self.pol_color, **self.kw)
        
    def get_params(self):
        '''Return a list of dicts with all panels' indices.'''
        d_params = [i.Params.get() for i in self.panels]
        name = '{} to {}'.format(*[i.strftime(format='%d-%m-%Y') \
                       for i in self.__str__()]) if              \
                       self.way == 'by_time' else self.__str__()
        
        [i.update(dict(Group=name)) for i in d_params]
        return d_params
    
    def get_all_points(self):
        ''' Get all points of the Piper diagram. Columns: xy coordinates and 'ID'
        group or Date column.
        Return a pandas.DataFrame'''
        column = dict(by_time='Time', by_groups='Group_ID')
        return pd.concat([i.get_column(column[self.way]) for i in
                          self.panels], axis=0, ignore_index=True)

class Panel:

    def __init__(self, df):
        self.df = df

    def get_column(self, column):
        ''' Get points of this panel with their XY coordinates and ANOTHER
        column (introduced by the 'keep_column' arg). Created to represent points
        according to a third variable inside the Diagram class'''        
        a = pd.DataFrame(self.points, columns=['X','Y'])
        b = self.df.copy().reset_index(drop=True)[[column]]
        return pd.concat([a,b], axis=1)
    
class Cation(Panel):
    '''Contains all cation panel data and methods'''
     
    def __init__(self, df, params=True, **kw):
        '''params: boolean Default True. Calculate all the indices with the
        Indixes class'''
        super().__init__(df)
        self.panel = 'cation'
        self.components = ['NaK_epm','Mg_epm']
        self.points = self.transform(df.filter(items=self.components).to_numpy())
        self.Params = Indixes(self.points, df, self.panel, **kw)

    def transform(self, points):
        ''' Make all the required euclidean transformation to the panel's
        points via matrixes-product (@)'''
        T = Transform()
        return points @ T.scale() @ T.t_shear()
    
    def info(self):
        ''' Contains the info to represent on the diagram'''
        x,y = 0.72, 0.81
        text = '\n'.join([
            "{:^34}".format(r'$\bf{}$ ({})'.format('Cation', self.Params.num_points)),
            "{:<}: {:>17.1f}".format('Area (%)',      self.Params.area),
            "{:<}: {:>8.1f}".format('Shape Idx. (%)', self.Params.shape_idx),
            "{:<}: {:>8.2f}".format('Orientation ({})'.format('$^{o}$'),
                                                      self.Params.angle),
            "{:<}: {:>17.1f}".format('Blau (%)',      self.Params.blau_idx[0]),
            "{:<}: {:>20.2f}".format('SD (%)',      self.Params.sd)
                        ])
        return x,y,text

class Anion(Panel):
    '''Contains all anion panel data and methods'''
    
    def __init__(self, df, params=True, **kw):
        '''params: boolean Default True. Calculate all the indices with the
        Indixes class'''
        super().__init__(df)
        self.panel = 'anion'
        self.components = ['Cl_epm','SO4_epm']
        self.points = self.transform(df.filter(items=self.components).to_numpy())
        self.Params = Indixes(self.points, df, self.panel, **kw)

    def transform(self, points):
        ''' Make all the required euclidean transformation to the panel's
        points'''
        T = Transform()
        return (points @ T.scale() @ T.t_shear()) + T.a_translation()

    def info(self):
        ''' Contains the info to represent on the diagram'''
        x,y = 0.72, 0.58
        text = '\n'.join([
            "{:^36}".format(r'$\bf{}$ ({})'.format('Anion', self.Params.num_points)),
            "{:<}: {:>17.1f}".format('Area (%)',      self.Params.area),
            "{:<}: {:>8.1f}".format('Shape Idx. (%)', self.Params.shape_idx),
            "{:<}: {:>8.2f}".format('Orientation ({})'.format('$^{o}$'),
                                                      self.Params.angle),
            "{:<}: {:>17.1f}".format('Blau (%)',      self.Params.blau_idx[0]),
            "{:<}: {:>20.2f}".format('SD (%)',      self.Params.sd)
                        ])
        return x,y,text

class Diamond(Panel):
    '''Contains all diamond panel data and methods'''
    def __init__(self, df, params=True, **kw):
        '''params: boolean Default True. Calculate all the indices with the
        Indixes class'''
        super().__init__(df)
        self.panel = 'diamond'
        self.components = ['NaK_epm','HCO3CO3_epm']
        points = df.filter(items=self.components)
        points['HCO3CO3_epm']= points['HCO3CO3_epm'].map(lambda x: 100 - x)
        self.points = self.transform(points.to_numpy())
        self.Params = Indixes(self.points, df, self.panel, **kw) if params else None
        
    def transform(self, points):
        ''' Make all the required euclidean transformation to the panel's
        points'''
        T = Transform()
        p_trans = (
        points @ T.scale() @ T.d_shear() @ T.rotation() ) + T.d_translation()
        return p_trans
   
    def info(self):
        ''' Contains the info to represent on the diagram'''
        x,y = 0.07, 0.68
        text = '\n'.join([
            "{:^34}".format(r'$\bf{}$ ({})'.format('Diamond', self.Params.num_points)),
            "{:<}: {:>17.1f}".format('Area (%)',      self.Params.area),
            "{:<}: {:>8.1f}".format('Shape Idx. (%)', self.Params.shape_idx),
            "{:<}: {:>8.2f}".format('Orientation ({})'.format('$^{o}$'),
                                                      self.Params.angle),
            "{:<}: {:>17.1f}".format('Blau (%)',      self.Params.blau_idx[0]),
            "{:<}: {:>20.2f}".format('SD (%)',      self.Params.sd)
                        ])
        return x,y,text
