# -*- coding: utf-8 -*-
"""
===============================================================================
======================== Module of Qfacies package ============================
===============================================================================
Contains all the visual representation methods of 
Qfacies.

Classes that need the Piper diagram for plotting the data are subclassed from
Skeleton (the Piper diagram skeleton class). Thesea are Overlap, Plot_Diagram
and PLot_all

The rest of the classes (Evolution_time and Evolution_groups), that allow to
represent all the extracted information in a three-graphs plot, are subclassed
from Evolution class
===============================================================================
===============================================================================
"""
import pandas as pd
import numpy as np
import warnings

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap as lsc
from matplotlib.colors import Normalize
import matplotlib.cm as cm 

_extension_graph = 'png' 
_dpi = 400

class Skeleton:
    '''Create basic structure of Piper diagram.'''

    def __init__(self):
        global _extension_graph, _dpi
        self._extension_graph, self._dpi = _extension_graph, _dpi
        # --------------------------------------------------------------------
        self.fig, self.ax = plt.subplots()
        self.offset,self.l_offset, self.margin = 22, 10, 20
        self.sin, self.cos = np.sin(np.pi/3), np.cos(np.pi/3)
        self.tan = np.tan(np.pi/3)
        self.tri_params = {'facecolor':'None', 'lw':0.5, 'ec':'black', 'zorder':3}
        # --------------------------------------------------------------------
        [self.ax.add_patch(i) for i in (self._cation(), self._anion(), self._diamond())]
        self._grid()
        self._labels()
        self._skeleton_layout()

    def _cation(self):
        # ---------------- DRAW CATION EQUILATERAL TRIANGLE -------------------
        cat_TRI = np.array([( 0,   0),
                            (50, self.sin * 100),
                            (100,  0)])
        self.cation = plt.Polygon(cat_TRI, **self.tri_params)
        return self.cation
    
    def _anion(self):
        # ---------------- DRAW ANION EQUILATERAL TRIANGLE --------------------
        an_TRI = np.array([(100 + self.offset,  0),
                           (150 + self.offset, self.sin * 100),
                           (200 + self.offset,  0)])
         
        self.anion = plt.Polygon(an_TRI, **self.tri_params)
        return self.anion
    
    def _diamond(self):
        offset = self.offset
        # --------------- DRAW DIAMOND EQUILATERAL TRIANGLE -------------------
        # --------------- Conformed by fours point (H, J, I, G)----------------
        H = ((100 + offset) * 0.5        , self.sin * (100 + offset))   
        J = (((100 + offset) * 0.5 + 100), self.sin  * (100 + offset)) 
        I = ((200 + offset) * 0.5        , self.sin  * (200 + offset))
        G = ((200 + offset) * 0.5        , self.sin  * (200 + offset) - 
                                           2 * (self.sin  * 100))
        self.diamond = plt.Polygon(np.array([H,I,J,G]), **self.tri_params)
        return self.diamond
   
    def _grid(self):
        '''Reference lines within each panel. Created for every 20 units'''
        sin, cos, offset = self.sin, self.cos, self.offset
        #######################################################################
        #######################------- Visual Grid ----------##################
        #######################################################################
        gl_params = {'lw': 0.2, 'color':'black'}
        # Grid lines set first leftward dipping, then rightwards and then horizontal
        # --------------------- Cation Triangle -------------------------------
        [self.ax.plot([espaciado, (espaciado/2) + 50], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [self.ax.plot([100 - espaciado, -(espaciado/2) + 50], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [self.ax.plot([cos * espaciado, cos * espaciado + (100 - espaciado)], 
                  [sin * espaciado, sin * espaciado], **gl_params)
                  for espaciado in range(20,100,20)]
        
        # --------------------- Anion Triangle --------------------------------
        [self.ax.plot([espaciado + 100 + offset, (espaciado/2) + 50 + 100 + offset], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [self.ax.plot([100 - espaciado + 100 + offset, -(espaciado/2) + 50 + 100 + offset], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [self.ax.plot([cos * espaciado + 100 + offset, cos * espaciado + 200 - espaciado + offset], 
                  [sin * espaciado, sin * espaciado],**gl_params)
                  for espaciado in range(20,100,20)]
        
        # --------------------- Diamond Panel ---------------------------------
        [self.ax.plot([((100 + offset) + espaciado) * cos, (((100 + offset) + espaciado) + 100) * cos], 
                  [sin * (offset + 100 - espaciado),
                   sin * (offset + 100 - espaciado) + sin * 100], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [self.ax.plot([((100 + offset) + espaciado) * cos,(((100 + offset) + espaciado) + 100) * cos], 
                  [(sin * (100 + offset + espaciado)), sin * (offset + espaciado)],
                  **gl_params)
                  for espaciado in range(20,100,20)]
    
    def _labels(self):
        '''LABELLING '''
        sin, cos, tan, offset = self.sin, self.cos, self.tan, self.offset
        l_offset = self.l_offset
        # Labels set by line equations.
        
        text_params = {'va':'center', 'ha':'center', 'style':'italic', 'fontsize':5}
        arrow_params = {'linewidth':0.2, 'head_width':3.5, 'head_length':4,
                        'facecolor':'black', 'edgecolor':'black'}
        x_array = np.linspace(0,50,9)
        
        # r - line:         y = tan(60º)*x + l_offset 
        ######################## -------- Mg ------ ###########################
        x, y = x_array[2],  tan * x_array[2] + l_offset
        self.ax.arrow(x, y, cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = x_array[2:4].mean(), sin * 22 + y
        self.ax.text(x_txt, y_txt, '$Mg^{2+}$', rotation = 60, **text_params)
        ######################## -------- SO4 + Cl ------ #####################
        x = x + 50 + (cos * offset)
        y = y + sin * 100 + (sin * offset)
        self.ax.arrow(x, y, cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = x_txt + 50 + (cos * offset), y + sin * 22
        self.ax.text(x_txt, y_txt, '$SO_{4}^{2-}$ + $Cl^-$', rotation = 60, **text_params)
        # r' - line:       y = tan(60º) * (x - 100 - offset)  l_offset 
        ###################### -------- CO3 + HCO3 ------ #####################
        x = x_array[4:5].mean() + 100 + offset
        y = tan * (x - 100 - offset) + l_offset
        self.ax.arrow(x, y, -cos * 25, -sin * 25, **arrow_params)
        x_txt, y_txt = 98 + offset + x_array[3], y - 3
        self.ax.text(x_txt, y_txt, '$CO_{3}^{2-}$ + $HCO_{3}^{-}$', rotation = 60, **text_params)
        # p - line:        y = tan(60º) * (-x + 100 + l_offset)  
        ######################### -------- Na + K ------ ######################
        x = np.linspace(50,100,9)[4]
        y = tan * (-x + 100 + (l_offset/tan)) 
        self.ax.arrow(x, y, cos * 25, -sin * 25, **arrow_params)
        x_txt,y_txt = 5 + np.linspace(50,100,9)[5:6].mean(), y -  10
        self.ax.text(x_txt, y_txt, '$Na^{+} + K^{+}$', rotation = 300, **text_params)
        # p' - line:       y = tan(60º) * (-x + 200 + offset + l_offset/tan) 
        ####################### -------- Ca + Mg ------ #######################
        x = 100 + offset/2 + np.linspace(0,50,9)[6]
        y = tan * x_array[2] + l_offset + sin * 100 + (sin * offset)
        self.ax.arrow(x, y, -cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = 99 + offset/2 + np.linspace(0,50,9)[6], y + sin * 19
        self.ax.text(x_txt, y_txt, '$Ca^{2+} + Mg^{2+}$', rotation = 300, **text_params)
        ######################### -------- SO4 ------ #########################
        x = 100 + offset + np.linspace(50,100,9)[6]
        y = tan * x_array[2] + l_offset
        self.ax.arrow(x, y, -cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = 149 + offset + x_array[6], sin * 18 + y
        self.ax.text(x_txt, y_txt, '$SO_{4}^{2-}$', rotation = 300, **text_params)
        # h:                 y = 0 
        ######################### -------- Ca ------ ##########################
        x, y = np.linspace(0,100,9)[5], -(cos * l_offset)
        self.ax.arrow(x, y, -25, 0, **arrow_params)
        x_txt, y_txt =  np.linspace(0,100,9)[4], y - 8
        self.ax.text(x_txt, y_txt, '$Ca^{2+}$', **text_params)
        ######################### -------- Cl ------ ##########################
        x, y = 100 + offset + np.linspace(0,100,9)[3], - (cos * l_offset)
        self.ax.arrow(x, y, 25, 0, **arrow_params)
        x_txt, y_txt = 100 + offset + np.linspace(0,100,9)[4], y - 8
        self.ax.text(x_txt, y_txt, '$Cl^{-}$', **text_params)
    
    def _skeleton_layout(self):
        '''LAYOUT PROPERTIES'''
        margin, offset, sin = self.margin, self.offset, self.sin
        # Boundaries Box:
        box = np.array([(-margin, -margin),
                        (200 + offset + margin, -margin),
                        (200 + offset + margin, margin + sin * (offset + 200)),
                        (-margin, margin + sin * (offset + 200))])
        self.ax.add_patch(plt.Polygon(box, lw = 1, ec = 'black', fc = 'None'))
        self.ax.set_aspect('equal')
        self.ax.axis('off')
        self.fig.tight_layout()
        
    def get_clip_path(self):
        '''Create a plt.Compound path and return it as a Patch for clipping.'''
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path
        
        all_polys = [self.cation, self.anion, self.diamond]
        vertices = np.concatenate([i.get_path().vertices for i in all_polys])
        codes = np.concatenate([i.get_path().codes for i in all_polys])
        return PathPatch(Path(vertices, codes), transform=self.ax.transData)
         
    def get_canvas(self):
        '''Get fig and ax objects'''
        return self.fig, self.ax

class Overlap(Skeleton):
    
    def __init__(self, diagrams, way, df_name, overlap_var='Area', **kw):
        ''' Create a figure with all the groups information. Groups will be
        drawn by their convex hull or centroids according to kw['overfig_way'].
        Groups will be graduated by the values of one of the following variables:
        'Area','Shape','Blau','Angle','Time', according to kw['overfig_var']
        Note that 'Time' only available if way=='by_time'.
        '''
                
        super().__init__()
        self.kw_format = kw
        if (way == 'by_groups') & (overlap_var == 'Time'):
            warnings.warn(f"Cannot create Overlap figure by {overlap_var}.\n\
            when the analysis method is 'by_groups'. Automatically changed \n\
            to by 'Area'.")
            overlap_var = 'Area'
        
        self.cmap = self._cmap(kw['cmap_name'])
        self.cmap = self.cmap.reversed() if kw['cmap_reversed'] else self.cmap
        self.df_name = df_name
        self.diagrams, self.overlap_var = diagrams, overlap_var
        self.diagrams.sort(key=lambda x: x.cation.Params.area, reverse=True)
        
        # Draw groups by their convex hulls or centroids:
        if kw['overlap_way'] == 'centroids':
            self.extract_centroid(); self.colorbar(kw['overlap_way'])
        elif kw['overlap_way'] == 'polygons':
            self.extract_pols();     self.colorbar(kw['overlap_way'])
        else:
            raise KeyError("Only possible by 'centroid' or 'polygons")
                           
        self.layout()
        
    def extract_pols(self):
        '''Extract polygons and its values and keep it in a plt.Collection'''
        pols = [plt.Polygon(j.Params.ch_points) for
                i in self.diagrams for j in i.panels]
        values = np.array([j.Params.get()[self.overlap_var] for i in
            self.diagrams for j in i.panels if self.overlap_var != 'Time'])
        
        if self.overlap_var == 'Time':
            values = np.array([mdates.date2num(i.__str__()).mean()
                               for i in self.diagrams])
            values = np.repeat(values, 3) # Repeat values for each panel
            
        '''Creation of a Matplotlib collection for representation'''
        kw_col = {'match_original':False, 'lw': 1,'facecolor':'none',
                  'cmap':self.cmap, 'alpha':0.8, 'zorder':3}
        
        col = PatchCollection(pols, **kw_col)
        col.set_array(values)
        col.set_clip_path(self.get_clip_path())
        self.ax.add_collection(col)
        self.col, self.values = col, values
   
    def extract_centroid(self):
        '''Extract polygons and its values and keep it in a plt.Collection'''
        points = [j.Params.centroid() for i in self.diagrams for j in i.panels]
        df_pts = pd.DataFrame(points, columns=['X','Y'])
            
        if self.overlap_var == 'Time':
            values = [mdates.date2num(i.__str__()).mean() for i in self.diagrams]
            df_pts['Values'] = np.repeat( np.array(values), 3)
        else:
            df_pts['Values'] = [j.Params.get()[self.overlap_var] for i in \
            self.diagrams for j in i.panels if self.overlap_var != 'Time']
        
        kw = {'s':8, 'alpha':0.8, 'cmap':self.cmap, 'lw':0, 'zorder':3}
        
        self.ax.scatter(df_pts['X'], df_pts['Y'],
                c=df_pts['Values'], **kw).set_clip_path(self.get_clip_path())

        self.values = df_pts['Values'].to_numpy()
        
    def colorbar(self, overfig_way):
        '''Colorbar'''
        mappable = cm.ScalarMappable(norm=Normalize(self.values.min(),
                                                    self.values.max()),
                   cmap=self.cmap) if overfig_way == 'centroids' else self.col
                
        cbaxes = self.fig.add_axes([0.21, 0.5, 0.025,0.3])
        cbar = self.fig.colorbar(mappable, cax=cbaxes, spacing = 'uniform',
          extendrect = True, extendfrac = 'auto', drawedges = False)
        cbar.outline.set_linewidth(0.3) # Contour colorbar lines.
        cbar.ax.tick_params(axis='y',width=0.3,length=2, labelsize=5)    
        cbar.solids.set_edgecolor("face")  # Do not draw edge lines
        
        if self.overlap_var != 'Time':
            cbar.set_ticks([min(self.values), max(self.values)])
        else:
            ticks = np.linspace(min(self.values), max(self.values), 5)
            ticklabels = [mdates.num2date(i).strftime("%d-%m-%Y") for i in ticks]
            cbar.set_ticks(ticks)
            cbar.ax.set_yticklabels(np.array(ticklabels), fontsize=4)
            
    def _cmap(self, cmap_name):
        '''Edit cmap'''
        cmap = cm.get_cmap(cmap_name)
        # Shortening colormap range not to display white colors:
        color_list = [cmap(i) for i in np.linspace(40, cmap.N, 3).astype(int)]
        return lsc.from_list('Custom map', color_list, N=cmap.N)
    
    def layout(self):
        '''Figure layout'''
        units = {'Area':' (%)','Shape':' (%)','Angle':' ({})'.format('$^{o}$'),
                'Blau':' (%)','Dispersion':'', 'Time':''}
        
        text= "Graduation of all groups' convex-hulls \naccording to {}{}. {}"\
               .format(  r'$\bf{}$'.format(self.overlap_var),
                       units[self.overlap_var],
                       '${}$'.format(self.kw_format['overlap_way'].capitalize()))
        self.ax.text(-15, 195,text,fontsize=5, linespacing=1.5)
        self.fig.savefig(f'Graphics/Overlapped_{self.df_name}.{self._extension_graph}',
                         bbox_inches = 'tight', pad_inches=0, dpi=400)

class Plot_Diagram(Skeleton):
    ''' Class that allows to plot all Piper diagram information'''
    
    def __init__(self, panels, name, way, df_name,
                 labelling=False, df=None, **kw):
        self.kw_format = kw
        self.name, self.panels = name, panels
        self.way, self.df_name =  way, df_name
        self.name = '{} to {}'.format(*[i.strftime(format='%d-%m-%Y')
                    for i in name]) if way == 'by_time' else name

        super().__init__()
        self.points(); self.ch_points(); self.facies()
        self.lines(clip=kw['clip_lines']); self.info()
        self.outliers() if (kw['plot_outliers'] & kw['lof']) else None
        self.layout()
        
    def points(self):
        '''Plot all group's points for each panel and clip them'''
        kw = {'s':self.kw_format['point_size'],
              'alpha':self.kw_format['point_transparency'],
              'color':self.kw_format['color_inner_points'],
              'lw':0, 'zorder':3}
        points = np.concatenate(tuple(i.Params.points for i in self.panels), axis=0)
        scatter = self.ax.scatter(points[:,0], points[:,1], **kw)
        scatter.set_clip_path(self.get_clip_path())
        
    def outliers(self):
        '''Plot outlier points for each panel and clip them'''
        kw = {'s':self.kw_format['point_size'],
              'alpha':self.kw_format['point_transparency'],
              'color':self.kw_format['color_outliers'],
              'lw':0, 'zorder':3}
        points = np.concatenate(tuple(i.Params.outliers for i in self.panels), axis=0)
        scatter = self.ax.scatter(points[:,0], points[:,1], **kw)
        scatter.set_clip_path(self.get_clip_path())
    
    def ch_points(self):
        '''Plot only convex-hull's points'''
        kw = {'s':self.kw_format['point_size'],
              'alpha':self.kw_format['point_transparency'],
              'color':self.kw_format['color_outter_points'],
              'lw':0, 'zorder':4}
        points = np.concatenate(tuple(i.Params.ch_points for i in self.panels), axis=0)
        scatter = self.ax.scatter(points[:,0], points[:,1], **kw)
        scatter.set_clip_path(self.get_clip_path())
        
    def facies(self):
        ''' Plot convex hull polygons'''
        kw = {'lw':self.kw_format['pol_lw'], 
              'ec':self.kw_format['pol_ec'],
              'fc':self.kw_format['pol_color'],
              'alpha':self.kw_format['pol_transparency'],
              'zorder':1}
        pols = [plt.Polygon(i.Params.ch_points) for i in self.panels]
        col = PatchCollection(pols, **kw)
        col.set_clip_path(self.get_clip_path())
        self.ax.add_collection(col)
         
    def lines(self, clip=False):
        '''Plot regression lines'''
        def regresion(array):
             x,y = array[:,0], array[:,1]
             coeffs = np.polyfit(x,y,1)
             poly = np.poly1d(coeffs)
             y_fit= poly(x)
             return np.stack((x, y_fit)).T
         
        lines = [regresion(i.Params.points) for i in self.panels]
        line_col = LineCollection(lines, linestyle='solid',  zorder=4,
              color=self.kw_format['line_color'],
              lw=self.kw_format['line_lw'])
        line_col.set_clip_path(self.get_clip_path()) if clip else None
        self.ax.add_collection(line_col)
        
    def info(self):
        '''Plot the facies' indexes of each group as a text box.'''
        fancy_kw = dict(lw=0.2, pad=0.5)
        box_kw   = dict(boxstyle='round', facecolor='#F3F3F3',**fancy_kw)
        kw_text  = dict(transform=self.ax.transAxes, va='center', bbox=box_kw, 
                        fontfamily='sans-serif', fontsize=5.5, linespacing=1.8)
        info = [i.info() for i in self.panels]
        [self.ax.text(i[0], i[1], i[2], **kw_text) for i in info]
        
        # -----------------       Title Info    ------------------------------
        text = '\n'.join([
                "Dataset: ${}$".format(self.df_name.replace('_','\_')),
                "{}: {:<}".format('Group' if self.way=='by_groups' else 'Window',
                                   self.name)])
                # "Points: {:<}".format(  len(self.panels[0].Params.points)  )])
        self.ax.text(0.06, 0.875, text, **dict(transform=self.ax.transAxes,
                     fontsize=7), linespacing=1.65)
        
        # Outliers legend
        self.ax.text(0.06, 0.84, "Outliers (neighbours={})".format(
                                 self.kw_format['lof_neighbours']),
                     fontsize=5.5, color=self.kw_format['color_outliers'],
                     **dict(transform=self.ax.transAxes, ha='left')) \
                     if (self._check_outlier() & self.kw_format['lof']) else None
    
    def layout(self):
        '''Layout propertires: saving plot and closing figure'''
        self.fig.savefig(f"Graphics/{self.name}_{self.df_name}.{self._extension_graph}",
                         bbox_inches = 'tight', pad_inches=0, dpi=self._dpi)
        plt.close()
        
    def _check_outlier(self):
        '''Check if any outlier has been detected'''
        try:
            outlier_list = [True if self.panels[i].Params.outliers.size != 0 \
                       else False for i in range(3)]
            return any(outlier_list)
        except:
            return False

class Plot_all(Skeleton):
    
    def __init__(self, way, df_name, df=None,
                        diagrams=None, **kw):
        ''' Plot all the dataset points and graduate them by time.
        Only available if way == 'by_time'. '''
        super().__init__()
        self.df_name = df_name
        # Points features:
        cmap = cm.get_cmap(kw['cmap_name'])
        cmap = cmap.reversed() if kw['cmap_reversed'] else cmap
        kw = {'s':8, 'alpha':0.8, 'cmap':cmap, 'lw':0, 'zorder':3}

        if way=='by_time':
            self.df = df
            self.df['Time_fmt'] = mdates.date2num(self.df['Time'])       

            self.scatt = self.ax.scatter(df['X'], df['Y'], c=df['Time_fmt'], **kw)
            self.scatt.set_clip_path(self.get_clip_path())
            self.colorbar();  self.layout(way)

        elif way=='by_groups':
            dfs = [i.get_all_points() for i in diagrams]

            [self.ax.scatter(i['X'], i['Y'],
                    label=i['Group_ID'].iloc[0], **kw
                    ).set_clip_path(self.get_clip_path())
                    for i in dfs]
            self.layout(way)

    def colorbar(self):
        '''Colorbar'''
        cbaxes = self.fig.add_axes([0.21, 0.5, 0.025,0.3])
        ticks = np.linspace(self.df['Time_fmt'].iloc[0],
                             self.df['Time_fmt'].iloc[-1],  5)
        ticklabels = [mdates.num2date(i).strftime("%d-%m-%Y") for i in ticks]
        
        cbar = self.fig.colorbar(self.scatt, cax=cbaxes, spacing = 'uniform',
          extendrect = True, extendfrac = 'auto', drawedges = False,
          ticks=ticks)
        cbar.outline.set_linewidth(0.3) # Contour colorbar lines.
        cbar.ax.tick_params(axis='y',width=0.3,length=2, labelsize=4)    
        cbar.solids.set_edgecolor("face")  # Do not draw edge lines
        cbar.ax.set_yticklabels(ticklabels)
               
    def layout(self, way):
        '''Figure layout'''
        # Information plot on the figure:

        if way=='by_time':
            text= "Graduation of all data-file\npoints according to {}.\
            ".format(r'$\bf{Time}$')

        else:
            text= "Display of all dataset {}".format(r'$\bf{groups}$')
            #-------------------------- Legend ------------------------------------
            legend = self.ax.legend(loc = 'upper right', bbox_to_anchor = (0.95, 0.86),
                      frameon = False, ncol = 1, fontsize = 6, handletextpad = 0.6)


        self.ax.text(-15, 196,text,fontsize=5, linespacing=1.5)
        self.ax.text(155, 200, f"Dataset: ${self.df_name}$",
                 fontsize=6, linespacing=1.5)
        # Saving figure at 'Graphics'
        self.fig.savefig(f'Graphics/All points_{self.df_name}.{self._extension_graph}',
                         bbox_inches = 'tight', pad_inches=0, dpi=400)

class Evolution:
    '''Create a three-graphs plot with all the indexes info.'''
    
    def __init__(self, diagrams, df, df_name, way, mapa_dict):
        ''' Take all groups' params and plot them throughout time or by groups
            depending on the classifier method ('way') '''
        self.fig, self.ax = plt.subplots(3, figsize=(20,17), sharex=True, sharey=True)
        self.ax = self.ax.tolist()
        self.axt = [i.twinx() for i in self.ax]
        self.axis = (self.ax, self.axt)  #Tuple of axis for representing
        self.df, self.df_name = df, df_name
        self.mapa_dict = mapa_dict
        
        # Colors
        self.mapping = dict(Area='#DC2828',  # Red
                            Shape='#4178D4', # Blue
                            Blau='#32A941',  # Green
                            Angle='orange')  # Orange
        
    def legend(self, way, zoom10=False, zoom100=False):
        '''Format the legend of the figure. Extract labels, add units and
        plot legend. Extra info is plotted if kwargs are True'''
        
        handle1, labels1 = self.ax[0].get_legend_handles_labels()
        handle2, labels2 = self.axt[0].get_legend_handles_labels()
        handles, labels = [*handle1,*handle2], [*labels1,*labels2]
        
        if way == 'by_groups':
            labels = [i.split(' ')[1].replace(')','') for i in labels]
            
        labels = map(self.mapa_dict.get, labels)
        self.fig.legend(handles, labels,bbox_to_anchor=(0.42,0.76,0.2,0.2),
                   loc= 'upper center',fontsize= 20, ncol=4, borderpad=0.5)
        
        # Plot aditional info if some ts is zoomed:
        self.fig.text(0.3,0.14, 'Dashed line (---):   zoomed x10', ha='center',
             fontsize=18) if zoom10 else None
        self.fig.text(0.7,0.14, 'Dashed-dotted line (-.-):   zoomed x100',
             ha='center', fontsize=18) if zoom100 else None

    def layout(self):
        global _extension_graph
        # Erase tick from first and second plots.
        [i.tick_params(axis='x', bottom=False) for i in (self.ax[0], self.ax[1])]
        # Labelsize of all labels
        [j.tick_params(axis='both', labelsize=20) for i in self.axis for j in i]

        # Axis titles
        [i.set_title(r'$\bf{}$'.format(j), fontdict=dict(fontsize=18, ha='center'),
        pad=10) for i,j in zip(self.ax, ('Cation','Anion','Diamond'))]
        
        # Set a grid
        [i.grid('both',which='major') for i in self.ax]
        
        # Remove spines from all plots
        [j.spines[l].set_visible(False) for i in self.axis for j in i
                                      for l in ['top','bottom','left','right']]
        
        self.fig.subplots_adjust(hspace=0.1)  # Vertical space between subplots 
        self.fig.savefig(f'Graphics/Stability Evolution_{self.df_name}.{_extension_graph}',
            dpi=400, facecolor='#E1E1E1', bbox_inches = 'tight', pad_inches=1)
        
class Evolution_time(Evolution):
    '''Represent groups' indexes as temporal series.'''
    
    def __init__(self, diagrams, df, df_name, way, mapa_dict):
        
        super().__init__(diagrams, df, df_name, way, mapa_dict)
        
        # Preparing legend information:
        self.zooms10 = []
        self.zooms100 = []

        # MAIN:
        self._format_for_time()
        [self._lines(i,**dict(color=j)) for i,j in self.mapping.items()]
        self.legend(way, zoom10=any(self.zooms10), zoom100=any(self.zooms100))
        self._ticks(way); self.layout()
            
    def _format_for_time(self):
        ''' Create a field with the midpoint of the the group-time interval.'''
        self.df['Onset'] = self.df.Group.map(lambda x: pd.to_datetime(x.split(' to ')[0]))
        self.df['End'] = self.df.Group.map(lambda x: pd.to_datetime(x.split(' to ')[1]))
        self.df['Time'] = self.df['Onset'] + (self.df['End'] - self.df['Onset'])/2
        self.df = self.df.sort_values('Time').drop(['Onset','End','Group'],
                                                    axis=1).set_index('Time')
    def _lines(self, var, **kw):
        '''Plot the variable time-series on each panel iteratively'''
        # Define pandas plotting keywargs:
        kw_com = dict(lw=2, label=var, legend=False, xlabel='', **kw)
        kw_com = dict(ylabel='%', ylim=(0,100), **kw_com) if var!='Angle' else \
                 dict(ylabel='$^{o}$', ylim=(0,180), **kw_com)
        kw_com_x10 = dict(style='--', **kw_com) # For plotting zoom ts (x10).
        kw_com_x100 = dict(style='-.', **kw_com) # For plotting zoom ts (x100).
        
        # Pivot the df by certain variable and sort groups.
        ts = pd.pivot_table(self.df,values=var, index='Time', columns=['Panel'])
        side = self.axis[0] if var!='Angle' else self.axis[1]
        
        def check(ts, column):
            '''Check temporal series and zoom by 10 if max(ts) < 10.'''
            if (ts[column].max() < 10) & (ts[column].mean() != 0):
                ts[column] = ts[column]*10
                zoom = True
            else:
                zoom = False
            return ts, zoom
        
        # Check each ts' variable and plot it depending on the zoom made:            
        for column, axis in zip(ts.columns, side):
            # First Zoom (x10)
            ts, zoomx10 = check(ts, column)
            # Second Zoom (x100) if max. still lower than 10.
            ts, zoomx100 = check(ts, column)
            
            if zoomx100:
                ts.plot(y=column, ax=axis, **kw_com_x100)
            elif zoomx10:
                ts.plot(y=column, ax=axis, **kw_com_x10)
            else:
                ts.plot(y=column, ax=axis, **kw_com)
                
            self.zooms10.append(zoomx10)
            self.zooms100.append(zoomx100)
        
    def _ticks(self, way):
        '''Format ticks of axis figure depending on 'way' param chosen.'''
        ax, axt = self.ax, self.axt
        # Y-axis format:
        [i.set_yticks(range(0,101,20)) for i in ax]
        [i.set_yticks(range(0,181,30)) for i in axt]
        [i.set_ylabel('%', fontsize=20, rotation=180, ha='left',va='center') for i in ax]
        [i.set_ylabel('Degrees', fontsize=20, rotation=270,ha='center',va='baseline') for i in axt]
        
        # X-axis format:
        years_fmt = mdates.DateFormatter('%Y')
        ax[2].xaxis.set_major_locator(mdates.YearLocator())
        ax[2].xaxis.set_major_formatter(years_fmt)
        
        months_fmt = mdates.DateFormatter('')
        ax[2].xaxis.set_minor_locator(mdates.MonthLocator())
        ax[2].xaxis.set_minor_formatter(months_fmt)

        ax[2].set_xlim(self.df.index.min(),self.df.index.max())
        ax[2].tick_params(axis='x',labelrotation=0, pad=12)
        [i.set_horizontalalignment('center') for i in ax[2].get_xticklabels()]
        
        # Erase month ticks set wrongly on the first and second plot
        [ax[i].xaxis.set_ticks_position('none')  for i in (0,1)]
        
class Evolution_groups(Evolution):
    '''Represent groups' indexes as stacked bar histograms'''
    
    def __init__(self, diagrams, df, df_name, way, mapa_dict):
        super().__init__(diagrams, df, df_name, way, mapa_dict)
        
        self.hist(); self.legend(way); self.ticks(way); self.layout();
        
    def hist(self):
        '''Plot the variable time-series on each panel iteratively'''
        d = pd.pivot_table(self.df,values=['Area','Shape','Blau','Angle'],
                            index='Group', columns=['Panel'])
        d.sort_index(inplace=True) # Sort the index on the way pandas does.

        for i,j in zip(('cation','anion','diamond'), self.ax):
            bars = d.swaplevel(0,1,axis=1).sort_index(axis=1)[[i]]
            bars.plot(kind="bar", stacked=True, ax=j, color=self.mapping.values(),
                      legend=False, xlabel='')
        
    def ticks(self, way):
        '''Format ticks of axis figure depending on 'way' param chosen.'''
        [i.remove() for i in self.axt]

        # Axis labels format
        self.ax[1].set_ylabel('Stacked Values', fontsize=20, rotation=90,
                              ha='center',va='center', labelpad=15)
        
        [j.set_horizontalalignment('right') for i in self.ax
                                            for j in i.get_xticklabels()]
        self.ax[2].tick_params(axis='x',labelrotation=30, pad=12)
        
    
# =============================================================================
# ========================= CUSTOM WARNING =====================================
# =============================================================================

def format_warning(msg, *args, **kwargs):
    return f'\nWarning: {msg}\n'

warnings.formatwarning = format_warning