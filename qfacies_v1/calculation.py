# -*- coding: utf-8 -*-
"""
===============================================================================
======================== Module of Q-Facies package ============================
===============================================================================
Contains two classes that allow to apply all the euclidean transformations 
requiered, calculate all of the indeces for each panel, and to identify outliers.
===============================================================================
===============================================================================
"""
import numpy as np

_lof = False  # Whether to consider outlier analyses or not.

def import_skl():
    try:
        from sklearn.neighbors import LocalOutlierFactor
        from matplotlib.path import Path
        global _lof
        _lof = True
        global LocalOutlierFactor, Path
        
    except ModuleNotFoundError:
        raise ('Cannot find sklearn package.')
        
class Indixes:
    ''' Class that includes all calculation methods of Q-Facies indices.'''
    
    def __init__(self, points, df, panel, **kw):
        
        global _lof
        self.panel = panel
        self.df = df
        #self.df.dropna(inplace=True)
        self.points = points
        
        # Calculate convex hull taking into account outliers or not.
        self.LOF(neighbours=kw['lof_neighbours']) if _lof else None
        self.ch_points = self.ConvexHull(self.points)
        # Create indixes
        self.main()                        
                      
    def main(self):
        self.num_points = self.df.shape[0]
        self.area       = self.gauss()
        self.perimeter  = self.perim()
        self.shape_idx  = self.shape(self.area, self.perimeter)
        self.blau_idx   = self.blau()  # Tuple (%, values)
        self.angle      = self.orientation()
        self.sd         = self.dispersion()
        
    def get(self):
        '''Get all Group parameters.'''
        return dict(Area=self.area, Shape=self.shape_idx,
                   Angle=self.angle, Panel=self.panel,
                   Blau=self.blau_idx[0], Dispersion = self.sd,
                   points=self.df.shape[0], Dominant = self.blau_idx[2],
                   **self.blau_idx[1])
         
    def ConvexHull(self, points):
        '''
        Create the convex hull of a given array and return the resulting coordinates
        array. The convex hull is calculated by the 'Graham scan' method [1]
        and implemented from RodolfoFerro's code [2].
        
        [1] Graham, R. L. (1972). An efficient algorithm for determining the
            convex hull of a finite planar set. Info. Pro. Lett., 1, 132-133.
            
        [2] RodolfoFerro, ConvexHull.GrahamScan (2015), Github.
            https://github.com/RodolfoFerro/ConvexHull
        '''
        
        def RightTurn(p1, p2, p3):
            ''' Function to determine if we have a counterclock-wise turn.'''
            if (p3[1]-p1[1])*(p2[0]-p1[0]) >= (p2[1]-p1[1])*(p3[0]-p1[0]):
                return False
            return True

        def GrahamScan(P):
        	P.sort()			# Sort the x-coordintate of the set of points
        	L_upper = [P[0], P[1]]		# Initialize with the left-most point
        	# Compute the upper part of the hull
        	for i in range(2,len(P)):
        		L_upper.append(P[i])
        		while len(L_upper) > 2 and not RightTurn(L_upper[-1],L_upper[-2],L_upper[-3]):
        			del L_upper[-2]
        	L_lower = [P[-1], P[-2]]	# Initialize the lower part
        	# Compute the lower part of the hull
        	for i in range(len(P)-3,-1,-1):
        		L_lower.append(P[i])
        		while len(L_lower) > 2 and not RightTurn(L_lower[-1],L_lower[-2],L_lower[-3]):
        			del L_lower[-2]
        	del L_lower[0], L_lower[-1]
        	L = L_upper + L_lower		# Build the full hull
            
        	return np.array(L)
        
        def main(points):
            ''' Execute the Convex-Hull. First, from array to list of tuples.'''
            points = [(x,y) for x,y in zip(points[:,0], points[:,1])]
            L = GrahamScan(points)
            
            return L
       
        return main(points)
    
    def centroid(self):
        ''' Centroid coordinates calculation: axis means'''
        return self.points[:,0].mean(), self.points[:,1].mean()
    
    def LOF(self, neighbours=50):
        '''
        Two-steps outlier detection method:
            1_ Unsupervised Outlier Detection using Local Outlier Factor (LOF),
               from sklearn.neighbors.LocalOutlierFactor.
            2_ Exclude from analysis only those outliers outside the
                convex hull polygon.

        Return None. It just removes outlier points by modifying self.points
        attribute and create a new variable to contain them.'''
        
        mask_1 = LocalOutlierFactor(neighbours).fit_predict(self.points)
        mask_1 = np.where(mask_1 == 1, True, False).astype(bool)
        
        # Check if preliminary outliers lies within the ch polygon
        initial_inliers = np.compress(mask_1, self.points, axis=0)
        ch_path = Path(self.ConvexHull(initial_inliers))
        isin_pol = np.apply_along_axis(lambda x: ch_path.contains_point(x),
                                        1, self.points).astype(bool)
        # Do not erase outliers if they fall within ch polygon
        mask = np.logical_or(mask_1, isin_pol)
        lof_points = np.compress(mask, self.points, axis=0)
        self.outliers = np.compress(~mask, self.points, axis=0)
        self.points = lof_points # Redefine points to erase Outliers.
        
        # Delete selected outliers from df:
        self.df = self.df[mask]
    
    def gauss(self):
        """
        Calculation of polygon area using Gauss theorem.
        This is computed with the points sorted in a non-clockwise sense.
        
        a) first diagonal of the determinant b) second one c) Absolute value
        
        return: polygon area normalized to the polygon extent
        """
        # Last coordinate must be repeted in order to close the polygon.
        points = np.concatenate(([self.ch_points[-1]], self.ch_points))

        a = np.array([i * j for i, j in zip(points[:-1,0], points[1:,1])]).sum()
        b = np.array([i * j for i, j in zip(points[1:,0], points[:-1,1])]).sum()
        c = abs(a - b)/2
        
        if self.panel == 'anion' or self.panel == 'cation':  # Area as panel percentage
             return c/(50*np.sin(np.pi/3)) # Normalized area to the total triangle area 
        elif self.panel == 'diamond':
             return c/(100*np.sin(np.pi/3)) # Normalized area to the total diamond area 

    def perim(self):
        ''' Given an array with points sorted in a non-clockwise sense,
        it calculates the polygon perimiter'''
        # Last coordinate must be repeated in order to close the polygon.
        points = np.concatenate(([self.ch_points[-1]], self.ch_points))
        # Elements differences
        aristas = abs(np.diff(points,axis=0))
        # Pythagorean theorem for the calculation of the hypotenuse
        h = sum([(x**2 + y**2)**0.5 for x,y in zip(aristas[:,0], aristas[:,1])])
        return h
    
    def shape(self, area, perimeter):
        ''' Shape Index, introduced by Richardson in 1961 (see [3])
        Whilst values close to 100% indicate a circular-like form, 
        values close to zero mean the contrary (linear-like form).
        
        [3] Haggett, P., Cliff, A. D., & Frey, A. (1977). Locational analysis in
                 human geography (2nd ed.). London: Edward Arnold Ltd.
        
        '''
        area = area*50 if self.area == 'anion' or 'cation' else area*100
        return ( (4 * np.pi * area)/(perimeter**2) ) * 100
   
    def orientation(self):
        ''' Return the angle of the linear regression model.'''
        x,y = self.points[:,0], self.points[:,1]
        coeffs = np.polyfit(x,y,1)
        angle = (np.arctan(coeffs[0])) * 180/np.pi
        return angle + 180 if angle<0 else angle

    def blau(self):
        '''
        Calculation of the Blau index [4] that ranges between 0.25 and 1. Values are 
        normalized to the interval [0,100].
        
        This index gives information about the data dispersion among the four
        possible facies (fs) that conform every panel, being calculated for
        each one of them.
        
        [4] Blau, P.M. (1977) Inequality and Heterogeneity: A Primitive Theory
                  of Social, Nev York: The Free Press.
        '''
        df = self.df
        total = df.shape[0]     # Total number of samples
        
        # Facies maps:
        cation_map = dict(A='Magnesium',B='Calcium',C='Sodium-Potassium',D='Mixed')
        anion_map  = dict(A='Sulphate',B='Bicarbonate',C='Chloride',D='Mixed')
        panel_map = dict(cation=cation_map, anion=anion_map)
        
        def get_dominant(facies, panel_map):
            '''Extract the dominant facies on cation and anion panels. Could be
            more than one '''
            max_value = max(facies.values())
            dominants = [key for key,value in facies.items() if value==max_value]
            dominants = [panel_map[x] for x in dominants]
            return '-'.join(dominants) if len(dominants)>1 else dominants[0]
        
        # CATION FACIES (facies' points relative frequency).
        if self.panel == 'cation':
            facies = dict(
            A = df.Mg_epm.loc[df.Mg_epm >= 50].count()/total,  
            B = df.Ca_epm.loc[df.Ca_epm >= 50].count()/total,   
            C = df.NaK_epm.loc[df.NaK_epm >= 50].count()/total )
            facies['D'] = round( 1 - sum(facies.values()), 2)
            
        # ANION FACIES (facies' points relative frequency).
        elif self.panel == 'anion':
            facies = dict(
            A = df.SO4_epm.loc[df.SO4_epm >= 50].count()/total,
            B = df.HCO3CO3_epm.loc[df.HCO3CO3_epm >= 50].count()/total,
            C = df.Cl_epm.loc[df.Cl_epm >= 50].count()/total   )
            facies['D'] = round( 1 - sum(facies.values()), 2)
            
        # DIAMOND FACIES (facies' points relative frequency).
        elif self.panel == 'diamond':
            up_carb, down_carb = (df.HCO3CO3_epm >= 50), (df.HCO3CO3_epm < 50)
            up_NaK, down_NaK = (df.NaK_epm >= 50), (df.NaK_epm < 50)
            
            facies = dict(
            A = df.loc[down_carb & down_NaK].shape[0]/total,
            B = df.loc[down_carb & up_NaK].shape[0]/total,
            C = df.loc[up_carb & up_NaK].shape[0]/total,
            D = df.loc[down_NaK & up_carb].shape[0]/total      )

        ''' Blau formula and normalize values to 0 - 100 range'''
        index = 1 - sum([i**2 for i in facies.values()])
        normalized_index = round(     (index/0.75) * 100, 2   )
        dominant = get_dominant(facies, panel_map[self.panel]) if self.panel != \
            'diamond' else None
        
        return [normalized_index, facies, dominant]
 
    def dispersion(self):
        '''Standard distance index. Calculation based on the 'typical distance'
        concept by Roberto Bachi [5] *. Values ares returned normalized to the
        maximum dispersion value of each panel:
            For cation and anion panels: 57.73 distance units
            For diamond panel: 70.71 distance units
        
        [5] Bachi, R. (1963). Standard distance measures and related methods for
            spatial analysis. Papers of the Regional Science Association 10,
            83–132 (1963). https://doi.org/10.1007/BF01934680
        
        * Well explained at: 
        https://volaya.github.io/libro-sig/chapters/Estadistica_espacial.html
        '''  
        x, y = self.points[:,0], self.points[:,1]
        sd_distance =   np.sqrt( ((np.sum(x**2)/len(x)) - np.mean(x)**2) +
                                 ((np.sum(y**2)/len(y)) - np.mean(y)**2) )
        if self.panel == 'anion' or self.panel == 'cation':
            return sd_distance * 100/57.735026918962575
        elif self.panel == 'diamond':
            return sd_distance * 100/70.71067811865476
        
class Transform:
    '''
    Affine 2D transformation for plotting points into a Piper diagram.
    Transformation consists in a set of scale, shear, rotation and
    translocation operations that are applied to an xy-array of points.
    '''
    
    def __init__(self):
        self.offset = 22 # This must match with Plot.Skeleton.offset value.
    
    def rotation(self):
        '''Rotation transformation of 300 degrees. Aply for diamond panel'''
        return np.array([(np.cos(np.radians(300)), np.sin(np.radians(300))),
                            (-np.sin(np.radians(300)), np.cos(np.radians(300)))])
     
    def scale(self):
        '''Scale transformation. For cation and anion panels. '''
        return np.array([(1,      0        ),
                          (0,np.sin(np.pi/3))])
     
    def t_shear(self):
        '''Shear transformation. For cation and anion panels. '''
        return np.array([(       1       ,0),
                          (np.tan(np.pi/6),1)])
     
    def d_shear(self):
        '''Shear transformation for diamond panel'''
        return np.array([(1,0),
                          (-np.tan(np.pi/6),1)])
     
    def d_translation(self):
        '''Translation of diamond points'''
        Ax = 50 + self.offset/2
        Ay = np.sin(np.pi/3) * (100 + self.offset)
        return np.array([(Ax),
                         (Ay)])
    
    def a_translation(self):
        '''Translation of anion points'''
        return np.array([(self.offset + 100),
                    (0)])
