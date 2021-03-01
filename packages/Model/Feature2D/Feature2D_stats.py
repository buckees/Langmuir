"""Diagnose the feature model."""

import pandas as pd

class STATS(object):
    """Mesh object."""

    def __init__(self, fname='xxx'):
        """
        Init the MESH2D.
        
        fname: str, filename of xxx_Species.csv.
        """
        # self.col = ['Launched', 'Init Posn', 'Init Vel', 'Init Erg', 'Init Ang',
        #             'Escaped', 'Reflection', 'Etch', 'Terminated']
        self.df_sp = pd.read_csv(fname + '_Species.csv', header=0,
                                 index_col=0)
        self.df_sp['Launched'] = 0
        self.df_sp['Init Posn'] = self.df_sp.apply(lambda x: [], axis=1)
        self.df_sp['Init Vel'] = self.df_sp.apply(lambda x: [], axis=1)
        self.df_sp['Init Erg'] = self.df_sp.apply(lambda x: [], axis=1)
        self.df_sp['Init Ang'] = self.df_sp.apply(lambda x: [], axis=1)
        self.df_sp['Escaped'] = 0
        self.df_sp['Reflection'] = self.df_sp.apply(lambda x: [], axis=1)
        self.df_sp['Etch'] = 0
        self.df_sp['Terminated'] = 0
        
        self.df_mat = pd.read_csv(fname + '_Material.csv', header=0, 
                                  index_col=0)
        self.df_mat['Etch'] = 0