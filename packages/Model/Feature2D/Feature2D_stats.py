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
    
    def save2csv(self, fcsv='Feature2D'):
        """Save stats to csv file."""        
        self.df_sp.to_csv(fcsv + '_Stats.csv', mode='w', 
                          columns=['Launched', 'Escaped', 'Escaped Pct',
                                'Etch', 'Etch Pct',
                                'Terminated', 'Terminated Pct'],
                          index=True, header=True, 
                          float_format='%.2f', na_rep='NA')
        
        with open(fcsv + '_Stats.csv', 'a') as f:
            f.write('\n')
    
        stats.df_mat.to_csv(fcsv + '_Stats.csv', mode='a', 
                            index=True, header=True,
                            float_format='%.2f', na_rep='NA')