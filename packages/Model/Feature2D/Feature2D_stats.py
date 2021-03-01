"""Diagnose the feature model."""

import pandas as pd

class STATS(object):
    """Mesh object."""

    def __init__(self, fname='Species'):
        """
        Init the MESH2D.
        
        fname: str, filename of xxx_Species.csv.
        """
        self.col = ['Launched', 'Init Posn', 'Init Vel', 'Init Erg', 'Init Ang',
                    'Escaped', 'Reflection', 'Etch', 'Terminated']
        df_sp = pd.read_csv(fname + '_Species.csv', header=0)
        self.df = pd.DataFrame(columns=self.col, index=df_sp['Species'])
        self.df['Launched'] = 0
        self.df['Init Posn'] = self.df.apply(lambda x: [], axis=1)
        self.df['Init Vel'] = self.df.apply(lambda x: [], axis=1)
        self.df['Init Erg'] = self.df.apply(lambda x: [], axis=1)
        self.df['Init Ang'] = self.df.apply(lambda x: [], axis=1)
        self.df['Escaped'] = 0
        self.df['Reflection'] = self.df.apply(lambda x: [], axis=1)
        self.df['Etch'] = 0
        self.df['Terminated'] = 0