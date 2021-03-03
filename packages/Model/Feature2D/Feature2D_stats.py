"""Diagnose the feature model."""

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

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
    
        self.df_mat.to_csv(fcsv + '_Stats.csv', mode='a', 
                            index=True, header=True,
                            float_format='%.2f', na_rep='NA')
        
    def plot(self, fpng='Feature2D'):
        """Plot stats to png file."""
        for sp in self.df_sp.index:
            fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=600,
                                       constrained_layout=True)
            
            ax = axes[0]
            ax.hist(self.df_sp.loc[sp, 'Init Erg'], bins=100, density=False)
            ax.set_title('Energy Distribution')
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('a.u.')
            
            ax = axes[1]
            ax.hist(self.df_sp.loc[sp, 'Init Ang'], bins=100, density=False)
            ax.set_title('Velocity Distribution')
            ax.set_xlabel('Angle (degree)')
            ax.set_ylabel('a.u.')
            
            fig.suptitle(sp, fontsize=20)
            fig.savefig('Init_' + sp + '.png', dpi=600)
            plt.close()
        
        n_sp = len(self.df_sp)
        fig, axes = plt.subplots(1, n_sp, figsize=(6*n_sp, 6), dpi=600,
                                   constrained_layout=True)
        for i, sp in enumerate(self.df_sp.index):
            ax = axes[i]
            temp = self.df_sp.loc[sp, 'Reflection']
            cout = Counter(temp)
            ax.bar(cout.keys(), cout.values())
            ax.set_title(sp)
            ax.set_xlabel('Number of Reflections')
            ax.set_ylabel('Counts')
        fig.suptitle('Stats of Reflection', fontsize=20)
        fig.savefig('Reflection.png', dpi=600)
        plt.close()    