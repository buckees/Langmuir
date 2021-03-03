"""Diagnose the sheath model."""

import numpy as np
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
        self.col = ['Init_Vel_x', 'Init_Vel_z', 'Init_Vel_y', 
                    'Init_Erg', 'Init_Ang',
                    'End_Vel_x', 'End_Vel_z', 'End_Vel_y', 
                    'End_Erg', 'End_Ang',
                    'Collision', 'Lifetime'
                    'hitWafer']
        self.df = pd.DataFrame(columns=self.col)

    def save2csv(self, fcsv='Sheath2D'):
        """Save stats to csv file."""        
        self.df.to_csv(fcsv + '_Stats.csv', mode='w', 
                       index=True, header=True,
                       na_rep='NA')
    
    def plot(self, fpng='Sheath2D'):
        """Plot stats to png file."""
        fig, axes = plt.subplots(2, 2, figsize=(8, 6), dpi=600,
                                 constrained_layout=True)
        
        ax = axes[0, 0]
        self.df['Init_Erg'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Ion Energy Distribution')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Count')
        
        ax = axes[0, 1]
        self.df['Init_Ang'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Ion Angular Distribution')
        ax.set_xlabel('Angle (degree)')
        ax.set_ylabel('Count')
        
        ax = axes[1, 0]
        self.df['End_Erg'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Ion Energy Distribution')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Count')
        
        ax = axes[1, 1]
        self.df['End_Ang'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Ion Angular Distribution')
        ax.set_xlabel('Angle (degree)')
        ax.set_ylabel('Count')
        
        fig.savefig(fpng + '_IAEDF.png', dpi=600)
        plt.close()
        
        fig, axes = plt.subplots(2, 2, figsize=(8, 6), dpi=600,
                                 constrained_layout=True)
        
        ax = axes[0, 0]
        
        temp = self.df['Collision']
        cout = Counter(temp)
        ax.bar(cout.keys(), cout.values())
        ax.set_title('Collision')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Count')
        
        ax = axes[0, 1]
        self.df['Lifetime'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Lifetime Distribution')
        ax.set_xlabel('Lifetime (s)')
        ax.set_ylabel('Count')
        
        ax = axes[1, 0]
        self.df['End_Erg'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Ion Energy Distribution')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Count')
        
        ax = axes[1, 1]
        self.df['End_Ang'].hist(bins=100, density=False, ax=ax)
        ax.set_title('Ion Angular Distribution')
        ax.set_xlabel('Angle (degree)')
        ax.set_ylabel('Count')
        
        fig.savefig(fpng + '_Stats.png', dpi=600)
        plt.close()

