"""Process chemistry input files."""

import pandas as pd
import numpy as np
# import random

class CHEMISTRY(object):
    """Mesh object."""

    def __init__(self, name='Chemistry'):
        """
        Init the MESH2D.
        
        name: str, var, name of the CHEMISTRY.
        """
        self.name = name
        self.col = ['Species', 'Material', 'Reaction_Type',
                    'Product1', 'Product2', 'Product3',
                    'Removal1', 'Removal2', 'Removal3',
                    'Energy_Type', 'Special_Number',
                    'A', 'n', 'erg_th', 'erg_ref',
                    'Angle_Type']
        
    def load(self, fname):
        """
        Load chemistry info from 3 files.
        
        fname: str, var, basename for chemistry input files.
        """
        self.df_sp = pd.read_csv(fname + '_Species.csv', header=0)
        self.df_mat = pd.read_csv(fname + '_Material.csv', header=0)
        df_rct = pd.DataFrame(columns=self.col)
        # create default reflections with assigned default values 
        for sp in self.df_sp['Species']:
            for mat in self.df_mat['Material']:
                #       Sp, Mat, Rct_Type, 
                temp = [sp, mat, 'Reflect',
                #       Prod1, Prod2, Prod3, Rm1, Rm2, Rm3, 
                        sp, mat, 'None', 'None', 'None', 'None',
                #       Erg_Type, Sp_Num,
                        'Constant', 0,
                #       A, n, erg_th, erg_ref, Angle_Type
                        1.0, 0.0, 0.0, 0.0, 0]
                df_rct = df_rct.append(pd.Series(temp, index=self.col), 
                                       ignore_index=True)
        # readin reaction info
        df_rct2 = pd.read_csv(fname + '_Reaction.csv', header=0)
        # append new reactions or overwrite reflection values
        for i in range(len(df_rct2)):
            s_rct2 = df_rct2.loc[i]
            if s_rct2['Reaction_Type'] == 'Reflect':
                for j in range(len(df_rct)):
                    s_rct = df_rct.loc[j]
                    if s_rct['Reaction_Type'] == 'Reflect':
                        if (s_rct['Species'] == 
                            s_rct2['Species']) and (s_rct['Material'] == 
                                                    s_rct2['Material']):
                            df_rct.loc[j, 'A'] = s_rct2['A']
                            
            else:
                df_rct = df_rct.append(s_rct2, ignore_index=True)
        # add df_rct as attribute
        self.df_rct = df_rct

    def save(self, fname='df_rct'):
        """Save reaction info to file."""        
        self.df_rct.to_csv(fname + '.csv', index=True)
    
    def choose_react(self, species, material, energy=0.0, angle=0.0, 
                     idiag=False):
        """
        Return a Series format chosen reaction/reflection with given parameters.
        
        species: str, var
        material: str, var
        energy: float, var
        angle: float, var
        """
        df_sub = self.df_rct[(self.df_rct['Species'] == species) & 
                             (self.df_rct['Material'] == material)].copy()
        df_sub['Yield'] = 0.0
          
        for idx, row in df_sub.iterrows():
            # calc energy-dependent yield
            Yield = 0.0
            if row['Energy_Type'] == 'Constant':
                Yield = row['A']
            elif row['Energy_Type'] == 'Power':
                A, n, erg_th, erg_ref = \
                    row['A'], row['n'], row['erg_th'], row['erg_ref']
                if energy >= erg_th:
                    Yield = A*(energy - erg_th)**n/(erg_ref - erg_th)**n
            elif row['Energy_Type'] == 'Function':
                inum = row['Special_Number']
                Yield = self.func_yield(inum=inum, energy=energy, angle=angle)
            
            # calc anlge-dependent yield
            if row['Angle_Type'] == 0:
                pass
            elif row['Angle_Typle'] == 1:
                pass
            
            df_sub.loc[idx, ['Yield']] = Yield
            
            df_sub = self.calc_prob(df_sub)
        
        n_row = len(df_sub)
        weight = df_sub['Prob'].tolist()
        idx = np.random.choice(n_row, 1, p=weight)[0]
        chosen_rct = df_sub.iloc[idx]
        # return 
        if idiag:
            return chosen_rct, df_sub
        else:
            return chosen_rct

    @staticmethod
    def func_yield(inum, energy, angle):
        """
        Return yield with defined function.
        
        inum: int, special number
        energy: float, eV, ion energy
        angle: float, degree, ion angle w.r.t. surface normal
        """
        pass
        return 1.00
    
    @staticmethod
    def calc_prob(df):
        df['Prob'] = 0.0
        yield_tot = sum(df['Yield'])
        for idx, row in df.iterrows():
            Yield = row['Yield']
            Prob = Yield/yield_tot
            df.loc[idx, ['Prob']] = Prob 
        return df
