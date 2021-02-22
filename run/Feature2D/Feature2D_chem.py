"""Process chemistry input files."""

import pandas as pd
import numpy as np

fsp = 'Ar_Cl2_v01_Species'
fmat = 'Ar_Cl2_v01_Material'
frct = 'Ar_Cl2_v01_Reaction'

col = ['Species', 'Material', 'Reaction_Type',
       'Product1', 'Product2', 'Product3',
       'Removal1', 'Removal2', 'Removal3',
       'Energy_Type', 'Special_Number',
       'A', 'n', 'erg_th', 'erg_ref',
       'Angle_Type']

df_sp = pd.read_csv(fsp + '.csv', header=0)
df_mat = pd.read_csv(fmat + '.csv', header=0)
df_rct = pd.DataFrame(columns=col)

for sp in df_sp['Species']:
    for mat in df_mat['Material']:
        #       Sp, Mat, Rct_Type, 
        temp = [sp, mat, 'Reflect',
        #       Prod1, Prod2, Prod3, Rm1, Rm2, Rm3, 
                sp, mat, 'None', 'None', 'None', 'None',
        #       Erg_Type, Sp_Num,
                'Constant', 0,
        #       A, n, erg_th, erg_ref, Angle_Type
                1.0, 0.0, 0.0, 0.0, 0]
        df_rct = df_rct.append(pd.Series(temp, index=col), ignore_index=True)

df_rct2 = pd.read_csv(frct + '.csv', header=0)

for i in range(len(df_rct2)):
    s_rct2 = df_rct2.loc[i]
    if s_rct2['Reaction_Type'] == 'Reflect':
        for j in range(len(df_rct)):
            s_rct = df_rct.loc[j]
            if s_rct['Reaction_Type'] == 'Reflect':
                if (s_rct['Species'] == 
                    s_rct2['Species']) and (s_rct['Material'] == 
                                            s_rct2['Material']):
                    print('Value Matched.',
                          df_rct.loc[j, 'A'], s_rct2['A'], '\n')
                    df_rct.loc[j, 'A'] = s_rct2['A']
                    
    else:
        df_rct = df_rct.append(s_rct2, ignore_index=True)

print(df_rct[['Species', 'Material', 'Reaction_Type', 'Energy_Type','A']])
df_rct.to_csv('df_rct.csv',index=True)

print('\nNow the reaction dataframe has been built!')

def func_yield(inum, energy, angle):
    """
    Return yield with defined function.
    
    inum: int, special number
    energy: float, eV, ion energy
    angle: float, degree, ion angle w.r.t. surface normal
    """
    pass
    return 1.00

def choose_react(species, material, energy=0.0, angle=0.0):
    """
    Return chosen reaction/reflection given parameters.
    
    species: str,
    material: str,
    energy: float,
    angle: float,
    """
    df_sub = df_rct[(df_rct['Species'] == species) & 
                    (df_rct['Material'] == material)].copy()
    df_sub['Yield'] = 0.0
    df_sub['Probability'] = 0.0
    print('\nThe sub dataframe is filtered.')
    print(df_sub[['Species', 'Material', 'Reaction_Type', 'Energy_Type',
                  'A', 'n', 'erg_th', 'erg_ref', 'Angle_Type']])
    
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
            Yield = func_yield(inum=inum, energy=energy, angle=angle)
        
        # calc anlge-dependent yield
        if row['Angle_Type'] == 0:
            pass
        elif row['Angle_Typle'] == 1:
            pass
        
        df_sub.loc[idx, ['Yield']] = Yield
    
    print('\n', df_sub[['Species', 'Material', 'Reaction_Type', 'Energy_Type',
                  'Yield', 'Probability']])

choose_react('Ar+', 'Si_', energy=10.0, angle=0)