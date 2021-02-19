"""Process chemistry input files."""

import pandas as pd
import numpy as np

fsp = 'Ar_Cl2_v01_Species'
fmat = 'Ar_Cl2_v01_Material'
frct = 'Ar_Cl2_v01_Reaction'

col = ['Species', 'Material',
       'Product1', 'Product2', 'Product3',
       'Removal1', 'Removal2', 'Removal3',
       'Reaction_Type',
       'Energy_Type', 'Special_Number',
       'A', 'n', 'erg_th', 'erg_ref',
       'Angle_Type']

df_sp = pd.read_csv(fsp + '.csv', header=0)
df_mat = pd.read_csv(fmat + '.csv', header=0)
df_rct = pd.read_csv(frct + '.csv', header=0)
# df_rct = pd.DataFrame(columns=col)

for sp in df_sp['Species']:
    for mat in df_mat['Material']:
        #       Sp, Mat, Prod1, Prod2, Prod3, Rm1, Rm2, Rm3, 
        temp = [sp, mat, sp, mat, None, None, None, None,
        #       Rct_Type, Erg_Type, Sp_Num,
                'Reflect', 'Constant', 0,
        #       A, n, erg_th, erg_ref, Angle_Type
                1.0, 0.0, 0.0, 0.0, 0]
        df_rct = df_rct.append(pd.Series(temp, index=col), ignore_index=True)

