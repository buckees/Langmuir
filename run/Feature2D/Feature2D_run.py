"""Feature Model 2D. Main program."""
import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Common.Particle import PARTICLE
from packages.Model.Feature2D.Feature2D_mesh import MESH2D
from packages.Model.Feature2D.Feature2D_chem import CHEMISTRY
from packages.Model.Feature2D.Feature2D_rflct import REFLECT
from packages.Model.Feature2D.Feature2D_rct import REACT
from packages.Model.Feature2D.Feature2D_stats import STATS
from packages.Model.Feature2D.Feature2D_main import MAIN

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 10000
oper.max_step = 1000
oper.step_fac = 0.5
oper.max_rflct = 5
oper.surf_norm_range = 3
oper.surf_norm_mode = 'Sum Vector'
oper.num_plot = 5
oper.prob_rflct = 0.5
oper.Tn = 0.025
oper.idiag = True
oper.fname = 'Ar_Cl2_v01'

# init mesh
mesh = MESH2D()
# readin mesh
fname = 'SiEtch_VerticalPR_Mesh'
mesh.readin_mesh(fname)

# init ptcl
ptcl = PARTICLE()
ptcl.load_database()
# ptcl.customize_ptcl('Ion', 40, 1)

def xFunc():
    return mesh.init_ptcl_posn()

vel_list = np.load('vel.npy')
num_vel_list = len(vel_list)
def vFunc():
    idx = random.randint(0, num_vel_list-1)
    return vel_list[idx]

ptcl.add_xFunc(xFunc)
ptcl.add_vFunc(vFunc)

# init chemistry
chem = CHEMISTRY()
chem.load(oper.fname)
chem.save(oper.fname)

# init reflection
rflct = REFLECT()

# init reaction
rct = REACT()

# init stats
if oper.idiag:
    stats = STATS(oper.fname)
    stats = MAIN(oper, ptcl, mesh, chem, rct, rflct, stats)
else:
    MAIN(oper, ptcl, mesh, chem, rct, rflct)

if oper.idiag:
    stats.df['Escaped Pct'] = stats.df['Escaped']/stats.df['Launched']
    stats.df['Etch Pct'] = stats.df['Etch']/stats.df['Launched']
    stats.df['Terminated Pct'] = stats.df['Terminated']/stats.df['Launched']
    stats.df.to_csv(fname + '_Stats.csv', index=True,
                   columns=['Launched', 'Escaped', 'Escaped Pct',
                            'Etch', 'Etch Pct',
                            'Terminated', 'Terminated Pct'],
                   float_format='%.2f',
                   na_rep='NA')
    
    for sp in stats.df.index:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), dpi=600,
                                   constrained_layout=True)
        
        ax = axes[0]
        ax.hist(stats.df.loc[sp, 'Init Erg'], bins=100, density=False)
        ax.set_title('Energy Distribution')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('a.u.')
        
        ax = axes[1]
        ax.hist(stats.df.loc[sp, 'Init Ang'], bins=100, density=False)
        ax.set_title('Velocity Distribution')
        ax.set_xlabel('Angle (degree)')
        ax.set_ylabel('a.u.')
        
        fig.suptitle(sp, fontsize=20)
        fig.savefig('init_' + sp + '.png', dpi=600)
        plt.close()
    
    n_sp = len(stats.df)
    fig, axes = plt.subplots(1, n_sp, figsize=(6*n_sp, 6), dpi=600,
                               constrained_layout=True)
    for i, sp in enumerate(stats.df.index):
        ax = axes[i]
        temp = stats.df.loc[sp, 'Reflection']
        # ax.hist(temp)
        # pd.Series(temp).value_counts().plot(kind='bar')
        cout = Counter(temp)
        ax.bar(cout.keys(), cout.values())
        ax.set_title(sp)
        ax.set_xlabel('Number of Reflections')
        ax.set_ylabel('Counts')
    fig.suptitle('Stats of Reflection', fontsize=20)
    fig.savefig('Reflection.png', dpi=600)
    plt.close()    
