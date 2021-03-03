"""Diagnose the sheath model."""

import numpy as np
import pandas as pd

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
                    'Collision', 'hitWafer']
        self.df = pd.DataFrame(columns=self.col)
    
    def append_row(self):
        """Append an empty row."""
        row = [0.0 for i in range(10)]
        row.append([0, False])
        self.df.loc[len(self.df)] = row