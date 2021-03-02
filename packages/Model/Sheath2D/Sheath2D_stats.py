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
        self.col = ['Init_Vel', 'Init_Erg', 'Init_Ang',
                    'End_Vel', 'End_Erg', 'End_Ang',
                    'Collision', 'hitWafer']
        self.df = pd.DataFrame(columns=self.col)
    
    def append_row(self):
        """Append an empty row."""
        row = [np.zeros(3), 0.0, 0.0,
               np.zeros(3), 0.0, 0.0,
               0, False]
        self.df.append(row, columns=self.col)