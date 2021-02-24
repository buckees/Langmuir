"""Feature Model 2D. Main program."""


from packages.Model.Feature2D.Feature2D_chem import CHEMISTRY

chem = CHEMISTRY('chemistry')
chem.load('Ar_Cl2_v01')
chem.save('Ar_Cl2_v01')
chosen_rct = chem.choose_react('Cl+', 'Si_', 100.0, 0.0)
print(chosen_rct)