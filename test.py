
import numpy as np
import PyPsi

cart = np.array([[9, 0.0, 0.0, 0.0], 
                 [1, 0.0, 0.0, 1.0]])
basisSet = "6-31g*"

path = "./PyPsi/"


pypsi = PyPsi.PyPsi(cart, basisSet)

print(pypsi.Molecule_NumAtoms())
print(pypsi.BasisSet_FuncToAngular())
print(pypsi.Integrals_Overlap())
print(pypsi.Integrals_Dipole())

