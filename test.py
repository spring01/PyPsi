
import numpy as np
import PyPsi

cart = np.array([[9, 0.0, 0.0, 0.0], 
                 [1, 0.0, 0.0, 1.0]])

#~ cart = np.array([[9, 0.0, 0.0, 0.0]])

basisSet = "6-31g*"

path = "./PyPsi/"


pypsi = PyPsi.PyPsi(cart, basisSet)

occOrb = np.random.rand(17, 5) - 0.5
occOrb = occOrb / np.tile(np.sqrt(sum(occOrb**2)), [np.shape(occOrb)[0], 1])

print(pypsi.Molecule_NumAtoms())
print(pypsi.BasisSet_FuncToAngular())

print(pypsi.Integrals_Dipole())
print(pypsi.Integrals_Overlap())
print(pypsi.Integrals_Potential())
print(pypsi.Integrals_PotentialEachCore())


print(occOrb)
print(sum(occOrb**2))


