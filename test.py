
import numpy as np
import PyPsi

cart = np.array([[9, 1.0, 2.0, 3.0], 
                 [1, 4.0, 5.0, 6.0]])
basisSet = "6-31g*"

path = "/home/haichen/working/PyPsi/lib"

dims = cart.shape

print(cart)
print(type(cart))
print(type(dims))
print(cart.shape[0])


pypsi = PyPsi.PyPsi(cart, basisSet, 0, 1, path)
print(pypsi.Molecule_NumAtoms())

