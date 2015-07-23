
import numpy as np
import PyPsi

cart = np.array([[9, 0.0, 0.1, 0.2], 
                 [1, 0.3, 0.4, 1.0]])


basisSet = "6-31g*"

path = "./PyPsi/"


pypsi = PyPsi.PyPsi(cart, basisSet)

occOrb = np.random.rand(17, 5) - 0.5
occOrb = occOrb / np.tile(np.sqrt(sum(occOrb**2)), [np.shape(occOrb)[0], 1])

print("Settings_MaxNumCPUCores", pypsi.Settings_MaxNumCPUCores())
print("Settings_SetMaxNumCPUCores", pypsi.Settings_SetMaxNumCPUCores(2), pypsi.Settings_MaxNumCPUCores())
print("Settings_MaxMemoryInGB", pypsi.Settings_MaxMemoryInGB())
print("Settings_SetMaxMemory", pypsi.Settings_SetMaxMemory("2gb"), pypsi.Settings_MaxMemoryInGB())
print("Settings_PsiDataDir", pypsi.Settings_PsiDataDir())
print("Settings_SetPsiDataDir", pypsi.Settings_SetPsiDataDir(path), pypsi.Settings_PsiDataDir())
print("Settings_TempDir", pypsi.Settings_TempDir())
pypsi.SCF_RunSCF()
print("Settings_TempDir after SCF_RunSCF", pypsi.Settings_TempDir())

print("Molecule_NumAtoms", pypsi.Molecule_NumAtoms())
print("Molecule_NumElectrons", pypsi.Molecule_NumElectrons())
print("Molecule_Geometry", pypsi.Molecule_Geometry())
geom = pypsi.Molecule_Geometry()
print("Molecule_SetGeometry", pypsi.Molecule_SetGeometry(geom + 1), pypsi.Molecule_Geometry())
print("Molecule_AtomicNumbers", pypsi.Molecule_AtomicNumbers())
print("Molecule_ChargeMult", pypsi.Molecule_ChargeMult())
print("Molecule_SetChargeMult", pypsi.Molecule_SetChargeMult(0, 3), pypsi.Molecule_ChargeMult())
print("Molecule_Fix", pypsi.Molecule_Fix())
print("Molecule_Free", pypsi.Molecule_Free())

#~ print(pypsi.Molecule_NumAtoms())
#~ print(pypsi.Molecule_NumAtoms())
#~ print(pypsi.Molecule_NumAtoms())
#~ print(pypsi.Molecule_NumAtoms())
#~ 
#~ print(pypsi.Molecule_NumAtoms())
#~ print(pypsi.BasisSet_FuncToAngular())
#~ 
#~ print(pypsi.Integrals_Dipole())
#~ print(pypsi.Integrals_Overlap())
#~ print(pypsi.Integrals_Potential())
#~ print(pypsi.Integrals_PotentialEachCore())
#~ 
#~ 
#~ print(occOrb)
#~ print(sum(occOrb**2))


