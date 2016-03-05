
import numpy as np
import PyPsi

xyz = np.array([[9, 0.0, 0.1, 0.2], 
                 [1, 0.3, 0.4, 1.0]])

basis = "6-31g*"
path = "./PyPsi/"

pypsi = PyPsi.PyPsi(xyz, basis)
pypsi3 = PyPsi.PyPsi(xyz, basis, 0, 3)

occOrb = np.random.rand(17, 5) - 0.5
occOrb = occOrb / np.tile(np.sqrt(sum(occOrb**2)), [np.shape(occOrb)[0], 1])
dens = np.dot(occOrb, occOrb.transpose())

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
print("Molecule_AtomicNumbers", pypsi.Molecule_AtomicNumbers())
atomNums = pypsi.Molecule_AtomicNumbers()
print("Molecule_ChargeMult", pypsi.Molecule_ChargeMult())

print("BasisSet_Name", pypsi.BasisSet_Name())
print("BasisSet_IsSpherical", pypsi.BasisSet_IsSpherical())
print("BasisSet_NumShells", pypsi.BasisSet_NumShells())
print("BasisSet_NumFunctions", pypsi.BasisSet_NumFunctions())
print("BasisSet_ShellTypes", pypsi.BasisSet_ShellTypes())
print("BasisSet_ShellNumPrimitives", pypsi.BasisSet_ShellNumPrimitives())
print("BasisSet_ShellNumFunctions", pypsi.BasisSet_ShellNumFunctions())
print("BasisSet_ShellToCenter", pypsi.BasisSet_ShellToCenter())
print("BasisSet_FuncToCenter", pypsi.BasisSet_FuncToCenter())
print("BasisSet_FuncToShell", pypsi.BasisSet_FuncToShell())
print("BasisSet_FuncToAngular", pypsi.BasisSet_FuncToAngular())
print("BasisSet_PrimExp", pypsi.BasisSet_PrimExp())
print("BasisSet_PrimCoeffUnnorm", pypsi.BasisSet_PrimCoeffUnnorm())

print("Integrals_Overlap", pypsi.Integrals_Overlap())
print("Integrals_Kinetic", pypsi.Integrals_Kinetic())
print("Integrals_Potential", pypsi.Integrals_Potential())
print("Integrals_PotentialEachCore", pypsi.Integrals_PotentialEachCore())
print("Integrals_PotentialPtQ", pypsi.Integrals_PotentialPtQ(np.concatenate((np.reshape(atomNums, [np.size(atomNums), -1]), geom), axis=1)))
print("Integrals_Dipole", pypsi.Integrals_Dipole())
print("Integrals_ijkl", pypsi.Integrals_ijkl(0, 1, 2, 2))

print("JK_Initialize", pypsi.JK_Initialize("pkjk"))
print("JK_Type", pypsi.JK_Type())
print("JK_CalcAllFromOccOrb", pypsi.JK_CalcAllFromOccOrb([occOrb]))
print("JK_CalcAllFromOccOrb", pypsi.JK_CalcAllFromOccOrb([occOrb, occOrb]))
print("JK_CalcAllFromDens", pypsi.JK_CalcAllFromDens([dens]))
print("JK_CalcAllFromDens", pypsi.JK_CalcAllFromDens([dens, dens]))
print("JK_GetJ", pypsi.JK_GetJ())
print("JK_GetK", pypsi.JK_GetK())
print("JK_OccOrbToJ", pypsi.JK_OccOrbToJ([occOrb]))
print("JK_OccOrbToK", pypsi.JK_OccOrbToK([occOrb, occOrb]))
print("JK_DensToJ", pypsi.JK_DensToJ([dens]))
print("JK_DensToK", pypsi.JK_DensToK([dens, dens]))
pypsi.JK_Initialize("dfjk")
print("JK_DFTensor_AuxPriPairs", pypsi.JK_DFTensor_AuxPriPairs())
print("JK_DFMetric_InvJHalf", pypsi.JK_DFMetric_InvJHalf())
pypsi.JK_Initialize("pkjk")

print("DFT_Initialize", pypsi.DFT_Initialize("svwn5"))
print("DFT_DensToV", pypsi.DFT_DensToV([dens]))
print("DFT_DensToV", pypsi3.DFT_DensToV([dens, dens]))
print("DFT_OccOrbToV", pypsi.DFT_OccOrbToV([occOrb]))
print("DFT_OccOrbToV", pypsi3.DFT_OccOrbToV([occOrb, occOrb]))
print("DFT_EnergyXC", pypsi.DFT_EnergyXC())

allOrb = pypsi.SCF_OrbitalAlpha()
print("SCF_SetSCFType", pypsi.SCF_SetSCFType("rks"))
print("SCF_SetGuessOrb", pypsi.SCF_SetGuessOrb([allOrb]))
print("SCF_RunSCF", pypsi.SCF_RunSCF())
print("SCF_SetGuessType", pypsi.SCF_SetGuessType("SAD"))
print("SCF_SetGuessType", pypsi.SCF_SetGuessType("CORE"))
print("SCF_TotalEnergy", pypsi.SCF_TotalEnergy())
print("SCF_OrbitalAlpha", pypsi.SCF_OrbitalAlpha())
print("SCF_OrbitalBeta", pypsi.SCF_OrbitalBeta())
print("SCF_OrbEigValAlpha", pypsi.SCF_OrbEigValAlpha())
print("SCF_OrbEigValBeta", pypsi.SCF_OrbEigValBeta())
print("SCF_DensityAlpha", pypsi.SCF_DensityAlpha())
print("SCF_DensityBeta", pypsi.SCF_DensityBeta())
print("SCF_FockAlpha", pypsi.SCF_FockAlpha())
print("SCF_FockBeta", pypsi.SCF_FockBeta())
print("SCF_GuessDensity", pypsi.SCF_GuessDensity())
print("SCF_RunSCF", pypsi.SCF_RunSCF())

print("SCF_Gradient", pypsi.SCF_Gradient())


