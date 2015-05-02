% A simple test script creates a MatPsi object named matpsi 

cart = [ ...
    8 0.0 0.0 0.0;
    1 0.0 0.0 1.0;
    1 0.0 1.0 0.0;];
matpsi = MatPsi2(cart, '6-31g*');

% Settings
matpsi.Settings_SetMaxNumCPUCores(4);
matpsi.Settings_MaxNumCPUCores();
matpsi.Settings_SetMaxMemory('4gb');
matpsi.Settings_MaxMemoryInGB();
matpsi.Settings_SetPsiDataDir('./@MatPsi2');
matpsi.Settings_PsiDataDir();
matpsi.Settings_TempDir();

% Molecule
matpsi.Molecule_Fix();
matpsi.Molecule_Free();
matpsi.Molecule_NumAtoms();
matpsi.Molecule_NumElectrons();
geom = matpsi.Molecule_Geometry();
matpsi.Molecule_SetGeometry(geom + 1);
atomNums = matpsi.Molecule_AtomicNumbers();
matpsi.Molecule_NucRepEnergy();
matpsi.Molecule_SetChargeMult(0,1);

% BasisSet
matpsi.BasisSet_Name();
matpsi.BasisSet_SetBasisSet('6-31g')
matpsi.BasisSet_IsSpherical();
matpsi.BasisSet_NumFunctions();
matpsi.BasisSet_NumShells();
matpsi.BasisSet_ShellTypes();
matpsi.BasisSet_ShellNumPrimitives();
matpsi.BasisSet_ShellNumFunctions();
matpsi.BasisSet_ShellToCenter();
matpsi.BasisSet_FuncToCenter();
matpsi.BasisSet_FuncToShell();
matpsi.BasisSet_FuncToAngular();
matpsi.BasisSet_PrimExp();
matpsi.BasisSet_PrimCoeffUnnorm();

% Integrals
testMat = matpsi.Integrals_Overlap();
matpsi.Integrals_Kinetic();
matpsi.Integrals_Potential();
matpsi.Integrals_PotentialEachCore();
matpsi.Integrals_PotentialPtQ([atomNums' geom]);
matpsi.Integrals_Dipole();
matpsi.Integrals_ijkl(3,4,5,5);
matpsi.Integrals_NumUniqueTEIs();
matpsi.Integrals_AllUniqueTEIs();
matpsi.Integrals_AllTEIs();
matpsi.Integrals_IndicesForK();

% JK
matpsi.JK_Initialize('PKJK');
matpsi.JK_Type();
matpsi.JK_DensToJ(testMat);
matpsi.JK_DensToK(testMat);
matpsi.JK_OccOrbToJ(testMat);
matpsi.JK_OccOrbToK(testMat);
matpsi.JK_Initialize('DFJK');
matpsi.JK_DFTensor_AuxPriPairs();
matpsi.JK_DFTensor_AuxPriPri();
matpsi.JK_DFMetric_InvJHalf();
matpsi.JK_Initialize('icjk');
matpsi.JK_Initialize('directjk');
matpsi.JK_Initialize('pkJK');

% SCF
matpsi.SCF_RunSCF();
orb = matpsi.SCF_OrbitalAlpha();
matpsi.SCF_SetSCFType('rks');
matpsi.SCF_RunSCF();
matpsi.SCF_SetSCFType('uhf');
matpsi.SCF_RunSCF();
matpsi.SCF_SetSCFType('uks');
matpsi.SCF_RunSCF();
matpsi.Molecule_SetChargeMult(1,2);
matpsi.Molecule_ChargeMult();
matpsi.SCF_SetSCFType('uhf');
matpsi.SCF_RunSCF();
matpsi.SCF_SetSCFType('uks');
matpsi.SCF_RunSCF();
matpsi.Molecule_SetChargeMult(0,1);
matpsi.SCF_SetSCFType('rhf');
matpsi.SCF_SetGuessOrb(orb);
matpsi.SCF_RunSCF();
matpsi.SCF_EnableMOM();
matpsi.SCF_DisableMOM();
matpsi.SCF_EnableDamping();
matpsi.SCF_DisableDamping();
matpsi.SCF_EnableDIIS();
matpsi.SCF_DisableDIIS();
matpsi.SCF_GuessSAD();
matpsi.SCF_GuessCore();
matpsi.SCF_TotalEnergy();
matpsi.SCF_OrbitalAlpha();
matpsi.SCF_OrbitalBeta();
matpsi.SCF_OrbEigValAlpha();
matpsi.SCF_OrbEigValBeta();
matpsi.SCF_DensityAlpha();
matpsi.SCF_DensityBeta();
matpsi.SCF_CoreHamiltonian();
matpsi.SCF_FockAlpha();
matpsi.SCF_FockBeta();
matpsi.SCF_Gradient();
matpsi.SCF_GuessDensity();
matpsi.SCF_RHF_J();
matpsi.SCF_RHF_K();


