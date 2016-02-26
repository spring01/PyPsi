#include <string>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/numeric.hpp>
#include "PyPsi.hh"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(JK_Initialize_overloads,
                                       JK_Initialize, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SCF_EnableMOM_overloads,
                                       SCF_EnableMOM, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SCF_EnableDamping_overloads,
                                       SCF_EnableDamping, 0, 1)

using namespace psi;

BOOST_PYTHON_MODULE(PyPsi)
{
    using namespace boost::python;
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    
    class_<PyPsi>("PyPsi", init<boost::python::numeric::array&, std::string>())
        .def(init<boost::python::numeric::array&, std::string, int, int>())
        .def(init<boost::python::numeric::array&, std::string, int, int, std::string>())
        .def("Settings_MaxNumCPUCores", &PyPsi::Settings_MaxNumCPUCores)
        .def("Settings_SetMaxNumCPUCores", &PyPsi::Settings_SetMaxNumCPUCores)
        .def("Settings_MaxMemoryInGB", &PyPsi::Settings_MaxMemoryInGB)
        .def("Settings_SetMaxMemory", &PyPsi::Settings_SetMaxMemory)
        .def("Settings_PsiDataDir", &PyPsi::Settings_PsiDataDir)
        .def("Settings_SetPsiDataDir", &PyPsi::Settings_SetPsiDataDir)
        .def("Settings_TempDir", &PyPsi::Settings_TempDir)
        
        .def("Molecule_NumAtoms", &PyPsi::Molecule_NumAtoms)
        .def("Molecule_NumElectrons", &PyPsi::Molecule_NumElectrons)
        .def("Molecule_Geometry", &PyPsi::Molecule_Geometry)
        .def("Molecule_NucRepEnergy", &PyPsi::Molecule_NucRepEnergy)
        .def("Molecule_AtomicNumbers", &PyPsi::Molecule_AtomicNumbers)
        .def("Molecule_ChargeMult", &PyPsi::Molecule_ChargeMult)
        
        .def("BasisSet_Name", &PyPsi::BasisSet_Name)
        .def("BasisSet_IsSpherical", &PyPsi::BasisSet_IsSpherical)
        .def("BasisSet_NumFunctions", &PyPsi::BasisSet_NumFunctions)
        .def("BasisSet_NumShells", &PyPsi::BasisSet_NumShells)
        .def("BasisSet_ShellTypes", &PyPsi::BasisSet_ShellTypes)
        .def("BasisSet_ShellNumPrimitives", &PyPsi::BasisSet_ShellNumPrimitives)
        .def("BasisSet_ShellNumFunctions", &PyPsi::BasisSet_ShellNumFunctions)
        .def("BasisSet_ShellToCenter", &PyPsi::BasisSet_ShellToCenter)
        .def("BasisSet_FuncToCenter", &PyPsi::BasisSet_FuncToCenter)
        .def("BasisSet_FuncToShell", &PyPsi::BasisSet_FuncToShell)
        .def("BasisSet_FuncToAngular", &PyPsi::BasisSet_FuncToAngular)
        .def("BasisSet_PrimExp", &PyPsi::BasisSet_PrimExp)
        .def("BasisSet_PrimCoeffUnnorm", &PyPsi::BasisSet_PrimCoeffUnnorm)
        
        .def("Integrals_Overlap", &PyPsi::Integrals_Overlap)
        .def("Integrals_Kinetic", &PyPsi::Integrals_Kinetic)
        .def("Integrals_Potential", &PyPsi::Integrals_Potential)
        .def("Integrals_PotentialEachCore", &PyPsi::Integrals_PotentialEachCore)
        .def("Integrals_PotentialPtQ", &PyPsi::Integrals_PotentialPtQ)
        .def("Integrals_Dipole", &PyPsi::Integrals_Dipole)
        .def("Integrals_ijkl", &PyPsi::Integrals_ijkl)
        
        .def("JK_Initialize", &PyPsi::JK_Initialize, JK_Initialize_overloads())
        .def("JK_Type", &PyPsi::JK_Type)
        .def("JK_DensToJ", &PyPsi::JK_DensToJ)
        .def("JK_DensToK", &PyPsi::JK_DensToK)
        .def("JK_OccOrbToJ", &PyPsi::JK_OccOrbToJ)
        .def("JK_OccOrbToK", &PyPsi::JK_OccOrbToK)
        .def("JK_CalcAllFromDens", &PyPsi::JK_CalcAllFromDens)
        .def("JK_CalcAllFromOccOrb", &PyPsi::JK_CalcAllFromOccOrb)
        .def("JK_RetrieveJ", &PyPsi::JK_RetrieveJ)
        .def("JK_RetrieveK", &PyPsi::JK_RetrieveK)
        .def("JK_DFTensor_AuxPriPairs", &PyPsi::JK_DFTensor_AuxPriPairs)
        .def("JK_DFTensor_AuxPriPri", &PyPsi::JK_DFTensor_AuxPriPri)
        .def("JK_DFMetric_InvJHalf", &PyPsi::JK_DFMetric_InvJHalf)
        
        .def("DFT_Initialize", &PyPsi::DFT_Initialize)
        .def("DFT_DensToV", &PyPsi::DFT_DensToV)
        .def("DFT_OccOrbToV", &PyPsi::DFT_OccOrbToV)
        .def("DFT_EnergyXC", &PyPsi::DFT_EnergyXC)
        
        .def("SCF_SetSCFType", &PyPsi::SCF_SetSCFType)
        .def("SCF_SetGuessOrb", &PyPsi::SCF_SetGuessOrb)
        .def("SCF_RunSCF", &PyPsi::SCF_RunSCF)
        .def("SCF_EnableMOM", &PyPsi::SCF_EnableMOM, SCF_EnableMOM_overloads())
        .def("SCF_EnableDamping", &PyPsi::SCF_EnableDamping, SCF_EnableDamping_overloads())
        .def("SCF_EnableDIIS", &PyPsi::SCF_EnableDIIS)
        .def("SCF_DisableDIIS", &PyPsi::SCF_DisableDIIS)
        .def("SCF_SetGuessType", &PyPsi::SCF_SetGuessType)
        .def("SCF_TotalEnergy", &PyPsi::SCF_TotalEnergy)
        .def("SCF_OrbitalAlpha", &PyPsi::SCF_OrbitalAlpha)
        .def("SCF_OrbitalBeta", &PyPsi::SCF_OrbitalBeta)
        .def("SCF_OrbEigValAlpha", &PyPsi::SCF_OrbEigValAlpha)
        .def("SCF_OrbEigValBeta", &PyPsi::SCF_OrbEigValBeta)
        .def("SCF_DensityAlpha", &PyPsi::SCF_DensityAlpha)
        .def("SCF_DensityBeta", &PyPsi::SCF_DensityBeta)
        .def("SCF_FockAlpha", &PyPsi::SCF_FockAlpha)
        .def("SCF_FockBeta", &PyPsi::SCF_FockBeta)
        .def("SCF_RHF_J", &PyPsi::SCF_RHF_J)
        .def("SCF_RHF_K", &PyPsi::SCF_RHF_K)
        .def("SCF_GuessDensity", &PyPsi::SCF_GuessDensity)
        
        .def("SCF_Gradient", &PyPsi::SCF_Gradient)
        ;
}
