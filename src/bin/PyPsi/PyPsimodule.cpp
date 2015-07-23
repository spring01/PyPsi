
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/numeric.hpp>

#include <string>

#include "PyPsi.hh"

using namespace psi;

BOOST_PYTHON_MODULE(PyPsi)
{
    using namespace boost::python;
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
    
    class_<PyPsi>("PyPsi", init<boost::python::numeric::array&, std::string, int, int, std::string>())
        .def(init<boost::python::numeric::array&, std::string, optional<int, int> >())
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
        .def("Molecule_SetGeometry", &PyPsi::Molecule_SetGeometry)
        .def("Molecule_NucRepEnergy", &PyPsi::Molecule_NucRepEnergy)
        .def("Molecule_AtomicNumbers", &PyPsi::Molecule_AtomicNumbers)
        .def("Molecule_ChargeMult", &PyPsi::Molecule_ChargeMult)
        .def("Molecule_SetChargeMult", &PyPsi::Molecule_SetChargeMult)
        .def("Molecule_Fix", &PyPsi::Molecule_Fix)
        .def("Molecule_Free", &PyPsi::Molecule_Free)
        
        .def("BasisSet_Name", &PyPsi::BasisSet_Name)
        .def("BasisSet_SetBasisSet", &PyPsi::BasisSet_SetBasisSet)
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
        
        .def("SCF_RunSCF", &PyPsi::SCF_RunSCF)
        ;
}
