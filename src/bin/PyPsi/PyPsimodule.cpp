
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
        .def("Molecule_NumAtoms", &PyPsi::Molecule_NumAtoms)
        .def("BasisSet_FuncToAngular", &PyPsi::BasisSet_FuncToAngular)
        .def("BasisSet_Name", &PyPsi::BasisSet_Name)
        .def("BasisSet_NumFunctions", &PyPsi::BasisSet_NumFunctions)
        .def("Integrals_Dipole", &PyPsi::Integrals_Dipole)
        .def("Integrals_Overlap", &PyPsi::Integrals_Overlap)
        .def("Integrals_Potential", &PyPsi::Integrals_Potential)
        .def("Integrals_PotentialEachCore", &PyPsi::Integrals_PotentialEachCore)
        .def("SCF_RunSCF", &PyPsi::SCF_RunSCF)
        ;
}
