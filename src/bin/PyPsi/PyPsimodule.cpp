
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
        .def("SCF_RunSCF", &PyPsi::SCF_RunSCF)
        ;
}
