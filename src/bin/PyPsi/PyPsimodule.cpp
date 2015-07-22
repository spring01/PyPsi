
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <iostream>
#include <string>

//~ #include "PyPsi.cc"
#include "PyPsi.hh"

using namespace psi;

BOOST_PYTHON_MODULE(PyPsi)
{
    using namespace boost::python;
    class_<PyPsi>("PyPsi", init<SharedMatrix, std::string, int, int, std::string>())
        // Add a regular member function.
        .def("Molecule_NumAtoms", &PyPsi::Molecule_NumAtoms)
        ;
}
