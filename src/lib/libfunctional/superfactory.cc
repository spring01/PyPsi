/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <boost/shared_ptr.hpp>
//~ #include <boost/python.hpp>
//~ #include <boost/python/object.hpp>
#include <liboptions/liboptions.h>
//~ #include <liboptions/liboptions_python.h>
#include <psi4-dec.h>
#include "superfunctional.h"
#include "functional.h"

//~ #define PY_TRY(ptr, command)  \
     //~ if(!(ptr = command)){    \
         //~ PyErr_Print();       \
         //~ exit(1);             \
     //~ }

//~ using namespace boost::python;

namespace psi {

boost::shared_ptr<SuperFunctional> SuperFunctional::current(Options& options, int npoints, int deriv)
{
    //~ throw PSIEXCEPTION("SuperFunctional::current: Python elimination causes a problem.");
    if (npoints == -1) {
        npoints = options.get_int("DFT_BLOCK_MAX_POINTS");
    }

    boost::shared_ptr<SuperFunctional> super;
    //~ if (options.get_str("DFT_FUNCTIONAL") == "GEN" || options.get_str("DFT_FUNCTIONAL") == "") {
        //~ boost::python::object pySuper = dynamic_cast<PythonDataType*>(options["DFT_CUSTOM_FUNCTIONAL"].get())->to_python();
        //~ super = boost::python::extract<boost::shared_ptr<SuperFunctional> >(pySuper);
        //~ if (!super) {
            //~ throw PSIEXCEPTION("Custom Functional requested, but nothing provided in DFT_CUSTOM_FUNCTIONAL");
        //~ }
    //~ } else {
        super = SuperFunctional::build(options.get_str("DFT_FUNCTIONAL"), npoints, deriv);
        if (options["DFT_OMEGA"].has_changed() && super->is_x_lrc()) {
            super->set_x_omega(options.get_double("DFT_OMEGA"));
        }
        if (options["DFT_ALPHA"].has_changed()) {
            super->set_x_alpha(options.get_double("DFT_ALPHA"));
        }
    //~ }

    if (npoints != super->max_points())
        super->set_max_points(npoints);
    if (deriv != super->deriv())
        super->set_deriv(deriv);

    return super;
}

boost::shared_ptr<SuperFunctional> build_b3lyp_superfunctional(int max_points, int deriv);

boost::shared_ptr<SuperFunctional> SuperFunctional::build(const std::string& alias, int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;

    //~ if (Py_IsInitialized()) {
        //~ try {
            //~ // Grab the SuperFunctional off of the Python plane
            //~ PyObject *functional;
            //~ PY_TRY(functional, PyImport_ImportModule("functional") );
            //~ PyObject *function;
            //~ PY_TRY(function, PyObject_GetAttrString(functional, "build_superfunctional"));
            //~ PyObject *pargs;
            //~ PY_TRY(pargs, Py_BuildValue("(sii)", alias.c_str(), max_points, deriv));
            //~ PyObject *ret;
            //~ PY_TRY(ret, PyEval_CallObject(function, pargs));
//~ 
            //~ // Extract the SuperFunctional
            //~ super = boost::python::extract<boost::shared_ptr<SuperFunctional> >(ret);
//~ 
            //~ // Decref Python env pointers
            //~ Py_DECREF(ret);
            //~ Py_DECREF(pargs);
            //~ Py_DECREF(function);
            //~ Py_DECREF(functional);
        //~ }
        //~ catch (error_already_set const& e)
        //~ {
            //~ PyErr_Print();
            //~ exit(1);
        //~ }
    //~ }
    //~ else {
        //~ throw PSIEXCEPTION("Unable to parse superfunctional.\n");
    //~ }
    
    // add a B3LYP test functional; added by spring
    std::string name = alias;
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    if(!name.compare("b3lyp")) {
        super = build_b3lyp_superfunctional(max_points, deriv);
    } else {
        throw PSIEXCEPTION("Superfunctional not supported yet.");
    }

    return super;
}

boost::shared_ptr<SuperFunctional> build_b3lyp_superfunctional(int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // add a B3LYP test functional; added by spring
    super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);
    super->set_name("B3LYP");
    super->set_description("    B3LYP Hybrid-GGA Exchange-Correlation Functional\n");
    super->set_citation("    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n");
    
    // B3_X
    boost::shared_ptr<Functional> b3 = Functional::build_base("B88_X");
    b3->set_name("B3_X");
    b3->set_description("    Becke88 GGA Exchange (B3LYP weighting)\n");
    b3->set_citation("    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n");
    b3->set_gga(true);
    b3->set_meta(false);
    b3->set_alpha(0.8);
    b3->set_omega(0.0);
    b3->set_parameter("B88_d", 0.0042);
    b3->set_parameter("B88_a", 0.9000);
    b3->set_alpha(1.0);
    super->add_x_functional(b3);
    
    // LYP_C
    boost::shared_ptr<Functional> lyp = Functional::build_base("LYP_C");
    lyp->set_name("LYP_C");
    lyp->set_alpha(0.81);
    
    // VWN3RPA_C
    boost::shared_ptr<Functional> vwn = Functional::build_base("VWN3_C");
    vwn->set_name("VWN3RPA_C");
    vwn->set_description("    VWN3 (RPA) LSDA Correlation\n");
    vwn->set_citation("    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n");
    vwn->set_gga(false);
    vwn->set_meta(false);
    vwn->set_alpha(1.0);
    vwn->set_omega(0.0);
    vwn->set_parameter("EcP_2", -0.409286);
    vwn->set_parameter("EcP_3", 13.0720);
    vwn->set_parameter("EcP_4", 42.7198);
    vwn->set_parameter("EcF_2", -0.743294);
    vwn->set_parameter("EcF_3", 20.1231);
    vwn->set_parameter("EcF_4", 101.578);
    vwn->set_alpha(0.19);
    super->add_c_functional(vwn);
    super->add_c_functional(lyp);
    
    // Set GKS up after adding functionals
    super->set_x_omega(0.0);
    super->set_c_omega(0.0);
    super->set_x_alpha(0.2);
    super->set_c_alpha(0.0);
    super->allocate();

    return super;
}


}
