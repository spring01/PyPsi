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
#include <liboptions/liboptions.h>
#include <psi4-dec.h>
#include "superfunctional.h"
#include "functional.h"
#include <libdisp/dispersion.h>


namespace psi {

boost::shared_ptr<SuperFunctional> SuperFunctional::current(Options& options, int npoints, int deriv)
{
    if (npoints == -1)
        npoints = options.get_int("DFT_BLOCK_MAX_POINTS");

    boost::shared_ptr<SuperFunctional> super;
    super = SuperFunctional::build(options, npoints, deriv);
    if (options["DFT_OMEGA"].has_changed() && super->is_x_lrc())
        super->set_x_omega(options.get_double("DFT_OMEGA"));
    if (options["DFT_ALPHA"].has_changed())
        super->set_x_alpha(options.get_double("DFT_ALPHA"));

    if (npoints != super->max_points())
        super->set_max_points(npoints);
    if (deriv != super->deriv())
        super->set_deriv(deriv);

    return super;
}

boost::shared_ptr<Functional> build_b3_x_functional() {
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
    return b3;
}

boost::shared_ptr<Functional> build_s_x_functional() {
	boost::shared_ptr<Functional> s_x = Functional::build_base("S_X");
    s_x->set_name("S_X");
    s_x->set_description("    Slater LSDA Exchange\n");
    s_x->set_citation("    J.C. Slater, Phys. Rev., 81(3):385-390, 1951\n");
    s_x->set_gga(false);
    s_x->set_meta(false);
    s_x->set_alpha(1.0);
    s_x->set_omega(0.0);
    return s_x;
}

boost::shared_ptr<Functional> build_vwn3rpa_c_functional() {
	boost::shared_ptr<Functional> vwn3rpa = Functional::build_base("VWN3_C");
    vwn3rpa->set_name("VWN3RPA_C");
    vwn3rpa->set_description("    VWN3 (RPA) LSDA Correlation\n");
    vwn3rpa->set_citation("    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n");
    vwn3rpa->set_gga(false);
    vwn3rpa->set_meta(false);
    vwn3rpa->set_alpha(1.0);
    vwn3rpa->set_omega(0.0);
    vwn3rpa->set_parameter("EcP_2", -0.409286);
    vwn3rpa->set_parameter("EcP_3", 13.0720);
    vwn3rpa->set_parameter("EcP_4", 42.7198);
    vwn3rpa->set_parameter("EcF_2", -0.743294);
    vwn3rpa->set_parameter("EcF_3", 20.1231);
    vwn3rpa->set_parameter("EcF_4", 101.578);
    return vwn3rpa;
}

boost::shared_ptr<Functional> build_vwn5_c_functional() {
    boost::shared_ptr<Functional> vwn5 = Functional::build_base("VWN5_C");
    vwn5->set_name("VWN5_C");
    vwn5->set_description("    VWN5 LSDA Correlation\n");
    vwn5->set_citation("    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n");
    vwn5->set_gga(false);
    vwn5->set_meta(false);
    vwn5->set_alpha(1.0);
    vwn5->set_omega(0.0);
    vwn5->set_parameter("EcP_2", -0.10498);
    vwn5->set_parameter("EcP_3", 3.72744);
    vwn5->set_parameter("EcP_4", 12.9352);
    vwn5->set_parameter("EcF_2", -0.32500);
    vwn5->set_parameter("EcF_3", 7.06042);
    vwn5->set_parameter("EcF_4", 18.0578);
    vwn5->set_parameter("Ac_2", -0.00475840);
    vwn5->set_parameter("Ac_3", 1.13107);
    vwn5->set_parameter("Ac_4", 13.0045);
    return vwn5;
}

boost::shared_ptr<Functional> build_vwn5rpa_c_functional() {
    boost::shared_ptr<Functional> vwn5rpa = Functional::build_base("VWN5_C");
    vwn5rpa->set_name("VWN5RPA_C");
    vwn5rpa->set_description("    VWN5 (RPA) LSDA Correlation\n");
    vwn5rpa->set_citation("    S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200-1211, 1980\n");
    vwn5rpa->set_gga(false);
    vwn5rpa->set_meta(false);
    vwn5rpa->set_alpha(1.0);
    vwn5rpa->set_omega(0.0);
    vwn5rpa->set_parameter("EcP_2", -0.409286);
    vwn5rpa->set_parameter("EcP_3", 13.0720);
    vwn5rpa->set_parameter("EcP_4", 42.7198);
    vwn5rpa->set_parameter("EcF_2", -0.743294);
    vwn5rpa->set_parameter("EcF_3", 20.1231);
    vwn5rpa->set_parameter("EcF_4", 101.578);
    vwn5rpa->set_parameter("Ac_2", -0.228344);
    vwn5rpa->set_parameter("Ac_3", 1.06835);
    vwn5rpa->set_parameter("Ac_4", 11.4813);
    return vwn5rpa;
}

boost::shared_ptr<SuperFunctional> build_b3lyp_superfunctional(int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // add a B3LYP functional; added by spring
    super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);
    super->set_name("B3LYP");
    super->set_description("    B3LYP Hybrid-GGA Exchange-Correlation Functional\n");
    super->set_citation("    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n");
    
    // B3_X
    super->add_x_functional(build_b3_x_functional());
    
    // LYP_C
    boost::shared_ptr<Functional> lyp = Functional::build_base("LYP_C");
    lyp->set_name("LYP_C");
    lyp->set_alpha(0.81);
    
    // VWN3RPA_C
    boost::shared_ptr<Functional> vwn = build_vwn3rpa_c_functional();
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

boost::shared_ptr<SuperFunctional> build_b3lypv5_superfunctional(int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // add a B3LYP functional; added by spring
    super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);
    super->set_name("B3LYP");
    super->set_description("    B3LYP Hybrid-GGA Exchange-Correlation Functional\n");
    super->set_citation("    P.J. Stephens et. al., J. Phys. Chem., 98, 11623-11627, 1994\n");
    
    // B3_X
    super->add_x_functional(build_b3_x_functional());
    
    // LYP_C
    boost::shared_ptr<Functional> lyp = Functional::build_base("LYP_C");
    lyp->set_name("LYP_C");
    lyp->set_alpha(0.81);
    
    // VWN5
    boost::shared_ptr<Functional> vwn = build_vwn5_c_functional();
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

boost::shared_ptr<SuperFunctional> build_lsda_superfunctional(int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // add a SVWN3 functional; added by spring
    super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);
    super->set_name("SVWN3");
    super->set_description("    SVWN3 (RPA) LSDA Functional\n");
    super->set_citation("    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n");
    
    // Add member functionals
    super->add_x_functional(build_s_x_functional());
    super->add_c_functional(build_vwn3rpa_c_functional());

    // Set GKS up after adding functionals
    super->set_x_omega(0.0);
    super->set_c_omega(0.0);
    super->set_x_alpha(0.0);
    super->set_c_alpha(0.0);

    // => End User-Customization <= #

    // Call this last
    super->allocate();

    return super;
}

boost::shared_ptr<SuperFunctional> build_svwn5_superfunctional(int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // add a SVWN5 functional; added by spring
    super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);
    super->set_name("SVWN5");
    super->set_description("    SVWN5 Functional\n");
    super->set_citation("    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n");
    
    // Add member functionals
    super->add_x_functional(build_s_x_functional());
    super->add_c_functional(build_vwn5_c_functional());

    // Set GKS up after adding functionals
    super->set_x_omega(0.0);
    super->set_c_omega(0.0);
    super->set_x_alpha(0.0);
    super->set_c_alpha(0.0);

    // => End User-Customization <= #

    // Call this last
    super->allocate();

    return super;
}

boost::shared_ptr<SuperFunctional> build_wb97xd_superfunctional(Options& options, int max_points, int deriv)
{

    // Call this first
    boost::shared_ptr<SuperFunctional> super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);

    // => User-Customization <= //

    // No spaces, keep it short and according to convention
    super->set_name("wB97X-D");
    // Tab in, trailing newlines
    super->set_description("    Parameterized Hybrid LRC B97 GGA XC Functional with Dispersion\n");
    // Tab in, trailing newlines
    super->set_citation("    J.-D. Chai and M. Head-Gordon, Phys. Chem. Chem. Phys., 10, 6615-6620, 2008\n");

    // Add member functionals
    double alpha = 2.22036E-1;
    double omega = 0.2;
    boost::shared_ptr<Functional> wb97x_x = Functional::build_base("wB97_X");
    wb97x_x->set_name("wB97X_X");
    wb97x_x->set_alpha(1.0 / (1.0 - alpha));

    wb97x_x->set_parameter("B97_gamma", 0.004);
    wb97x_x->set_parameter("B97_a0", 7.77964E-1);     // Table 1: c_{x\sigma,0}
    wb97x_x->set_parameter("B97_a1", 6.61160E-1);     // Table 1: c_{x\sigma,1}
    wb97x_x->set_parameter("B97_a2", 5.74541E-1);     // Table 1: c_{x\sigma,2}
    wb97x_x->set_parameter("B97_a3", -5.25671E0);     // Table 1: c_{x\sigma,3}
    wb97x_x->set_parameter("B97_a4", 1.16386E1);      // Table 1: c_{x\sigma,4}

    boost::shared_ptr<Functional> wb97x_c = Functional::build_base("B_C");
    wb97x_c->set_name("wB97X_C");

    wb97x_c->set_parameter("B97_os_gamma", 0.006);
    wb97x_c->set_parameter("B97_os_a0", 1.0);         // Table 1: c_{c\alpha\beta,0}
    wb97x_c->set_parameter("B97_os_a1", 1.79413E0);   // Table 1: c_{c\alpha\beta,1}
    wb97x_c->set_parameter("B97_os_a2", -1.20477E1);  // Table 1: c_{c\alpha\beta,2}
    wb97x_c->set_parameter("B97_os_a3", 1.40847E1);   // Table 1: c_{c\alpha\beta,3}
    wb97x_c->set_parameter("B97_os_a4", -8.50809E0);  // Table 1: c_{c\alpha\beta,4}

    wb97x_c->set_parameter("B97_ss_gamma", 0.2);
    wb97x_c->set_parameter("B97_ss_a0", 1.0);         // Table 1: c_{c\sigma\sigma,0}
    wb97x_c->set_parameter("B97_ss_a1", -6.90539E0);  // Table 1: c_{c\sigma\sigma,1}
    wb97x_c->set_parameter("B97_ss_a2", 3.13343E1);   // Table 1: c_{c\sigma\sigma,2}
    wb97x_c->set_parameter("B97_ss_a3", -5.10533E1);  // Table 1: c_{c\sigma\sigma,3}
    wb97x_c->set_parameter("B97_ss_a4", 2.64423E1);   // Table 1: c_{c\sigma\sigma,4}

    super->add_x_functional(wb97x_x);
    super->add_c_functional(wb97x_c);

    // Set GKS up after adding functionals
    super->set_x_omega(omega);   // Table 1: omega
    super->set_c_omega(0.0);
    super->set_x_alpha(alpha);   // Table 1: c_x
    super->set_c_alpha(0.0);

    // => -D2 (CHG Damping Function) <= //
    super->set_dispersion(Dispersion::build(options, "-CHG", 1.0, 0.0, 0.0, 0.0));

    // => End User-Customization <= //

    // Call this last
    super->allocate();
    return super;
}

/// building blocks
boost::shared_ptr<SuperFunctional> build_s_x_superfunctional(Options& options, int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // add a S_X functional; added by spring
    super = SuperFunctional::blank();
    super->set_max_points(max_points);
    super->set_deriv(deriv);
    super->set_name("S_X");
    super->set_description("    S_X Functional\n");
    super->set_citation("    Adamson et. al., J. Comput. Chem., 20(9), 921-927, 1999\n");
    
    // Add member functionals
    super->add_x_functional(build_s_x_functional());

    // Call this last
    super->allocate();

    return super;
}

// superfunctional builder
boost::shared_ptr<SuperFunctional> SuperFunctional::build(Options& options, int max_points, int deriv)
{
    boost::shared_ptr<SuperFunctional> super;
    
    // build a superfunctional; added by spring
    std::string name = options.get_str("DFT_FUNCTIONAL");
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    if(!name.compare("b3lyp")) {
        super = build_b3lyp_superfunctional(max_points, deriv);
    } else if(!name.compare("b3lypv5")) {
        super = build_b3lypv5_superfunctional(max_points, deriv);
    } else if(!name.compare("lsda") || !name.compare("svwn") || !name.compare("svwn3")) {
        super = build_lsda_superfunctional(max_points, deriv);
	} else if(!name.compare("svwn5")) {
        super = build_svwn5_superfunctional(max_points, deriv);
	} else if(!name.compare("wb97xd")) {
        super = build_wb97xd_superfunctional(options, max_points, deriv);
	} else {
        throw PSIEXCEPTION("Superfunctional not supported yet.");
    }

    return super;
}

}
