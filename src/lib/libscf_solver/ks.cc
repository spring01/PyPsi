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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <libfock/jk.h>
#include <libfock/v.h>
#include <libfunctional/superfunctional.h>
#include <libdisp/dispersion.h>
#include <lib3index/3index.h>
#include "ks.h"

using namespace std;
using namespace psi;
using namespace boost;

namespace psi { namespace scf {

// added by spring
KS::KS(Process::Environment& process_environment_in) :
    process_environment_(process_environment_in), 
    options_(process_environment_in.options), 
    psio_(process_environment_in.psio())
{
    common_init();
}

KS::KS(Process::Environment& process_environment_in, Options & options, boost::shared_ptr<PSIO> psio) :
    process_environment_(process_environment_in), options_(options), psio_(psio)
{
    common_init();
}
KS::~KS()
{
}
void KS::common_init()
{
    // Take the molecule from the environment
    molecule_ = process_environment_.molecule();

    // Load in the basis set
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basisset_ = BasisSet::construct(process_environment_, parser, molecule_, "BASIS");
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    sobasisset_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, fact));

    potential_ = VBase::build_V(process_environment_, KS::options_,(options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));
    potential_->initialize();
    functional_ = potential_->functional();

    // Print the KS-specific stuff
    potential_->print_header();

}

// added by spring
RKS::RKS(Process::Environment& process_environment_in, boost::shared_ptr<JK> jk_in) :
    RHF(process_environment_in, jk_in), KS(process_environment_in)
{
    common_init();
}

RKS::RKS(Process::Environment& process_environment_in, Options & options, boost::shared_ptr<JK> jk_in, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    RHF(process_environment_in, options, jk_in, psio, chkpt), KS(process_environment_in, options,psio)
{
    common_init();
}
RKS::RKS(Process::Environment& process_environment_in, Options & options, boost::shared_ptr<JK> jk_in, boost::shared_ptr<PSIO> psio) :
    RHF(process_environment_in, options, jk_in, psio), KS(process_environment_in, options,psio)
{
    common_init();
}
void RKS::common_init()
{
    wK_ = factory_->create_shared_matrix("wKa (Long-Range Hartree-Fock Exchange)");
}
RKS::~RKS()
{
}
void RKS::integrals()
{
    HF::integrals();

    if (!functional_->is_x_lrc()) return;
           
    if (KS::options_.get_str("SCF_TYPE") == "DIRECT") {
    } else if (KS::options_.get_str("SCF_TYPE") == "DF") {
    } else if (KS::options_.get_str("SCF_TYPE") == "OUT_OF_CORE") {
    } else if (KS::options_.get_str("SCF_TYPE") == "PK") {
    } else {
        throw PSIEXCEPTION("SCF_TYPE is not supported by RC functionals");
    }
}
void RKS::form_V()
{
    // Push the C matrix on
    std::vector<SharedMatrix> & C = potential_->C();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    
    // Run the potential object
    potential_->compute();

    // Pull the V matrices off 
    const std::vector<SharedMatrix> & V = potential_->V();
    V_ = V[0];
}
void RKS::form_G()
{
    //~ timer_on("Form V");
    form_V();
    //~ timer_off("Form V");

    // Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    
    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();
    const std::vector<SharedMatrix> & wK = jk_->wK();
    J_ = J[0];
    J_->scale(2.0);
    if (functional_->is_x_hybrid()) {
        K_ = K[0];
    }
    if (functional_->is_x_lrc()) {
        wK_ = wK[0];
    }
    
    G_->copy(J_);
        
    G_->add(V_);

    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;

    if (alpha != 0.0) {
        K_->scale(alpha);
        G_->subtract(K_);
        K_->scale(1.0/alpha);
    } else {
        K_->zero();
    }

    if (functional_->is_x_lrc()) {
        wK_->scale(beta);
        G_->subtract(wK_);
        wK_->scale(1.0/beta);
    } else {
        wK_->zero();
    }

    if (debug_ > 2) {
        J_->print();
        K_->print();
        wK_->print();
        V_->print();
    }
}
double RKS::compute_E()
{
    // E_DFT = 2.0 D*H + 2.0 D*J - \alpha D*K + E_xc
    double one_electron_E = 2.0*D_->vector_dot(H_);
    double coulomb_E = D_->vector_dot(J_);
    
    std::map<std::string, double>& quad = potential_->quadrature_values();  
    double XC_E = quad["FUNCTIONAL"];
    double exchange_E = 0.0;
    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;
    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha*Da_->vector_dot(K_);
    }
    if (functional_->is_x_lrc()) {
        exchange_E -=  beta*Da_->vector_dot(wK_);
    }

    double dashD_E = 0.0;
    boost::shared_ptr<Dispersion> disp = functional_->dispersion();
    if (disp) {
        dashD_E = disp->compute_energy(HF::molecule_);
    }

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += coulomb_E;
    Etotal += exchange_E;
    Etotal += XC_E; 
    Etotal += dashD_E;

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = coulomb_E + exchange_E;
    energies_["XC"] = XC_E;
    energies_["-D"] = dashD_E;

    if (debug_) {
        fprintf(outfile, "   => Energetics <=\n\n");
        fprintf(outfile, "    Nuclear Repulsion Energy = %24.14f\n", nuclearrep_);
        fprintf(outfile, "    One-Electron Energy =      %24.14f\n", one_electron_E);
        fprintf(outfile, "    Coulomb Energy =           %24.14f\n", coulomb_E);
        fprintf(outfile, "    Hybrid Exchange Energy =   %24.14f\n", exchange_E);
        fprintf(outfile, "    XC Functional Energy =     %24.14f\n", XC_E); 
        fprintf(outfile, "    -D Energy =                %24.14f\n\n", dashD_E);
    }

    return Etotal;
}
void RKS::finalize()
{
    RHF::finalize();
}
void RKS::stability_analysis()
{
    throw PSIEXCEPTION("DFT stabilty analysis has not been implemented yet.  Sorry :(");
}

// added by spring
UKS::UKS(Process::Environment& process_environment_in, boost::shared_ptr<JK> jk_in) :
    UHF(process_environment_in, jk_in), KS(process_environment_in)
{
    common_init();
}

UKS::UKS(Process::Environment& process_environment_in, Options & options, boost::shared_ptr<JK> jk_in, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    UHF(process_environment_in, options, jk_in, psio, chkpt), KS(process_environment_in, options,psio)
{
    common_init();
}
UKS::UKS(Process::Environment& process_environment_in, Options & options, boost::shared_ptr<JK> jk_in, boost::shared_ptr<PSIO> psio) :
    UHF(process_environment_in, options, jk_in, psio), KS(process_environment_in, options,psio)
{
    common_init();
}
void UKS::common_init()
{
    wKa_ = factory_->create_shared_matrix("wKa (Long-Range Hartree-Fock Exchange)");
    wKb_ = factory_->create_shared_matrix("wKb (Long-Range Hartree-Fock Exchange)");
}
UKS::~UKS()
{
}
void UKS::integrals()
{
    HF::integrals();

    if (!functional_->is_x_lrc()) return;
           
    if (KS::options_.get_str("SCF_TYPE") == "DIRECT") {
    } else if (KS::options_.get_str("SCF_TYPE") == "DF") {
    } else if (KS::options_.get_str("SCF_TYPE") == "OUT_OF_CORE") {
    } else if (KS::options_.get_str("SCF_TYPE") == "PK") {
    } else {
        throw PSIEXCEPTION("SCF_TYPE is not supported by RC functionals");
    }    
}
void UKS::form_V()
{
    // Push the C matrix on
    std::vector<SharedMatrix> & C = potential_->C();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));
    
    // Run the potential object
    potential_->compute();

    // Pull the V matrices off 
    const std::vector<SharedMatrix> & V = potential_->V();
    Va_ = V[0];
    Vb_ = V[1];
}
void UKS::form_G()
{
    //~ timer_on("Form V");
    form_V();
    //~ timer_off("Form V");

    // Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));
    
    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();
    const std::vector<SharedMatrix> & wK = jk_->wK();
    J_->copy(J[0]);
    J_->add(J[1]);
    if (functional_->is_x_hybrid()) {
        Ka_ = K[0];
        Kb_ = K[1];
    }
    if (functional_->is_x_lrc()) {
        wKa_ = wK[0];
        wKb_ = wK[1];
    }
    Ga_->copy(J_);
    Gb_->copy(J_);
        
    Ga_->add(Va_);
    Gb_->add(Vb_);

    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;
    if (alpha != 0.0) {
        Ka_->scale(alpha);
        Kb_->scale(alpha);
        Ga_->subtract(Ka_);
        Gb_->subtract(Kb_);
        Ka_->scale(1.0/alpha);
        Kb_->scale(1.0/alpha);
    } else {
        Ka_->zero();
        Kb_->zero();
    }

    if (functional_->is_x_lrc()) {
        wKa_->scale(beta);
        wKb_->scale(beta);
        Ga_->subtract(wKa_);
        Gb_->subtract(wKb_);
        wKa_->scale(1.0/beta);
        wKb_->scale(1.0/beta);
    } else {
        wKa_->zero();
        wKb_->zero();
    }

    if (debug_ > 2) {
        J_->print(outfile);
        Ka_->print(outfile);
        Kb_->print(outfile);
        wKa_->print(outfile);
        wKb_->print(outfile);
        Va_->print();
        Vb_->print();
    }
}
double UKS::compute_E()
{
    // E_DFT = 2.0 D*H + D*J - \alpha D*K + E_xc
    double one_electron_E = Da_->vector_dot(H_);
    one_electron_E += Db_->vector_dot(H_);
    double coulomb_E = Da_->vector_dot(J_);
    coulomb_E += Db_->vector_dot(J_);

    std::map<std::string, double>& quad = potential_->quadrature_values();  
    double XC_E = quad["FUNCTIONAL"];
    double exchange_E = 0.0;
    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;
    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha*Da_->vector_dot(Ka_);
        exchange_E -= alpha*Db_->vector_dot(Kb_);
    }
    if (functional_->is_x_lrc()) {
        exchange_E -=  beta*Da_->vector_dot(wKa_);
        exchange_E -=  beta*Db_->vector_dot(wKb_);
    }

    double dashD_E = 0.0;
    boost::shared_ptr<Dispersion> disp = functional_->dispersion();
    if (disp) {
        dashD_E = disp->compute_energy(HF::molecule_);
    }

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5 * coulomb_E;
    Etotal += 0.5 * exchange_E;
    Etotal += XC_E; 
    Etotal += dashD_E;

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = 0.5 * (coulomb_E + exchange_E);
    energies_["XC"] = XC_E;
    energies_["-D"] = dashD_E;

    if (debug_) {
        fprintf(outfile, "   => Energetics <=\n\n");
        fprintf(outfile, "    Nuclear Repulsion Energy = %24.14f\n", nuclearrep_);
        fprintf(outfile, "    One-Electron Energy =      %24.14f\n", one_electron_E);
        fprintf(outfile, "    Coulomb Energy =           %24.14f\n", 0.5 * coulomb_E);
        fprintf(outfile, "    Hybrid Exchange Energy =   %24.14f\n", 0.5 * exchange_E);
        fprintf(outfile, "    XC Functional Energy =     %24.14f\n", XC_E); 
        fprintf(outfile, "    -D Energy =                %24.14f\n", dashD_E);
    }

    return Etotal;
}
void UKS::finalize()
{
    UHF::finalize();
}
void UKS::stability_analysis()
{
    throw PSIEXCEPTION("DFT stabilty analysis has not been implemented yet.  Sorry :(");
}

}}
