
#ifndef _psi_src_bin_pypsi_pypsi_hh_
#define _psi_src_bin_pypsi_pypsi_hh_

#include <libmints/mints.h>
#include <libfock/jk.h>
#include <libfock/v.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <boost/shared_array.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libscf_solver/rhf.h>
#include <libscf_solver/ks.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <libscfgrad/scf_grad.h>


#include <boost/python/numeric.hpp>
#include <boost/python.hpp>


using namespace psi;

typedef boost::python::numeric::array NPArray;
typedef boost::shared_ptr<NPArray> SharedNPArray;

typedef boost::python::list PyList;
typedef boost::shared_ptr<PyList> SharedPyList;

typedef std::vector<SharedMatrix> VecSharedMat;

class PyPsi {
protected:

    Process::Environment process_environment_;
    boost::shared_ptr<LocalCommWrapper> worldcomm_;
    
    boost::shared_ptr<PSIO> psio_;
    
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
    boost::shared_ptr<MatrixFactory> matfac_;
    boost::shared_ptr<JK> jk_;
    boost::shared_ptr<DFJK> dfjk_;
    boost::shared_ptr<VBase> dftPotential_;
    boost::shared_ptr<scf::HF> wfn_;
    
    std::vector<SharedMatrix> guessOrbital_;
    
    // create psio object 
    void create_psio();
    
    // create basis object and one & two electron integral factories 
    void create_basis_and_integral_factories();
    
    // initialize wfn_ 
    void create_wfn();
    
    // exception function for DFJK utilities 
    void JK_DFException(std::string functionName);
    
    void Construct(const NPArray& xyz, const std::string& basis,
                   const int charge, const int multiplicity,
                   const std::string& path);
    
public:
    // constructors
    PyPsi(const NPArray& xyz, const std::string& basis);
    PyPsi(const NPArray& xyz, const std::string& basis, 
          const int charge, const int multiplicity);
    PyPsi(const NPArray& xyz, const std::string& basis, 
          const int charge, const int multiplicity, const std::string& path);
    
    // destructor 
    virtual ~PyPsi();
    
    // CPU and memory controll 
    int Settings_MaxNumCPUCores()
    {
        return process_environment_.get_n_threads();
    }
    double Settings_MaxMemoryInGB()
    {
        return (double)process_environment_.get_memory() / 1E+9;
    }
    std::string Settings_PsiDataDir()
    {
        return process_environment_("PSIDATADIR");
    }
    std::string Settings_TempDir()
    {
        return psio_->_psio_manager_->get_default_path();
    }
    void Settings_SetMaxNumCPUCores(const int ncores);
    void Settings_SetMaxMemory(const std::string&);
    void Settings_SetPsiDataDir(const std::string& path)
    {
        process_environment_.set("PSIDATADIR", path);
    }
    
    
    //*** Molecule properties 
    int Molecule_NumAtoms() // number of atoms 
    {
        return molecule_->natom();
    }
    int Molecule_NumElectrons(); // number of electrons 
    NPArray Molecule_Geometry(); // geometry in Bohr 
    double Molecule_NucRepEnergy() // nuclear repulsion energy 
    {
        return molecule_->nuclear_repulsion_energy();
    }
    NPArray Molecule_AtomicNumbers(); // atomic number list vector 
    NPArray Molecule_ChargeMult();
    
    
    //*** Basis set properties 
    const std::string BasisSet_Name() // basis set name string 
    {
        return basis_->name();
    }
    bool BasisSet_IsSpherical()
    {
        return basis_->has_puream();
    }
    int BasisSet_NumShells()
    {
        return basis_->nshell();
    }
    int BasisSet_NumFunctions() // number of basis functions 
    {
        return basis_->nbf();
    }
    NPArray BasisSet_ShellTypes();
    NPArray BasisSet_ShellNumPrimitives();
    NPArray BasisSet_ShellNumFunctions();
    NPArray BasisSet_ShellToCenter();
    NPArray BasisSet_FuncToCenter(); // basis function to its center atom
    NPArray BasisSet_FuncToShell();
    NPArray BasisSet_FuncToAngular(); // basis function to its angular momentum
    NPArray BasisSet_PrimExp();
    NPArray BasisSet_PrimCoeffUnnorm();
    
    
    //*** Integral package
    NPArray Integrals_Overlap(); // <i|j>
    NPArray Integrals_Kinetic(); // <i|T|j>
    NPArray Integrals_Potential(); // total, <i|sum(1/R)|j>
    PyList Integrals_PotentialEachCore(); // atom-separated, <i|1/R|j>
    NPArray Integrals_PotentialPtQ(NPArray&); // <i|1/R_p|j>
    PyList Integrals_Dipole(); // dipole matrices <i|x|j>, <i|y|j>, <i|z|j>
    int Integrals_NumUniqueTEIs(); // number of unique TEIs 
    double Integrals_ijkl(int, int, int, int); // (ij|kl)
    // ## HIGH MEMORY COST METHODS ## 
    void Integrals_AllUniqueTEIs(double*); // all unique TEIs in a vector
    void Integrals_AllTEIs(double*); // all (repetitive) TEIs in a 4D-array
    void Integrals_IndicesForK(double*, double*); // sorted TEI vecs for K
    // ## HIGH MEMORY COST METHODS ## 
    
    
    //*** JK related
    // use different types of JK 
    void JK_Initialize(std::string jktype,
                       std::string auxBasis = "CC-PVDZ-JKFIT");
    const std::string JK_Type();
    
    // methods computing J/K 
    PyList JK_DensToJ(PyList&);
    PyList JK_DensToK(PyList&);
    PyList JK_OccOrbToJ(PyList&);
    PyList JK_OccOrbToK(PyList&);
    void JK_CalcAllFromDens(PyList&);
    void JK_CalcAllFromOccOrb(PyList&);
    PyList JK_RetrieveJ();
    PyList JK_RetrieveK();
    
    // specially for density-fitting JK
    NPArray JK_DFTensor_AuxPriPairs();
    NPArray JK_DFMetric_InvJHalf();
    
    
    //*** DFT related
    void DFT_Initialize(std::string);
    PyList DFT_DensToV(PyList&);
    PyList DFT_OccOrbToV(PyList&);
    double DFT_EnergyXC();
    
    
    //*** SCF related
    // method of doing SCF calculations 
    void SCF_SetSCFType(std::string);
    void SCF_SetGuessOrb(PyList&);
    double SCF_RunSCF();
    
    // methods controlling SCF algorithm 
    void SCF_SetGuessType(const std::string&);
    
    // methods extracting SCF results
    double SCF_TotalEnergy();
    NPArray SCF_OrbitalAlpha();
    NPArray SCF_OrbitalBeta();
    NPArray SCF_OrbEigValAlpha();
    NPArray SCF_OrbEigValBeta();
    NPArray SCF_DensityAlpha();
    NPArray SCF_DensityBeta();
    NPArray SCF_FockAlpha();
    NPArray SCF_FockBeta();
    NPArray SCF_GuessDensity();
    
    NPArray SCF_Gradient();
    
    
    
};

#endif // _psi_src_bin_pypsi_pypsi_hh_
