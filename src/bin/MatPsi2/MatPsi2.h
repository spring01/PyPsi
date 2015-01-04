#include <libmints/mints.h>
#include <libfock/jk.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <boost/shared_array.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libscf_solver/rhf.h>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>



using namespace std;
using namespace psi;
using namespace boost;

class MatPsi2 {
protected:
    
    std::string molstring_;
    std::string basisname_;

    Process::Environment process_environment_;
    boost::shared_ptr<worldcomm> worldcomm_;
    
    boost::shared_ptr<PSIO> psio_;
    
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
    boost::shared_ptr<MatrixFactory> matfac_;
    boost::shared_ptr<JK> jk_;
    boost::shared_ptr<scf::RHF> rhf_;
    
    // create basis object 
    void create_basis();
    
    // create basis object and one & two electron integral factories 
    void create_basis_and_integral_factories();
    
    // exception function for DFJK utilities
    void DFJKException(std::string functionName);
    
public:
    // constructor; takes in 2 strings and parse them 
    MatPsi2(const std::string& molstring, const std::string& basisname, 
        int charge, int multiplicity, const std::string& path);
    
    // destructor 
    virtual ~MatPsi2();
    
    // Construcing properties 
    std::string& InputInfo_MoleculeString() { return molstring_; } // the string describing the molecule 
    std::string& InputInfo_BasisSet() { return basisname_; } // basis set name string 
    
    // CPU and memory controll 
    int Settings_MaxNumCPUCores() { return process_environment_.get_n_threads(); }
    double Settings_MaxMemoryInGB() { return (double)process_environment_.get_memory() / 1E+9; }
    std::string Settings_PsiDataDir() { return process_environment_("PSIDATADIR"); }
    std::string Settings_TempDir() { return psio_->_psio_manager_->get_file_path(0); }
    void Settings_SetMaxNumCPUCores(int ncores);
    void Settings_SetMaxMemory(std::string);
    void Settings_SetPsiDataDir(std::string path) { process_environment_.set("PSIDATADIR", path); }
    
    
    
    //*** Molecule properties 
    int Molecule_NumAtoms() { return molecule_->natom(); } // number of atoms 
    int Molecule_NumElectrons(); // number of electrons 
    SharedMatrix Molecule_Geometry() { return molecule_->geometry().clone(); } // geometry in Bohr 
    void Molecule_SetGeometry(SharedMatrix newGeom); // set a new geometry in Bohr 
    double Molecule_NuclearRepulsionEnergy() { return molecule_->nuclear_repulsion_energy(); } // nuclear repulsion energy 
    SharedVector Molecule_AtomicNumbers(); // atomic number list vector 
    void Molecule_SetCharge(int charge) { molecule_->set_molecular_charge(charge); }
    
    //*** Molecule operations 
    void Molecule_Fix();
    void Molecule_Free();
    
    
    //*** Basis set properties 
    void BasisSet_SetBasisSet(const std::string& basisname); // set a new basis set 
    bool BasisSet_IsSpherical() { return basis_->has_puream(); }
    int BasisSet_NumFunctions() { return basis_->nbf(); } // number of basis functions 
    int BasisSet_NumShells() { return basis_->nshell(); }
    SharedVector BasisSet_ShellTypes();
    SharedVector BasisSet_ShellNumPrimitives();
    SharedVector BasisSet_ShellNumFunctions();
    SharedVector BasisSet_ShellToCenter();
    SharedVector BasisSet_FunctionToCenter(); // map basis function number to the number of atom it is centred on 
    SharedVector BasisSet_FunctionToAngularMomentum(); // map basis function number to its angular momentum 
    SharedVector BasisSet_PrimitiveExponents();
    SharedVector BasisSet_PrimitiveCoefficients();
    
    
    //*** Integral package
    SharedMatrix Integrals_Overlap(); // overlap matrix S <i|j>
    SharedMatrix Integrals_Kinetic(); // kinetic energy matrix KE 
    SharedMatrix Integrals_Potential(); // total potential energy matrix EN <i|sum(1/R)|j>
    std::vector<SharedMatrix> Integrals_Dipole(); // dipole matrices <i|x|j>, <i|y|j>, <i|z|j>
    std::vector<SharedMatrix> Integrals_PotentialEachCore(); // atom-separated EN 
    SharedMatrix Integrals_PotentialPointCharges(SharedMatrix Zxyz_list); // compute from a given point charge list the environment potential energy matrix ENVI 
    int Integrals_NumUniqueTEIs(); // number of unique TEIs 
    double Integrals_ijkl(int i, int j, int k, int l); // (ij|kl), chemist's notation 
    // ## HIGH MEMORY COST METHODS ## 
    void Integrals_AllUniqueTEIs(double*); // all unique TEIs in a vector 
    void Integrals_AllTEIs(double*); // all (repetitive) TEIs in a 4D-array 
    void Integrals_IndicesForExchange(double*, double*); // pre-arrange TEI vectors for forming J/K 
    // ## HIGH MEMORY COST METHODS ## 
    
    
    //*** JK related
    // use different types of JK 
    void JK_Initialize(std::string jktype, std::string auxiliaryBasisSetName = "CC-PVDZ-JKFIT");
    const std::string& JK_Type();
    
    // methods computing J/K/G 
    SharedMatrix JK_DensityToJ(SharedMatrix);
    SharedMatrix JK_DensityToK(SharedMatrix);
    SharedMatrix JK_OrbitalToJ(SharedMatrix);
    SharedMatrix JK_OrbitalToK(SharedMatrix);
    SharedMatrix JK_OccupiedOrbitalToJ(SharedMatrix);
    SharedMatrix JK_OccupiedOrbitalToK(SharedMatrix);
    
    // specially for density-fitting JK
    SharedMatrix DFJK_mnQMatrixUnique();
    std::vector<SharedMatrix> DFJK_mnQTensorFull();
    SharedMatrix DFJK_mnAMatrixUnique();
    std::vector<SharedMatrix> DFJK_mnATensorFull();
    SharedMatrix DFJK_InverseJHalfMetric();
    
    
    //*** SCF related
    // create/reset RHF object 
    void RHF_Reset();
    
    // method of doing RHF calculations 
    double RHF_DoSCF(); // "regular" restricted Hartree-Fock 
    
    // methods controlling RHF algorithm 
    void RHF_EnableMOM(int mom_start);
    void RHF_EnableDamping(double damping_percentage);
    void RHF_EnableDIIS();
    void RHF_DisableDIIS();
    void RHF_GuessSAD();
    void RHF_GuessCore();
    
    // methods extracting restricted Hartree-Fock results
    double RHF_TotalEnergy();       // restricted Hartree-Fock energy 
    SharedMatrix RHF_Orbital();   // molecular orbital coefficients
    SharedVector RHF_OrbitalEnergies(); // molecular orbital engenvalues 
    SharedMatrix RHF_Density();   // density matrix 
    SharedMatrix RHF_CoreHamiltonian();  // core Hamiltonian 
    SharedMatrix RHF_JMatrix();   // Coulomb interaction matrix J 
    SharedMatrix RHF_KMatrix();   // exchange interaction matrix K
    SharedMatrix RHF_FockMatrix();   // entire Fock matrix 
    
};
