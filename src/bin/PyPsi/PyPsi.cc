
#include "PyPsi.hh"
#include <read_options.cc>
#include "numpy/noprefix.h"

namespace psi {
#ifdef PSIDEBUG
    FILE *outfile = stdout;
#else
    FILE *outfile = fopen("/dev/null", "w");
#endif
    char *psi_file_prefix = (char*)"pypsi";
    std::string outfile_name = "";
    extern int read_options(const std::string& name, Options& options, 
                            bool suppress_printing = false);
}

static unsigned long ParseMemoryStr(const std::string& mem_str)
{
    using namespace boost;
    std::string mem_str_ = mem_str;
    algorithm::to_lower(mem_str_);
    cmatch cm;
    regex_search(mem_str_.c_str(), cm, regex("([+]?[0-9]*.?[0-9]+)"));
    double mem = lexical_cast<double>(std::string(cm[1].first, cm[1].second));
    regex_search(mem_str_.c_str(), cm, regex("([a-z])"));
    std::string unit_str(cm[1].first, cm[1].second);
    unsigned long unit;
    if (iequals(unit_str, "b"))
        unit = 1ul;
    else if (iequals(unit_str, "k"))
        unit = 100ul;
    else if (iequals(unit_str, "m"))
        unit = 1000000ul;
    else
        unit = 1000000000ul;
    if ((unsigned long)(mem * unit) < 100000000ul) // less than 100 mb 
        return 100000000ul; // if less than 100mb then return 100 mb 
    else
        return (unsigned long)(mem * unit);
}

static void CheckMatDim(const NPArray& npArray,
                        const int numRow, const int numCol)
{
    using namespace boost::python;
    if (extract<int>(npArray.attr("ndim")) != 2)
        throw PSIEXCEPTION("CheckMatDim: Input must be a 2d array");
    if (numRow != -1 && extract<int>(npArray.attr("shape")[0]) != numRow)
        throw PSIEXCEPTION("CheckMatDim: Number of rows does not agree.");
    if (numCol != -1 && extract<int>(npArray.attr("shape")[1]) != numCol)
        throw PSIEXCEPTION("CheckMatDim: Number of columns does not agree.");
}

static SharedMatrix NPArrayToSharedMatrix(const NPArray& npArray)
{
    using namespace boost::python;
    int dim0 = extract<int>(npArray.attr("shape")[0]);
    int dim1 = extract<int>(npArray.attr("shape")[1]);
    SharedMatrix sharedMatrix(new Matrix(dim0, dim1));
    double **matPtr = sharedMatrix->pointer();
    for (int i = 0; i < dim0; i++)
        for (int j = 0; j < dim1; j++)
            matPtr[i][j] = extract<double>(npArray[i][j]);
    return sharedMatrix;
}

static SharedNPArray NewSharedNPArray(const int numDim, const int *dims,
                                      const int dataType = PyArray_DOUBLE)
{
    boost::scoped_array<Py_intptr_t> dims_(new Py_intptr_t[numDim]);
    for (int i = 0; i < numDim; i++)
        dims_[i] = dims[i];
    return SharedNPArray(new NPArray(boost::python::detail::new_reference(
                         PyArray_ZEROS(numDim, dims_.get(), dataType, 0))));
}

static SharedNPArray NewSharedNPArrayIntVector(const int length)
{
    const int dims[1] = {length};
    return NewSharedNPArray(1, dims, PyArray_INT64);
}

static SharedNPArray NewSharedNPArrayDoubleVector(const int length)
{
    const int dims[1] = {length};
    return NewSharedNPArray(1, dims);
}

static SharedNPArray SharedVectorToSharedNPArray(SharedVector sharedVec)
{
    npy_intp size = sharedVec->dim();
    return SharedNPArray(new NPArray(boost::python::handle<>(
                         PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE,
                                                   sharedVec->pointer()))));
}

static SharedNPArray SharedMatrixToSharedNPArray(SharedMatrix sharedMat)
{
    npy_intp dims[2] = {sharedMat->nrow(), sharedMat->ncol()};
    return SharedNPArray(new NPArray(boost::python::handle<>(
                         PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE,
                                                   sharedMat->get_pointer()))));
}

static SharedPyList VecSharedMatToSharedPyList(const VecSharedMat& vecSharedMat)
{
    SharedPyList sharedPyList(new PyList);
    for (int i = 0; i < vecSharedMat.size(); i++)
        sharedPyList->append(*SharedMatrixToSharedNPArray(vecSharedMat[i]));
    return sharedPyList;
}

static std::string FindDefaultPath()
{
    using namespace boost::python;
    object module = extract<object>(import("imp").attr("find_module"));
    tuple info = extract<tuple>(module("PyPsi"));
    return std::string(extract<char const*>(extract<str>(info[1])));
}

static int FindNumElectrons(const NPArray& xyz, const int charge)
{
    CheckMatDim(xyz, -1, 4);
    SharedVector atomNums = NPArrayToSharedMatrix(xyz)->get_column(0, 0);
    int nelectron = 0;
    for (int i = 0; i < atomNums->dim(); i++)
        nelectron += (int)atomNums->get(i);
    return nelectron - charge;
}

// Constructors
PyPsi::PyPsi(const NPArray& xyz, const std::string& basis)
{
    Construct(xyz, basis,
              0, 1 + FindNumElectrons(xyz, 0) % 2, FindDefaultPath());
}

PyPsi::PyPsi(const NPArray& xyz, const std::string& basis,
             const int charge, const int multiplicity)
{
    Construct(xyz, basis, charge, multiplicity, FindDefaultPath());
}

PyPsi::PyPsi(const NPArray& xyz, const std::string& basis,
             const int charge, const int multiplicity, const std::string& path)
{
    Construct(xyz, basis, charge, multiplicity, path);
}

void PyPsi::Construct(const NPArray& xyz, const std::string& basis,
                      const int charge, const int multiplicity,
                      const std::string& path)
{
    using namespace boost;
    
    CheckMatDim(xyz, -1, 4);
    
    // to output NumPy Array
    import_array();
    
    // some necessary initializations
    process_environment_.initialize();
    
    // set cores and memory 
    process_environment_.set_n_threads(1);
    process_environment_.set_memory(ParseMemoryStr("1000mb"));
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    process_environment_.set_worldcomm(worldcomm_);
    
    // read in options 
    process_environment_.options.set_read_globals(true);
    read_options("", process_environment_.options, true);
    process_environment_.options.set_read_globals(false);
    process_environment_.set("PSIDATADIR", path);
    process_environment_.options.set_global_int("MAXITER", 100);
    
    Wavefunction::initialize_singletons();
    
    // initialize psio 
    create_psio();
    
    // create molecule object and set its basis set name 
    int nelectron = FindNumElectrons(xyz, charge);
    if (multiplicity > nelectron + 1 || multiplicity % 2 == nelectron % 2)
        throw PSIEXCEPTION("PyPsi::common_init: charge/mult not compatible.");
    molecule_ = Molecule::create(NPArrayToSharedMatrix(xyz),
                                 charge, multiplicity);
    molecule_->set_reinterpret_coordentry(false);
    molecule_->set_basis_all_atoms(basis);
    molecule_->set_point_group(shared_ptr<PointGroup>(new PointGroup("C1"))); 
    process_environment_.set_molecule(molecule_);
    
    // create basis object and one & two electron integral factories & rhf 
    create_basis_and_integral_factories();
    
    // create matrix factory object 
    int nbf[] = {basis_->nbf()};
    matfac_ = shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
    
    // set default SCF reference according to multiplicity
    if (molecule_->multiplicity() > 1)
        process_environment_.options.set_global_str("REFERENCE", "UHF");
    else
        process_environment_.options.set_global_str("REFERENCE", "RHF");
    
    // set default DFT functional to B3LYP
    process_environment_.options.set_global_str("DFT_FUNCTIONAL", "B3LYP");
}

void PyPsi::create_psio()
{
    psio_ = boost::shared_ptr<PSIO>(new PSIO);
    process_environment_.set_psio(psio_);
}

void PyPsi::create_basis_and_integral_factories()
{
    using namespace boost;
    // create basis object 
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basis_ = BasisSet::construct(process_environment_, parser,
                                 molecule_, "BASIS");  
    
    // create integral factory object 
    intfac_ = shared_ptr<IntegralFactory>
        (new IntegralFactory(basis_, basis_, basis_, basis_));
    
    // create two electron integral generator
    eri_ = shared_ptr<TwoBodyAOInt>(intfac_->eri());
}

// destructor 
PyPsi::~PyPsi()
{
    if (wfn_ != NULL)
        wfn_->extern_finalize();
    if (jk_ != NULL)
        jk_->finalize();
}

void PyPsi::Settings_SetMaxNumCPUCores(const int ncores)
{
    // set cores and update worldcomm_ 
    process_environment_.set_n_threads(ncores);
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    process_environment_.set_worldcomm(worldcomm_);
}

void PyPsi::Settings_SetMaxMemory(const std::string& memory_str)
{
    // set memory and update worldcomm_ 
    process_environment_.set_memory(ParseMemoryStr(memory_str));
    worldcomm_ = initialize_communicator(0, NULL, process_environment_);
    process_environment_.set_worldcomm(worldcomm_);
}

NPArray PyPsi::Molecule_Geometry()
{
    return *SharedMatrixToSharedNPArray(molecule_->geometry().clone());
}

NPArray PyPsi::Molecule_AtomicNumbers()
{
    SharedNPArray zlistvec = NewSharedNPArrayIntVector(molecule_->natom());
    for (int i = 0; i < molecule_->natom(); i++)
        (*zlistvec)[i] = molecule_->Z(i);
    return *zlistvec;
}

int PyPsi::Molecule_NumElectrons()
{
    int nelectron  = 0;
    for (int i = 0; i < molecule_->natom(); i++)
        nelectron += (int)molecule_->Z(i);
    return nelectron - molecule_->molecular_charge();
}

NPArray PyPsi::Molecule_ChargeMult()
{
    SharedNPArray charge_mult = NewSharedNPArrayIntVector(2);
    (*charge_mult)[0] = molecule_->molecular_charge();
    (*charge_mult)[1] = molecule_->multiplicity();
    return *charge_mult;
}

NPArray PyPsi::BasisSet_ShellTypes()
{
    const int numShell = basis_->nshell();
    SharedNPArray shellTypesVec = NewSharedNPArrayIntVector(numShell);
    for (int i = 0; i < numShell; i++)
        (*shellTypesVec)[i] = basis_->shell(i).am();
    return *shellTypesVec;
}

NPArray PyPsi::BasisSet_ShellNumFunctions()
{
    const int numShell = basis_->nshell();
    SharedNPArray shellNfuncsVec = NewSharedNPArrayIntVector(numShell);
    for (int i = 0; i < numShell; i++)
        (*shellNfuncsVec)[i] = basis_->shell(i).nfunction();
    return *shellNfuncsVec;
}

NPArray PyPsi::BasisSet_ShellNumPrimitives()
{
    const int numShell = basis_->nshell();
    SharedNPArray shellNprimsVec = NewSharedNPArrayIntVector(numShell);
    for (int i = 0; i < numShell; i++)
        (*shellNprimsVec)[i] = basis_->shell(i).nprimitive();
    return *shellNprimsVec;
}

NPArray PyPsi::BasisSet_ShellToCenter()
{
    const int numShell = basis_->nshell();
    SharedNPArray shell2centerVec = NewSharedNPArrayIntVector(numShell);
    for (int i = 0; i < numShell; i++)
        (*shell2centerVec)[i] = basis_->shell_to_center(i);
    return *shell2centerVec;
}

NPArray PyPsi::BasisSet_FuncToCenter()
{
    const int numFunc = basis_->nbf();
    SharedNPArray func2centerVec = NewSharedNPArrayIntVector(numFunc);
    for (int i = 0; i < numFunc; i++)
        (*func2centerVec)[i] = basis_->function_to_center(i);
    return *func2centerVec;
}

NPArray PyPsi::BasisSet_FuncToShell()
{
    const int numFunc = basis_->nbf();
    SharedNPArray func2shellVec = NewSharedNPArrayIntVector(numFunc);
    for (int i = 0; i < numFunc; i++)
        (*func2shellVec)[i] = basis_->function_to_shell(i);
    return *func2shellVec;
}

NPArray PyPsi::BasisSet_FuncToAngular()
{
    const int numFunc = basis_->nbf();
    SharedNPArray func2amVec = NewSharedNPArrayIntVector(numFunc);
    for (int i = 0; i < numFunc; i++)
        (*func2amVec)[i] = basis_->shell(basis_->function_to_shell(i)).am();
    return *func2amVec;
}

NPArray PyPsi::BasisSet_PrimExp()
{
    const int numShell = basis_->nshell();
    const int numPrim = basis_->nprimitive();
    SharedNPArray primExpsVec = NewSharedNPArrayDoubleVector(numPrim);
    std::vector<double> temp;
    for (int i = 0; i < numShell; i++) {
        std::vector<double> currExps = basis_->shell(i).exps();
        temp.insert(temp.end(), currExps.begin(), currExps.end());
    }
    for (int i = 0; i < numPrim; i++)
        (*primExpsVec)[i] = temp[i];
    return *primExpsVec;
}

NPArray PyPsi::BasisSet_PrimCoeffUnnorm()
{
    const int numShell = basis_->nshell();
    const int numPrim = basis_->nprimitive();
    SharedNPArray primCoefsVec = NewSharedNPArrayDoubleVector(numPrim);
    std::vector<double> temp;
    for (int i = 0; i < numShell; i++) {
        std::vector<double> currCoefs = basis_->shell(i).original_coefs();
        temp.insert(temp.end(), currCoefs.begin(), currCoefs.end());
    }
    for (int i = 0; i < numPrim; i++)
        (*primCoefsVec)[i] = temp[i];
    return *primCoefsVec;
}

NPArray PyPsi::Integrals_Overlap()
{
    SharedMatrix sMat(matfac_->create_matrix("Overlap"));
    boost::shared_ptr<OneBodyAOInt> sOBI(intfac_->ao_overlap());
    sOBI->compute(sMat);
    sMat->hermitivitize();
    return *SharedMatrixToSharedNPArray(sMat);
}

NPArray PyPsi::Integrals_Kinetic()
{
    SharedMatrix tMat(matfac_->create_matrix("Kinetic"));
    boost::shared_ptr<OneBodyAOInt> tOBI(intfac_->ao_kinetic());
    tOBI->compute(tMat);
    tMat->hermitivitize();
    return *SharedMatrixToSharedNPArray(tMat);
}

NPArray PyPsi::Integrals_Potential()
{
    SharedMatrix vMat(matfac_->create_matrix("Potential"));
    boost::shared_ptr<OneBodyAOInt> vOBI(intfac_->ao_potential());
    vOBI->compute(vMat);
    vMat->hermitivitize();
    return *SharedMatrixToSharedNPArray(vMat);
}

PyList PyPsi::Integrals_Dipole()
{
    VecSharedMat ao_dipole;
    ao_dipole.push_back(matfac_->create_shared_matrix("Dipole x"));
    ao_dipole.push_back(matfac_->create_shared_matrix("Dipole y"));
    ao_dipole.push_back(matfac_->create_shared_matrix("Dipole z"));
    
    boost::shared_ptr<OneBodyAOInt> dipoleOBI(intfac_->ao_dipole());
    dipoleOBI->compute(ao_dipole);
    ao_dipole[0]->hermitivitize();
    ao_dipole[1]->hermitivitize();
    ao_dipole[2]->hermitivitize();
    return *VecSharedMatToSharedPyList(ao_dipole);
}

PyList PyPsi::Integrals_PotentialEachCore()
{
    using namespace boost;
    VecSharedMat viMatVec;
    shared_ptr<PotentialInt> viPtI((PotentialInt*)intfac_->ao_potential());
    SharedMatrix Zxyz = viPtI->charge_field();
    SharedMatrix Zxyz_rowi(new Matrix(1, 4));
    const int natom = molecule_->natom();
    for (int i = 0; i < natom; i++) {
        SharedVector Zxyz_rowi_vec = Zxyz->get_row(0, i);
        Zxyz_rowi->set_row(0, 0, Zxyz_rowi_vec);
        viPtI->set_charge_field(Zxyz_rowi);
        viMatVec.push_back(matfac_->create_shared_matrix("PotentialEachCore"));
        viPtI->compute(viMatVec[i]);
        viMatVec[i]->hermitivitize();
    }
    return *VecSharedMatToSharedPyList(viMatVec);
}

NPArray PyPsi::Integrals_PotentialPtQ(NPArray& Zxyz_list)
{
    using namespace boost;
    CheckMatDim(Zxyz_list, -1, 4);
    shared_ptr<PotentialInt> viPtI((PotentialInt*)intfac_->ao_potential());
    viPtI->set_charge_field(NPArrayToSharedMatrix(Zxyz_list));
    SharedMatrix vZxyzListMat(matfac_->create_matrix("PotentialPointCharges"));
    viPtI->compute(vZxyzListMat);
    vZxyzListMat->hermitivitize();
    return *SharedMatrixToSharedNPArray(vZxyzListMat);
}

double PyPsi::Integrals_ijkl(int i, int j, int k, int l)
{
    int ish = basis_->function_to_shell(i);
    int jsh = basis_->function_to_shell(j);
    int ksh = basis_->function_to_shell(k);
    int lsh = basis_->function_to_shell(l);
    int ii = i - basis_->shell_to_basis_function(ish);
    int jj = j - basis_->shell_to_basis_function(jsh);
    int kk = k - basis_->shell_to_basis_function(ksh);
    int ll = l - basis_->shell_to_basis_function(lsh);
    int ni = basis_->shell(ish).nfunction();
    int nj = basis_->shell(jsh).nfunction();
    int nk = basis_->shell(ksh).nfunction();
    int nl = basis_->shell(lsh).nfunction();
    eri_->compute_shell(ish, jsh, ksh, lsh);
    const double *buffer = eri_->buffer();
    return buffer[ll + nl * (kk + nk * (jj + nj * ii))];
}

static inline int ij2I(int i, int j)
{
    if (i < j) {
        int tmp = i;
        i = j;
        j = tmp;
    }
    return i * ( i + 1 ) / 2 + j;
}

int PyPsi::Integrals_NumUniqueTEIs()
{
    int nbf = basis_->nbf();
    return (nbf * (nbf + 1) * (nbf * nbf + nbf + 2)) / 8;
}

void PyPsi::Integrals_AllUniqueTEIs(double* matpt)
{
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nuniq = Integrals_NumUniqueTEIs();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next())
            matpt[ij2I(ij2I(intIter.i(), intIter.j()),
                  ij2I(intIter.k(), intIter.l()))] = buffer[intIter.index()];
    }
}

void PyPsi::Integrals_AllTEIs(double* matpt)
{
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nbf = basis_->nbf();
    const double *buf = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            matpt[l + nbf * (k + nbf * (j + nbf * i))] = buf[intIter.index()];
            matpt[l + nbf * (k + nbf * (i + nbf * j))] = buf[intIter.index()];
            matpt[k + nbf * (l + nbf * (j + nbf * i))] = buf[intIter.index()];
            matpt[k + nbf * (l + nbf * (i + nbf * j))] = buf[intIter.index()];
            matpt[j + nbf * (i + nbf * (l + nbf * k))] = buf[intIter.index()];
            matpt[j + nbf * (i + nbf * (k + nbf * l))] = buf[intIter.index()];
            matpt[i + nbf * (j + nbf * (l + nbf * k))] = buf[intIter.index()];
            matpt[i + nbf * (j + nbf * (k + nbf * l))] = buf[intIter.index()];
        }
    }
}

void PyPsi::Integrals_IndicesForK(double* indices1, double* indices2)
{
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            indices1[ij2I(ij2I(i, j), ij2I(k, l))]
                = ij2I(ij2I(i, l), ij2I(k, j));
            indices2[ij2I(ij2I(i, j), ij2I(k, l))]
                = ij2I(ij2I(i, k), ij2I(j, l));
        }
    }
}

void PyPsi::JK_Initialize(std::string jktype, std::string auxBasis)
{
    using namespace boost;
    using namespace scf;
    std::transform(jktype.begin(), jktype.end(), jktype.begin(), ::toupper);
    if (jk_ != NULL)
        jk_->finalize();
    if (wfn_ == NULL) {
        std::string scfType = process_environment_.options.get_str("REFERENCE");
        if (scfType == "RHF" || scfType == "RKS") {
            if (molecule_->multiplicity() > 1)
                throw PSIEXCEPTION("JK_Initialize: "
                                   "RHF or RKS can handle singlets only.");
            wfn_ = shared_ptr<RHF>(new RHF(process_environment_, basis_));
        } else if (scfType == "UHF" || scfType == "UKS") {
            wfn_ = shared_ptr<UHF>(new UHF(process_environment_, basis_));
        } else {
            throw PSIEXCEPTION("JK_Initialize: "
                               "Reference SCF type not recognized.");
        }
        process_environment_.set_wavefunction(wfn_);
        wfn_->extern_finalize();
        wfn_.reset();
    }
    if (jktype == "PKJK") {
        jk_ = shared_ptr<JK>(new PKJK(process_environment_, basis_, psio_));
    } else if (jktype == "DFJK") {
        shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
        molecule_->set_basis_all_atoms(auxBasis, "DF_BASIS_SCF");
        shared_ptr<BasisSet> auxiliary
            = BasisSet::construct(process_environment_, parser,
                                  molecule_, "DF_BASIS_SCF");
        dfjk_ = shared_ptr<DFJK>(new DFJK(process_environment_, basis_,
                                          auxiliary, psio_));
        jk_ = shared_ptr<JK>(dfjk_);
        molecule_->set_basis_all_atoms(basis_->name());
    } else if (jktype == "ICJK") {
        jk_ = shared_ptr<JK>(new ICJK(process_environment_, basis_));
    } else if (jktype == "DIRECTJK") {
        jk_ = shared_ptr<JK>(new DirectJK(process_environment_, basis_));
    } else {
        throw PSIEXCEPTION("JK_Initialize: JK type not recognized.");
    }
    jk_->set_memory(process_environment_.get_memory());
    jk_->set_cutoff(0.0);
    jk_->initialize();
}

const std::string PyPsi::JK_Type()
{
    if (jk_ == NULL)
        throw PSIEXCEPTION("JK_Type: JK object has not been initialized.");
    return jk_->JKtype();
}

static SharedMatrix DensToEigVectors(SharedMatrix density)
{
    if (density == NULL)
        return density;
    int dim = density->ncol();
    SharedMatrix eigVectors(new Matrix(dim, dim));
    boost::shared_ptr<Vector> eigValues(new Vector(dim));
    density->diagonalize(eigVectors, eigValues);
    for (int i = 0; i < dim; i++)
        eigValues->set(i, sqrt(std::max(0.0, eigValues->get(i))));
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            eigVectors->set(i, j, eigVectors->get(i, j) * eigValues->get(j));
    return eigVectors;
}

static void TurnDensSetIntoFakeOccOrbSet(PyList& densSet, int nbf)
{
    int listLength = boost::python::len(densSet);
    for (int i = 0; i < listLength; i++) {
        NPArray dens = boost::python::extract<NPArray>(densSet[i]);
        CheckMatDim(dens, nbf, nbf);
        SharedMatrix fakeOccOrb = DensToEigVectors(NPArrayToSharedMatrix(dens));
        densSet[i] = *SharedMatrixToSharedNPArray(fakeOccOrb);
    }
}

PyList PyPsi::JK_DensToJ(PyList& densSet)
{
    TurnDensSetIntoFakeOccOrbSet(densSet, basis_->nbf());
	return JK_OccOrbToJ(densSet); 
}

PyList PyPsi::JK_DensToK(PyList& densSet)
{
    TurnDensSetIntoFakeOccOrbSet(densSet, basis_->nbf());
	return JK_OccOrbToK(densSet); 
}

PyList PyPsi::JK_OccOrbToJ(PyList& occOrbSet)
{
    if (jk_ == NULL)
        JK_Initialize("PKJK");
    jk_->set_do_K(false);
    JK_CalcAllFromOccOrb(occOrbSet);
    jk_->set_do_K(true);
    return *VecSharedMatToSharedPyList(jk_->J());
}

PyList PyPsi::JK_OccOrbToK(PyList& occOrbSet)
{
    if (jk_ == NULL)
        JK_Initialize("PKJK");
    jk_->set_do_J(false);
    JK_CalcAllFromOccOrb(occOrbSet);
    jk_->set_do_J(true);
    return *VecSharedMatToSharedPyList(jk_->K());
}

void PyPsi::JK_CalcAllFromDens(PyList& densSet)
{
    TurnDensSetIntoFakeOccOrbSet(densSet, basis_->nbf());
    JK_CalcAllFromOccOrb(densSet);
}

void PyPsi::JK_CalcAllFromOccOrb(PyList& occOrbSet)
{
    int listLength = boost::python::len(occOrbSet);
    if (listLength != 1 && listLength != 2)
        throw PSIEXCEPTION("JK_CalcAllFromOccOrb: "
                           "Input list length must be 1 or 2.");
    if (jk_ == NULL)
        JK_Initialize("PKJK");
    jk_->C_left().clear();
    for (int i = 0; i < listLength; i++) {
        NPArray occOrb = boost::python::extract<NPArray>(occOrbSet[i]);
        CheckMatDim(occOrb, basis_->nbf(), -1);
        jk_->C_left().push_back(NPArrayToSharedMatrix(occOrb));
    }
    jk_->compute();
}

PyList PyPsi::JK_RetrieveJ()
{
	if (jk_ == NULL)
        throw PSIEXCEPTION("JK_RetrieveJ: J/K calculation has not been done.");
    return *VecSharedMatToSharedPyList(jk_->J());
}

PyList PyPsi::JK_RetrieveK()
{
	if (jk_ == NULL)
        throw PSIEXCEPTION("JK_RetrieveK: J/K calculation has not been done.");
    return *VecSharedMatToSharedPyList(jk_->K());
}

void PyPsi::JK_DFException(std::string functionName)
{
    if (jk_ == NULL)
        JK_Initialize("DFJK");
    if (!boost::iequals(jk_->JKtype(), "DFJK") || dfjk_ == NULL)
        throw PSIEXCEPTION(functionName + ": Only for Density-fitting JK.");
}

NPArray PyPsi::JK_DFTensor_AuxPriPairs()
{
    JK_DFException("JK_DFTensor_AuxPriPairs");
    return *SharedMatrixToSharedNPArray(dfjk_->GetQmn());
}

NPArray PyPsi::JK_DFMetric_InvJHalf()
{
    JK_DFException("JK_DFMetric_InvJHalf");
    return *SharedMatrixToSharedNPArray(dfjk_->GetInvJHalf());
}

void PyPsi::DFT_Initialize(std::string functional)
{
    std::transform(functional.begin(), functional.end(),
                   functional.begin(), ::toupper);
    process_environment_.options.set_global_str("DFT_FUNCTIONAL", functional);
    std::string scfType = process_environment_.options.get_str("REFERENCE");
    std::string dftType;
    if (scfType == "RHF" || scfType == "RKS") {
        if (molecule_->multiplicity() > 1)
            throw PSIEXCEPTION("DFT_Initialize: "
                               "RHF or RKS can handle singlets only.");
        dftType = "RV";
    } else if (scfType == "UHF" || scfType == "UKS") {
        dftType = "UV";
    } else {
        throw PSIEXCEPTION("DFT_Initialize: "
                           "Reference SCF type not recognized.");
    }
    dftPotential_ = VBase::build_V(process_environment_,
                                   process_environment_.options, dftType);
    dftPotential_->initialize();
}

PyList PyPsi::DFT_DensToV(PyList& densSet)
{
    TurnDensSetIntoFakeOccOrbSet(densSet, basis_->nbf());
    return DFT_OccOrbToV(densSet);
}

PyList PyPsi::DFT_OccOrbToV(PyList& occOrbSet)
{
    if (dftPotential_ == NULL)
        DFT_Initialize(process_environment_.options.get_str("DFT_FUNCTIONAL"));
    int listLength = boost::python::len(occOrbSet);
    if (dynamic_cast<UV*>(dftPotential_.get()) != NULL && listLength != 2)
        throw PSIEXCEPTION("DFT_OccOrbToV: UV needs 2 orbital sets.");
    else if (dynamic_cast<RV*>(dftPotential_.get()) != NULL && listLength != 1)
        throw PSIEXCEPTION("DFT_OccOrbToV: RV needs only 1 orbital set.");
    dftPotential_->C().clear();
    for (int i = 0; i < listLength; i++) {
        NPArray occOrb = boost::python::extract<NPArray>(occOrbSet[i]);
        CheckMatDim(occOrb, basis_->nbf(), -1);
        dftPotential_->C().push_back(NPArrayToSharedMatrix(occOrb));
    }
    dftPotential_->compute();
    return *VecSharedMatToSharedPyList(dftPotential_->V());
}

double PyPsi::DFT_EnergyXC()
{
    return dftPotential_->quadrature_values()["FUNCTIONAL"];
}

void PyPsi::SCF_SetSCFType(std::string scfType)
{
    std::transform(scfType.begin(), scfType.end(), scfType.begin(), ::toupper);
    process_environment_.options.set_global_str("REFERENCE", scfType);
    process_environment_.options.set_global_str("GUESS", "CORE");
}

void PyPsi::SCF_SetGuessOrb(PyList& guessOrbSet)
{
    using namespace boost::python;
    int listLength = len(guessOrbSet);
    if (listLength != 1 && listLength != 2)
        throw PSIEXCEPTION("Input list length must be 1 or 2.");
    process_environment_.options.set_global_str("GUESS", "ORBITAL");
    NPArray guessOrbAlpha = extract<NPArray>(guessOrbSet[0]);
    CheckMatDim(guessOrbAlpha, basis_->nbf(), -1);
    guessOrbital_.clear();
    guessOrbital_.push_back(NPArrayToSharedMatrix(guessOrbAlpha));
    if (listLength == 2) {
        NPArray guessOrbBeta = extract<NPArray>(guessOrbSet[1]);
        CheckMatDim(guessOrbBeta, basis_->nbf(), -1);
        guessOrbital_.push_back(NPArrayToSharedMatrix(guessOrbBeta));
    } else {
        guessOrbital_.push_back(NPArrayToSharedMatrix(guessOrbAlpha));
    }
}

void PyPsi::create_wfn()
{
    using namespace boost;
    using namespace scf;
    if (jk_ == NULL)
        JK_Initialize("PKJK");
    std::string scfType = process_environment_.options.get_str("REFERENCE");
    if (scfType == "RHF") {
        if (molecule_->multiplicity() > 1)
            throw PSIEXCEPTION("RHF can handle singlets only.");
        wfn_ = shared_ptr<RHF>(new RHF(process_environment_, jk_));
    } else if (scfType == "UHF") {
        wfn_ = shared_ptr<UHF>(new UHF(process_environment_, jk_));
    } else if (scfType == "RKS") {
        if (molecule_->multiplicity() > 1)
            throw PSIEXCEPTION("RKS can handle singlets only.");
        wfn_ = shared_ptr<RKS>(new RKS(process_environment_, jk_));
    } else if (scfType == "UKS") {
        wfn_ = shared_ptr<UKS>(new UKS(process_environment_, jk_));
    } else {
        throw PSIEXCEPTION("SCF type not recognized.");
    }
}

double PyPsi::SCF_RunSCF()
{
    create_wfn();
    if (process_environment_.options.get_str("GUESS") == "ORBITAL")
        wfn_->SetGuessOrbital(guessOrbital_);
    process_environment_.set_wavefunction(wfn_);
    return wfn_->compute_energy();
}

void PyPsi::SCF_SetGuessType(const std::string& guessType)
{
    process_environment_.options.set_global_str("GUESS", guessType);
}

double PyPsi::SCF_TotalEnergy()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return wfn_->EHF(); 
}

NPArray PyPsi::SCF_OrbitalAlpha()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray(wfn_->Ca()); 
}

NPArray PyPsi::SCF_OrbitalBeta()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray(wfn_->Cb()); 
}

NPArray PyPsi::SCF_OrbEigValAlpha()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedVectorToSharedNPArray(wfn_->epsilon_a()); 
}

NPArray PyPsi::SCF_OrbEigValBeta()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedVectorToSharedNPArray(wfn_->epsilon_b()); 
}

NPArray PyPsi::SCF_DensityAlpha()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray(wfn_->Da()); 
}

NPArray PyPsi::SCF_DensityBeta()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray(wfn_->Db()); 
}

NPArray PyPsi::SCF_FockAlpha()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray(wfn_->Fa()); 
}

NPArray PyPsi::SCF_FockBeta()
{ 
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray(wfn_->Fb()); 
}

NPArray PyPsi::SCF_Gradient()
{
    if (wfn_ == NULL)
        SCF_RunSCF();
    return *SharedMatrixToSharedNPArray
        (scfgrad::SCFGrad(process_environment_).compute_gradient());
}

NPArray PyPsi::SCF_GuessDensity()
{
    create_wfn();
    process_environment_.set_wavefunction(wfn_);
    return *SharedMatrixToSharedNPArray(wfn_->GuessDensity());
}


