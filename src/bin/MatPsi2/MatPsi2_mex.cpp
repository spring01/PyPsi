#include "mex.h"
#include "class_handle.hpp"
#include "MatPsi2.h"

using namespace std;
using namespace psi;
using namespace boost;

SharedMatrix InputMatrix(const mxArray*& Mat_m) {
    int nrow = mxGetM(Mat_m);
    int ncol = mxGetN(Mat_m);
    SharedMatrix Mat_c(new Matrix(nrow, ncol));
    double* Mat_m_pt = mxGetPr(Mat_m);
    double** Mat_c_pt = Mat_c->pointer();
    for(int i = 0; i < ncol; i++)
        for(int j = 0; j < nrow; j++)
            Mat_c_pt[j][i] = *Mat_m_pt++; // Matlab loops over a column first, but C++ loops over a row first 
    return Mat_c;
}

double InputScalar(const mxArray*& Mat_m) {
    double* Mat_m_pt = mxGetPr(Mat_m);
    return *Mat_m_pt;
}

void OutputMatrix(mxArray*& Mat_m, SharedMatrix Mat_c) {
    int nrow = Mat_c->nrow();
    int ncol = Mat_c->ncol();
    Mat_m = mxCreateDoubleMatrix( nrow, ncol, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    double** Mat_c_pt = Mat_c->pointer();
    for(int i = 0; i < ncol; i++)
        for(int j = 0; j < nrow; j++)
            *Mat_m_pt++ = Mat_c_pt[j][i];
}

void OutputVector(mxArray*& Mat_m, SharedVector Vec_c) {
    int dim = Vec_c->dim();
    Mat_m = mxCreateDoubleMatrix( 1, dim, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    double* Vec_c_pt = Vec_c->pointer();
    for(int i = 0; i < dim; i++) {
        *Mat_m_pt++ = *Vec_c_pt++;
    }
}

void OutputScalar(mxArray*& Mat_m, double scalar) {
    Mat_m = mxCreateDoubleMatrix( 1, 1, mxREAL);
    double* Mat_m_pt = mxGetPr(Mat_m);
    *Mat_m_pt = scalar;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    
    // Constructor
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("MatPsi2 Constructor: One output expected.");
        if ( nrhs!=6 || !mxIsDouble(prhs[1]) || mxGetN(prhs[1]) != 4 || !mxIsChar(prhs[2]) || !mxIsDouble(prhs[3]) ||!mxIsDouble(prhs[4]) || !mxIsChar(prhs[5]) )
            mexErrMsgTxt("MatPsi2 Constructor: MatPsi(cartesian, basis_name, charge, multiplicity) input expected.");
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<MatPsi2>(new MatPsi2(
            InputMatrix(prhs[1]), 
            (std::string)mxArrayToString(prhs[2]), 
            (int)InputScalar(prhs[3]), 
            (int)InputScalar(prhs[4]),
            (std::string)mxArrayToString(prhs[5]) + "/"));
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<MatPsi2>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    MatPsi2* MatPsi_obj = convertMat2Ptr<MatPsi2>(prhs[1]);
    
    int nbf = MatPsi_obj->BasisSet_NumFunctions();
    
    //*** Call the various class methods 
    //*** Settings 
    if (!strcmp("Settings_MaxNumCPUCores", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Settings_MaxNumCPUCores());
        return;
    }
    if (!strcmp("Settings_MaxMemoryInGB", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->Settings_MaxMemoryInGB());
        return;
    }
    if (!strcmp("Settings_PsiDataDir", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->Settings_PsiDataDir()).c_str());
        return;
    }
    if (!strcmp("Settings_TempDir", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->Settings_TempDir()).c_str());
        return;
    }
    if (!strcmp("Settings_SetMaxNumCPUCores", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->Settings_SetMaxNumCPUCores(1);
            return;
        }
        if (nrhs!=3 ||  mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("Settings_SetMaxNumCPUCores(ncores): Integer input expected.");
        MatPsi_obj->Settings_SetMaxNumCPUCores((int)InputScalar(prhs[2]));
        return;
    }
    if (!strcmp("Settings_SetMaxMemory", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->Settings_SetMaxMemory("1000mb");
            return;
        }
        if (nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("Settings_SetMaxMemory(\"memory\"): String input expected.");
        MatPsi_obj->Settings_SetMaxMemory((std::string)mxArrayToString(prhs[2]));
        return;
    }
    if (!strcmp("Settings_SetPsiDataDir", cmd)) {
        if (nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("Settings_SetPsiDataDir(\"psiDataDir\"): String input expected.");
        MatPsi_obj->Settings_SetPsiDataDir((std::string)mxArrayToString(prhs[2]));
        return;
    }
    
    //*** Molecule 
    if (!strcmp("Molecule_Fix", cmd)) {
        MatPsi_obj->Molecule_Fix();
        return;
    }
    if (!strcmp("Molecule_Free", cmd)) {
        MatPsi_obj->Molecule_Free();
        return;
    }
    if (!strcmp("Molecule_NumAtoms", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Molecule_NumAtoms());
        return;
    }
    if (!strcmp("Molecule_NumElectrons", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Molecule_NumElectrons());
        return;
    }
    if (!strcmp("Molecule_Geometry", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Molecule_Geometry());
        return;
    }
    if (!strcmp("Molecule_SetGeometry", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != MatPsi_obj->Molecule_NumAtoms() || mxGetN(prhs[2]) != 3)
            mexErrMsgTxt("Molecule_SetGeometry(newGeom): NumAtoms by 3 matrix input expected.");
        // Call the method
        MatPsi_obj->Molecule_SetGeometry(InputMatrix(prhs[2]));
        return;
    }
    if (!strcmp("Molecule_AtomicNumbers", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->Molecule_AtomicNumbers());
        return;
    }
    if (!strcmp("Molecule_NuclearRepulsionEnergy", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->Molecule_NuclearRepulsionEnergy());
        return;
    }
    if (!strcmp("Molecule_SetCharge", cmd)) {
        if ( nrhs!=3 || !mxIsDouble(prhs[2]))
            mexErrMsgTxt("Molecule_SetCharge(charge): Integer input expected.");
        MatPsi_obj->Molecule_SetCharge((int)InputScalar(prhs[2]));
        return;
    }
    
    //*** BasisSet 
    if (!strcmp("BasisSet_Name", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->BasisSet_Name()).c_str());
        return;
    }
    if (!strcmp("BasisSet_SetBasisSet", cmd)) {
        if ( nrhs!=3 || !mxIsChar(prhs[2]))
            mexErrMsgTxt("BasisSet_SetBasisSet(\"basis\"): String input expected.");
        MatPsi_obj->BasisSet_SetBasisSet((std::string)mxArrayToString(prhs[2]));
        return;
    }
    if (!strcmp("BasisSet_IsSpherical", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->BasisSet_IsSpherical());
        return;
    }
    if (!strcmp("BasisSet_NumFunctions", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->BasisSet_NumFunctions());
        return;
    }
    if (!strcmp("BasisSet_NumShells", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->BasisSet_NumShells());
        return;
    }
    if (!strcmp("BasisSet_ShellTypes", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_ShellTypes());
        return;
    }
    if (!strcmp("BasisSet_ShellNumPrimitives", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_ShellNumPrimitives());
        return;
    }
    if (!strcmp("BasisSet_ShellNumFunctions", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_ShellNumFunctions());
        return;
    }
    if (!strcmp("BasisSet_ShellToCenter", cmd)) {
        SharedVector shell2centerVec = MatPsi_obj->BasisSet_ShellToCenter();
        for(int i = 0; i < shell2centerVec->dim(); i++)
            shell2centerVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], shell2centerVec);
        return;
    }
    if (!strcmp("BasisSet_FunctionToCenter", cmd)) {
        SharedVector func2centerVec = MatPsi_obj->BasisSet_FunctionToCenter();
        for(int i = 0; i < func2centerVec->dim(); i++)
            func2centerVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], func2centerVec);
        return;
    }
    if (!strcmp("BasisSet_FunctionToShell", cmd)) {
        SharedVector func2shellVec = MatPsi_obj->BasisSet_FunctionToShell();
        for(int i = 0; i < func2shellVec->dim(); i++)
            func2shellVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], func2shellVec);
        return;
    }
    if (!strcmp("BasisSet_FunctionToAngularMomentum", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_FunctionToAngularMomentum());
        return;
    }
    if (!strcmp("BasisSet_PrimitiveExponents", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_PrimitiveExponents());
        return;
    }
    if (!strcmp("BasisSet_PrimitiveCoefficients", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->BasisSet_PrimitiveCoefficients());
        return;
    }
    
    
    //*** Integral 
    if (!strcmp("Integrals_Overlap", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_Overlap());
        return;
    }
    if (!strcmp("Integrals_Kinetic", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_Kinetic());
        return;
    }
    if (!strcmp("Integrals_Potential", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_Potential());
        return;
    }
    if (!strcmp("Integrals_PotentialEachCore", cmd)) {
        std::vector<SharedMatrix> viMatArray = MatPsi_obj->Integrals_PotentialEachCore();
        int ncol = viMatArray[0]->ncol();
        int nrow = viMatArray[0]->nrow();
        int natom = MatPsi_obj->Molecule_NumAtoms();
        mwSize dims[3] = {ncol, nrow, natom};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int iatom = 0; iatom < natom; iatom++) {
            double* tmp_pt = viMatArray[iatom]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[iatom*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    if (!strcmp("Integrals_PotentialPointCharges", cmd)) {
        // Check parameters
        if (nrhs!=3)
            mexErrMsgTxt("Integrals_PotentialPointCharges(Zxyz_mat): (number of point charges) by 4 matrix input expected.");
        if (mxGetN(prhs[2]) != 4)
            mexErrMsgTxt("Integrals_PotentialPointCharges: Zxyz list matrix dimension does not agree.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->Integrals_PotentialPointCharges(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("Integrals_Dipole", cmd)) {
        std::vector<SharedMatrix> dipole = MatPsi_obj->Integrals_Dipole();
        if(nlhs == 3) {
            OutputMatrix(plhs[0], dipole[0]);
            OutputMatrix(plhs[1], dipole[1]);
            OutputMatrix(plhs[2], dipole[2]);
            return;
        }
        int ncol = dipole[0]->ncol();
        int nrow = dipole[0]->nrow();
        mwSize dims[3] = {ncol, nrow, 3};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int idim = 0; idim < 3; idim++) {
            double* tmp_pt = dipole[idim]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[idim*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    if (!strcmp("Integrals_ijkl", cmd)) {
        // Check parameters
        if (nrhs!=6 || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]) || !mxIsDouble(prhs[5]))
            mexErrMsgTxt("Integrals_ijkl(i, j, k, l): 4 integers input expected.");
        // Call the method
        int ind[4];
        for(int i = 0; i < 4; i++) {
            ind[i] = (int)mxGetScalar(prhs[2+i]) - 1; // -1 convert Matlab convention to C++ convention
            if(ind[i] < 0 || ind[i] >= nbf)
                mexErrMsgTxt("Integrals_ijkl: Required index not within scale.");
        }
        OutputScalar(plhs[0], MatPsi_obj->Integrals_ijkl(ind[0], ind[1], ind[2], ind[3]));
        return;
    }
    if (!strcmp("Integrals_NumUniqueTEIs", cmd)) {
        OutputScalar(plhs[0], (double)MatPsi_obj->Integrals_NumUniqueTEIs());
        return;
    }
    if (!strcmp("Integrals_AllUniqueTEIs", cmd)) {
        plhs[0] = mxCreateDoubleMatrix( 1, MatPsi_obj->Integrals_NumUniqueTEIs(), mxREAL);
        double* matpt = mxGetPr(plhs[0]);
        MatPsi_obj->Integrals_AllUniqueTEIs(matpt);
        return;
    }
    if (!strcmp("Integrals_AllTEIs", cmd)) {
        mwSize dims[4] = {nbf, nbf, nbf, nbf};
        plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
        double* matpt = mxGetPr(plhs[0]);
        MatPsi_obj->Integrals_AllTEIs(matpt);
        return;
    }
    if (!strcmp("Integrals_IndicesForExchange", cmd)) {
        plhs[0] = mxCreateDoubleMatrix( 1, MatPsi_obj->Integrals_NumUniqueTEIs(), mxREAL);
        double* matpt1 = mxGetPr(plhs[0]);
        plhs[1] = mxCreateDoubleMatrix( 1, MatPsi_obj->Integrals_NumUniqueTEIs(), mxREAL);
        double* matpt2 = mxGetPr(plhs[1]);
        MatPsi_obj->Integrals_IndicesForExchange(matpt1, matpt2);
        return;
    }
    
    //*** JK related  
    if (!strcmp("JK_Initialize", cmd)) {
        if ( (nrhs!=3 && nrhs!=4) || !mxIsChar(prhs[2]) )
            mexErrMsgTxt("JK_Initialize(\"jktype\"): String input expected.");
        std::string auxiliaryBasisSetName;
        if (nrhs==3) {
            MatPsi_obj->JK_Initialize((std::string)mxArrayToString(prhs[2]));
        }else {
            MatPsi_obj->JK_Initialize((std::string)mxArrayToString(prhs[2]), (std::string)mxArrayToString(prhs[3]));
        }
        return;
    }
    if (!strcmp("JK_Type", cmd)) {
        plhs[0] = mxCreateString((MatPsi_obj->JK_Type()).c_str());
        return;
    }
    if (!strcmp("JK_DensityToJ", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != nbf || mxGetN(prhs[2]) != nbf)
            mexErrMsgTxt("JK_DensityToJ(density): nbf by nbf matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_DensityToJ(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("JK_DensityToK", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != nbf || mxGetN(prhs[2]) != nbf)
            mexErrMsgTxt("JK_DensityToK(density): nbf by nbf matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_DensityToK(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("JK_OrbitalToJ", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != nbf || mxGetN(prhs[2]) < MatPsi_obj->Molecule_NumElectrons()/2)
            mexErrMsgTxt("JK_OrbitalToJ(orbital): nbf by (at least numElec/2) matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_OrbitalToJ(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("JK_OrbitalToK", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != nbf || mxGetN(prhs[2]) < MatPsi_obj->Molecule_NumElectrons()/2)
            mexErrMsgTxt("JK_OrbitalToK(orbital): nbf by (at least numElec/2) matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_OrbitalToK(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("JK_OccupiedOrbitalToJ", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != nbf)
            mexErrMsgTxt("JK_OrbitalToJ(occupiedOrbital): nbf by any matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_OccupiedOrbitalToJ(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("JK_OccupiedOrbitalToK", cmd)) {
        // Check parameters
        if (nrhs!=3 || mxGetM(prhs[2]) != nbf)
            mexErrMsgTxt("JK_OrbitalToK(occupiedOrbital): nbf by any matrix input expected.");
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->JK_OccupiedOrbitalToK(InputMatrix(prhs[2])));
        return;
    }
    if (!strcmp("DFJK_mnQMatrixUnique", cmd)) {
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->DFJK_mnQMatrixUnique());
        return;
    }
    if (!strcmp("DFJK_mnQTensorFull", cmd)) {
        std::vector<SharedMatrix> mnQFull = MatPsi_obj->DFJK_mnQTensorFull();
        int ncol = mnQFull[0]->ncol();
        int nrow = mnQFull[0]->nrow();
        int naux = mnQFull.size();
        mwSize dims[3] = {ncol, nrow, naux};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int iaux = 0; iaux < naux; iaux++) {
            double* tmp_pt = mnQFull[iaux]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[iaux*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    if (!strcmp("DFJK_mnAMatrixUnique", cmd)) {
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->DFJK_mnAMatrixUnique());
        return;
    }
    if (!strcmp("DFJK_mnATensorFull", cmd)) {
        std::vector<SharedMatrix> mnAFull = MatPsi_obj->DFJK_mnATensorFull();
        int ncol = mnAFull[0]->ncol();
        int nrow = mnAFull[0]->nrow();
        int naux = mnAFull.size();
        mwSize dims[3] = {ncol, nrow, naux};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int iaux = 0; iaux < naux; iaux++) {
            double* tmp_pt = mnAFull[iaux]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[iaux*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    if (!strcmp("DFJK_InverseJHalfMetric", cmd)) {
        // Call the method
        OutputMatrix(plhs[0], MatPsi_obj->DFJK_InverseJHalfMetric());
        return;
    }
    
    //*** SCF related 
    if (!strcmp("SCF_RunRHF", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->SCF_RunRHF());
        return;
    }
    if (!strcmp("SCF_RunUHF", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->SCF_RunUHF());
        return;
    }
    if (!strcmp("SCF_RunRKS", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->SCF_RunRKS());
        return;
    }
    if (!strcmp("SCF_EnableMOM", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->SCF_EnableMOM(20);
            return;
        }
        if (nrhs!=3 || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("SCF_EnableMOM(mom_start): Integer input expected.");
        MatPsi_obj->SCF_EnableMOM((int)InputScalar(prhs[2]));
        return;
    }
    if (!strcmp("SCF_DisableMOM", cmd)) {
        MatPsi_obj->SCF_EnableMOM(0);
        return;
    }
    if (!strcmp("SCF_EnableDamping", cmd)) {
        if (nrhs == 2) {
            MatPsi_obj->SCF_EnableDamping(20.0);
            return;
        }
        if (nrhs!=3 || mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
            mexErrMsgTxt("SCF_EnableDamping(damping_percentage): 1 double input expected.");
        MatPsi_obj->SCF_EnableDamping(InputScalar(prhs[2]));
        return;
    }
    if (!strcmp("SCF_DisableDamping", cmd)) {
        MatPsi_obj->SCF_EnableDamping(0.0);
        return;
    }
    if (!strcmp("SCF_EnableDIIS", cmd)) {
        MatPsi_obj->SCF_EnableDIIS();
        return;
    }
    if (!strcmp("SCF_DisableDIIS", cmd)) {
        MatPsi_obj->SCF_DisableDIIS();
        return;
    }
    if (!strcmp("SCF_GuessSAD", cmd)) {
        MatPsi_obj->SCF_GuessSAD();
        return;
    }
    if (!strcmp("SCF_GuessCore", cmd)) {
        MatPsi_obj->SCF_GuessCore();
        return;
    }
    if (!strcmp("SCF_TotalEnergy", cmd)) {
        OutputScalar(plhs[0], MatPsi_obj->SCF_TotalEnergy());
        return;
    }
    if (!strcmp("SCF_OrbitalAlpha", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_OrbitalAlpha());
        return;
    }
    if (!strcmp("SCF_OrbitalBeta", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_OrbitalBeta());
        return;
    }
    if (!strcmp("SCF_OrbitalEnergiesAlpha", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->SCF_OrbitalEnergiesAlpha());
        return;
    }
    if (!strcmp("SCF_OrbitalEnergiesBeta", cmd)) {
        OutputVector(plhs[0], MatPsi_obj->SCF_OrbitalEnergiesBeta());
        return;
    }
    if (!strcmp("SCF_DensityAlpha", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_DensityAlpha());
        return;
    }
    if (!strcmp("SCF_DensityBeta", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_DensityBeta());
        return;
    }
    if (!strcmp("SCF_CoreHamiltonian", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_CoreHamiltonian());
        return;
    }
    if (!strcmp("SCF_FockAlpha", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_FockAlpha());
        return;
    }
    if (!strcmp("SCF_FockBeta", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_FockBeta());
        return;
    }
    if (!strcmp("SCF_Gradient", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_Gradient());
        return;
    }
    if (!strcmp("SCF_InitialGuessDensity", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_InitialGuessDensity());
        return;
    }
    if (!strcmp("SCF_RHF_Coulomb", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_RHF_Coulomb());
        return;
    }
    if (!strcmp("SCF_RHF_Exchange", cmd)) {
        OutputMatrix(plhs[0], MatPsi_obj->SCF_RHF_Exchange());
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

