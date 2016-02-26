#include "PyPsi.hpp"

int main()
{
    SharedMatrix geom(new Matrix(3, 4));
    geom->set(0, 0, 8.0); geom->set(0, 1, 0.000000); geom->set(0, 2,  0.000000); geom->set(0, 3,  0.110200);
    geom->set(1, 0, 1.0); geom->set(1, 1, 0.000000); geom->set(1, 2,  0.711600); geom->set(1, 3, -0.440800);
    geom->set(2, 0, 1.0); geom->set(2, 1, 0.000000); geom->set(2, 2, -0.711600); geom->set(2, 3, -0.440800);
    
    std::string basisSet = "sto-3g";
    PyPsi pypsi(geom, basisSet);
    
    pypsi.JK_Initialize("dfjk");
    pypsi.Integrals_Overlap();
    pypsi.SCF_SetSCFType("rks");
    pypsi.SCF_RunSCF();
    
    printf("%f\n", pypsi.SCF_TotalEnergy());
    
    
    return 0;
}


