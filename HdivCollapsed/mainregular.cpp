#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "pzcmesh.h"
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "TPZAnalyticSolution.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "mixedpoisson.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZVecL2.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "TPZCompElHDivCollapsed.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
TLaplaceExample1 LaplaceExact;
#endif

inline void Laplace_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    LaplaceExact.Solution(xv, val, deriv);
}

TPZGeoMesh * SetupSquareMesh(int nelx, bool hdivcollapsed);
TPZCompMesh * flux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);
TPZCompMesh * fluxhdivcollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder);
TPZCompMesh * pressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder, bool hdivcollapsed);
TPZMultiphysicsCompMesh * multiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder, bool hdivcollapsed);

enum MMATID {Emat1, Ebc1, Ebc2, Ebc3, Ebc4, Eintmat1, Eintbc1, Eintbc2, Eintbc3, Eintbc4, ESkeleton};

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif

#ifndef _AUTODIFF
    std::cout << "This program needs FAD to run \n";
    DebugStop();
#endif

    // Initial data
    int maxnelxcount = 5, maxporder = 6;

    bool useexact = true;
#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
#endif
    // If hdivcollapsed is true, the multiphysics simulation will be performed using the TPZCompElHDivCollapsed with 1D elements,
    // if not, it will use 2D elements, and traditional Hdiv elements too.
    bool hdivcollapsed = true; 
    int numthreads = 4;

    bool printgeomesh = true;
    bool printcmesh = true;
    bool printstifmatrix = false;

    for (int POrder = 1; POrder <= maxporder; POrder ++)
    {
        int nelx = 3;
        auto gmesh = SetupSquareMesh(nelx, hdivcollapsed);
        
        TPZCheckGeom gcheck(gmesh);
        gcheck.CheckIds();
        gcheck.CheckUniqueId();

        if(printgeomesh)
        {
            std::ofstream out("Geometry.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(gmesh, out, true);
        }

        std::cout << "Building computational mesh...\n";
        TPZCompMesh * cmeshf;
        auto cmeshp = pressure(gmesh, POrder, hdivcollapsed);
        if (hdivcollapsed)
        {
            cmeshf = fluxhdivcollapsed(gmesh, POrder);
        }
        else
        {
            cmeshf = flux(gmesh, POrder);
        }
        
        if(printcmesh)
        {
            std::ofstream outf("CMeshFlux.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintCMeshVTK(cmeshf, outf, true);
            std::ofstream outp("CMeshPressure.vtk");
            vtk.PrintCMeshVTK(cmeshp, outp, true);
        }
        std::cout << "Creating multiphysics mesh\n";
        TPZManVector< TPZCompMesh *, 2> meshvector(2);
        meshvector[0] = cmeshf;
        meshvector[1] = cmeshp;
        auto cmeshm = multiphysics(gmesh, POrder, hdivcollapsed);
        TPZManVector<int> active(2,1);
        cmeshm->BuildMultiphysicsSpace(active, meshvector);
        cmeshm->LoadReferences();

        if(printcmesh)
        {
            std::ofstream out("CMeshMultiphysics.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintCMeshVTK(cmeshm, out, true);
            std::ofstream outcmesh("CMeshMultiphysics.txt");
            cmeshm->Print(outcmesh);
        }

        std::cout << "Analysis...\n";
        std::cout << "neq = " << cmeshm->NEquations() << std::endl;
        
        // I'm running without optimizing bandwidth to help to find where is the bug
        // If it is false, the code breaks in the ComputeRequiredData of the TPZCompElHDivCollapsed - I believe the problem is there.
        bool optimizeBandwidth = false;
        TPZAnalysis an(cmeshm, optimizeBandwidth);
        an.SetStep(POrder);
        an.SetExact(Laplace_exact);

        TPZSymetricSpStructMatrix matskl(cmeshm);
        an.SetStructuralMatrix(matskl);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();

        if(printstifmatrix)
        {
            std::ofstream output("stiffness.txt");
            an.Solver().Matrix()->Print("K = ",output,EMathematicaInput);
        }
        an.Solve();
        
        std::cout << "Post Processing...\n";
        std::string plotfile("DarcyP.vtk");
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Pressure");  
        vecnames.Push("Flux");
        int postProcessResolution = 3;

        int dim = gmesh->Dimension();
        an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
        an.PostProcess(postProcessResolution,dim);

        TPZManVector<REAL,3> error;
        ofstream outerror("error.txt");
        an.PostProcessError(error, false, outerror);
        
    }
        
    std::cout << "Check:: Calculation finished successfully\n";
    return EXIT_SUCCESS;
}

TPZGeoMesh * SetupSquareMesh(int nelx, bool hdivcollapsed)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh, Emat1);

    // Setting boundary conditions
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);

    gmesh->BuildConnectivity();

    if (hdivcollapsed)
    {
        TPZGeoMesh *gmeshcollapsed = new TPZGeoMesh();
        gmeshcollapsed->SetDimension(1);
        gmeshcollapsed->NodeVec() = gmesh->NodeVec();

        for (auto gel : gmesh->ElementVec())
        {
            if (!gel) continue;
            if (gel->MaterialId() != Emat1) continue;
            for (auto isides = 0; isides < gel->NSides(); isides++)
            {
                TPZGeoElSide gelside = gel->Neighbour(isides);
                if (!gelside || gelside.Dimension() != 1) continue;
                if (gelside.Element()->Index() < gel->Index()) continue; // if the neighbour is a previous element, continue (the interface element was already created)
                
                TPZManVector<int64_t> nodes(2,0);
                nodes[0] = gelside.SideNodeIndex(0);
                nodes[1] = gelside.SideNodeIndex(1);
                int64_t index;
                gmeshcollapsed->CreateGeoElement(EOned, nodes, ESkeleton, index);
            }
        }
        
        gmeshcollapsed->BuildConnectivity();
        
        for (auto gel : gmeshcollapsed->ElementVec())
        {
            if(!gel) continue;
            TPZManVector<REAL,3> xco1(3), xco2(3);
            gel->Node(0).GetCoordinates(xco1);
            gel->Node(1).GetCoordinates(xco2);
            if (IsZero(xco1[1] + 1) && IsZero(xco2[1] + 1))
            {
                TPZGeoElBC(gel,2,Ebc1);
            }
            else if(IsZero(xco1[0] - 1) && IsZero(xco2[0] - 1))
            {
                TPZGeoElBC(gel,2,Ebc2);
            }
            else if(IsZero(xco1[1] - 1) && IsZero(xco2[1] - 1))
            {
                TPZGeoElBC(gel,2,Ebc3);
            }
            else if(IsZero(xco1[0] + 1) && IsZero(xco2[0] + 1))
            {
                TPZGeoElBC(gel,2,Ebc4);
            }
        }
        gmesh = gmeshcollapsed;
        gmesh->ResetConnectivities();
    }
    
    return gmesh;
}

TPZCompMesh * flux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    
    // TPZVecL2 *mat = new TPZVecL2(ESkeleton);
    TPZMixedPoisson *mat = new TPZMixedPoisson(Emat1,dim);
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    
    cmesh->InsertMaterialObject(mat); //Insere material na malha
    cmesh->SetDimModel(dim);

    int nstate = 1;
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

TPZCompMesh * fluxhdivcollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = 1;
    TPZCompMesh * cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->CleanUp();
    cmeshcollapsed->SetDimModel(dim);

    TPZMixedPoisson *mat = new TPZMixedPoisson(ESkeleton,dim);
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    
    cmeshcollapsed->InsertMaterialObject(mat); //Insere material na malha
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->Dimension() != 1) continue;
        if (gel->MaterialId() != ESkeleton) continue;
        int64_t index;
        auto celhdivc = new TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>(*cmeshcollapsed, gel, index);
    }
    cmeshcollapsed->AutoBuild();

    int nstate = 1;
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    cmeshcollapsed->AutoBuild();
    cmeshcollapsed->AdjustBoundaryElements();
    cmeshcollapsed->CleanUpUnconnectedNodes();
    
    return cmeshcollapsed;
}

TPZCompMesh * pressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder, bool hdivcollapsed)
{        
    auto fmat = Emat1;
    auto dim = 2;
    if (hdivcollapsed)
    {
        fmat = ESkeleton;
        dim = 1;
    }

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    
    TPZMixedPoisson *mat = new TPZMixedPoisson(fmat,dim);
    cmesh->InsertMaterialObject(mat); //Insere material na malha

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZMultiphysicsCompMesh *  multiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder, bool hdivcollapsed)
{
    // gmesh->ResetReference();
    auto dim = gmesh->Dimension();
    auto fmat = Emat1;
    if (hdivcollapsed)
    {
        fmat = ESkeleton;
        dim = 1;
    }

    TPZMultiphysicsCompMesh * cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    
    TPZMixedPoisson *mat = new TPZMixedPoisson(fmat, dim);
    
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}