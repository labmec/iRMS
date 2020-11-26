#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "SBFemHdiv.h"
#include "TPZHybridizeHDiv.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif

//-----------------------------------------
    // CONFIGURATIONS:

    // If hdivcollapsed is true, the multiphysics simulation will be performed using the TPZCompElHDivCollapsed with 1D elements,
    // if not, it will use 2D elements, and traditional Hdiv elements too.
    gHdivcollapsed = false;

    // If sbfemhdiv is true, the class TPZCompElHDivSBFem will be used.
    gSbfemhdiv = true;

    // If hybrid is true, the SBFem-Hdiv mesh will create interface elements
    gHybrid = true;

    // Initial data
    int maxporder = 6;
    bool useexact = true;
    int numthreads = 8;

#ifdef PZDEBUG
    // Outputs
    bool printvtk = true;
    bool printcmesh = true;
    bool printstifmatrix = true;
#endif
//-----------------------------------------

#ifdef _AUTODIFF
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
#endif
    
    for (int POrder = 1; POrder <= maxporder; POrder ++)
    {
        std::cout << "Building geometric mesh...\n";
        int nelx = 1;
        auto gmesh = new TPZGeoMesh;
        if (gSbfemhdiv)
        {
            gmesh = SetupCollapsedMesh(nelx);
        }
        else
        {
            gmesh = SetupSquareMesh(nelx);
        }
        
        TPZCheckGeom gcheck(gmesh);
        gcheck.CheckIds();
        gcheck.CheckUniqueId();

#ifdef PZDEBUG
        if(printvtk)
        {
            std::ofstream out("Geometry.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintGMeshVTK(gmesh, out, true);
        }
#endif

        std::cout << "Building computational mesh...\n";
        TPZCompMesh * cmeshf;
        auto cmeshp = pressure(gmesh, POrder);
        if (gSbfemhdiv)
        {
            cmeshf = fluxhdivsbfem(gmesh, POrder);
#ifdef PZDEBUG
            if (printcmesh)
            {
                std::ofstream out("cmeshsbfemhdiv.txt");
                cmeshf->Print(out);
            }
#endif
        }
        else if(gHdivcollapsed)
        {
            cmeshf = fluxhdivcollapsed(gmesh, POrder);
        }
        else
        {
            cmeshf = flux(gmesh, POrder);
        }

#ifdef PZDEBUG        
        if(printvtk)
        {
            std::ofstream outf("CMeshFlux.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintCMeshVTK(cmeshf, outf, true);
            std::ofstream outp("CMeshPressure.vtk");
            vtk.PrintCMeshVTK(cmeshp, outp, true);
        }
#endif

        auto cmeshm = new TPZMultiphysicsCompMesh;
        if (gHdivcollapsed || gSbfemhdiv)
        {
            cmeshm = multiphysicscollapsed(gmesh, cmeshp, cmeshf, POrder);
        }
        else
        {
            cmeshm = multiphysics(gmesh, cmeshp, cmeshf, POrder);
        }
        
#ifdef PZDEBUG
        if(printcmesh)
        {
            std::ofstream out("CMeshMultiphysics.vtk");
            TPZVTKGeoMesh vtk;
            vtk.PrintCMeshVTK(cmeshm, out, true);
            std::ofstream outcmesh("CMeshMultiphysics.txt");
            cmeshm->Print(outcmesh);
        }
#endif

        if (gHybrid)
        {
            TPZHybridizeHDiv hybrid;
            auto hybridsbfem = hybrid.Hybridize(cmeshm,false);
            // hybrid.HybridizeGivenMesh(*cmeshm, false);
            hybridsbfem->CleanUpUnconnectedNodes();
            hybridsbfem->AdjustBoundaryElements();
            cmeshm->CleanUp();
            cmeshm = hybridsbfem;
        }
        
        std::cout << "Analysis...\n";
        std::cout << "neq = " << cmeshm->NEquations() << std::endl;
        
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
            TPZElementMatrix ek, ef;
            cmeshm->Element(0)->CalcStiff(ek,ef);
            ek.fMat.Print("ek = ", std::cout, EMathematicaInput);
            std::ofstream output("stiffness.txt");
            an.Solver().Matrix()->Print("K = ", output, EMathematicaInput);
        }
        an.Solve();
        
        if (!gHdivcollapsed && !gSbfemhdiv)
        {
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
        
    }
        
    std::cout << "Check:: Calculation finished successfully\n";
    return EXIT_SUCCESS;
}