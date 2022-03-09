// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(string& filenameFine, string& filenameCoarse);
TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EFracture, EIntersection, EIntersectionEnd, EPLossAtIntersect,
    /*14 HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
};

auto exactSol = [](const TPZVec<REAL>& loc,
                   TPZVec<STATE>& u,
                   TPZFMatrix<STATE>& gradU){
    const auto& x = loc[0];
    const auto& y = loc[1];
    const auto& z = loc[2];
    u[0] = 1. - z;
    gradU(0,0) = 0.; // not used
};

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");    
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    string basemeshpath(FRACMESHES);
    string filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/twoElCoarse.msh";
    std::string filenameFine = basemeshpath + "/verificationMHMNoHybrid/twoElFine.msh";
    RunProblem(filenameFine,filenameCoarse);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
void RunProblem(string& filenamefine, string& filenamecoarse)
{
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    ReadMeshes(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    TMRSDataTransfer sim_data;
    FillDataTransfer(sim_data);
    
    // ----- Printing gmesh -----
#ifdef PZDEBUG
    const bool printgmesh = true;
    if (printgmesh) {
        if(gmeshfine){
            std::ofstream name("GeoMesh_Fine_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
        }
        if(gmeshcoarse){
            std::ofstream name("GeoMesh_Coarse_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, name);
        }
    }
#endif
    
    // ----- Approximation space -----
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    
    // ----- Setting gmesh -----
    // Code takes a fine and a coarse mesh to generate MHM data structure
    TMRSApproxSpaceGenerator aspace;
    aspace.InitMatIdForMergeMeshes() = EInitVolumeMatForMHM;
    aspace.FractureMatId() = EFracture;
    aspace.SetGeometry(gmeshfine,gmeshcoarse);
//    aspace.SetGeometry(gmeshfine);
    
    // ----- Setting the global data transfer -----
    aspace.SetDataTransfer(sim_data);
    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = false;
    int n_threads = 0;
    bool UsingPzSparse = false;
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // ----- Post processing -----
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    const int dimToPost = 2;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    
    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, unit_pressure = 1.0, zero_pressure = 0.;
    
    // Domain material
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Volume"] = EVolume;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,unit_pressure);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,zero_pressure);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(EFaceBCPressure,D_Type,unit_pressure);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = EFracture;
    
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(1);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);

    
    // Simulation properties
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    //FracAndReservoirProperties
    REAL kappa=1.0;
    int  id1=EVolume;
    std::vector<std::pair<int, REAL>> idPerm(1);
    idPerm[0]= std::make_pair(id1,kappa);
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    // ----- Use function as BC -----
//    aspace.SetForcingFunctionBC(exactSol);
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    return sim_data;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse){
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c1";
    for (int ivol = 1; ivol <= 2; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 1; ivol <= 2; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = EFracture;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux;
            
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine){
            
    // Creating gmsh reader
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);

    return gmeshFine;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

