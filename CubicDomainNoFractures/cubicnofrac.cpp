// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>

// ----- Functions -----

void RunTest(const int caseToSim);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TPZGeoMesh*& gmeshcoarse,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);

TPZGeoMesh *ReadFractureMeshCase0();
void ReadFractureMeshCase1(TPZGeoMesh*& gmesh,TPZGeoMesh*& gmeshcoarse);


enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EDONTUSEGLOBALFRACID, EIntersection, EIntersectionEnd, EPLossAtIntersect,
    /*14 HAS TO BE LAST*/EInitVolumeMatForMHM
};

int globFracID = 10;
// ----- End of Functions -----

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----


//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------

#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");    
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    // 0: 2x1x1 domain
    // 1: 2x2x2 domain with a coarse and a fine mesh to generate MHM structure
    const int caseToSim = 1;
    RunTest(caseToSim);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void RunTest(const int caseToSim)
{
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    TMRSDataTransfer sim_data;
    CreateGMeshAndDataTransfer(gmeshfine,gmeshcoarse,sim_data,caseToSim);
    
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
    
    // ----- Approximation space -----
    TMRSApproxSpaceGenerator aspace;
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    
    // ----- Setting gmesh -----
    aspace.SetGeometry(gmeshfine,gmeshcoarse);
    
    // ----- Setting the global data transfer -----
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
    aspace.SetDataTransfer(sim_data);

    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsingPzSparse = false;
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    mixAnalisys->Assemble();
//    const char* name = "RHS";
//    TPZMatrixSolver<STATE>* matsol = dynamic_cast<TPZMatrixSolver<STATE>*>(mixAnalisys->Solver());
//    std::ofstream oo("matbad.txt");
//    matsol->Matrix()->Print(oo);
//    mixAnalisys->Rhs().Print(name);
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // ----- Post processing -----
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    const int dimToPost = 3;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    
    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TPZGeoMesh*& gmeshcoarse, TMRSDataTransfer &sim_data, const int caseToSim) {
    switch (caseToSim) {
        case 0: {
            gmesh = ReadFractureMeshCase0();
        }
            break;
        case 1: {
            ReadFractureMeshCase1(gmesh,gmeshcoarse);
        }
            break;
        default:
            DebugStop();
    }
            
    sim_data  = SettingFracturesSimple(caseToSim);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer SettingFracturesSimple(const int caseToSim){
    
    // Fracture material
    TMRSDataTransfer sim_data;
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL pressure_in = 1.0 ;
    REAL zero_flux = 0.0;
    
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = globFracID;

    // NS: What are these?
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mIterface_material_idFracBound = 104;
    

    
    // Boundary conditions
    if (caseToSim < 2) {
//        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
//        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EFaceBCPressure,D_Type,pressure_in);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,2.);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,0.);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(EFaceBCPressure,D_Type,1.);
        

//        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(1);
//        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EFracPressure,N_Type,zero_flux);
    
    }
    else {
        DebugStop();
    }
    
    // Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.1;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 1;
    sim_data.mTNumerics.m_max_iter_sfi = 1;
    
    // BorderElementOrder
    sim_data.mTNumerics.m_BorderElementPresOrder=1;
    sim_data.mTNumerics.m_BorderElementFluxOrder=1;
    
    // Other properties?
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0; //*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    //FracAndReservoirProperties
    sim_data.mTFracProperties.m_Permeability = 1.;
    REAL kappa=1.0;
    int  id1=EVolume;
    std::vector<std::pair<int, REAL>> idPerm(1);
    idPerm[0]= std::make_pair(id1,kappa);
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
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
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    return sim_data;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase0(){
            
    TPZGeoMesh* gmesh = nullptr;

    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,1,2};
    const MMeshType elType = MMeshType::EHexahedral;
    
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);
    
    // ----- Fracture element -----
//    int64_t index;
//    TPZManVector<int64_t,2> nodesId = {4,5};
//    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
//    nodesId = {5,7};
//    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
//    nodesId = {7,6};
//    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
//    nodesId = {6,4};
//    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
//    TPZManVector<int64_t,4> nodesIdVec = {4,6,7,5};
//    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
//    gmesh = gen3d.BuildBoundaryElements(EInlet, ENoflux, ENoflux, ENoflux, ENoflux, EOutlet);
    gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure,EFaceBCPressure,EFaceBCPressure,EFaceBCPressure,EFaceBCPressure,EFaceBCPressure);
    
    gmesh->BuildConnectivity();
        
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void ReadFractureMeshCase1(TPZGeoMesh*& gmeshfine,TPZGeoMesh*& gmeshcoarse){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    string basemeshpath(FRACMESHES);
    
    // ===================> Coarse mesh
    
    // Gmsh file
    std::string filenamecoarse = basemeshpath + "/verificationMHMNoHybrid/twoElCoarse.msh";
    
    // Domain
    std::string volbase = "c1";
    for (int ivol = 1; ivol <= 2; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenamecoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
    
    // Gmsh file
    std::string filenamefine = basemeshpath + "/verificationMHMNoHybrid/twoElFine.msh";
        
    // Domain
    for (int ivol = 1; ivol <= 2; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
//        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
//    // Fractures
//    dim_name_and_physical_tagFine[2]["Fracture12"] = globFracID;
//
//    // Fractures BCs
//    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracPressure;
            
    gmeshfine = generateGMeshWithPhysTagVec(filenamefine,dim_name_and_physical_tagFine);
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

