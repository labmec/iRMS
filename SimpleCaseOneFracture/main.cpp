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
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void ChangeBCsToNoFlux(TPZGeoMesh* gmesh);

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase1(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase2(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase3();
TPZGeoMesh *ReadFractureMeshCase4();
TPZGeoMesh *ReadFractureMeshCase5();

enum EMatid {ENone, EVolume, EInlet, EOutlet, ENoflux, EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure, EDONTUSEGLOBALFRACID, EIntersection, EIntersectionEnd, EPLossAtIntersect};
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
static TPZLogger mainlogger("onefrac");
#endif

int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");    
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for SimpleCaseOneFracture target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    const int caseToSim = 4;
    // 0: 1 frac cte pressure
    // 1: 1 frac linear pressure variation
    // 2: 1 frac cte pressure | frac domain touching frac bnd
    // 3: 2 frac cte pressure | frac domain touching frac bnd | w/ frac intersection
    // 4: 2 frac linear pressure variation | frac domain touching frac bnd | w/ frac intersection
    // 5: 1 frac linear pressure variation | frac domain touching frac bnd
    RunTest(caseToSim);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void RunTest(const int caseToSim)
{
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmesh = nullptr;
    TMRSDataTransfer sim_data;
    CreateGMeshAndDataTransfer(gmesh,sim_data,caseToSim);
    
    const bool printgmesh = true;
    if (printgmesh) {
        std::ofstream name("GeoMeshFractures.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name);
    }
    
    // ----- Approximation space -----
    TMRSApproxSpaceGenerator aspace;
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
    
    // ----- Setting gmesh -----
    aspace.SetGeometry(gmesh);
    
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
    RenumType renumtype = RenumType::EDefault;
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, renumtype);
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
    delete gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim) {
    string basemeshpath(FRACMESHES);
    std::string filename = basemeshpath + "/2DMeshes/1fracNoBnd.msh";
    switch (caseToSim) {
        case 0: {
            gmesh = ReadFractureMeshCase0(filename);
            
        }
            break;
        case 1: {
            gmesh = ReadFractureMeshCase1(filename);
        }
            break;
        case 2: {
            std::string filename2 = basemeshpath + "/verifications/1frac2el.msh";
            gmesh = ReadFractureMeshCase2(filename2);
        }
            break;
        case 3: {
            gmesh = ReadFractureMeshCase3();
        }
            break;
        case 4: {
            gmesh = ReadFractureMeshCase4();
        }
            break;
        case 5: {
            gmesh = ReadFractureMeshCase5();
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
    
    
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fractures"] = globFracID;

    // NS: What are these?
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    

    
    // Boundary conditions
    if (caseToSim < 2) {
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,1.);
        
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracPressure] = std::make_pair(N_Type, zero_flux);
    }
    else if (caseToSim < 4){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EFaceBCPressure] = std::make_pair(D_Type,1.);

		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracPressure] = std::make_pair(D_Type, 1.);
    }
    else if (caseToSim == 4){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
		
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracInlet] = std::make_pair(D_Type, 2.);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracOutlet] = std::make_pair(D_Type, 0.);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracNoFlux] = std::make_pair(N_Type, zero_flux);
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EPLossAtIntersect] = std::make_pair(Mixed_Type, 0.5);
        
        sim_data.mTFracIntersectProperties.m_IntersectionPressureLossId = EPLossAtIntersect;
    }
    else if (caseToSim == 5){
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EInlet] = std::make_pair(D_Type,2.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[EOutlet] = std::make_pair(D_Type,0.);
		sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[ENoflux] = std::make_pair(N_Type,zero_flux);
        
		sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[EFracNoFlux] = std::make_pair(N_Type, zero_flux);
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
    sim_data.mTNumerics.m_MortarBorderElementPresOrder=1;
    sim_data.mTNumerics.m_MortarBorderElementFluxOrder=1;
    
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
    
    //FracAndReservoirProperties
	DebugStop(); // the next few lines were added without testing. Please delete DebugStop and check if everything is fine
	TMRSDataTransfer::TFracProperties::FracProp fracprop;
	fracprop.m_perm = 1.;
	fracprop.m_width = 1.;
	fracprop.m_fracbc.insert(EFracNoFlux);
	fracprop.m_fracIntersectMatID = EIntersection;
	sim_data.mTFracProperties.m_fracprops[globFracID] = fracprop;
	
	REAL kappa=1.0;
    int  id1=EVolume;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa;
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
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    return sim_data;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // domain
    std::string volbase = "c";

    // Volumes
    for (int ivol = 1; ivol < 28; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracPressure;
            
    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase1(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // domain
    std::string volbase = "c";

    // Volumes
    for (int ivol = 1; ivol < 28; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracPressure;
        
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    ChangeBCsToNoFlux(gmesh);
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase2(std::string &filename){
    
    const bool fromgmesh = 0;
    TPZGeoMesh* gmesh = nullptr;
    if (fromgmesh) {
        TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
        
        // domain
        std::string volbase = "c";
        
        // Volumes
        for (int ivol = 1; ivol < 3; ivol++) {
            std::string ivolstring = volbase + to_string(ivol);
            dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
        }
        dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
        
        // Fractures
        dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;
        
        // Fractures BCs
        dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracPressure;
        
        gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    }
    else{
        // ----- Create Geo Mesh -----
        const TPZVec<REAL> minX = {-1.,-1.,-1.};
        const TPZVec<REAL> maxX = {1.,1.,1.};
        const TPZVec<int> nelDiv = {1,1,2};
        const MMeshType elType = MMeshType::EHexahedral;

        TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
        gmesh = gen3d.BuildVolumetricElements(EVolume);

        // ----- Fracture element -----
        int64_t index;
        TPZManVector<int64_t,2> nodesId = {4,5};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
        nodesId = {5,7};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
        nodesId = {7,6};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
        nodesId = {6,4};
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
        TPZManVector<int64_t,4> nodesIdVec = {4,6,7,5};
        new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);

        // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
        gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
        
        gmesh->BuildConnectivity();
    }
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase3(){
    
    TPZGeoMesh* gmesh = nullptr;
    
    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,2,2};
    const MMeshType elType = MMeshType::EHexahedral;

    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);

    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {6,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {7,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {9,11};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {10,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {8,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);

    nodesId = {9,15};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {14,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {8,2};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {2,3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);
    nodesId = {3,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracPressure,*gmesh,index);

    nodesId = {8,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EIntersection,*gmesh,index);
    
    TPZManVector<int64_t,4> nodesIdVec = {6,7,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {2,3,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
    gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
    
    gmesh->BuildConnectivity();

    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase4(){
    
    TPZGeoMesh* gmesh = nullptr;
    
    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,2,2};
    const MMeshType elType = MMeshType::EHexahedral;

    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);

    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {6,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracInlet,*gmesh,index);
    nodesId = {7,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {9,11};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracOutlet,*gmesh,index);
    nodesId = {10,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {8,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);

    nodesId = {9,15};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {14,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {8,2};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {2,3};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {3,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);

    nodesId = {8,9};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EIntersection,*gmesh,index);
    
    TPZManVector<int64_t,4> nodesIdVec = {6,7,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,11,10};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {2,3,9,8};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    nodesIdVec = {8,9,15,14};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
    gmesh = gen3d.BuildBoundaryElements(ENoflux, ENoflux, EInlet, ENoflux, EOutlet, ENoflux);
    
    gmesh->BuildConnectivity();

    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase5(){
    
    TPZGeoMesh* gmesh = nullptr;

    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,1,2};
    const MMeshType elType = MMeshType::EHexahedral;
    
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    gmesh = gen3d.BuildVolumetricElements(EVolume);
    
    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {4,5};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {5,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {7,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    nodesId = {6,4};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EFracNoFlux,*gmesh,index);
    TPZManVector<int64_t,4> nodesIdVec = {4,6,7,5};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    
    // OBS: For some reason, the code leads to wrong results if these bcs are created before the fracture
    gmesh = gen3d.BuildBoundaryElements(EInlet, ENoflux, ENoflux, ENoflux, ENoflux, EOutlet);
    
    gmesh->BuildConnectivity();
    
    
    return gmesh;
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

void ChangeBCsToNoFlux(TPZGeoMesh* gmesh) {
    for(auto gel : gmesh->ElementVec()){
        if(gel->Dimension() != 2) continue;
        if(gel->MaterialId() != EFaceBCPressure) continue;

        const int pos = 2;
        TPZManVector<REAL,2> qsicent(2,-1.);
        TPZManVector<REAL,3> cent(3,-1.);
        gel->CenterPoint(gel->NSides()-1, qsicent);
        gel->X(qsicent, cent);
        const REAL distM1 = fabs(cent[pos] + 1.);
        const REAL distP1 = fabs(cent[pos] - 1.);
        const REAL tol = 1.e-8;
        if (distM1 < tol) {
//            gel->SetMaterialId(EFaceBCPressure);
            gel->SetMaterialId(EInlet);
            continue;
        }
        if (distP1 < tol) {
//            gel->SetMaterialId(EFaceBCPressure);
            gel->SetMaterialId(EOutlet);
            continue;
        }
//        gel->SetMaterialId(EFaceBCPressure);
        gel->SetMaterialId(ENoflux);
    }
}
