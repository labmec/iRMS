// C++ includes
// nothing

// PZ includes
#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TPZRefPatternTools.h"
#include "TPZReservoirTools.h"
#include "pzlog.h"
#include "imrs_config.h"

// ----- Functions -----

void RunTest(const int caseToSim);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void ChangeBCsToNoFlux(TPZGeoMesh* gmesh);

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase1(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase2(std::string &filename);

enum EMatid {ENone, EDomain, EInlet, EOutlet, ENoflux, EPressure, EIntersection, EIntersectionEnd, EVolume, EFaceBCPressure};
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

static TPZLogger mainlogger("onefrac");

int main(){
    TPZLogger::InitializePZLOG("log4cxx.cfg");    
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for SimpleCaseOneFracture target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
    
    const int caseToSim = 2;
    // 0: 1 frac cte pressure
    // 1: 1 frac linear pressure variation
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
    aspace.SetDataTransfer(sim_data);
    aspace.MatIDFracIntesect() = EIntersection;
    
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
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // ----- Post processing -----
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    const int dimToPost = 2;
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
        

        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EPressure,N_Type,zero_flux);
    
    }
    else if (caseToSim < 4){
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EFaceBCPressure,D_Type,1.);
        

        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EPressure,D_Type,1.);
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
    sim_data.mTFracProperties.m_Permeability = 1.e4;
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
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
            
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
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
        
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    ChangeBCsToNoFlux(gmesh);
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase2(std::string &filename){
    
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
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
        
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
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
