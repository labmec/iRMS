#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// C++ includes

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
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "pzlog.h"

// ----- Functions -----

void CaseOnlyFractures(const int caseToSim);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
    
TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase1(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase2(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase3(std::string &filename);

enum EMatid {ENone, EDomain, EInlet, EOutlet, ENoflux, EPressure, EIntersection};
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

const bool GLOBALctePressureCase = 0;

int main(){
    TPZLogger::InitializePZLOG();
    const int caseToSim = 3;
    // 0: 2 perpendicular fractures, cte pressure
    // 1: 2 perpendicular fractures, 1D flow
    // 2: Flemisch example 1
    // 3: Flemisch example 2
    // 4: Flemisch example 3
    // 5: Flemisch example 4
    CaseOnlyFractures(caseToSim);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CaseOnlyFractures(const int caseToSim)
{
    // Reading ONLY the fractures of mesh from DFN
    TPZGeoMesh *gmesh = nullptr;
    TMRSDataTransfer sim_data;
    
    CreateGMeshAndDataTransfer(gmesh,sim_data,caseToSim);
    
    const bool printgmesh = true;
    if (printgmesh) {
        std::ofstream name("GeoMeshFractures.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name);
    }
        
    TMRSApproxSpaceGenerator aspace;

    // Approximation space
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E2Space;
    
    // Setting gmesh
    aspace.SetGeometry(gmesh);
    
    // Setting the global data transfer
    aspace.SetDataTransfer(sim_data);
    aspace.MatIDFracIntesect() = EIntersection;
    
    // Creates de multiphysics compmesh
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    // Analysis parameters
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsingPzSparse = false;
    bool UsePardiso_Q = true; // This parameter should be always "true"
    
    // Setting analysis
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    // -------------- End of running problem --------------
    
    // Post processing
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    mixAnalisys->PostProcessTimeStep();
    
    // Cleaning up
    delete gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim) {
    string basemeshpath(FRACMESHES);
    switch (caseToSim) {
        case 0: {
            std::string filename = basemeshpath + "/2DMeshes/2fracFromDfn.msh";
            gmesh = ReadFractureMeshCase0(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 1: {
//            std::string filename = basemeshpath + "/2DMeshes/2frac1DfluxFromDfn.msh";
            std::string filename = basemeshpath + "/2DMeshes/2fracRef.msh";            
            gmesh = ReadFractureMeshCase1(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 2: {
            std::string filename = basemeshpath + "/2DMeshes/flemisch1.msh";
            gmesh = ReadFractureMeshCase2(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 3: {
            std::string filename = basemeshpath + "/2DMeshes/flemisch2.msh";
            gmesh = ReadFractureMeshCase3(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 4: {
            DebugStop();
        }
            break;
        case 5: {
            DebugStop();
        }
            break;
        default:
            DebugStop();
    }
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer SettingFracturesSimple(const int caseToSim){
    
    // Fracture material
    TMRSDataTransfer sim_data;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["Fractures"] = EDomain;

    // NS: What are these?
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mIterface_material_idFracBound = 104;
    
    int D_Type = 0;
    int N_Type = 1;
    REAL pressure_in = 1.0 ;
    
    // Boundary conditions
    if (caseToSim == 0) {
        int bcfracid = EPressure;
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(bcfracid,D_Type,pressure_in);
    }
    else if (caseToSim == 1) {
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(3);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,pressure_in);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,0.);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,0.);
    }
    else if (caseToSim == 2){
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(3);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,pressure_in);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,0.);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,0.);
    }
    else if (caseToSim == 3) {
        int bcfracid = EPressure;
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(bcfracid,D_Type,pressure_in);
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
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    //FracAndReservoirProperties
    sim_data.mTFracProperties.m_Permeability = 0.00001;
    REAL kappa=1.0;
    int  id1=1;
    int  id2=2;
    std::vector<std::pair<int, REAL>> idPerm(2);
    idPerm[0]= std::make_pair(id1,kappa);
    idPerm[1]= std::make_pair(id2,kappa);
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

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = EDomain;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
        
    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase1(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = EDomain;

    // Fractures BCs
    const bool fluxThroughIntersect = false;
    if (fluxThroughIntersect){
        dim_name_and_physical_tagFine[1]["BCfrac0_0"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac0_1"] = EInlet;
        dim_name_and_physical_tagFine[1]["BCfrac0_2"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac0_3"] = EOutlet;
        dim_name_and_physical_tagFine[1]["BCfrac1_0"] = EInlet;
        dim_name_and_physical_tagFine[1]["BCfrac1_1"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac1_2"] = EOutlet;
        dim_name_and_physical_tagFine[1]["BCfrac1_3"] = ENoflux;
    }
    else {
        dim_name_and_physical_tagFine[1]["BCfrac0_0"] = EInlet;
        dim_name_and_physical_tagFine[1]["BCfrac0_1"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac0_2"] = EOutlet;
        dim_name_and_physical_tagFine[1]["BCfrac0_3"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac1_0"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac1_1"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac1_2"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac1_3"] = ENoflux;
    }
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;

    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase2(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0_0"] = EInlet;
    dim_name_and_physical_tagFine[1]["BCfrac0_1"] = ENoflux;
    dim_name_and_physical_tagFine[1]["BCfrac0_2"] = EOutlet;
    dim_name_and_physical_tagFine[1]["BCfrac0_3"] = ENoflux;
    
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    // No intersection in Flemisch example 1!

    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase3(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac2"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac3"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac4"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac5"] = EPressure;
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_2"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_4"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_8"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_1_2"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_8"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_4"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_3_4"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_8"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_4_5"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_8"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_5_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_6_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_8"] = EIntersection;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_7_8"] = EIntersection;

    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine){
            
    // Creating gmsh reader
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    GeometryFine.SetFormatVersion("4.1");
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh4(filename,nullptr,false);

    return gmeshFine;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
