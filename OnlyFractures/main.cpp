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
void ChangeMeshToImposePressureOnIntersection(TMRSApproxSpaceGenerator& aspace,TPZMultiphysicsCompMesh *mixed_operator);
void FixGMeshForIntersectionOverBCs(TPZGeoMesh* gmesh);

auto exactSol = [](const TPZVec<REAL> &loc,
                   TPZVec<STATE>&u,
                   TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    u[0]= 1-x;
    gradU(0,0) = 0.;
//    gradU(1,0) = 0.;
};

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase1(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase2(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase3(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase4(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase5(std::string &filename);

enum EMatid {ENone, EDomain, EInlet, EOutlet, ENoflux, EPressure, EIntersection, EIntersectionEnd};
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

static TPZLogger mainlogger("onlyfractures");

int main(){
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for OnlyFractures target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
        
    const int caseToSim = 0;
    // 0: 2 perpendicular fractures, cte pressure
    // 1: 2 perpendicular fractures, 1D flow
    // 2: Flemisch example 1
    // 3: Flemisch example 2
    // 4: Flemisch example 3
    // 5: Flemisch example 4
    // 6: inclined fracture
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
    
    const bool isCtePressVariation = true;
    
    const bool printgmesh = true;
    if (printgmesh) {
        std::ofstream name("GeoMeshFractures.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name);
    }
        
    TMRSApproxSpaceGenerator aspace;
    if (isCtePressVariation) {
        aspace.SetForcingFunctionBC(exactSol);
    }
    
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
    
    if (true) {
        ChangeMeshToImposePressureOnIntersection(aspace, mixed_operator);
    }
        
    // Analysis parameters
    bool must_opt_band_width_Q = false;
    int n_threads = 0;
    bool UsingPzSparse = false;
    bool UsePardiso_Q = true;
    
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
    std::ofstream outP("outp.txt");
    mixed_operator->MeshVector()[1]->Print(outP);
    std::ofstream outF("outf.txt");
    mixed_operator->MeshVector()[0]->Print(outF);
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
//            std::string filename = basemeshpath + "/2DMeshes/flemisch2.msh";
            std::string filename = basemeshpath + "/2DMeshes/flemisch2_new.msh";
            gmesh = ReadFractureMeshCase3(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 4: {
//            std::string filename = basemeshpath + "/2DMeshes/flemisch3_4frac_2.msh";
            std::string filename = basemeshpath + "/2DMeshes/flemisch3_new.msh";
            gmesh = ReadFractureMeshCase4(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
            
        }
            break;
        case 5: {
            std::string filename = basemeshpath + "/2DMeshes/flemisch4.msh";
            gmesh = ReadFractureMeshCase5(filename); // same names
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 6: {
            std::string filename = basemeshpath + "/2DMeshes/inclinedxy.msh";
            gmesh = ReadFractureMeshCase0(filename); // same names
            sim_data  = SettingFracturesSimple(caseToSim);
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
    else if (caseToSim == 3 || caseToSim == 4 || caseToSim == 5 || caseToSim == 6) {
        int bcfracid = EPressure;
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(bcfracid,D_Type,pressure_in);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EInlet,D_Type,pressure_in);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(EOutlet,D_Type,0.);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(ENoflux,N_Type,0.);
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
    
    const bool isCtePressure = true;
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture12"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture13"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture14"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture15"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture16"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture17"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture18"] = EDomain;
    
    
    // Fractures BCs
    if (isCtePressure) {
        dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac2"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac3"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac4"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac5"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac6"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac7"] = EPressure;
        dim_name_and_physical_tagFine[1]["BCfrac8"] = EPressure;
    }
    else{
        dim_name_and_physical_tagFine[1]["BCfrac0"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac1"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac2"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac3"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac4"] = ENoflux;
        dim_name_and_physical_tagFine[1]["BCfrac5"] = ENoflux;
    }
    
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
    
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    /*
    if (!isCtePressure) {
        const int icoor = 0;
        const REAL tol = 1.e-10, coor0 = 0., coor1 = 1.;
        for (auto* gel : gmesh->ElementVec()){
            const int geldim = gel->Dimension();
            const int gelmatid = gel->MaterialId();
            if (geldim == 1 and gelmatid == ENoflux) {
                int count0 = 0, count1 = 0;
                for (int i = 0; i < gel->NNodes(); i++) {
                    TPZGeoNode* gnod = gel->NodePtr(i);
                    const REAL dif0 = fabs(gnod->Coord(icoor) - coor0);
                    const REAL dif1 = fabs(gnod->Coord(icoor) - coor1);
                    if (dif0 < tol) count0++;
                    if (dif1 < tol) count1++;
                }
                if (count0 == 2)
                    gel->SetMaterialId(EInlet);
                if (count1 == 2)
                    gel->SetMaterialId(EOutlet);
            }
        }
    }
     */
//    for (auto* gel : gmesh->ElementVec()){
//        if (gel->MaterialId() != EIntersection) {
//            continue;
//        }
//        TPZGeoElBC *gbc = new TPZGeoElBC(gel,gel->NSides()-1,EPressure);
//    }
//    gmesh->SetDimension(2);
//    gmesh->BuildConnectivity();
    
    FixGMeshForIntersectionOverBCs(gmesh);
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase4(std::string &filename){
    
    const bool isCtePressure = true;
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture12"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture13"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture14"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture15"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture16"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture17"] = EDomain;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac2"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac3"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac4"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac5"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac6"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac7"] = EPressure;
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_7"] = EIntersection;
    
    // for isolated fracture cases
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_2"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    
    FixGMeshForIntersectionOverBCs(gmesh);
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase5(std::string &filename){

    // base names
    std::string fracbase = "Fracture";
    const int fracidskip = 10;
    std::string bcfracbase = "BCfrac";
    std::string fracintersecbase = "fracIntersection_";
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    for (int ifrac = 0; ifrac < 52; ifrac++) {
        std::string ifracstring = fracbase + to_string(ifrac+fracidskip);
        dim_name_and_physical_tagFine[2][ifracstring] = EDomain;
    }
//    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;

    // Fractures BCs

    for (int ifrac = 0; ifrac < 52; ifrac++) {
        std::string ifracbcstring = bcfracbase + to_string(ifrac);
        dim_name_and_physical_tagFine[1][ifracbcstring] = EPressure;
    }
    
    for (int ifrac = 0; ifrac < 52; ifrac++) {
        for (int jfrac = 0; jfrac < 52; jfrac++) {
            std::string ifracintersecstring = fracintersecbase + to_string(ifrac) + "_" + to_string(jfrac);
            dim_name_and_physical_tagFine[1][ifracintersecstring] = EIntersection;
        }
    }
    
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    
    FixGMeshForIntersectionOverBCs(gmesh);
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


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

void ChangeMeshToImposePressureOnIntersection(TMRSApproxSpaceGenerator& aspace,TPZMultiphysicsCompMesh *mixed_operator) {
    
    if (!aspace.Hybridizer()) {
        return;
    }
    
  {
    const int lagrangematidend = aspace.Hybridizer()->lagrangeInterfaceEndMatId();
    
    TPZMaterial *mat = mixed_operator->FindMaterial(lagrangematidend);
    if (!mat) {
        DebugStop();
    }
    
    mixed_operator->DeleteMaterial(lagrangematidend);
    
    TPZMaterial *matdomain = mixed_operator->FindMaterial(EDomain);
    if (!matdomain)
        DebugStop();
    
    TMRSDarcyFlowWithMem<TMRSMemory> *matdf = dynamic_cast<TMRSDarcyFlowWithMem<TMRSMemory>* >(matdomain);
    if (!matdf) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,1.);
    TPZBndCondT<REAL>* bnd = matdf->CreateBC(matdf, lagrangematidend, 0, val1, val2);
    if (aspace.HasForcingFunctionBC()){
        bnd->SetForcingFunctionBC(exactSol);
    }
    mixed_operator->InsertMaterialObject(bnd);
    
    mixed_operator->LoadReferences();
  }
   
    
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FixGMeshForIntersectionOverBCs(TPZGeoMesh* gmesh) {
    for (auto gel : gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != EIntersection) {
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.HasNeighbour(EPressure);
        while (neigh) {
            gel->SetMaterialId(EIntersectionEnd);
            neigh.Element()->SetMaterialId(ENone);
            neigh = neigh.HasNeighbour(EPressure);
        }
    }
    gmesh->BuildConnectivity();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


TPZGeoMesh *ReadFractureMeshCase6(std::string &filename){
    
    DebugStop(); // Fix names!!!
    
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture12"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture13"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture14"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture15"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture16"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture17"] = EDomain;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac2"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac3"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac4"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac5"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac6"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac7"] = EPressure;
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_6"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_7"] = EIntersection;
    
    // for isolated fracture cases
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_2"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    
    FixGMeshForIntersectionOverBCs(gmesh);
    
    return gmesh;
}

