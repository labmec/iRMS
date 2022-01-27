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

void VerificationCases(const int caseToSim);
TMRSDataTransfer SettingFracturesSimple(const int caseToSim);
void CreateGMeshAndDataTransfer(TPZGeoMesh*& gmesh,TMRSDataTransfer &sim_data, const int caseToSim);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void ChangeMeshToImposePressureOnIntersection(TMRSApproxSpaceGenerator& aspace,TPZMultiphysicsCompMesh *mixed_operator);
void FixGMeshForIntersectionOverBCs(TPZGeoMesh* gmesh);
void ChangeBCsToNoFlux(TPZGeoMesh* gmesh);

auto exactSol = [](const TPZVec<REAL> &loc,
                   TPZVec<STATE>&u,
                   TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
//    u[0]= 1-x;
    u[0]= 1-z;
    gradU(0,0) = 0.;
//    gradU(1,0) = 0.;
};

TPZGeoMesh *ReadFractureMeshCase0(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase1(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase2(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase3(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase4(std::string &filename);
TPZGeoMesh *ReadFractureMeshCase5(std::string &filename);

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

#ifdef PZ_LOG
static TPZLogger mainlogger("onlyfractures");
#endif

int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for OnlyFractures target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
        
    const int caseToSim = 0;
    // 0: 2 perpendicular fractures, cte pressure
    // 1: 2 fractures, aligned, 1D flow
    // 2: 2 fractures, slightly not aligned, 1D flow
    // 3: 2 fractures, slightly not aligned, 1D flow, no matching frac-domain boundary
    // 4: Flemisch example 3
    // 5: Flemisch example 4
    // 6: inclined fracture
    VerificationCases(caseToSim);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void VerificationCases(const int caseToSim)
{
    // Reading ONLY the fractures of mesh from DFN
    TPZGeoMesh *gmesh = nullptr;
    TMRSDataTransfer sim_data;
    
    CreateGMeshAndDataTransfer(gmesh,sim_data,caseToSim);
    
    const bool isCtePressVariation = false;
    
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
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
//    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E2Space;
    
    // Setting gmesh
    aspace.SetGeometry(gmesh);
    
    // Setting the global data transfer
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
    aspace.SetDataTransfer(sim_data);
        
    // Creates de multiphysics compmesh
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    if (false) {
        ChangeMeshToImposePressureOnIntersection(aspace, mixed_operator);
    }
        
    // Analysis parameters
    bool must_opt_band_width_Q = true;
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
//    std::ofstream outP("outp.txt");
//    mixed_operator->MeshVector()[1]->Print(outP);
//    std::ofstream outF("outf.txt");
//    mixed_operator->MeshVector()[0]->Print(outF);
    
    const int dimToPost = 3;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    
    const std::string varname = "state";
    std::set<int> matids;
    matids.insert(EInlet);
//    matids.insert(EOutlet);
    
    TPZCompMesh *fluxmesh = mixed_operator->MeshVector()[0];
    fluxmesh->Reference()->ResetReference();
    fluxmesh->LoadReferences();
    TPZVec<STATE> vecint = fluxmesh->Integrate(varname, matids);
    if (vecint.size())
        std::cout << "\nint inlet = " << vecint[0] << std::endl;
    matids.clear();
    matids.insert(EOutlet);
    vecint = fluxmesh->Integrate(varname, matids);
    if (vecint.size())
        std::cout << "\nint outlet = " << vecint[0] << std::endl;
    
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
            std::string filename = basemeshpath + "/2DMeshes/2parallel.msh";
            gmesh = ReadFractureMeshCase1(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 2: {
            std::string filename = basemeshpath + "/verifications/2fracNotAlignedcut.msh";
            gmesh = ReadFractureMeshCase2(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 3: {
//            std::string filename = basemeshpath + "/2DMeshes/flemisch2.msh";
//            std::string filename = basemeshpath + "/verifications/verif0.msh";
            std::string filename = basemeshpath + "/verifications/2blocking.msh";
            gmesh = ReadFractureMeshCase3(filename);
            sim_data  = SettingFracturesSimple(caseToSim);
        }
            break;
        case 4: {
//            std::string filename = basemeshpath + "/2DMeshes/flemisch3_4frac_2.msh";
//            std::string filename = basemeshpath + "/2DMeshes/flemisch3_new.msh";
            std::string filename = basemeshpath + "/verifications/flcase3.msh";
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
    if (caseToSim < 4) {
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
    else if (caseToSim < 5) {
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EFaceBCPressure,D_Type,pressure_in);
        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(1);
        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EPressure,D_Type,pressure_in);
    }
    else if (caseToSim == 10) {
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
//    int  id2=2;
    std::vector<std::pair<int, REAL>> idPerm(1);
    idPerm[0]= std::make_pair(id1,kappa);
//    idPerm[1]= std::make_pair(id2,kappa);
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
    dim_name_and_physical_tagFine[3]["c1"] = EVolume;
    dim_name_and_physical_tagFine[3]["c2"] = EVolume;
    dim_name_and_physical_tagFine[3]["c3"] = EVolume;
    dim_name_and_physical_tagFine[3]["c4"] = EVolume;
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
        
    return generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase1(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[3]["c1"] = EVolume;
    dim_name_and_physical_tagFine[3]["c2"] = EVolume;
    dim_name_and_physical_tagFine[3]["c3"] = EVolume;
    dim_name_and_physical_tagFine[3]["c4"] = EVolume;
    dim_name_and_physical_tagFine[3]["c5"] = EVolume;
    dim_name_and_physical_tagFine[3]["c6"] = EVolume;
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture11"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    
    
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    ChangeBCsToNoFlux(gmesh);
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase2(std::string &filename){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Fractures
    dim_name_and_physical_tagFine[3]["c1"] = EVolume;
    dim_name_and_physical_tagFine[3]["c2"] = EVolume;
    dim_name_and_physical_tagFine[3]["c3"] = EVolume;
    dim_name_and_physical_tagFine[3]["c4"] = EVolume;
    dim_name_and_physical_tagFine[3]["c5"] = EVolume;
    dim_name_and_physical_tagFine[3]["c6"] = EVolume;
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture11"] = globFracID;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    
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

    std::string volbase = "c";

    // Volumes
    for (int ivol = 1; ivol < 28; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // BC domain
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture11"] = globFracID;
    
    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EPressure;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EPressure;
    
    
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    ChangeBCsToNoFlux(gmesh);
//    FixGMeshForIntersectionOverBCs(gmesh);
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadFractureMeshCase4(std::string &filename){
    
    const bool isCtePressure = true;
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D

    // Volumes
    std::string volbase = "c";
    for (int ivol = 229; ivol < 445; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // BC domain
    dim_name_and_physical_tagFine[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc6"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture11"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture12"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture13"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture14"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture15"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture16"] = globFracID;
    dim_name_and_physical_tagFine[2]["Fracture17"] = globFracID;

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
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_5_6"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_6_7"] = EIntersection;
    
    // for isolated fracture cases
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_3"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_1_2"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_1_3"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    
    TPZGeoMesh* gmesh = generateGMeshWithPhysTagVec(filename,dim_name_and_physical_tagFine);
    
//    FixGMeshForIntersectionOverBCs(gmesh);
    
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
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);

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


void ChangeBCsToNoFlux(TPZGeoMesh* gmesh) {
    for(auto gel : gmesh->ElementVec()){
        if(gel->Dimension() != 2) continue;
        if(gel->MaterialId() != EFaceBCPressure) continue;

        const int pos = 1;
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
