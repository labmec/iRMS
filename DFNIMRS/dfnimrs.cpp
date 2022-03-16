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
void RunProblem(string& filenameFine, string& filenameCoarse, const int simcase);
TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase1(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase3(string& filenameFine, string& filenameCoarse,
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
    string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
    string logpath = basemeshpath + "/../DFNIMRS/log4cxx.cfg";
    TPZLogger::InitializePZLOG(logpath);
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    
    // 0: two elements, 1 frac
    // 1: 4 elements, 2 frac, w/ intersection
    // 2: Flemisch case 1
    // 3: Flemisch case 2
    // 4: Flemisch case 3
    int simcase = 0;
    string filenameCoarse, filenameFine;
    switch (simcase) {
        case 0:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/twoElCoarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/twoElFine.msh";
            break;
        case 1:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/intersectCoarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/intersectFine.msh";
            break;
        case 2:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case1_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case1_fine.msh";
            break;
        case 3:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case2_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case2_fine.msh";
            break;
        case 4:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case3_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case3_fine.msh";
            break;
        default:
            break;
    }
    RunProblem(filenameFine,filenameCoarse,simcase);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
void RunProblem(string& filenamefine, string& filenamecoarse, const int simcase)
{
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    if (simcase < 2) {
        ReadMeshes(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    }
    else if (simcase == 2){
        ReadMeshesFlemischCase1(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    }
    else if (simcase == 3){
        ReadMeshesFlemischCase2(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    }
    else if (simcase == 4){
        ReadMeshesFlemischCase3(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    }
    else {
        DebugStop();
    }
    
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
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
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
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    cout << "\n--------------------- Assembling ---------------------\n" << endl;
    cout << "Number of equations: " << mixed_operator->NEquations() << endl;
    mixAnalisys->Assemble();

    cout << "\n--------------------- Solving ---------------------\n" << endl;
    mixAnalisys->Solve();
    
    // The system is solve as non linear, so have to multiply by -1
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false; // So we can postprocess variables correctly
    mixed_operator->LoadSolution(mixed_operator->Solution());
    
    // ----- Post processing -----
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    int dimToPost = 2;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    dimToPost = 3;
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
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    
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
    std::string volbase = "c";
    for (int ivol = 1; ivol <= 4; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 1; ivol <= 4; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
    
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = EFracture;
    
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EFracture;

    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EFracNoFlux;
    
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
            
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

void ReadMeshesFlemischCase1(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    for (int ivol = 33; ivol <= 44; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc6"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
    
    
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 33; ivol <= 44; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-33);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc6"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux;
    
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
            
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    for (int ivol = 1; ivol <= 512; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
        
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = 1; ivol <= 512; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-1);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 8; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
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
            
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesFlemischCase3(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    const int initc = 229, endc = 444;
    for (int ivol = initc; ivol <= endc; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagCoarse[2]["bc6"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
        
    // ===================> Fine mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
            
    // Domain
    for (int ivol = initc; ivol <= endc; ivol++) {
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagFine[3][ivolstring] = EInitVolumeMatForMHM + (ivol-initc);
//        dim_name_and_physical_tagFine[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagFine[2]["bc4"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc5"] = EFaceBCPressure;
    dim_name_and_physical_tagFine[2]["bc6"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture11"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture12"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture13"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture14"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture15"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture16"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture17"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 7; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;

    dim_name_and_physical_tagFine[1]["fracIntersection_0_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_3"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_7"] = EIntersection;
    
    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
