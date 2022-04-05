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
TMRSDataTransfer FillDataTransferCase1(TMRSDataTransfer& sim_data);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void fixPossibleMissingIntersections(TPZGeoMesh* gmesh);

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase1(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase3(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase4LF(string& filenameFine, string& filenameCoarse,
                               TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase4Debug(string& filenameFine, string& filenameCoarse,
                                  TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesWell(string& filenameFine, string& filenameCoarse,
                    TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesIP3D(string& filenameFine, string& filenameCoarse,
                    TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);


enum EMatid {/*0*/ENone, EVolume, EVolume2, EInlet, EOutlet, ENoflux,
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
    // 5: Flemisch case 4
    // 6: Flemisch case 4 with less fractures (and no overlap)
    // 7: Flemisch case 4 with much less fractures (for debugging)
    // 8: Well mesh (Initially idealized just for generating a beautiful mesh)
    // 9: IP3D mesh (Initially idealized just for generating a beautiful mesh)
    int simcase = 9;
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
        case 5:
            DebugStop(); // Need to generate mesh without overlap or need to treat overlap
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine.msh";
            break;
        case 6:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse_lf.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine_lf.msh";
            break;
        case 7:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse_debug.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine_debug.msh";
            break;
        case 8:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/wellmesh_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/wellmesh_fine.msh";
            break;
        case 9:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/ip3dmesh_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/ip3dmesh_fine.msh";
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
    auto start_time = std::chrono::steady_clock::now();
    
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    if (simcase < 2)
        ReadMeshes(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 2)
        ReadMeshesFlemischCase1(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 3)
        ReadMeshesFlemischCase2(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 4)
        ReadMeshesFlemischCase3(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 5)
        DebugStop();
    else if (simcase == 6)
        ReadMeshesFlemischCase4LF(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 7)
        ReadMeshesFlemischCase4Debug(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 8)
        ReadMeshesWell(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 9)
        ReadMeshesIP3D(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else
        DebugStop();
        
    fixPossibleMissingIntersections(gmeshfine); // read about this in the function
    
    
    TMRSDataTransfer sim_data;
    
    if (simcase == 2)
        FillDataTransferCase1(sim_data);
    else
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
    bool must_opt_band_width_Q = true; // This makes solving very fast!
    int n_threads = 8;
    bool UsingPzSparse = true; // Necessary to use multithread for now...
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    
    // -------------- Running problem --------------
    cout << "\n--------------------- Assembling ---------------------\n" << endl;
    cout << "Number of elements: " << mixed_operator->NElements() << endl;
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
    
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()/1000.;
    cout << "\n\n\t--------- Total time of simulation = " << total_time << " seconds -------\n" << endl;

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

TMRSDataTransfer FillDataTransferCase1(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, inlet_pressure = 4.0, outlet_pressure = 1.0;
    
    // Domain material
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Volume2"] = EVolume2;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,inlet_pressure);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,outlet_pressure);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
            
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
    sim_data.mTFracProperties.m_Permeability = 1.0e-3;
    
    REAL kappa1=1.0e-5;
    REAL kappa2=1.0e-6;
    int id1 = EVolume2;
    int id2 = EVolume;
    std::vector<std::pair<int, REAL>> idPerm(2);
    idPerm[0]= std::make_pair(id1,kappa1);
    idPerm[1]= std::make_pair(id2,kappa2);
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
    for (int ivol = 33; ivol <= 36; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume2;
    }
    for (int ivol = 37; ivol <= 44; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }

    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc4"] = EInlet;
    dim_name_and_physical_tagCoarse[2]["bc5"] = EOutlet;
    dim_name_and_physical_tagCoarse[2]["bc6"] = ENoflux;
        
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
    dim_name_and_physical_tagFine[2]["bc4"] = EInlet;
    dim_name_and_physical_tagFine[2]["bc5"] = EOutlet;
    dim_name_and_physical_tagFine[2]["bc6"] = ENoflux;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EFracNoFlux;
    
    // Intersections
    // No intersections in this case!
            
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

void ReadMeshesFlemischCase4LF(string& filenameFine, string& filenameCoarse,
                               TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    const int initc = 1, endc = 1000;
    for (int ivol = initc; ivol <= endc; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
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
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 36; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
//    dim_name_and_physical_tagFine[1]["fracIntersection_1_25"] = 10025;
//    dim_name_and_physical_tagFine[1]["fracIntersection_1_31"] = 10031;
//    dim_name_and_physical_tagFine[1]["fracIntersection_2_25"] = 20025;
//    dim_name_and_physical_tagFine[1]["fracIntersection_3_5" ] = 30005;
//    dim_name_and_physical_tagFine[1]["fracIntersection_4_25"] = 40025;
//    dim_name_and_physical_tagFine[1]["fracIntersection_5_25"] = 50025;
//    dim_name_and_physical_tagFine[1]["fracIntersection_6_19"] = 60019;
//    dim_name_and_physical_tagFine[1]["fracIntersection_6_29"] = 60029;
//    dim_name_and_physical_tagFine[1]["fracIntersection_7_11"] = 70011;
//    dim_name_and_physical_tagFine[1]["fracIntersection_8_34"] = 80034;
//    dim_name_and_physical_tagFine[1]["fracIntersection_17_22"] =170022;
//    dim_name_and_physical_tagFine[1]["fracIntersection_25_36"] =250036;
//    dim_name_and_physical_tagFine[1]["fracIntersection_26_27"] =260027;
//    dim_name_and_physical_tagFine[1]["fracIntersection_29_32"] =290032;
//    dim_name_and_physical_tagFine[1]["fracIntersection_30_33"] =300033;
//    dim_name_and_physical_tagFine[1]["fracIntersection_30_36"] =300036;
//    dim_name_and_physical_tagFine[1]["fracIntersection_31_33"] =310033;
//    dim_name_and_physical_tagFine[1]["fracIntersection_31_36"] =310036;
//    dim_name_and_physical_tagFine[1]["fracIntersection_34_35"] =340035;
    
    dim_name_and_physical_tagFine[1]["fracIntersection_1_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_1_31"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_2_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_3_5" ] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_25"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_19"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_29"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_7_11"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_8_34"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_17_22"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_25_36"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_26_27"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_29_32"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_30_33"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_30_36"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_31_33"] = -121321313;
    dim_name_and_physical_tagFine[1]["fracIntersection_31_36"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_34_35"] = EIntersection;

    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesFlemischCase4Debug(string& filenameFine, string& filenameCoarse,
                                  TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    const int initc = 1, endc = 1000;
    for (int ivol = initc; ivol <= endc; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
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
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 19; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_4_7"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_8"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_5_11"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_8"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_6_11"] =  EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_9_10"] =  EIntersection;


    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesWell(string& filenameFine, string& filenameCoarse,
                    TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    
    // Domain
    std::string volbase = "c";
    const int initc = 1, endc = 6;
    for (int ivol = initc; ivol <= endc; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
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
    dim_name_and_physical_tagFine[2]["bc1"] = EFaceBCPressure;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture12"] = EFracture;
    dim_name_and_physical_tagFine[2]["Fracture13"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 1; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
    // No intersections in this mesh

    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesIP3D(string& filenameFine, string& filenameCoarse,
                    TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmesh) {
    
    string basemeshpath(FRACMESHES);
    // ===================> Coarse mesh <=======================
    gmesh = new TPZGeoMesh;
    const int nnodes = 32 + 1;
    gmesh->NodeVec().Resize(nnodes);
    int id = 1;
    std::ifstream read(basemeshpath + "/verificationMHMNoHybrid/ip3dmesh_nodes.txt");
    string str;
    REAL c1, c2, c3;
    TPZManVector<REAL,3> coord(3,0.);
    TPZGeoNode nodebog(id++,coord,*gmesh);
    gmesh->NodeVec()[0] = nodebog;
    int inod = 1;
    read >> coord[0];
    while (read) {
        cout << "inod = " << inod << endl;
        read >> coord[1];
        read >> coord[2];
        TPZGeoNode node(id++,coord,*gmesh);
        gmesh->NodeVec()[inod] = node;
        inod++;
        read >> coord[0];
    }
    
    int64_t index;
    // ========> Outer cubes
    const int outercubematid = 1;
    TPZManVector<int64_t,8> cubeind = {17,1,5,20,18,3,7,19};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, outercubematid, index);
    cubeind = {2,21,24,6,4,22,23,8};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, outercubematid, index);
    cubeind = {3,7,31,32,4,8,30,29};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, outercubematid, index);
    cubeind = {1,25,26,5,2,28,27,6};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, outercubematid, index);
    
    // ========> Inside cubes
    const int innercubematid = 2;
    cubeind = {5,7,3,1,13,15,11,9};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, innercubematid, index);
    cubeind = {5,6,8,7,13,14,16,15};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, innercubematid, index);
    cubeind = {6,8,4,2,14,16,12,10};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, innercubematid, index);
    cubeind = {4,3,1,2,12,11,9,10};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, innercubematid, index);
    cubeind = {6,5,1,2,14,13,9,10};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, innercubematid, index);
    cubeind = {8,7,3,4,16,15,11,12};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, innercubematid, index);

    // ========> Cube well
    const int wellcubematid = 10;
    cubeind = {13,15,11,9,14,16,12,10};
    gmesh->CreateGeoElement(MElementType::ECube, cubeind, wellcubematid, index);
    
    // ========> Prisms
    const int prismmatid = 3;
    TPZManVector<int64_t,6> prismind = {17,25,1,20,26,5};
    gmesh->CreateGeoElement(MElementType::EPrisma, prismind, prismmatid, index);
    prismind = {27,24,6,28,21,2};
    gmesh->CreateGeoElement(MElementType::EPrisma, prismind, prismmatid, index);
    prismind = {4,22,29,8,23,30};
    gmesh->CreateGeoElement(MElementType::EPrisma, prismind, prismmatid, index);
    prismind = {7,31,19,3,32,18};
    gmesh->CreateGeoElement(MElementType::EPrisma, prismind, prismmatid, index);
    
    // just so it works...
    gmeshfine = gmesh;
    
    // ===================> Fine mesh <=======================
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void fixPossibleMissingIntersections(TPZGeoMesh* gmesh){
    // it may happen that the supplied mesh from DFN has missing intersections where 2 fractures intersect.
    // We treat this here. But it would be ideal to make DFN robust enough so we could remove this function whatsoever
    // The main idea is to check if a fracture element has more than 1 fracture neightbor. If so, an intersection element is
    // needed to hybridzie the region.
    cout << "\n\t==================== Searching for problematic intersection elements ====================\n" << endl;
    
    int nInterCreated = 0;
    for (auto gel : gmesh->ElementVec()) {
        if (!gel || gel->MaterialId() != EFracture) {
            continue;
        }
        if (gel->Dimension() != 2) {
            DebugStop(); // Should it work with 2D problems? (1d fractures)
        }
        
        const int firstedge = gel->FirstSide(1);
        const int lastedge = gel->FirstSide(2);
        for (int iside = firstedge; iside < lastedge; iside++) {
            TPZStack<TPZGeoEl*> interEls;
            int nFracElsForSide = 0;
            TPZGeoElSide gside(gel,iside);
            TPZGeoElSide neig = gside.Neighbour();
            for (; gside != neig; neig++) {
                TPZGeoEl* neigel = neig.Element();
                if (neigel->MaterialId() == EIntersection) {
                    interEls.push_back(neigel);
                }
                if (neigel->MaterialId() == EFracture) {
                    nFracElsForSide++;
                }
            }
            if (interEls.size() > 1) {
                if (interEls.size() > 2) {
                    DebugStop(); // there are 3 intersection elements in the same place! Please check why...
                }
                cout << "Found two intersection elements at the same place!" << endl;
                cout << "Manually deleting intersection geoel..." << endl;
                TPZGeoEl* gelduplicate = interEls[1];
                const int64_t duplicateIndex = gelduplicate->Index();
                gelduplicate->RemoveConnectivities();
                delete gelduplicate;
                gmesh->ElementVec()[duplicateIndex] = nullptr;
            }
            if (nFracElsForSide > 1 && !interEls.size()) {
                cout << "nfracs for this side: " << nFracElsForSide << endl;
                cout << "Manually creating intersection geoel..." << endl;
                TPZGeoElBC(gel, iside, EIntersection);
                nInterCreated++;
            }
            
        }
    }
    if (nInterCreated) {
        cout << "\n\t\t==================== MESSAGE =================" << endl;
        cout << "\n\t- Imported mesh has " << nInterCreated << " missing intersection elements" << endl;
        cout << "\t- These were created here and this shouldn't be a problem...\n" << endl;
    }
        
    gmesh->BuildConnectivity();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
