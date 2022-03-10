// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>

// Unit test includes
#include <catch2/catch.hpp>

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(string& filenameFine, string& filenameCoarse, const bool isLinPVar);
TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void modifyBCsForLinPVar(TPZGeoMesh* gmesh);

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse, const bool isLinPVar);

const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname);

enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EFracture, EIntersection, EIntersectionEnd, EPLossAtIntersect,
    /*15 HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
};

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("cubicdomain");
#endif

// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("constant_pressure","[test_intersection_3D]"){
    string basemeshpath(FRACMESHES);
    string filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/twoElCoarse.msh";
    string filenameFine = basemeshpath + "/verificationMHMNoHybrid/twoElFine.msh";
    const bool isLinPVar = false;
    RunProblem(filenameFine,filenameCoarse,isLinPVar);
}
// ---- Test 0 ----
TEST_CASE("linear_pressure","[test_intersection_3D]"){
    string basemeshpath(FRACMESHES);
    string filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/intersectCoarse.msh";
    string filenameFine = basemeshpath + "/verificationMHMNoHybrid/intersectFine.msh";
    const bool isLinPVar = true;
    RunProblem(filenameFine,filenameCoarse,isLinPVar);
}

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
void RunProblem(string& filenamefine, string& filenamecoarse, const bool isLinPVar)
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    ReadMeshes(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse,isLinPVar);
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
    const int dimToPost = 3;
    mixAnalisys->PostProcessTimeStep(dimToPost);
    
    // ----- Compute integral of pressure and flux over domain and compare with analytical solution -----
    const std::string pvarname = "Pressure";
    const STATE integratedpressure = ComputeIntegralOverDomain(mixed_operator,pvarname);
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    const std::string qvarname = "Flux";
    STATE integratedflux = ComputeIntegralOverDomain(mixed_operator,qvarname);
    if (fabs(integratedflux) < 1.e-14 ) integratedflux = 0.; // to make Approx(0.) work
    std::cout << "\nintegral of flux  = " << integratedflux << std::endl;
    
    // ----- Comparing with analytical solution -----
    REQUIRE( integratedpressure == Approx( 8.0 ) ); // Approx is from catch2 lib
    if (isLinPVar) // linear pressure variation
        REQUIRE( integratedflux == Approx( 8./3. ) ); // Approx is from catch2 lib
    else // cte ppressyre
        REQUIRE( integratedflux == Approx( 0.) ); // Approx is from catch2 lib
    
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
    REAL zero_flux = 0.0, unit_pressure = 1.0, zero_pressure = 0., pressure_two = 2.;
    
    // Domain material
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Volume"] = EVolume;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,pressure_two);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,zero_pressure);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(EFaceBCPressure,D_Type,unit_pressure);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = EFracture;
    
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[1] = std::make_tuple(EFracInlet,D_Type,pressure_two);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[2] = std::make_tuple(EFracOutlet,D_Type,zero_pressure);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[3] = std::make_tuple(EPLossAtIntersect,Mixed_Type,0.5);
    
    sim_data.mTFracIntersectProperties.m_IntersectionPressureLossId = EPLossAtIntersect;
    
    // Simulation properties
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    // ReservoirProperties
    REAL kappa=1.0;
    int  id1=EVolume;
    std::vector<std::pair<int, REAL>> idPerm(1);
    idPerm[0]= std::make_pair(id1,kappa);
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    // Frac
    sim_data.mTFracProperties.m_Permeability = 1.;
    
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
                TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, const bool isLinPVar){
    
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
    
    if(isLinPVar){
        modifyBCsForLinPVar(gmeshfine);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void modifyBCsForLinPVar(TPZGeoMesh* gmesh) {
    
    const REAL cmin = -1, cmax = 1.;
    const int pos = 1;
    for (auto gel : gmesh->ElementVec()) {
        const int gelmatid = gel->MaterialId();
        if (gelmatid != EFaceBCPressure && gelmatid != EFracNoFlux)
            continue;
        
        TPZManVector<REAL,3> centqsi(2,0.), cent(3,0.);
        if (gelmatid == EFracNoFlux)
            centqsi.Resize(1, 0.);
        
        gel->CenterPoint(gel->NSides()-1, centqsi);
        gel->X(centqsi, cent);
        const REAL minsub = fabs(cent[pos] - cmin);
        const REAL maxsub = fabs(cent[pos] - cmax);
        if (minsub < ZeroTolerance()) {
            if(gelmatid == EFracNoFlux){
                gel->SetMaterialId(EFracInlet);
            }
            else{
                gel->SetMaterialId(EInlet);
            }
            
        }
        else if (maxsub < ZeroTolerance()) {
            if(gelmatid == EFracNoFlux){
                gel->SetMaterialId(EFracOutlet);
            }
            else {
                gel->SetMaterialId(EOutlet);
            }
            
        }
        else {
            if(gelmatid == EFracNoFlux){
                gel->SetMaterialId(EFracNoFlux);
            }
            else{
                gel->SetMaterialId(ENoflux);
            }
            
        }
    }
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

const STATE ComputeIntegralOverDomain(TPZCompMesh* cmesh, const std::string& varname) {
    std::set<int> matids;
    matids.insert(EVolume);
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences(); // compute integral in the multiphysics mesh
    TPZVec<STATE> vecint = cmesh->Integrate(varname, matids);
    if ((varname == "Pressure" && vecint.size() != 1) ||
        (varname == "Flux" && vecint.size() != 3)){
        DebugStop();
    }
    if (varname == "Pressure")
        return vecint[0];
    else if (varname == "Flux")
        return vecint[1];
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
