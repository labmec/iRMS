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
void RunProblem(const int& caseToSim);
TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data, int const simcase);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void modifyBCsForInletOutlet(TPZGeoMesh* gmesh);

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse, int const  caseToSim);
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
TEST_CASE("Two_El_One_Frac_cte_pressure","[test_frac_4spaces_3D]"){
    RunProblem(0);
}
// ---- Test 1 ----
TEST_CASE("Two_El_One_Frac_lin_pressure","[test_frac_4spaces_3D]"){
    RunProblem(1);
}
// ---- Test 2 ----
TEST_CASE("Four_El_Two_Frac_cte_pressure","[test_frac_4spaces_3D]"){
    RunProblem(2);
}
// ---- Test 3 ----
TEST_CASE("Four_El_Two_Frac_lin_pressure","[test_frac_4spaces_3D]"){
    RunProblem(3);
}
// ---- Test 4 ----
TEST_CASE("Two_El_One_Frac_cte_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(4);
}
// ---- Test 5 ----
TEST_CASE("Two_El_One_Frac_lin_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(5);
}
// ---- Test 6 ----
TEST_CASE("Four_El_Two_Frac_cte_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(6);
}
// ---- Test 7 ----
TEST_CASE("Four_El_Two_Frac_lin_pressure_Transp","[test_frac_4spaces_3D]"){
    RunProblem(7);
}

// Case0: TwoElements, one fracture constant pressure
// Case1: TwoElements, one fracture linear pressure
// Case2: FourElements, two intersecting fractures constant pressure
// Case3: FourElements, two intersecting fractures linear pressure

// Case4: TwoElements, one fracture constant pressure with tranport
// Case5: TwoElements, one fracture linear pressure with tranport
// Case6: FourElements, two intersecting fractures constant pressure with tranport
// Case7: FourElements, two intersecting fractures linear pressure with tranport

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------

//int main(){
//    int CaseToSim = 7;
//    RunProblem(CaseToSim);
//    return 0;
//}
void RunProblem( const int &caseToSim){
    string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for Cubic Domain problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    string filenameCoarse,filenameFine;
    if (caseToSim==0 || caseToSim==1 || caseToSim==4 || caseToSim==5) {
        //OneFractures
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/TwoElsUnitTestCoarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/TwoElsUnitTestFine.msh";
    }else if (caseToSim==2 || caseToSim==3 || caseToSim==6 || caseToSim==7){
        filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/intersectCoarse.msh";
        filenameFine = basemeshpath + "/verificationMHMNoHybrid/intersectFine.msh";
    }
    else{
        DebugStop();
    }
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    ReadMeshes(filenameFine,filenameCoarse,gmeshfine,gmeshcoarse,caseToSim);
    TMRSDataTransfer sim_data;
    FillDataTransfer(sim_data, caseToSim);
    
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
//  aspace.SetGeometry(gmeshfine);
    
    // ----- Setting the global data transfer -----
    aspace.SetDataTransfer(sim_data);
    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    std::ofstream name("MixedOperator.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(mixed_operator, name);
    
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = false;
    int n_threads = 0;
    bool UsingPzSparse = true;
    bool UsePardiso_Q = true;
    
    // ----- Setting analysis -----
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    if(caseToSim == 4 || caseToSim == 5 || caseToSim == 6 || caseToSim == 7){
        aspace.BuildAuxTransportCmesh();
        TPZCompMesh * transport_operator = aspace.GetTransportOperator();
        TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
        sfi_analysis->SetDataTransfer(&sim_data);
        sfi_analysis->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
        int n_steps = sim_data.mTNumerics.m_n_steps;
        REAL dt = sim_data.mTNumerics.m_dt;
        
        TPZStack<REAL,100> reporting_times;
        reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
        REAL sim_time = 0.0;
        int pos =0;
        REAL current_report_time = reporting_times[pos];
        int npos = reporting_times.size();
        
        sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
        REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        std::cout << "Mass report at time : " << 0.0 << std::endl;
        std::cout << "Mass integral :  " << initial_mass << std::endl;
        std::ofstream fileCilamce("IntegratedSat.txt");
        TPZFastCondensedElement::fSkipLoadSolution = false;
        bool first=true;
        for (int it = 1; it <= n_steps; it++) {
            sim_time = it*dt;
            sfi_analysis->m_transport_module->SetCurrentTime(dt);
            sfi_analysis->PostProcessTimeStep(2);
            sfi_analysis->RunTimeStep(); // runs mixed and transport problems
            mixed_operator->LoadSolution(mixed_operator->Solution());
            if (sim_time >=  current_report_time) {
                std::cout << "Time step number:  " << it << std::endl;
                std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
                mixed_operator->UpdatePreviousState(-1.);
                sfi_analysis->PostProcessTimeStep();
                pos++;
                current_report_time =reporting_times[pos];
                
                // computes integral of saturation for elements of certain matid
                REAL InntMassFrac=sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(aspace.FractureMatId());
                fileCilamce<<current_report_time/(86400*365)<<", "<<InntMassFrac<<std::endl;
               
                REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
                std::cout << "Mass report at time : " << sim_time << std::endl;
                std::cout << "Mass integral :  " << mass << std::endl;
            }
        }
    }
    else{
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
    }
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
    if (caseToSim == 1 || caseToSim == 3 || caseToSim == 5 || caseToSim == 7) // linear pressure variation
        REQUIRE( integratedflux == Approx( 8./3. ) ); // Approx is from catch2 lib
    if (caseToSim == 0 || caseToSim == 2 || caseToSim == 4 || caseToSim == 6)  // cte ppressyre
        REQUIRE( integratedflux == Approx( 0.) ); // Approx is from catch2 lib
    
    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TMRSDataTransfer FillDataTransfer(TMRSDataTransfer& sim_data, int const simcase){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, unit_pressure = 1.0, zero_pressure = 0., pressure_two = 2.;
    
    // Domain material
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = 10;
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[1]["EIntersection"] = 11;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(EInlet,D_Type,pressure_two);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,zero_pressure);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(EFaceBCPressure,D_Type,unit_pressure);
    
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[1] = std::make_tuple(EFracInlet,D_Type,pressure_two);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[2] = std::make_tuple(EFracOutlet,D_Type,zero_pressure);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[3] = std::make_tuple(EPLossAtIntersect,Mixed_Type,0.5);
    sim_data.mTFracIntersectProperties.m_IntersectionPressureLossId = EPLossAtIntersect;
    
    if (simcase==0 || simcase==2 || simcase==4 || simcase==6 ) {
        sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(EOutlet,D_Type,pressure_two);
        sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[2] = std::make_tuple(EFracOutlet,D_Type,pressure_two);
    }
    

    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mIterface_material_idFracBound = 104;
    
//    sim_data.mTFracIntersectProperties.m_IntersectionId = EPLossAtIntersect;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(5);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(EOutlet,D_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(EInlet,N_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(ENoflux,N_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(EFracOutlet,D_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[4] = std::make_tuple(EFracInlet,D_Type,1.);
    sim_data.mTGeometry.mIterface_material_idFracBound = EFracNoFlux;
    
 
    
    
    // Other properties
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
    sim_data.mTNumerics.m_n_steps = 10;
    sim_data.mTNumerics.m_dt      = 0.5;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
    //FracAndReservoirProperties
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
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    
    return sim_data;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int const caseToSim){
    
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
    modifyBCsForInletOutlet(gmeshfine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void modifyBCsForInletOutlet(TPZGeoMesh* gmesh) {
    
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
