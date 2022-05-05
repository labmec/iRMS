// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>
#include "json.hpp"

// ----- Namespaces -----
using namespace std;
// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(string& filenameFine, string& filenameCoarse, const string &filenamejson, const int simcase);
void FillDataTransfer(TMRSDataTransfer& sim_data);
void FillDataTransferCase1(TMRSDataTransfer& sim_data);
void FillDataTransferCase2(TMRSDataTransfer& sim_data);
void FillDataTransferCase3(TMRSDataTransfer& sim_data);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void fixPossibleMissingIntersections(TPZGeoMesh* gmesh);

void ReadMeshesDFN(string& filenameBase, TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int& initVolForMergeMeshes);

void ReadMeshes(string& filenameFine, string& filenameCoarse,
                TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase1(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse, const string &filenamejson,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);
void ModifyPermeabilityForCase2(TPZGeoMesh* gmesh);
void ModifyBCsForCase2(TPZGeoMesh* gmesh);

void ReadMeshesFlemischCase3(string& filenameFine, string& filenameCoarse,
                             TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);
void ModifyBCsForCase3(TPZGeoMesh* gmesh);

void ReadMeshesFlemischCase4LF(string& filenameFine, string& filenameCoarse,
                               TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);
void ModifyBCsForCase4(TPZGeoMesh* gmesh);

void ReadMeshesFlemischCase4Debug(string& filenameFine, string& filenameCoarse,
                                  TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);

void ReadMeshesWell(string& filenameFine, string& filenameCoarse,
                    TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);
void ReadMeshesIP3D(string& filenameFine, string& filenameCoarse,
                    TPZGeoMesh*& gmesh, TPZGeoMesh*& gmeshcoarse);


enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EFracture, EIntersection, EIntersectionEnd, EPLossAtIntersect, EVolume2,
    /*15 HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
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
    // 1: Flemisch case 1
    // 2: Flemisch case 2
    // 3: Flemisch case 3
    // 4: Flemisch case 4
    // 5: Flemisch case 4 with less fractures (and no overlap)
    // 6: Flemisch case 4 with much less fractures (for debugging)
    // 7: Well mesh (Initially idealized just for generating a beautiful mesh)
    // 8: IP3D mesh (Initially idealized just for generating a beautiful mesh)
	// 9: 4 elements, 2 frac, w/ intersection
	// 10: Automated case 1
    int simcase = 10;
    string filenameCoarse, filenameFine,filenamejson;
    // @TODO define a root name and extend it with _coarse, _fine, and also .json
    switch (simcase) {
        case 0:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/twoElCoarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/twoElFine.msh";
            break;
        case 1:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case1_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case1_fine_NotSoFine.msh";
//            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case1_fine.msh";
            break;
        case 2:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case2_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case2_fine.msh";
            filenamejson = basemeshpath + "/verificationMHMNoHybrid/fl_case2.json";
            break;
        case 3:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case3_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case3_fine.msh";
            break;
        case 4:
            DebugStop(); // Need to generate mesh without overlap or need to treat overlap
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine.msh";
            break;
        case 5:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse_lf.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine_lf.msh";
            break;
        case 6:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/fl_case4_coarse_debug.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/fl_case4_fine_debug.msh";
            break;
        case 7:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/wellmesh_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/wellmesh_fine.msh";
            break;
        case 8:
            filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/ip3dmesh_coarse.msh";
            filenameFine = basemeshpath + "/verificationMHMNoHybrid/ip3dmesh_fine.msh";
            break;
		case 9:
			filenameCoarse = basemeshpath + "/verificationMHMNoHybrid/intersectCoarse.msh";
			filenameFine = basemeshpath + "/verificationMHMNoHybrid/intersectFine.msh";
			break;
		case 10:
			filenameCoarse = basemeshpath + "/dfnimrs/fl_case1";
			break;
		default:
            break;
    }
    RunProblem(filenameFine,filenameCoarse,filenamejson,simcase);
    return 0;
}
// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
void RunProblem(string& filenamefine, string& filenamecoarse, const string &filenamejson, const int simcase)
{
    auto start_time = std::chrono::steady_clock::now();
    
    bool isRefineMesh = false;
    const bool isPostProc = true;
	const bool isRunWithTranport = false;
	int initVolForMergeMeshes = EInitVolumeMatForMHM;
    
    // ----- Creating gmesh and data transfer -----
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
    if (simcase == 0 || simcase == 9)
        ReadMeshes(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 1)
        ReadMeshesFlemischCase1(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 2){
        ReadMeshesFlemischCase2(filenamefine,filenamecoarse,filenamejson, gmeshfine,gmeshcoarse);
    }
    else if (simcase == 3)
        ReadMeshesFlemischCase3(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 4)
        DebugStop();
    else if (simcase == 5)
        ReadMeshesFlemischCase4LF(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 6)
        ReadMeshesFlemischCase4Debug(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 7)
        ReadMeshesWell(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
    else if (simcase == 8)
        ReadMeshesIP3D(filenamefine,filenamecoarse,gmeshfine,gmeshcoarse);
	else if (simcase == 10)
		ReadMeshesDFN(filenamecoarse, gmeshfine, gmeshcoarse, initVolForMergeMeshes);
    else
        DebugStop();
        
    fixPossibleMissingIntersections(gmeshfine); // read about this in the function
    
    
    TMRSDataTransfer sim_data;
	// ----- Approximation space -----
	sim_data.mTNumerics.m_four_approx_spaces_Q = true;
	sim_data.mTNumerics.m_mhm_mixed_Q = true;
	sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;

    if (simcase == 1)
        FillDataTransferCase1(sim_data);
    else if (simcase == 2)
        FillDataTransferCase2(sim_data);
    else if (simcase == 3)
        FillDataTransferCase3(sim_data);
    else if (simcase == 0 || simcase == 9)
        FillDataTransfer(sim_data);
	else
		DebugStop();
	
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
    
    // ----- Setting gmesh -----
    // Code takes a fine and a coarse mesh to generate MHM data structure
    TMRSApproxSpaceGenerator aspace;
	aspace.InitMatIdForMergeMeshes() = initVolForMergeMeshes;
    sim_data.mTFracProperties.m_matid = EFracture;
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
	
	// ----- Setting the global data transfer -----
	aspace.SetDataTransfer(sim_data);

	if(isRefineMesh){
		cout << "\n---------------------- Uniformly refining geomesh ----------------------" << endl;
		gRefDBase.InitializeUniformRefPattern(ECube);
		gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
		gRefDBase.InitializeUniformRefPattern(EOned);
		for (auto gel : gmeshfine->ElementVec()){
			TPZManVector<TPZGeoEl*,10> children;
			gel->Divide(children);
		}
	}
	
	aspace.SetGeometry(gmeshfine,gmeshcoarse);
	{
		std::ofstream name("GeoMesh_Fine_AfterMergeMeshes.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
        std::ofstream name2("GeoMesh_MHM_domain.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name2,aspace.mSubdomainIndexGel);
	}
	

    
//    aspace.SetGeometry(gmeshfine);
    
    // ----- Creates the multiphysics compmesh -----
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true; // This makes solving very fast!
    int n_threads = 8;
    bool UsingPzSparse = true; // Necessary to use multithread for now...
    bool UsePardiso_Q = true;
    
    cout << "\n---------------------- Creating Analysis (Might optimize bandwidth) ----------------------" << endl;

    if((simcase == 1 ||simcase == 2 || simcase == 3) && isRunWithTranport){

        aspace.BuildAuxTransportCmesh();
        TPZCompMesh * transport_operator = aspace.GetTransportOperator();
        std::string name("mesh");
        aspace.PrintGeometry(name);
        std::ofstream name2("TransportOperator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(transport_operator, name2);

        TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
        sfi_analysis->SetDataTransferAndBuildAlgDatStruct(&sim_data);
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
		const int typeToPPinit = 1; // 0: both, 1: p/flux, 2: saturation
		const int typeToPPsteps = 2; // 0: both, 1: p/flux, 2: saturation
        for (int it = 1; it <= n_steps; it++) {
            sim_time = it*dt;
            sfi_analysis->m_transport_module->SetCurrentTime(dt);
            sfi_analysis->RunTimeStep();
            if(it==1){
                sfi_analysis->PostProcessTimeStep(typeToPPinit);
            }
            mixed_operator->LoadSolution(mixed_operator->Solution());
            if (sim_time >=  current_report_time) {
				cout << "\n---------------------- SFI Step " << it << " ----------------------" << endl;
                std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
                mixed_operator->UpdatePreviousState(-1.);
                sfi_analysis->PostProcessTimeStep(typeToPPsteps);
                pos++;
                current_report_time =reporting_times[pos];
                REAL InntMassFrac=sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(10);
                fileCilamce<<current_report_time/(86400*365)<<", "<<InntMassFrac<<std::endl;

                REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
                std::cout << "Mass report at time : " << sim_time << std::endl;
                std::cout << "Mass integral :  " << mass << std::endl;
            }
        }
    }
    else{
        // -------------- Running problem --------------
        TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
        mixAnalisys->SetDataTransfer(&sim_data);
        n_threads = 8;
        mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
        {
            std::ofstream out("mixedCMesh.txt");
            mixed_operator->Print(out);
        }
        mixAnalisys->Assemble();
        if(0)
        {
            int64_t nc = mixed_operator->NConnects();
            TPZFMatrix<STATE> &sol = mixed_operator->Solution();
            for (int64_t ic = 0; ic<nc; ic++) {
                TPZConnect &c = mixed_operator->ConnectVec()[ic];
                int64_t seqnum = c.SequenceNumber();
                if(seqnum < 0) continue;
                unsigned char lagrange = c.LagrangeMultiplier();
                STATE fill = 0.;
                if(lagrange == 0 || lagrange == 2 || lagrange == 6)
                {
                }
                else {
                    fill = 4.;
                }
                int ndof = c.NShape();
                for (int idf = 0; idf < ndof ; idf++) {
                    int64_t index = mixed_operator->Block().Index(seqnum, idf);
                    sol(index,0) = fill;
                }

            }
        }
        else
        {
            mixAnalisys->Solve();

        }
        TPZFastCondensedElement::fSkipLoadSolution = false;
        mixed_operator->LoadSolution(mixed_operator->Solution());
		
		// The problem is linear, and therefore, we can just call assemble and solve once.
		// However, the system is assembled in a "nonlinear" fashion, thus, the solution represents
		// -DeltaU. So, to obtain the correct solution, we multiply it by -1.
		// Note: This has to be done after LoadSolution()!
		mixed_operator->UpdatePreviousState(-1.);
        
        // ----- Post processing -----
        if (isPostProc) {
            mixAnalisys->fsoltransfer.TransferFromMultiphysics();
            int dimToPost = 3;
            mixAnalisys->PostProcessTimeStep(dimToPost);
            dimToPost = 2;
            mixAnalisys->PostProcessTimeStep(dimToPost);
        }
    }
    
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()/1000.;
    cout << "\n\n\t--------- Total time of simulation = " << total_time << " seconds -------\n" << endl;

    // ----- Cleaning up -----
    delete gmeshfine;
    delete gmeshcoarse;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransfer(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, unit_pressure = 1.0, zero_pressure = 0.;
    
    // Domain material
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainNameAndMatId["Volume2"] = EVolume2;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[0] = std::make_tuple(EInlet,D_Type,unit_pressure);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[1] = std::make_tuple(EOutlet,D_Type,zero_pressure);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[3] = std::make_tuple(EFaceBCPressure,D_Type,unit_pressure);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracDimNameAndMatId[2]["Fractures"] = EFracture;
    
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue.Resize(1);
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);

    
    // Simulation properties
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    
    //FracAndReservoirProperties
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
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransferCase1(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, inlet_pressure = 4.0, outlet_pressure = 1.0;
    
    // Domain material
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainNameAndMatId["Volume2"] = EVolume2;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[0] = std::make_tuple(EInlet,D_Type,inlet_pressure);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[1] = std::make_tuple(EOutlet,D_Type,outlet_pressure);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracDimNameAndMatId[2]["Fractures"] = EFracture;
    
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue.Resize(1);
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);

    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 102;
    sim_data.mTGeometry.mInterface_material_idFracSup = 101;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
//    sim_data.mTFracIntersectProperties.m_IntersectionId = EPLossAtIntersect;
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue.Resize(5);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[0] = std::make_tuple(EOutlet,N_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[1] = std::make_tuple(EInlet,D_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[2] = std::make_tuple(ENoflux,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[3] = std::make_tuple(EFracOutlet,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[4] = std::make_tuple(EFracInlet,5,1.);
    sim_data.mTGeometry.mInterface_material_idFracBound = EFracNoFlux;
    
    // Other properties
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0; //*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    sim_data.mTNumerics.m_n_steps = 100;
    sim_data.mTNumerics.m_dt      = 1.0e7;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
    //FracAndReservoirProperties
//    sim_data.mTFracProperties.m_Permeability = 1.0;
//    REAL kappa1=1.0;
//    REAL kappa2=1.0;
    sim_data.mTFracProperties.m_Permeability = 1.0e-3;
    REAL kappa1=1.0e-5;
	REAL kappa2=1.0e-6;
    int id1 = EVolume2;
    int id2 = EVolume;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa1;
    idPerm[id2]= kappa2;
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    REAL resPorosity1 = 0.2;
    REAL resPorosity2 = 0.25;
    REAL FracPorosity = 0.4;
    REAL fracFactor = 0.1;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(3);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[0] = std::make_tuple(EVolume,resPorosity1, 1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[1] = std::make_tuple(EVolume2, resPorosity2,1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[2] = std::make_tuple(EFracture, FracPorosity, fracFactor);
    
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames, scalnamesTransport;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    scalnamesTransport.Push("Sw");
    scalnamesTransport.Push("So");
    
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    sim_data.mTPostProcess.m_scalnamesTransport =scalnamesTransport;
    
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
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransferCase2(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, inlet_pressure = 1.0, outlet_flux = -1.0;
    
    // Domain material
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    sim_data.mTGeometry.mDomainNameAndMatId["Volume2"] = EVolume2;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[0] = std::make_tuple(EInlet,D_Type,inlet_pressure);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[1] = std::make_tuple(EOutlet,N_Type,outlet_flux);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracDimNameAndMatId[2]["Fractures"] = EFracture;
    sim_data.mTGeometry.mDomainFracDimNameAndMatId[1]["Intersection"] = EIntersection;
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue.Resize(1);
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);

 
    
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
    sim_data.mTFracIntersectProperties.m_IntersectionId = EIntersection;
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue.Resize(7);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[0] = std::make_tuple(EOutlet,D_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[1] = std::make_tuple(EInlet,N_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[2] = std::make_tuple(ENoflux,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[3] = std::make_tuple(EFaceBCPressure,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[4] = std::make_tuple(EFracNoFlux,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[5] = std::make_tuple(EFracInlet,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[6] = std::make_tuple(EFracOutlet,5,1.);
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
    
    
    // Other properties
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0; //*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    sim_data.mTNumerics.m_n_steps = 100;
    sim_data.mTNumerics.m_dt      = 0.25;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
    //FracAndReservoirProperties
    // There are two cases, one with very condutive fracture (1e4) and other with very NON conductive fracture (1e-4)
    sim_data.mTFracProperties.m_Permeability = 1.e4;
//	sim_data.mTFracProperties.m_Permeability = 1.e-4;
    REAL kappa1=1.0;
    REAL kappa2=0.1;
    int id1 = EVolume;
    int id2 = EVolume2;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa1;
    idPerm[id2]= kappa2;
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
  
    REAL resPorosity = 0.2;
    REAL fracPorosity = 0.1;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(4);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[0] = std::make_tuple(EVolume, resPorosity,1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[1] = std::make_tuple(EVolume2, resPorosity,1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[2] = std::make_tuple(EFracture, fracPorosity, 0.2);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[3] = std::make_tuple(EIntersection, fracPorosity, 0.2*0.2);
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames, scalnamesTransport;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    scalnamesTransport.Push("Sw");
    scalnamesTransport.Push("So");
    
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    sim_data.mTPostProcess.m_scalnamesTransport =scalnamesTransport;
    
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
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransferCase3(TMRSDataTransfer& sim_data){
    
    int D_Type = 0;
    int N_Type = 1;
    int Mixed_Type = 2;
    REAL zero_flux = 0.0, outlet_pressure = 0.0, inlet_flux = -1.0;
    
    // Domain material
    sim_data.mTGeometry.mDomainNameAndMatId["Volume"] = EVolume;
    
    // Domain boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[0] = std::make_tuple(EInlet,N_Type,inlet_flux);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[1] = std::make_tuple(EOutlet,D_Type,outlet_pressure);
    sim_data.mTBoundaryConditions.mBCMixedMatIdTypeValue[2] = std::make_tuple(ENoflux,N_Type,zero_flux);
            
    // Fracture material
    sim_data.mTGeometry.mDomainFracDimNameAndMatId[2]["Fractures"] = EFracture;
    sim_data.mTGeometry.mDomainFracDimNameAndMatId[1]["Intersection"] = EIntersection;
    
    // Fracture boundary conditions
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue.Resize(1);
    sim_data.mTBoundaryConditions.mBCMixedFracMatIdTypeValue[0] = std::make_tuple(EFracNoFlux,N_Type,zero_flux);

    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
//    sim_data.mTFracIntersectProperties.m_IntersectionId = EPLossAtIntersect;
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue.Resize(5);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[0] = std::make_tuple(EOutlet,N_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[1] = std::make_tuple(EInlet,D_Type,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[2] = std::make_tuple(ENoflux,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[3] = std::make_tuple(EFracOutlet,5,1.);
    sim_data.mTBoundaryConditions.mBCTransportMatIdTypeValue[4] = std::make_tuple(EFracInlet,5,1.);
    sim_data.mTGeometry.mInterface_material_idFracBound = EFracNoFlux;
    
    // Other properties
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0; //*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    sim_data.mTNumerics.m_n_steps = 100;
    sim_data.mTNumerics.m_dt      = 0.01;//*day;
    sim_data.mTNumerics.m_max_iter_sfi=1;
    sim_data.mTNumerics.m_max_iter_mixed=1;
    sim_data.mTNumerics.m_max_iter_transport=1;
    
    //FracAndReservoirProperties
    // There are two cases, one with very condutive fracture (1e4) and other with very NON conductive fracture (1e-4)
    sim_data.mTFracProperties.m_Permeability = 1.e4;
//    sim_data.mTFracProperties.m_Permeability = 1.e-4;
//    REAL kappa1=1.0;
//    REAL kappa2=1.0;
    REAL kappa1=1.0;
    int id1 = EVolume;
    std::map<int, REAL> idPerm;
    idPerm[id1]= kappa1;
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    REAL resPorosity = 0.2;
    REAL fracPorosity = 0.2;
    REAL fracLength = 0.1;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(4);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[0] = std::make_tuple(EVolume, resPorosity,1.0);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[2] = std::make_tuple(EFracture, fracPorosity, fracLength);
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[3] = std::make_tuple(EIntersection, fracPorosity, fracLength*fracLength);
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames, scalnamesTransport;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    scalnamesTransport.Push("Sw");
    scalnamesTransport.Push("So");
    
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    sim_data.mTPostProcess.m_scalnamesTransport =scalnamesTransport;
    
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

void ReadMeshesDFN(string& filenameBase, TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int& initVolForMergeMeshes) {
	using json = nlohmann::json;
	std::string filenamejson = filenameBase + ".json";
	
	std::ifstream filejson(filenamejson);
	json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file

	// ===================> Coarse mesh <=======================
	TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
	TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
	
	std::set<int> allmatids; // use to check for repeated matids
	
	if(input.find("Domains") == input.end()) DebugStop();
	for(auto& domain : input["Domains"]){
		if(domain.find("matid") == domain.end()) DebugStop();
		if(domain.find("name") == domain.end()) DebugStop();
		const int matid = domain["matid"];
		const string name = domain["name"];
		dim_name_and_physical_tagCoarse[3][name] = matid;
		const bool is_in = allmatids.find(matid) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(matid);
	}
	
	if(input.find("Boundary") == input.end()) DebugStop();
	for(auto& bc : input["Boundary"]){
		if(bc.find("name") == bc.end()) DebugStop();
		if(bc.find("type") == bc.end()) DebugStop();
		if(bc.find("value") == bc.end()) DebugStop();
		const string name = bc["name"];
		const int matid = bc["matid"];
		dim_name_and_physical_tagCoarse[2][name] = matid;
		dim_name_and_physical_tagFine[2][name] = matid;
		const bool is_in = allmatids.find(matid) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(matid);
	}

	string filenameCoarse = filenameBase + "_coarse.msh";
	gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
	
	// ===================> Fine mesh <=======================
	// Note that boundary elements have been added previously to the set dim_name_and_physical_tagFine
		
	// Check if we should use the same PZ MatId for all the fractures.
	// 2022 May: that is the only possible case for now
	int fracUniqueMatId = -10000, bcFracUniqueMatId = -10000;
	if(input.find("FractureUniqueMatIDForIMRS") != input.end()) {
		fracUniqueMatId = input["FractureUniqueMatIDForIMRS"];
		const bool is_in = allmatids.find(fracUniqueMatId) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(fracUniqueMatId);
	} else{
		DebugStop(); // For now all fractures have to have same matid
	}
	if(input.find("FractureBCUniqueMatIDForIMRS") != input.end()) {
		bcFracUniqueMatId = input["FractureBCUniqueMatIDForIMRS"];
		const bool is_in = allmatids.find(bcFracUniqueMatId) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(bcFracUniqueMatId);
	} else{
		DebugStop(); // For now all fractures have to have same matid
	}
	
	if(input.find("FractureInitMatId") == input.end()) DebugStop();
	const int fracInitMatId = input["FractureInitMatId"];

	
	// Loop over fractures in Json
	int fracCounter = 0;
	if(input.find("Fractures") == input.end()) DebugStop();
	for(auto& frac : input["Fractures"]){
		const int matid = fracInitMatId + fracCounter*2;
		string fracname = "Fracture" + to_string(fracCounter);
		string bcfracname = "BCfrac" + to_string(fracCounter);
		if(fracUniqueMatId != -10000){
			dim_name_and_physical_tagFine[2][fracname] = fracUniqueMatId;
			dim_name_and_physical_tagFine[1][bcfracname] = bcFracUniqueMatId;
		}
		else{
			const bool is_in = allmatids.find(matid) != allmatids.end();
			if(is_in) DebugStop();
			allmatids.insert(matid);
			dim_name_and_physical_tagFine[2][fracname] = matid;
		}
		fracCounter++;
	}

	
	// Adding volume physical tags
	const int maxMatId = *allmatids.rbegin();
	initVolForMergeMeshes = maxMatId + 1;
	if(input.find("NCoarseGroups") == input.end()) DebugStop();
	const int nCoarseGroups = input["NCoarseGroups"];
	std::string volbase = "c";
	for (int ivol = 0; ivol <= nCoarseGroups; ivol++) {
		std::string ivolstring = volbase + to_string(ivol);
		dim_name_and_physical_tagFine[3][ivolstring] = initVolForMergeMeshes + ivol;
	}
	
	string filenameFine = filenameBase + "_fine.msh";
	gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);

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

void ReadMeshesFlemischCase2(string& filenameFine, string& filenameCoarse, const string &filenamejson,
                             TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse) {
    
    // ===================> Coarse mesh <=======================
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    TPZGmshReader readcoarse;
    std::ifstream filecoarse(filenameCoarse);
    readcoarse.ReadPhysicalProperties4(filecoarse);
    // Domain
    // @TODO  I dont understand, why cant the coarse element index be an indicator of the subdomain. Isnt that the purpose of the coarse mesh?
    // the material of the coarse mesh should be the material that will be applied to the fine mesh
    using json = nlohmann::json;

    json input;
    std::ifstream filejson(filenamejson);
    // file >> input;
    input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file

    auto &domains = input["Domains"];
    if(domains.find("vol1") != domains.end())
    {
        std::cout << "found vol1 \n";
    }
    else
    {
        std::cout << "did not find vol1\n";
    }
    std::string volbase = "c";
    for (int ivol = 1; ivol <= 512; ivol++) { // How to set this maximum?
        std::string ivolstring = volbase + to_string(ivol);
        dim_name_and_physical_tagCoarse[3][ivolstring] = EVolume;
    }
    
    // Domain BC
    dim_name_and_physical_tagCoarse[2]["bc1"] = EFaceBCPressure;
        
    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);
    ModifyPermeabilityForCase2(gmeshcoarse);
        
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
    // @TODO Nao podemos colocar essas condicoes de contorno no arquivo gmsh?
    ModifyBCsForCase2(gmeshfine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ModifyPermeabilityForCase2(TPZGeoMesh* gmesh) {
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != EVolume) continue; // only matrix
        
        TPZVec<REAL> masscent(3,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        
        const REAL x = xcenter[0], y = xcenter[1], z = xcenter[2];
        if (x > 0.5 && y < 0.5)
            gel->SetMaterialId(EVolume2);
        
        if (x > 0.75 && y > 0.5 && y < 0.75 && z > 0.5)
            gel->SetMaterialId(EVolume2);
        
        if (x > 0.625 && x < 0.75 && y > 0.5 && y < 0.625 && z > 0.5 && z < 0.75)
            gel->SetMaterialId(EVolume2);        
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ModifyBCsForCase2(TPZGeoMesh* gmesh) {
        
    const REAL inletDomain = 0.25, outletDomain = 0.875;
    const REAL zerotol = ZeroTolerance();
    
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != EFaceBCPressure) continue; // 2d faces on boundary only
        
        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        const bool isXzero = fabs(xcenter[0]) < zerotol, isYzero = fabs(xcenter[1]) < zerotol, isZzero = fabs(xcenter[2]) < zerotol;
        const bool isXone = fabs(xcenter[0]-1.) < zerotol, isYone = fabs(xcenter[1]-1.) < zerotol, isZone = fabs(xcenter[2]-1.) < zerotol;

        gel->SetMaterialId(ENoflux); // Default is no flux
        // Setting inlet BCs
        if(isXzero && (xcenter[1] < inletDomain && xcenter[2] < inletDomain))
                gel->SetMaterialId(EOutlet);
        if(isYzero && (xcenter[0] < inletDomain && xcenter[2] < inletDomain))
                gel->SetMaterialId(EOutlet);
        if(isZzero && (xcenter[0] < inletDomain && xcenter[1] < inletDomain))
                gel->SetMaterialId(EOutlet);

        // Setting outlet BCs
        if(isXone && (xcenter[1] > outletDomain && xcenter[2] > outletDomain))
                gel->SetMaterialId(EInlet);
        if(isYone && (xcenter[0] > outletDomain && xcenter[2] > outletDomain))
                gel->SetMaterialId(EInlet);
        if(isZone && (xcenter[0] > outletDomain && xcenter[1] > outletDomain))
                gel->SetMaterialId(EInlet);
    }
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
    // bc4: No flux
    // bc5: imposed flux = -1
    // bc6: imposed pressure h=0
    dim_name_and_physical_tagCoarse[2]["bc4"] = ENoflux;
    dim_name_and_physical_tagCoarse[2]["bc5"] = EInlet;
    dim_name_and_physical_tagCoarse[2]["bc6"] = EOutlet;
        
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
    dim_name_and_physical_tagFine[2]["bc4"] = ENoflux;
    dim_name_and_physical_tagFine[2]["bc5"] = ENoflux;
    dim_name_and_physical_tagFine[2]["bc6"] = ENoflux;
     
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
    ModifyBCsForCase3(gmeshfine);
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ModifyBCsForCase3(TPZGeoMesh* gmesh) {
    
    const REAL zerotol = ZeroTolerance();
    const REAL onethird = 1./3., twothirds = 2./3.;
    
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != ENoflux) continue; // 2d faces on boundary only
        
        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        const REAL x = xcenter[0], y = xcenter[1], z = xcenter[2];
        const bool isYzero = fabs(y) < zerotol;
        const bool isYend = fabs(y-2.25) < zerotol;
        
        // Setting inlet BCs
        if(isYzero && (z > onethird && z < twothirds))
                gel->SetMaterialId(EInlet);
    
        // Setting outlet BCs
        if(isYend && (z > 0. && z < onethird))
                gel->SetMaterialId(EOutlet);
        if(isYend && (z > twothirds && z < 1.))
                gel->SetMaterialId(EOutlet);
    }
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
    dim_name_and_physical_tagCoarse[2]["bc1"] = ENoflux;
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
    dim_name_and_physical_tagFine[2]["bc1"] = ENoflux;
     
    // Fractures
    dim_name_and_physical_tagFine[2]["Fracture10"] = EFracture;
    
    // Fractures BCs
    for (int i = 0; i <= 36; i++) {
        string bcfrac = "BCfrac" + to_string(i);
        dim_name_and_physical_tagFine[1][bcfrac] = EFracNoFlux;
    }
        
    // Intersections
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
    dim_name_and_physical_tagFine[1]["fracIntersection_31_33"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_31_36"] = EIntersection;
    dim_name_and_physical_tagFine[1]["fracIntersection_34_35"] = EIntersection;

    gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    ModifyBCsForCase4(gmeshfine);
}

void ModifyBCsForCase4(TPZGeoMesh* gmesh) {
    
    const REAL zerotol = ZeroTolerance();
    
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != ENoflux) continue; // 2d faces on boundary only
        
        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        const REAL x = xcenter[0], y = xcenter[1], z = xcenter[2];
        const bool isYend = fabs(y-1500.) < zerotol;
        const bool isXinit = fabs(x+500.) < zerotol;
        const bool isXend = fabs(x-350.) < zerotol;
        
        // Default is no flux already set previously
        
        // Setting inlet BCs
        if(isYend && (z > 300. && x < -200.))
                gel->SetMaterialId(EInlet);
        if(isXinit && (z > 300. && y > 1200.))
                gel->SetMaterialId(EInlet);

        // Setting outlet BCs
        if(isXinit && (y < 400. && z < 100.))
                gel->SetMaterialId(EOutlet);
        if(isXend && (y < 400. && z < 100.))
                gel->SetMaterialId(EOutlet);
    }
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
    cout << "\n---------------------- Searching for problematic intersection elements ----------------------" << endl;
    
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
            TPZStack<TPZGeoEl*> bcEls;
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
                if (neigel->MaterialId() == EFracNoFlux) {
                    bcEls.push_back(neigel);
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
            if ((interEls.size() || nFracElsForSide > 1) && bcEls.size()) {
                cout << "PLEASE CHECK CAREFULLY! An element may have a boundary even if it intersects\n";
                DebugStop();
                if(bcEls.size() > 1)
                    DebugStop(); // This is rather odd. Please check why there are two bcs at the same boundary
                for (auto bcgeoel : bcEls) {
                    bcgeoel->RemoveConnectivities();
                    delete bcgeoel;
                }
            }
            
            
        }
    }
    if (nInterCreated) {
        cout << "\n\t\t==================== MESSAGE =================" << endl;
        cout << "\n\t- Imported mesh has " << nInterCreated << " missing intersection elements" << endl;
        cout << "\t- These were created here and this shouldn't be a problem...\n" << endl;
    }
    else {
        cout << "\n====> Ok! No problematic intersections found" << endl;
    }
        
    gmesh->BuildConnectivity();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
