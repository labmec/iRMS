// C++ includes
// nothing

// PZ includes
#include "TPZGenGrid3D.h"
#include "TMRSApproxSpaceGenerator.h"
#include "imrs_config.h"
#include "pzlog.h"
#include <tpzgeoelrefpattern.h>
#include "json.hpp"
#include <filesystem>
#include <libInterpolate/Interpolate.hpp>
#include <libInterpolate/AnyInterpolator.hpp>
#include "TPZSimpleTimer.h"
#include "pzintel.h"
#include "pzsmanal.h"
// include dfn filereader
#include "filereader.h"

#include "DFNMesh.h"

// ----- Namespaces -----
using namespace std;
namespace fs = std::filesystem;

// ----- Global vars -----
const int glob_n_threads = 8;

// ----- End of namespaces -----

// ----- Functions -----
void RunProblem(string& filenameBase, const int simucase);
void ReadMeshesDFN(string& filenameBase, TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int& initVolForMergeMeshes, bool& isMHM, bool& needsMergeMeshes);
void CreateIntersectionElementForEachFrac(TPZGeoMesh* gmeshfine,
										  std::map<int,std::pair<int,int>>& matidtoFractures,
										  const int fracInitMatId, const int fracinc, const int FractureHybridPressureMatId);
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void fixPossibleMissingIntersections(TMRSDataTransfer& sim_data, TPZGeoMesh* gmesh);

void FillDataTransferDFN(string& filenameBase, string& outputFolder, TMRSDataTransfer& sim_data);

void FillPCteSol(TPZCompMesh* mpcmesh, const REAL pcte);

// Quick fix functions. Should be deleted in the future (May 2022)
void ModifyPermeabilityForCase2(TPZGeoMesh* gmesh);
void ModifyBCsForCase2(TPZGeoMesh* gmesh);
void ModifyBCsForCase3(TPZGeoMesh* gmesh);
void ModifyBCsForCase4(TPZGeoMesh* gmesh);
void ModifyBCsFor2ParallelFractures(TPZGeoMesh* gmesh);
bool fileExists(const fs::path& p, fs::file_status s = fs::file_status{});
void CreateOutputFolders(std::string& outputFolder);
void CopyInputFilesToOutputFolderAndFixFilename(std::string& filenameBase, std::string& outputFolder);

TPZGeoMesh * Transform2dMeshToUnisim3D(TPZGeoMesh* gmesh2d, int nLayers);
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers);
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
void FilterZeroNeumann(std::string& outputFolder, TMRSDataTransfer& sim_data, TPZAutoPointer<TPZStructMatrix> strmat, TPZCompMesh* cmesh);
void ComputeDiagnostics(std::string& outputFolder, TMRSDataTransfer& sim_data, std::set<int>& bcflux, TPZMultiphysicsCompMesh* mixed_operator);
void VerifyIfNeumannIsExactlyZero(const int matidNeumann, TPZMultiphysicsCompMesh* mixed_operator);

std::map<int,TPZVec<STATE>> computeIntegralOfNormalFlux(const std::set<int> &bcMatId, TPZMultiphysicsCompMesh *cmesh);

// TODO: Delete this enum? (May 2022)
enum EMatid {/*0*/ENone, EVolume, EInlet, EOutlet, ENoflux,
    /*5*/EFaceBCPressure, EFracInlet, EFracOutlet, EFracNoFlux, EFracPressure,
    /*10*/EFracture, EIntersection, EIntersectionEnd, EPLossAtIntersect, EVolume2,
    /*15 HAS TO BE LAST*/EInitVolumeMatForMHM /*Not sure if we will keep this structure */
};

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("imrs");
static TPZLogger fracIntersectLogger("imrs_fracIntersect");
#endif

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]){
    string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
    string logpath = basemeshpath + "/../DFNIMRS/log4cxx.cfg";
    TPZLogger::InitializePZLOG(logpath);
    if (mainlogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "\nLogger for DFNIMRS problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif
    
    string filenameBase;
    int simcase = 5;
    if (argc > 1) {
        std::cout << "\n===========> Running with provided argv path!" << std::endl;
        filenameBase = basemeshpath + argv[1];
    }
    else{
        // 0: two elements, 1 frac, no intersection
        // 1: Flemisch case 1
        // 2: Flemisch case 2
        // 3: Flemisch case 3
        // 4: Flemisch case 4
        // 5: 4 elements, 2 frac, w/ intersection
        // 6: Two parallel square fractures very close, but no overlap
        // 7: Two parallel square fractures very close, WITH overlap
        // 8: Flemisch case 3 with snapping of the bottom fracture to the middle fracture
        // 9: Study of snapping tolerances on case 3
        // 10: Case 3 snapping of middle fractures. NO snap of fractures to domain boundary
        // 11: Case 3 snapping of middle fractures. With snap of fractures to domain boundary
        // 12,13,14,15,16,17: Modified Case 3 where all fracs touch boundary. Snapping is 0.0001, 0.04, 0.05, 0.1, 0.01, 0.03 respectively
        // 18: joker path, edit at will
        // 19: Case4 mesh 2018
        // 20: Unisim
        // 21: Flemisch case 4 with constant pressure
        switch (simcase) {
            case 0:
                filenameBase = basemeshpath + "/dfnimrs/twoelCoarse";
                break;
            case 1:
                //            filenameBase = basemeshpath + "/dfnimrs/fl_case1";
                filenameBase = basemeshpath + "/dfnimrs/fl_case1/";
                break;
            case 2:
                filenameBase = basemeshpath + "/dfnimrs/fl_case2/";
                break;
            case 3:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3/";
                break;
            case 4:
                filenameBase = basemeshpath + "/dfnimrs/fl_case4_meshes/fl_case4_test2";
                break;
            case 5:
                filenameBase = basemeshpath + "/dfnimrs/intersect/";
                break;
            case 6:
                filenameBase = basemeshpath + "/dfnimrs/2parallel/";
                break;
            case 7:
                filenameBase = basemeshpath + "/dfnimrs/2paralleloverlap/";
                break;
            case 8:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_snap/";
                break;
            case 9:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/6x6x13/TestFunciona";
                break;
            case 10:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/";
                break;
            case 11:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/TestNoFunciona/";
                break;
            case 12:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/touchBound_s_0001/";
                break;
            case 13:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/touchBound_s_04/";
                break;
            case 14:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/touchBound_s_05/";
                break;
            case 15:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/touchBound_s_1/";
                break;
            case 16:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/touchBound_s_01/";
                break;
            case 17:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/touchBound_s_03/";
                break;
            case 18:
                filenameBase = basemeshpath + "/dfnimrs/fl_case3_meshes/testingNoSnapBound/";
                break;
            case 19:
                filenameBase = basemeshpath + "/dfnimrs/fl_case4_meshes/fl_case4_2018/";
                break;
            case 20:
                filenameBase = basemeshpath + "/dfnimrs/unisim_meshes/";
                break;
            case 21:
                filenameBase = basemeshpath + "/dfnimrs/fl_case4_meshes/fl_case4_lf";
                break;
            case 22:
                filenameBase = basemeshpath + "/dfnimrs/boxPerpFlux";
                break;
            case 23:
                filenameBase = basemeshpath + "/dfnimrs/boxPerpFlux6el";
                break;
            default:
                break;
        }
    }
    RunProblem(filenameBase,simcase);
    return 0;
}
// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

struct FractureQuantities {
    TPZMultiphysicsCompMesh *fCMesh = 0;
    std::set<int> fAllFracturesMatId;
    int fracMatid = 0;
    int fPressureIntersectMatid = 0;
    int fGlueMatid = 0;
    REAL fIntegrateFluxNorm = 0.;
    REAL fAveragePressure = 0.;
    REAL fdivintegral = 0.;
    REAL fdivpositive_integral = 0.;
    REAL fdivnegative_integral = 0.;
    REAL fBoundarySnapSize = 0.;
    REAL fFaceSnapSize = 0;
    REAL fFracBCFluxIntegral = 0.;
    REAL fFracIntersectFluxIntegral = 0.;
    // breakdown of fluid transmission by intersecting fracture set
    std::map<std::set<int>, std::pair<REAL,REAL> > fBoundaryTransmission;
    
    FractureQuantities(TPZMultiphysicsCompMesh *cmesh, const std::set<int> &allfracs, int pressureintersect, int matglue) : fCMesh(cmesh), fAllFracturesMatId(allfracs), fPressureIntersectMatid(pressureintersect),
    fGlueMatid(matglue) {
        
    }
    void ComputeFluxQuantities(int matid)
    {
        fracMatid = matid;
        IntegrateFractureFlux();
        AveragePressure();
        IntegrateDivFracture();
        ComputeTransmissionBreakdown();
    }
    
    void IntegrateFractureFlux()
    {
        std::set<int> matids = {fracMatid};
        auto flux = fCMesh->Integrate("FluxNorm", matids);
        fIntegrateFluxNorm = flux[0];
        TPZGeoMesh *gmesh = fCMesh->Reference();
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel ; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel) continue;
            if(gel->MaterialId() == fracMatid+1) {
                TPZGeoElSide gelside(gel);
                if(gelside.HasNeighbour(fPressureIntersectMatid)) {
                    fBoundarySnapSize += gel->Volume();
                }
            }
            if(gel->MaterialId() == fracMatid) {
                TPZGeoElSide gelside(gel);
                if(gelside.HasNeighbour(fGlueMatid)) {
                    fFaceSnapSize += gel->Volume();
                }
            }
        }
    }
    void AveragePressure()
    {
        std::set<int> matids = {fracMatid};
        auto pressure = fCMesh->Integrate("Pressure", matids);
        fAveragePressure = pressure[0];
        REAL area = 0.;
        TPZGeoMesh *gmesh  = fCMesh->Reference();
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel && gel->MaterialId() == fracMatid){
                area += gel->Volume();
            }
        }
        fAveragePressure /= area;
    }
    
    void IntegrateDivFracture()
    {
        std::set<int> matids = {fracMatid};
        auto divpos = fCMesh->Integrate("div_positive", matids);
        fdivpositive_integral = divpos[0];
        auto divneg = fCMesh->Integrate("div_negative", matids);
        fdivnegative_integral = divneg[0];
        auto divint = fCMesh->Integrate("div_q", matids);
        fdivintegral = divint[0];
        std::set<int> matidfracbound = {fracMatid+1};
        auto bcflux = fCMesh->Integrate("BCNormalFlux", matidfracbound);
        std::set<int> matidfracintersect = {fracMatid+2};
        auto intersectflux = fCMesh->Integrate("BCNormalFlux", matidfracintersect);
        fFracBCFluxIntegral = bcflux[0];
        fFracIntersectFluxIntegral = intersectflux[0];
    }
    
    void Print(std::ostream &out) {
        out << "\n------------ For fracture " << fracMatid << " ------------" << std::endl;
        out << "Average pressure " << fAveragePressure << std::endl;
        out << "Integrate fluxnorm " << fIntegrateFluxNorm << std::endl;
        out << "Integrated divergence " << fdivintegral << std::endl;
        out << "Integrated positive divergence " << fdivpositive_integral << std::endl;
        out << "Integrated negative divergence " << fdivnegative_integral << std::endl;
        out << "Integrated boundary flux " << fFracBCFluxIntegral << std::endl;
        out << "Integrated intersection flux " << fFracIntersectFluxIntegral << std::endl;
        out << "Size of snapped boundary " << fBoundarySnapSize << std::endl;
        out << "Area of overlapping faces " << fFaceSnapSize << std::endl;
        REAL sum = -fdivintegral;
        for(auto &it : fBoundaryTransmission) {
            out << "Flux transmitted to fractures ";
            for(auto id : it.first) out << id << " ";
            out << ": " << it.second.first << " and " << it.second.second << std::endl;
            sum += it.second.first + it.second.second;
        }
        out << "Conservation inside fracture " << sum << std::endl;;
    }
    
    void ComputeTransmissionBreakdown()
    {
        TPZCompMesh *fluxmesh = fCMesh->MeshVector()[0];
        TPZGeoMesh *gmesh = fluxmesh->Reference();
        gmesh->ResetReference();
        fluxmesh->LoadReferences();
        int fracintersect = fracMatid+2;
        // for each geometric element of fracintersect
        //      find all neighbours of matid in fractures (form a set)
        //      integrate the flux value
        //      add to the datastructure
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() != fracintersect) continue;
            TPZGeoElSide gelside(gel);
            std::set<int> fracmats;
            for(auto neigh = gelside.Neighbour(); neigh != gelside; neigh++)
            {
                int matid = neigh.Element()->MaterialId();
                if(fAllFracturesMatId.find(matid) != fAllFracturesMatId.end()) fracmats.insert(matid);
            }
            // remove the current fracture matid
            if(fracmats.find(fracMatid) == fracmats.end()) DebugStop();
            fracmats.erase(fracMatid);
            std::pair<REAL,REAL> previntegral = {0.,0.};
            if(fBoundaryTransmission.find(fracmats) == fBoundaryTransmission.end())
            {
                fBoundaryTransmission[fracmats] = previntegral;
            } else {
                previntegral = fBoundaryTransmission[fracmats];
            }
            TPZCompEl *cel = gel->Reference();
            TPZVec<REAL> fluxvec(1,0.);
            fluxmesh->SetDimModel(1);
            cel->Integrate(0, fluxvec);
            REAL flux = fluxvec[0];
            if(flux < 0.) previntegral.first += flux;
            else previntegral.second += flux;
            fBoundaryTransmission[fracmats] = previntegral;
        }
    }

};

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void RunProblem(string& filenameBase, const int simcase)
{
    auto start_time = std::chrono::steady_clock::now();
    
	// ----- Simulation and printing parameters -----
    const bool isRefineMesh = false;
    const bool isPostProc = true;
    const bool isPostProcessFracDiagnostics = true;
    const bool isFilterZeroNeumann = true;
	bool isMHM = true; // may be set in json
    bool needsMergeMeshes = true; // may be set in json
    
    // ----- output folder stuff -----
    if(filenameBase.back() != '/') filenameBase = filenameBase + "/";
    std::string outputFolder = "Output/"+filenameBase.substr(filenameBase.find("dfnimrs/") + 8);
    CreateOutputFolders(outputFolder);
    CopyInputFilesToOutputFolderAndFixFilename(filenameBase,outputFolder);
    outputFolder = outputFolder.substr(0,outputFolder.find_last_of("/"));
    outputFolder = outputFolder + "/";

    
    // ----- Creating gmesh and data transfer -----
	int initVolForMergeMeshes = -1000000;
    TPZGeoMesh *gmeshfine = nullptr, *gmeshcoarse = nullptr;
	ReadMeshesDFN(filenameBase, gmeshfine, gmeshcoarse, initVolForMergeMeshes,isMHM,needsMergeMeshes);
    
    if(simcase==20){
        gmeshfine = Transform2dMeshToUnisim3D(gmeshfine, 5);
    }
    // ----- Printing gmesh -----
#ifdef PZDEBUG
    if (1) {
        if(gmeshfine){
            gmeshfine->SetDimension(3);
            std::ofstream name(outputFolder + "GeoMesh_Fine_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
            std::ofstream name2(outputFolder + "GeoMesh_Fine_Initial.txt");
            gmeshfine->Print(name2);
        }
        if(gmeshcoarse){
            std::ofstream name(outputFolder + "GeoMesh_Coarse_Initial.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmeshcoarse, name);
        }
    }
#endif
	
	// ----- Approximation space -----
    TMRSDataTransfer sim_data;
	sim_data.mTNumerics.m_four_approx_spaces_Q = true;
	sim_data.mTNumerics.m_mhm_mixed_Q = isMHM;
	sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;

	// ----- Filling data transfer from json input file -----
	FillDataTransferDFN(filenameBase, outputFolder, sim_data);
	
	// ----- Creating intersection GeoElBCs -----
	fixPossibleMissingIntersections(sim_data,gmeshfine); // read about this in the function itself
    
	// ----- Creating approximation space creator and manager -----
	// Code takes a fine and a coarse mesh to generate MHM data structure
    TMRSApproxSpaceGenerator aspace;
	aspace.InitMatIdForMergeMeshes() = initVolForMergeMeshes;
	
	// ----- Setting the global data transfer in approximation space -----
	aspace.SetDataTransfer(sim_data);

	// ----- Refining mesh (always done in fine mesh) -----
	if(isRefineMesh){
		cout << "\n---------------------- Uniformly refining geomesh ----------------------" << endl;
		gRefDBase.InitializeUniformRefPattern(ECube);
		gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
		gRefDBase.InitializeUniformRefPattern(EOned);
		for (auto gel : gmeshfine->ElementVec()){
			if(!gel || gel->HasSubElement()) continue;
			TPZManVector<TPZGeoEl*,10> children;
			gel->Divide(children);
		}
	}
	
	// ----- Setting gmesh -----
	// Code takes a fine and a coarse mesh to generate MHM data structure
	aspace.SetGeometry(gmeshfine,gmeshcoarse);

	{
        {
            int dim = gmeshfine->Dimension();
            int64_t nel = gmeshfine->NElements();
            int error = 0;
            for (int64_t el = 0; el<nel; el++) {
                TPZGeoEl *gel = gmeshfine->Element(el);
                if(gel->Dimension() != dim) continue;
                int firstside = gel->FirstSide(2);
                for(int is = firstside; is< gel->NSides()-1; is++)
                {
                    TPZGeoElSide gelside(gel,is);
                    auto neigh = gelside.Neighbour();
                    if(neigh == gelside) {
                        std::cout << gelside << "Has no neighbour\n";
                        error = 1;
                        TPZGeoElBC gelbc(gelside,-1000);
                    }
                }
            }
            if(error == 1){
                // NOTE: This may happen if setting the boundary conditions wrongly. I suggest checking them
                std::cout << "3D element without neighbour \n";
                std::ofstream name(outputFolder + "GeoMesh_Fine_AfterMergeMeshes.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
                DebugStop();
            }
        }
//		std::ofstream name(outputFolder + "GeoMesh_Fine_AfterMergeMeshes.vtk");
//		TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name);
        std::ofstream name2(outputFolder + "GeoMesh_MHM_domain.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name2,aspace.mSubdomainIndexGel);
	}
    
	// ----- Changing BCs for some testing cases -----
    if(simcase == 6 || simcase == 7){
        //linear pressure...
        ModifyBCsFor2ParallelFractures(gmeshfine);
//        std::ofstream name3(outputFolder + "ModBCs.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmeshfine, name3);
    }
    
    // ----- Creates the multiphysics compmesh -----
	const int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
            
    // ----- Analysis parameters -----
    bool must_opt_band_width_Q = true; // This makes solving very fast!
    bool UsingPzSparse = true; // Necessary to use multithread for now...
    bool UsePardiso_Q = true; // lighting fast!
    
    cout << "\n---------------------- Creating Analysis (Might optimize bandwidth) ----------------------" << endl;
    if(sim_data.mTNumerics.m_run_with_transport){

		// Create transport mesh. TODO: Create transport data structure without the need for a mesh
        aspace.BuildAuxTransportCmesh();
        TPZCompMesh * transport_operator = aspace.GetTransportOperator();
#ifdef PZDEBUG
        std::string name(outputFolder + "mesh");
        aspace.PrintGeometry(name);
        std::ofstream name2(outputFolder + "TransportOperator.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(transport_operator, name2);
#endif

		// Creating coupled pressure/flow and transport analysis
        TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
        sfi_analysis->SetDataTransferAndBuildAlgDatStruct(&sim_data);
        sfi_analysis->Configure(glob_n_threads, UsePardiso_Q, UsingPzSparse);
        if(isFilterZeroNeumann) FilterZeroNeumann(outputFolder,sim_data,sfi_analysis->m_mixed_module->StructMatrix(),mixed_operator);
        const int n_steps = sim_data.mTNumerics.m_n_steps;
        const REAL dt = sim_data.mTNumerics.m_dt;

		// Times to report solution
        TPZStack<REAL,100> reporting_times;
        reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
        REAL sim_time = 0.0;
        int pos =0;
        REAL current_report_time = reporting_times[pos];
        int npos = reporting_times.size();

		// Initializing tranport solution
        sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
        REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        std::cout << "\nMass report at initial time : " << 0.0 << std::endl;
        std::cout << "Initial mass:  " << initial_mass << std::endl;
//        std::ofstream fileCilamce535(outputFolder + "IntegratedSatFrac365.txt");
//        std::ofstream fileCilamce515(outputFolder + "IntegratedSatFrac515.txt");
//        std::ofstream fileCilamce530(outputFolder + "IntegratedSatFrac530.txt");
        std::ofstream fileFracSat(outputFolder + "FracturesIntegratedSat.txt");
        fileFracSat << "t" << " ";
        for (int ifrac = 0 ; ifrac < sim_data.mTFracProperties.m_fracprops.size() ; ifrac++) fileFracSat << "frac" << ifrac << " ";
        fileFracSat << std::endl;

        TPZFastCondensedElement::fSkipLoadSolution = false;
		const int typeToPPinit = 1; // 0: both, 1: p/flux, 2: saturation
		const int typeToPPsteps = 2; // 0: both, 1: p/flux, 2: saturation
      
		// Looping over time steps
        for (int it = 1; it <= n_steps; it++) {
            sim_time = it*dt;
            sfi_analysis->m_transport_module->SetCurrentTime(dt);
            sfi_analysis->RunTimeStep();
            if(it == 1){
                sfi_analysis->PostProcessTimeStep(typeToPPinit);
                if(isPostProcessFracDiagnostics){
                    std::set<int> bcflux = {2,3,4}; // computes integral of quantity over these matids
                    ComputeDiagnostics(outputFolder, sim_data, bcflux, mixed_operator);
                }
                if(isFilterZeroNeumann) VerifyIfNeumannIsExactlyZero(4,mixed_operator);
            }
            mixed_operator->LoadSolution(mixed_operator->Solution());
			
			// Only post process based on reporting times
            if (sim_time >=  current_report_time) {
				cout << "\n---------------------- SFI Step " << it << " ----------------------" << endl;
                std::cout << "Simulation time:  " << sim_time << std::endl;
                mixed_operator->UpdatePreviousState(-1.);
                sfi_analysis->PostProcessTimeStep(typeToPPsteps);
                pos++;
                current_report_time = reporting_times[pos];
                REAL InntMassFrac = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(10);
                if(sim_data.mTNumerics.m_run_with_transport){
                    fileFracSat << current_report_time << " ";
                    for (auto& fprop : sim_data.mTFracProperties.m_fracprops) {
                        const int matid = fprop.first;
                        const REAL intMassThisFrac = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(matid);
                        fileFracSat << intMassThisFrac << " ";
                    }
                    fileFracSat << std::endl;
                }
                REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
                std::cout << "Mass report at time : " << sim_time << std::endl;
                std::cout << "Mass integral :  " << mass << std::endl;
            }
            sfi_analysis->m_transport_module->fAlgebraicTransport.VerifyConservation(it);
        }
    }
    else{
        // -------------- Running problem --------------
        TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
        mixAnalisys->SetDataTransfer(&sim_data);
        UsePardiso_Q = true;
        mixAnalisys->Configure(glob_n_threads, UsePardiso_Q, UsingPzSparse);
        if(isFilterZeroNeumann) FilterZeroNeumann(outputFolder,sim_data,mixAnalisys->StructMatrix(),mixed_operator);
        
        {
            std::ofstream out(outputFolder + "mixedCMesh.txt");
            mixed_operator->Print(out);
        }
        mixAnalisys->Assemble();
		
		// Testing if constant pressure leads to zero residual in cte pressure problem
		const bool testCtePressure = false;
		if(testCtePressure){
            const int neq = mixAnalisys->Mesh()->NEquations();
            TPZFMatrix<STATE> res(neq,1,0.);
            FillPCteSol(mixed_operator,1.);
            mixAnalisys->fsoltransfer.TransferFromMultiphysics();
//            mixAnalisys->PostProcessTimeStep(2);
//            mixAnalisys->PostProcessTimeStep(3);
            TPZMatrix<STATE>* mat = mixAnalisys->MatrixSolver<STATE>().Matrix().operator->();
            mat->Multiply(mixed_operator->Solution(), res);
            res = res + mixAnalisys->Rhs();
            std::ofstream out(outputFolder + "problematicElsGlob.txt");
            mixAnalisys->PrintVectorByElement(out, res, 1.e-6);
            if(sim_data.mTNumerics.m_mhm_mixed_Q){
                for(auto cel : mixed_operator->ElementVec()){
                    TPZSubCompMesh* subcmesh = dynamic_cast<TPZSubCompMesh*>(cel);
                    if(subcmesh){
                        FillPCteSol(subcmesh,1.);
//                        TPZSubMeshAnalysis* sanal = dynamic_cast<TPZSubMeshAnalysis*>(subcmesh->Analysis().operator->());
//                        TPZMatrix<STATE>* mat = sanal->Matrix().operator->();
//                        TPZMatRed<STATE, TPZFMatrix<STATE>>* matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE>>*>(mat);
//
//                        matred->ResetReduced();
//                        TPZFMatrix<STATE> res(subcmesh->TPZCompMesh::NEquations(),1,0.);
//                        TPZFMatrix<STATE> k = *mat;
//                        k.SetIsDecomposed(0);
//                        for (int i = 0; i < matred->Dim0(); i++) {
//                            for (int j = matred->Dim0(); j < matred->Dim0() + matred->Dim1(); j++) {
//                                k(i,j) = matred->fK01orig(i,j-matred->Dim0());
//                            }
//                        }
//
//                        TPZFMatrix<STATE> subcmeshsolint = subcmesh->TPZCompMesh::Solution();
//                        TPZEquationFilter& subcmeshfilter = sanal->StructMatrix()->EquationFilter();
//                        if(subcmeshfilter.IsActive()){
//                            TPZFMatrix<STATE> subcmeshsolintfilter(subcmeshfilter.NActiveEquations(),1,0.);
//                            subcmeshsolint.Resize(subcmesh->TPZCompMesh::NEquations(), 1);
//                            subcmeshfilter.Gather(subcmeshsolint, subcmeshsolintfilter);
//                            subcmeshsolint = subcmeshsolintfilter;
//                        }
//                        subcmeshsolint.Resize(k.Rows(), 1);
//                        subcmeshsolint.Print("solint=",std::cout,EMathematicaInput);
//
//                        k.Multiply(subcmeshsolint, res);
//                        res.Print("res=",std::cout,EMathematicaInput);
//                        k.Print("k=",std::cout,EMathematicaInput);
//                        TPZFMatrix<STATE> rhsint(matred->Dim0(),1,0.);
//                        for(int i = 0 ; i < matred->Dim0() ; i++) rhsint(i,0) = res(i,0);
//                        rhsint.Print("rhsint=",std::cout,EMathematicaInput);
//                        rhsint = rhsint + matred->fF0orig;
//                        rhsint.Resize(k.Rows(),1);
//                        matred->fF0orig.Print("f0=",std::cout,EMathematicaInput);
//                        std::string outname = outputFolder + "problematicEls_subcmesh_" + to_string(subcmesh->Index()) + ".txt";
//                        std::ofstream outs(outname);
//                        if(subcmeshfilter.IsActive()){
//                            TPZFMatrix<STATE> rhsexpand(subcmeshfilter.NEqExpand(),1,0.);
//                            subcmeshfilter.Scatter(rhsint, rhsexpand);
//                            rhsint = rhsexpand;
//                        }
//                        subcmesh->Analysis()->PrintVectorByElement(outs, rhsint);
                    }
                }
            }
		}
		
		// Solving problem
		mixAnalisys->Solve();
		mixAnalisys->VerifyElementFluxes();
		
        TPZFastCondensedElement::fSkipLoadSolution = false;

		
		// The problem is linear, and therefore, we can just call assemble and solve once.
		// However, the system is assembled in a "nonlinear" fashion, thus, the solution represents
		// -DeltaU. So, to obtain the correct solution, we multiply it by -1.
		// Note: This has to be done after LoadSolution()!

        REAL res_norm = Norm(mixAnalisys->Rhs());
        REAL normsol = Norm(mixed_operator->Solution());
//        mixed_operator->UpdatePreviousState(-1.);
        mixAnalisys->LoadSolution();
        mixAnalisys->fsoltransfer.TransferFromMultiphysics();
        mixAnalisys->LoadSolution();
        mixAnalisys-> Assemble();
//        mixAnalisys->fsoltransfer.TransferToMultiphysics();
      
        REAL res_norm1 = Norm(mixAnalisys->Rhs());
        REAL normsol1 = Norm(mixed_operator->Solution());
        
        mixAnalisys->fsoltransfer.TransferFromMultiphysics();
//        if(isFilterZeroNeumann) VerifyIfNeumannIsExactlyZero(4,mixed_operator);
        
        {
            // print the flux mesh with the solutions "loaded"
            TPZCompMesh *fluxmesh = mixed_operator->MeshVector()[0];
            std::ofstream flux(outputFolder + "fluxmesh.txt");
            fluxmesh->Print(flux);
        }

//        TPZCompMesh *pressure = mixed_operator->MeshVector()[1];
//        pressure->Solution().Print("pressure multipliers");
        // Computes the integral of the normal flux on the boundaries.
        if(isPostProcessFracDiagnostics){
            std::set<int> bcflux = {2,3,4}; // computes integral of quantity over these matids
            ComputeDiagnostics(outputFolder, sim_data, bcflux, mixed_operator);
        }
    
        
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

void FillDataTransferDFN(string& filenameBase, string& outputFolder, TMRSDataTransfer& sim_data) {
	
	using json = nlohmann::json;
	std::string filenamejson = filenameBase + ".json";
	
	std::ifstream filejson(filenamejson);
	json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file

	
	// ------------------------ Getting number of domains and fractures ------------------------
	if(input.find("Domains") == input.end()) DebugStop();
	const int ndom = input["Domains"].size();
	if(input.find("Fractures") == input.end());
	const int nfrac = input["Fractures"].size();
	sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(ndom+nfrac+1);
	int countPhi = 0;
	
	// ------------------------ Reading 3D Domain matids ------------------------
	std::map<int, REAL> idPerm;
	for(auto& domain : input["Domains"]){
		if(domain.find("matid") == domain.end()) DebugStop();
		if(domain.find("name") == domain.end()) DebugStop();
		if(domain.find("K") == domain.end()) DebugStop();
		if(domain.find("phi") == domain.end()) DebugStop();
		const int matid = domain["matid"];
		const string name = domain["name"];
		const REAL permeability = domain["K"];
		const REAL phi = domain["phi"];
		sim_data.mTGeometry.mDomainNameAndMatId[name] = matid;
		idPerm[matid]= permeability;
		sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(matid,phi, 1.0);
	}
    
	sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
	
	// ------------------------ Reading 3D Domain BC matids ------------------------
	if(input.find("Boundary") == input.end()) DebugStop();
	std::map<int,std::pair<int,REAL>>& BCFlowMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue;
    std::map<int,std::pair<int,REAL>>& BCTransportMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue;
	for(auto& bc : input["Boundary"]){
		if(bc.find("matid") == bc.end()) DebugStop();
		if(bc.find("type") == bc.end()) DebugStop();
		if(bc.find("value") == bc.end()) DebugStop();
		const int matid = bc["matid"];
		const int type = bc["type"];
		const REAL value = bc["value"];
		
		if(BCFlowMatIdToTypeValue.find(matid) != BCFlowMatIdToTypeValue.end()) DebugStop();
		BCFlowMatIdToTypeValue[matid] = std::make_pair(type, value);
        BCTransportMatIdToTypeValue[matid] = std::make_pair(type, value);
	}
	
	// ------------------------ Reading fractures and fracture bcs matids ------------------------
    int initfracmatid =input["FractureInitMatId"];
    int actualfracid = initfracmatid;
    REAL phiintersec=0;
    REAL inersecLenght =0;
    for(auto& fracture : input["Fractures"]){
        const int i = fracture["Index"];
        std::string name = "Fracture" + std::to_string(i);
        const int matid = actualfracid;
        const REAL permerm = fracture["K"];
		const REAL fracWidth = fracture["width"];
//		const REAL fracWidth = 0.1;
		TMRSDataTransfer::TFracProperties::FracProp fracprop;
        sim_data.mTGeometry.mDomainFracNameAndMatId[name] = matid;
        actualfracid +=5;
		
		// Fracture properties
		fracprop.m_perm = permerm;
		fracprop.m_width = fracWidth;
		fracprop.m_fracbc.insert(matid + 1);
        fracprop.m_fracIntersectMatID = matid + 2;
        
		// Fracture polygon
        int npoints = fracture["Nodes"].size();
        TPZFMatrix<REAL> polygonmatrix(3,npoints);
        for(int j=0; j<npoints; j++){
            for(int k=0; k<3; k++){
                polygonmatrix(k,j) = (REAL)fracture["Nodes"][j][k];
            }
        }
        fracprop.m_polydata.SetCornersX(polygonmatrix);
        
		// Setting data structure
        sim_data.mTFracProperties.m_fracprops[matid] = fracprop;
        
		// Fracture boundary conditions
        const REAL zero_flux = 0., zero_pressure = 0.;
        const int bcFracType_Dirichlet = 0, bcFracType_Neumann = 1;
        if(fracprop.m_fracbc.size() == 1){
            sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[*fracprop.m_fracbc.begin()] = std::make_pair(bcFracType_Neumann, zero_flux);
        }
        else{
            DebugStop(); // for now, dfnimrs only allows one bc which is zero flux
        }
        sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[fracprop.m_fracIntersectMatID] = std::make_pair(bcFracType_Dirichlet, zero_pressure);
		
		// Setting fracture porosity in reservoir properties data structure
		// NOTE: Width is stored in two place for now. In mPorosityAndVolumeScale and in fracprop.m_width. Dangerous...
		if(fracture.find("phi") == fracture.end()) DebugStop();
		const REAL phifrac = fracture["phi"];
		sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(matid, phifrac, fracWidth);
        
        phiintersec = phifrac;
        inersecLenght = fracWidth;

    }
//    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(299,0.2, 0.01);
	

    int FractureHybridPressureMatId = input["FractureHybridPressureMatId"];
    // this is the material id of the pressure between hybridized fluxes of intersecting fractures
    sim_data.mTGeometry.m_pressureMatId = FractureHybridPressureMatId;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(FractureHybridPressureMatId, phiintersec, inersecLenght*inersecLenght);
    
    
    // the material id of the conductivity between fractures
    if(input.find("FractureGlueMatId") != input.end())
    {
        int fractureGlueMatId = input["FractureGlueMatId"];
        sim_data.mTFracIntersectProperties.m_FractureGlueId = fractureGlueMatId;
        if(input.find("FractureGluePerm") == input.end()) DebugStop();
        REAL fractureGluePerm = input["FractureGluePerm"];
        sim_data.mTFracIntersectProperties.m_FractureGluePerm = fractureGluePerm;
    }
    
    // Transport properties
    if(input.find("RunWithTransport") != input.end()){
        sim_data.mTNumerics.m_run_with_transport = input["RunWithTransport"];
        if(sim_data.mTNumerics.m_run_with_transport){
            if(input.find("DeltaT") == input.end()) DebugStop();
            sim_data.mTNumerics.m_dt = input["DeltaT"];
            if(input.find("NSteps") == input.end()) DebugStop();
            sim_data.mTNumerics.m_n_steps = input["NSteps"];
        }
    }
    else{
        DebugStop(); // Please set run with tranport in json
//        sim_data.mTNumerics.m_n_steps = 100;
//        sim_data.mTNumerics.m_dt      = 1.e7;//*day;
    }
	
	// ------------------------ Setting extra stuff that is still not in JSON ------------------------
	const int D_Type = 0, N_Type = 1, Mixed_Type = 2;
	sim_data.mTGeometry.mInterface_material_id = 100;
	sim_data.mTGeometry.mInterface_material_idFracInf = 102;
	sim_data.mTGeometry.mInterface_material_idFracSup = 101;
	sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
	sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
	sim_data.mTGeometry.mSkeletonDiv = 0;
	sim_data.mTNumerics.m_sfi_tol = 0.0001;
	sim_data.mTNumerics.m_res_tol_transport = 0.0001;
	sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
	sim_data.mTNumerics.m_four_approx_spaces_Q = true;
	std::vector<REAL> grav(3,0.0);
	grav[1] = 0.0;//-9.8*(1.0e-6); // hor
	sim_data.mTNumerics.m_gravity = grav;
	sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = glob_n_threads;
	sim_data.mTNumerics.m_max_iter_sfi=1;
	sim_data.mTNumerics.m_max_iter_mixed=1;
	sim_data.mTNumerics.m_max_iter_transport=1;
	
	// PostProcess controls
//	std::string vtkfilename = filenameBase.substr(filenameBase.find("dfnimrs/") + 8);
//	std::replace( vtkfilename.begin(), vtkfilename.end(), '/', '_');
    std::string vtkfilename = filenameBase.substr(filenameBase.find_last_of("/") + 1);
	std::string pressurevtk = outputFolder + vtkfilename + ".vtk";
	std::string transportvtk = outputFolder + vtkfilename + "_transport.vtk";
	std::cout << "\n===> PostProcess file name: " << vtkfilename << std::endl;
	
	
	
	sim_data.mTPostProcess.m_file_name_mixed = pressurevtk;
	sim_data.mTPostProcess.m_file_name_transport = transportvtk;
	TPZStack<std::string,10> scalnames, vecnames, scalnamesTransport;
	vecnames.Push("Flux");
	scalnames.Push("Pressure");
    scalnames.Push("div_q");
	if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
		scalnames.Push("g_average");
		scalnames.Push("p_average");
	}
	scalnamesTransport.Push("Sw");
	scalnamesTransport.Push("So");
	
	sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
	sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
	sim_data.mTPostProcess.m_scalnamesTransport = scalnamesTransport;
	
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

void MapFractureIntersection(const std::string &filenameBase, std::map<std::string,int> &matmap, int firstfracintersect, std::map<int,std::pair<int,int>> &matidtoFractures)
{
    // This method loops over the physical tags for intersections created by DFNMesh
    // Creating gmsh reader
    TPZGmshReader  Geometry;
    std::string filename = filenameBase + "_fine.msh";
    std::ifstream input(filename);
    Geometry.ReadPhysicalProperties4(input); // There are some automatic associations that are not used. So dont worry if the code alerts about it
    auto physicaltags = Geometry.GetDimPhysicalTagName();
    // loop over all physical tags of dimension 1
    const std::string begins("fracIntersection");
    const int beglength = begins.size();
    // Here we select fracture intersection tags
    for(auto iter : physicaltags[1])
    {
        auto name = iter.second;
        if(name.compare(0,beglength,begins) == 0)
        {
            auto pos1 = name.find_first_of("_",0);
            auto pos2 = name.find_last_of("_");
            std::string strfac1 = name.substr(pos1+1,pos2-pos1-1);
            std::string strfac2 = name.substr(pos2+1);
            int frac1 = std::stoi(strfac1);
            int frac2 = std::stoi(strfac2);
            matmap[name] = firstfracintersect; // firstfracintersect is an id used that it is not any of the used ids so far (it is unique)
            
            // This map stores for each fracture with matid firstfracintersect the fracs frac1 and frac2. Each intersection has a unique id given by firstfracintersect
            matidtoFractures[firstfracintersect] = std::pair<int,int>(frac1,frac2);
            firstfracintersect++;
        }
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ReadMeshesDFN(string& filenameBase, TPZGeoMesh*& gmeshfine, TPZGeoMesh*& gmeshcoarse, int& initVolForMergeMeshes, bool& isMHM, bool& needsMergeMeshes) {
	using json = nlohmann::json;
	std::string filenamejson = filenameBase + ".json";
	
	std::ifstream filejson(filenamejson);
    if(!filejson) DebugStop();
	json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file

	// ===================> Coarse mesh <=======================
	TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
	TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
	
	std::set<int> allmatids; // used to check for repeated matids and highest matid
	
    // ------------------------ Check if isMHM is set on file ------------------------
    isMHM = true;
    if(input.find("useMHM") != input.end()) {
        isMHM = input["useMHM"];
    }
    
    if(input.find("needsMerge") != input.end()) {
        needsMergeMeshes = input["needsMerge"];
    }
    
    if(isMHM && !needsMergeMeshes) DebugStop();
    
    int SimulationDim= 3;
    if(input.find("SimulationDim") != input.end()) {
        SimulationDim = input["SimulationDim"];
    }
    
	// ------------------------ Get matids of 3D domain ------------------------
	if(input.find("Domains") == input.end()) DebugStop();
	for(auto& domain : input["Domains"]){
		if(domain.find("matid") == domain.end()) DebugStop();
		if(domain.find("name") == domain.end()) DebugStop();
		const int matid = domain["matid"];
		const string name = domain["name"];
		dim_name_and_physical_tagCoarse[SimulationDim][name] = matid;
		const bool is_in = allmatids.find(matid) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(matid);
	}
	
	// ------------------------ Get matids of BCs of 3D domain ------------------------
	// NOTE: This might not be needed for the coarse mesh. Please check...
	if(input.find("Boundary") == input.end()) DebugStop();
	for(auto& bc : input["Boundary"]){
		if(bc.find("name") == bc.end()) DebugStop();
		if(bc.find("matid") == bc.end()) DebugStop();
		const string name = bc["name"];
		const int matid = bc["matid"];
		dim_name_and_physical_tagCoarse[SimulationDim-1][name] = matid;
		dim_name_and_physical_tagFine[SimulationDim-1][name] = matid;
		const bool is_in = allmatids.find(matid) != allmatids.end();
		if(is_in) DebugStop();
		allmatids.insert(matid);
	}

	// ------------------------ Generate gmesh coarse ------------------------
    if(input.find("Mesh") == input.end()) DebugStop();
    std::string meshfile = input["Mesh"];
    {
        auto lastslash = filenamejson.find_last_of("/");
        auto meshdirname = filenamejson.substr(0,lastslash+1);
//        meshFile = meshFile.substr(meshFile.find("examples/") + 9,meshFile.length());
        meshfile = meshdirname + meshfile;
    }
    int ncoarse_vol = 0;
    if(needsMergeMeshes){
        gmeshcoarse = generateGMeshWithPhysTagVec(meshfile,dim_name_and_physical_tagCoarse);
        
        int64_t nelcoarse = gmeshcoarse->NElements();
        for(int64_t el = 0; el<nelcoarse; el++)
        {
            TPZGeoEl *gel = gmeshcoarse->Element(el);
            if(gel && gel->Dimension()==3) ncoarse_vol++;
        }
    }
	
//    string filenameCoarse = filenameBase + "_coarse.msh";
//    gmeshcoarse = generateGMeshWithPhysTagVec(filenameCoarse,dim_name_and_physical_tagCoarse);

   
	// ===================> Fine mesh <=======================
	// Note that boundary elements have been added previously to the set dim_name_and_physical_tagFine

	if(input.find("FractureInitMatId") == input.end()) DebugStop();
	const int fracInitMatId = input["FractureInitMatId"];
    const int fracinc = 5;
	
	// ------------------------ Loop over fractures in Json ------------------------
	int fracCounter = 0;
    bool isFracSim=false;
    if(input.find("Fractures") != input.end()){
	for(auto& frac : input["Fractures"]){
		const int matid = fracInitMatId + fracCounter*fracinc;
        const int bcmatid = fracInitMatId + fracCounter*fracinc+1;
		string fracname = "Fracture" + to_string(fracCounter);
		string bcfracname = "BCfrac" + to_string(fracCounter);
		const int currentFracId = fracInitMatId + fracinc * fracCounter;
        dim_name_and_physical_tagFine[2][fracname] = currentFracId;
        dim_name_and_physical_tagFine[1][bcfracname] = (currentFracId)+1;
        bool is_in = allmatids.find(currentFracId) != allmatids.end();
        if(is_in) DebugStop();
        allmatids.insert((currentFracId));
        is_in = allmatids.find((currentFracId)+1) != allmatids.end();
        if(is_in) DebugStop();
        allmatids.insert((currentFracId)+1);
		fracCounter++;
	}
        isFracSim=true;
    }
	
	// ------------------------ Adding volume physical tags------------------------
	const int maxMatId = *allmatids.rbegin();
	initVolForMergeMeshes = (1+maxMatId/100)*100;
	if(input.find("NCoarseGroups") == input.end()) DebugStop();
	const int nCoarseGroups = input["NCoarseGroups"];
    
    std::string volbase = "c";
    if(needsMergeMeshes){
        if(nCoarseGroups != ncoarse_vol) DebugStop();
        for (int ivol = 0; ivol < nCoarseGroups; ivol++) {
            std::string ivolstring = volbase + to_string(ivol);
            dim_name_and_physical_tagFine[SimulationDim][ivolstring] = initVolForMergeMeshes + ivol;
        }
    }
	// ------------------------ Generate gmesh fine ------------------------
	string filenameFine = filenameBase + "_fine.msh";
    
    /// add the intersection material ids to be read
    /// read only the header of the msh file
    /// identify the intersection groups by identify the substring
    if(input.find("FractureHybridPressureMatId") == input.end()) DebugStop();
    int FractureHybridPressureMatId = input["FractureHybridPressureMatId"];
    
    int firstfracintersect = initVolForMergeMeshes + ncoarse_vol;
    firstfracintersect = (1+firstfracintersect/100)*100;
    std::map<int,std::pair<int,int>> matidtoFractures;
    
    if(isFracSim){
        MapFractureIntersection(filenameBase, dim_name_and_physical_tagFine[1], firstfracintersect, matidtoFractures);
    }
    
    if (!isMHM && !needsMergeMeshes) {
        gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagCoarse);
    }
    else{
        gmeshfine = generateGMeshWithPhysTagVec(filenameFine,dim_name_and_physical_tagFine);
    }

    // for each intersection element create an intersection element for each fracture specifically
    if (isFracSim) {
        CreateIntersectionElementForEachFrac(gmeshfine,matidtoFractures,fracInitMatId,fracinc,FractureHybridPressureMatId);
    }
	
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateIntersectionElementForEachFrac(TPZGeoMesh* gmeshfine,
										  std::map<int,std::pair<int,int>>& matidtoFractures,
										  const int fracInitMatId, const int fracinc, const int FractureHybridPressureMatId) {
	int64_t nelem = gmeshfine->NElements();
	for(int64_t el = 0; el<nelem; el++)
	{
		TPZGeoEl *gel = gmeshfine->Element(el);
		if(!gel) continue;
		int matid = gel->MaterialId();
		auto it = matidtoFractures.find(matid);
		if(it != matidtoFractures.end())
		{
			if(gel->Dimension() > 1) DebugStop(); // has to be 1d (intersection)
			TPZGeoElSide gelside(gel);
			int frac1 = it->second.first;
            int matidfrac1 = fracInitMatId + frac1*fracinc;
			int matid1 = matidfrac1 + 2;
			int frac1bcid = matidfrac1 + 1;
			auto frac1bc = gelside.HasNeighbour(frac1bcid);
            auto frac1el = gelside.HasNeighbour(matidfrac1);
            if(!frac1el) DebugStop(); // Element has to have a fracture boundary
			if(frac1bc)
			{
#ifdef PZ_LOG
                if(fracIntersectLogger.isDebugEnabled()){
                    std::stringstream sout;
                    sout << "Changing matid of el with matid " << frac1bc.Element()->MaterialId() << " to matid " << matid1 << std::endl;
                    LOGPZ_DEBUG(fracIntersectLogger, sout.str())
                }
#endif
				frac1bc.Element()->SetMaterialId(matid1);
			}
			else
			{
#ifdef PZ_LOG
                if(fracIntersectLogger.isDebugEnabled()){
                    std::stringstream sout;
                    sout << "Creating boundary condition with matid " << matid1 << std::endl;
                    LOGPZ_DEBUG(fracIntersectLogger, sout.str())
                }
#endif
				TPZGeoElBC(gelside,matid1);
			}
			int frac2 = it->second.second;
            int matidfrac2 = fracInitMatId + frac2*fracinc;
			int matid2 = matidfrac2 + 2;
			int frac2bcid = matidfrac2 + 1;
			auto frac2bc = gelside.HasNeighbour(frac2bcid);
            auto frac2el = gelside.HasNeighbour(matidfrac2);
            if(!frac2el) DebugStop();
			if(frac2bc)
			{
#ifdef PZ_LOG
                if(fracIntersectLogger.isDebugEnabled()){
                    std::stringstream sout;
                    sout << "Changing matid of el with matid " << frac2bc.Element()->MaterialId() << " to matid " << matid2 << std::endl;
                    LOGPZ_DEBUG(fracIntersectLogger, sout.str())
                }
#endif

				frac2bc.Element()->SetMaterialId(matid2);
			}
			else
			{
#ifdef PZ_LOG
                if(fracIntersectLogger.isDebugEnabled()){
                    std::stringstream sout;
                    sout << "Creating boundary condition with matid " << matid2 << std::endl;
                    LOGPZ_DEBUG(fracIntersectLogger, sout.str())
                }
#endif
				TPZGeoElBC(gelside,matid2);
			}
			// create the geometric element for receiving the hybridized pressure element
			auto fracintersect = gelside.HasNeighbour(FractureHybridPressureMatId);
			if(!fracintersect)
			{
				TPZGeoElBC(gelside,FractureHybridPressureMatId);
			}
		}
	}
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ModifyPermeabilityForCase2(TPZGeoMesh* gmesh) {
//	DebugStop(); // fix me or generate the gmsh mesh correcty from the start and erase me
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != 1 && gel->MaterialId() != 2) continue; // only matrix

        TPZVec<REAL> masscent(3,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);

        const REAL x = xcenter[0], y = xcenter[1], z = xcenter[2];
		gel->SetMaterialId(1);
        if (x > 0.5 && y < 0.5)
            gel->SetMaterialId(2);

        if (x > 0.75 && y > 0.5 && y < 0.75 && z > 0.5)
            gel->SetMaterialId(2);

        if (x > 0.625 && x < 0.75 && y > 0.5 && y < 0.625 && z > 0.5 && z < 0.75)
            gel->SetMaterialId(2);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ModifyBCsForCase2(TPZGeoMesh* gmesh) {
	DebugStop(); // fix me or generate the gmsh mesh correcty from the start and erase me
        
//    const REAL inletDomain = 0.25, outletDomain = 0.875;
//    const REAL zerotol = ZeroTolerance();
//
//    for (auto gel: gmesh->ElementVec()) {
//        if (!gel) continue;
//        if (gel->MaterialId() != EFaceBCPressure) continue; // 2d faces on boundary only
//
//        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
//        gel->CenterPoint(gel->NSides()-1, masscent);
//        gel->X(masscent, xcenter);
//        const bool isXzero = fabs(xcenter[0]) < zerotol, isYzero = fabs(xcenter[1]) < zerotol, isZzero = fabs(xcenter[2]) < zerotol;
//        const bool isXone = fabs(xcenter[0]-1.) < zerotol, isYone = fabs(xcenter[1]-1.) < zerotol, isZone = fabs(xcenter[2]-1.) < zerotol;
//
//        gel->SetMaterialId(ENoflux); // Default is no flux
//        // Setting inlet BCs
//        if(isXzero && (xcenter[1] < inletDomain && xcenter[2] < inletDomain))
//                gel->SetMaterialId(EOutlet);
//        if(isYzero && (xcenter[0] < inletDomain && xcenter[2] < inletDomain))
//                gel->SetMaterialId(EOutlet);
//        if(isZzero && (xcenter[0] < inletDomain && xcenter[1] < inletDomain))
//                gel->SetMaterialId(EOutlet);
//
//        // Setting outlet BCs
//        if(isXone && (xcenter[1] > outletDomain && xcenter[2] > outletDomain))
//                gel->SetMaterialId(EInlet);
//        if(isYone && (xcenter[0] > outletDomain && xcenter[2] > outletDomain))
//                gel->SetMaterialId(EInlet);
//        if(isZone && (xcenter[0] > outletDomain && xcenter[1] > outletDomain))
//                gel->SetMaterialId(EInlet);
//    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ModifyBCsForCase3(TPZGeoMesh* gmesh) {
    
    const REAL zerotol = ZeroTolerance();
    const REAL onethird = 1./3., twothirds = 2./3.;
    
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
//        if (gel->MaterialId() != ENoflux) continue; // 2d faces on boundary only
        if(gel->Dimension() != 2){continue;}
        if((gel->MaterialId() != ENoflux) && (gel->MaterialId() != EInlet) && (gel->MaterialId() != EOutlet)){continue;}
        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        const REAL x = xcenter[0], y = xcenter[1], z = xcenter[2];
        const bool isYzero = fabs(y) < zerotol;
        const bool isYend = fabs(y-2.25) < zerotol;
        
        // Setting inlet BCs
        gel->SetMaterialId(ENoflux);
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

void ModifyBCsForCase4(TPZGeoMesh* gmesh) {
	DebugStop(); // fix me or generate the gmsh mesh correcty from the start and erase me
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

void ModifyBCsFor2ParallelFractures(TPZGeoMesh* gmesh) {
   
    const REAL zerotol = ZeroTolerance();
    for (auto gel: gmesh->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() != 2) continue; // 2d faces on boundary only
        
        TPZVec<REAL> masscent(2,0.0), xcenter(3,0.0);
        gel->CenterPoint(gel->NSides()-1, masscent);
        gel->X(masscent, xcenter);
        const REAL x = xcenter[0], y = xcenter[1], z = xcenter[2];
        const bool isXinit = fabs(x) < zerotol;
        const bool isXend = fabs(x-2) < zerotol;
        
        // Default is no flux already set previously
        // Setting inlet BCs
        if(isXinit)
                gel->SetMaterialId(3);
        if(isXend )
                gel->SetMaterialId(4);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void fixPossibleMissingIntersections(TMRSDataTransfer& sim_data, TPZGeoMesh* gmesh){
    // it may happen that the supplied mesh from DFN has missing intersections where 2 fractures intersect.
    // We treat this here. But it would be ideal to make DFN robust enough so we could remove this function whatsoever
    // The main idea is to check if a fracture element has more than 1 fracture neightbor. If so, an intersection element is
    // needed to hybridzie the region.
    cout << "\n---------------------- Searching for problematic intersection elements ----------------------" << endl;
    TPZSimpleTimer timer_fixinter;
    
    int nInterCreated = 0;
    for (auto gel : gmesh->ElementVec()) {
        if (!gel || !sim_data.mTFracProperties.isFracMatId(gel->MaterialId())) {
            continue;
        }
        int gelmatid = gel->MaterialId();
        if (gel->Dimension() != 2) {
            DebugStop(); // Should it work with 2D problems (1d fractures)?
        }
		
        auto fracprop = sim_data.mTFracProperties.m_fracprops[gelmatid];
//        const int matidintersect = sim_data.mTFracIntersectProperties.m_IntersectionId;
        const int matidintersect = fracprop.m_fracIntersectMatID;

        const int firstedge = gel->FirstSide(1);
        const int lastedge = gel->FirstSide(2);
        // loop over the sides of dimension 1
        for (int iside = firstedge; iside < lastedge; iside++) {
            TPZStack<TPZGeoEl*> interEls;
            TPZStack<TPZGeoEl*> bcEls;
            int nFracElsForSide = 0;
            // loop over the neighbours
            TPZGeoElSide gside(gel,iside);
            TPZGeoElSide neig = gside.Neighbour();
            for (; gside != neig; neig++) {
                TPZGeoEl* neigel = neig.Element();
				const int neighelmatid = neigel->MaterialId();
                if (neighelmatid == matidintersect) {
                    interEls.push_back(neigel);
                }
                // compute the number of fracture elements linked to this side
                if (sim_data.mTFracProperties.isFracMatId(neighelmatid)) {
                    nFracElsForSide++;
                }
                // compute the number of boundary condition elements linked to this side
                if (sim_data.mTFracProperties.isFracBCMatId(neighelmatid)) {
                    bcEls.push_back(neigel);
                }
            }
//            if(bcEls.size() > 1) DebugStop();
            if (interEls.size() > 1) {
                if (interEls.size() > 2) {
                    DebugStop(); // there are 3 intersection elements in the same place! Please check why...
                }
                cout << "Found two intersection elements at the same place! Is this a problem?" << endl;
                cout << "Manually deleting intersection geoel..." << endl;

                TPZGeoEl* gelduplicate = interEls[1];
                const int64_t duplicateIndex = gelduplicate->Index();
                gelduplicate->RemoveConnectivities();
                delete gelduplicate;
                gmesh->ElementVec()[duplicateIndex] = nullptr;
                
            }
            if (nFracElsForSide > 1 && !interEls.size() && !bcEls.size()) {
#ifdef PZDEBUG
                // Uncomment for debugging
//                cout << "Warning! Edge of el index " << gel->Index() << " has " << nFracElsForSide << " fracture elements as neighbors but does not have intersections nor bcs" << std::endl;
                // This can happen if the polygons do not intersect but the mesh generated by DFNMesh created small elements in this location.
                // However, the sideorient of the adjacent elements of a same fracture need to be opposite. This is checked later.
#endif
            }
            if (0 && (interEls.size() || nFracElsForSide > 1) && bcEls.size()) {
                cout << "PLEASE CHECK CAREFULLY! An element may have a boundary even if it intersects\n";
                cout << "Domains intersecting and boundary through snap\n Intersection matids ";
                for(auto it : interEls) cout << it->MaterialId() << " ";
                cout << "\n Boundary condition material ids ";
                for(auto it : bcEls) cout << it->MaterialId() << " ";
                cout << std::endl;
                cout << "matid of all neighbours " << gel->MaterialId();
                neig = gside.Neighbour();
                for (; gside != neig; neig++) {
                    TPZGeoEl* neigel = neig.Element();
                    const int neighelmatid = neigel->MaterialId();
                    cout << " " << neighelmatid;
                }
                cout << std::endl;
//                if(bcEls.size() > 1)
//                    DebugStop(); // This is rather odd. Please check why there are two bcs at the same boundary
//                for (auto bcgeoel : bcEls) {
//                    bcgeoel->RemoveConnectivities();
//                    delete bcgeoel;
//                }
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
    
    std::cout << "Total time: " << timer_fixinter.ReturnTimeDouble()/1000 << " seconds" << std::endl;

	// No need to buildConnectivity. It is already correct
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillPCteSol(TPZCompMesh* mpcmesh, const REAL pcte) {
	
	int64_t nc = mpcmesh->NConnects();
	TPZFMatrix<STATE> &sol = mpcmesh->Solution();
	for (int64_t ic = 0; ic<nc; ic++) {
		TPZConnect &c = mpcmesh->ConnectVec()[ic];
		int64_t seqnum = c.SequenceNumber();
		if(seqnum < 0) continue;
		unsigned char lagrange = c.LagrangeMultiplier();
		STATE fill = 0.;
		if(lagrange == 1 || lagrange == 3 || lagrange == 6 || lagrange == 7)
		{
			fill = pcte;
		}
		int ndof = c.NShape();
		for (int idf = 0; idf < ndof ; idf++) {
			int64_t index = mpcmesh->Block().Index(seqnum, idf);
			sol(index,0) = fill;
		}

	}
		
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

std::map<int,TPZVec<STATE>> computeIntegralOfNormalFlux(const std::set<int> &bcMatId, TPZMultiphysicsCompMesh *cmesh) {
	
	// modifying bc of top outlet of case 3 to compute integral separately
	TPZGeoMesh* gmesh = cmesh->Reference();
	const REAL zerotol = ZeroTolerance();
    std::map<int,TPZVec<STATE>> result;
	std::set<int> matidsInlet, matidsOutlet;
	std::string varname="BCNormalFlux";
	cmesh->Reference()->ResetReference();
	cmesh->LoadReferences(); // compute integral in the multiphysics mesh
	int nels = cmesh->NElements();
	for(auto it: bcMatId)
    {
        std::set<int> matid = {it};
        TPZVec<STATE> vecint = cmesh->Integrate(varname, matid);
        result[it] = vecint;
        std::cout << "Integral of normal flux for matid " << it << " = " << vecint << std::endl;
        
    }
    return result;
	
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

bool fileExists(const fs::path& p, fs::file_status s) {
    if(fs::status_known(s) ? fs::exists(s) : fs::exists(p))
        return true;
    else
        return false;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateOutputFolders(std::string& outputFolder) {

    std::string folders = outputFolder;
    char c = '/';
    std::string folderToCreate = "";
    int nfolders = 0;
    while (folders.find("/") != std::string::npos) {
        if(nfolders == 0) folderToCreate = folders.substr(0,folders.find("/"));
        else folderToCreate = folderToCreate + "/" + folders.substr(0,folders.find("/"));
        folders = folders.substr(folders.find("/")+1);
        if(!fileExists(folderToCreate)){
            if (!fs::create_directory(folderToCreate))
                DebugStop();
            else
                cout << "Directory created with name " << folderToCreate << endl;
        }
        nfolders++;
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CopyInputFilesToOutputFolderAndFixFilename(std::string& filenameBase, std::string& outputFolder){
    
    // Copying all files from input folder to output folder
    std::string onlyFolder = outputFolder.substr(0,outputFolder.find_last_of("/"));
    fs::copy(filenameBase , onlyFolder, fs::copy_options::update_existing | fs::copy_options::recursive);
    
    // Adding the json name to filenameBase
    int njson = 0;
    for (auto const& dir_entry : std::filesystem::directory_iterator{filenameBase}){
        std::string filename = dir_entry.path().string();
//        std::cout << dir_entry.path().string() << endl;
        
        if (filename.find(".") != std::string::npos){
            if(filename.substr(filename.find(".")) == ".json"){
                if(njson == 0){
                    filenameBase = filename.substr(0,filename.find("."));
                    njson++;
                }
                else {
                    cout << "\n\n=====> ERROR! There are two json files in the provided input folder" << endl;
                    DebugStop();
                }
            }
        }
    }

}

TPZGeoMesh * Transform2dMeshToUnisim3D(TPZGeoMesh* gmesh2d, int nLayers){
   
    REAL w = 200.0;
    std::string name2D("mesh2d.vtk");
    
    int topID= 5;
    int baseID = 5;
    TPZGeoMesh * returnedMesh = nullptr;
        TPZExtendGridDimension extend(gmesh2d, w);
        extend.SetElType(1);
        returnedMesh = extend.ExtendedMesh(nLayers,topID,baseID);
        ModifyTopeAndBase2(returnedMesh ,nLayers);
       
        std::ofstream file("unisim3d.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(returnedMesh,file );
    return returnedMesh;
    
}

void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers){
//    std::string filename2 = "Reservoir/base_unisimMOD.txt";
    
     std::string filename1 = "topeMOD.txt";
    std::string filename2 = "baseMOD.txt";
    std::vector<double> x, y, z, x1,y1,z1;
    ReadData(filename1, true, x, y, z);
    ReadData(filename2, true, x1, y1, z1);

    _2D::ThinPlateSplineInterpolator <double> interpTope;
    _2D::ThinPlateSplineInterpolator <double> interpBase;

    interpTope.setData(x,y,z);
    interpBase.setData(x1,y1,z1);

    int nCoordinates = gmesh->NodeVec().NElements();
    double sum=0.0;
    for (auto val:z) {
        sum += val;
    }
    double val_tope= sum / z.size();
    sum=0.0;
    for (auto val:z1) {
        sum += val;
    }
    double val_base= sum / z1.size();
//    val_base = 1000;
//    val_tope = 5000;
//
//    val_tope = 3000;
//    val_base = 3000;
    int npointsPerLayer = nCoordinates/(nlayers+1);
    double valinter=0.0;
    for (int ilay = 1; ilay <= nlayers+1; ilay++) {
        for (int ipoint = (ilay-1)*npointsPerLayer; ipoint<(ilay)*npointsPerLayer; ipoint++) {
            TPZGeoNode node = gmesh->NodeVec()[ipoint];
            TPZVec<REAL> co(3);
            node.GetCoordinates(co);
            double topeinterpol =interpTope(co[0],co[1]);
            double baseinterpol = interpBase(co[0],co[1]);
            if (topeinterpol==0) {
                topeinterpol = val_tope;
                if (co[0]>1000.00) {
                    topeinterpol -= 120;
                }
            }
            if (baseinterpol==0) {
                
                baseinterpol = val_base;
                if (co[0]>1000.00) {
                   baseinterpol = val_base-80;
                }

            }

            if (ilay==1) {
                valinter=topeinterpol;
//                valinter = 3500;
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
            if (ilay==nlayers+1) {
                valinter = baseinterpol;
//                valinter = 2850;
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
            if (ilay>1   && ilay < nlayers+1) {
                valinter = topeinterpol + (ilay-1)*(baseinterpol - topeinterpol)/(nlayers);
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
        }
    }
}
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z){
    
    bool modpoints = true;
    std::ifstream file;
    string basemeshpath(FRACMESHES);
    basemeshpath = basemeshpath + "/dfnimrs/unisim_meshes/Reservoir_props/" + name;
    file.open(basemeshpath);
    int i=1;
    
    
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        if(l != '/'){
            i=i+1;
            int val = i%15;
            if(val ==0){
                double a, b, c;
                if(iss >> a >> b >> c) ;
                if (modpoints) {
                    x.push_back(a - 350808.47);
                    y.push_back(b - 7.51376238e6);
                    z.push_back(c);
                }
                else{
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
                }
            };
        };
    };
    
    if(x.size() == 0){
        std::cout<<"No data read."<<std::endl;
        
        DebugStop();
    }
    if(print_table_Q){
        std::cout<<"*************************"<<std::endl;
        std::cout<<"Reading file... ok!"<<std::endl;
        std::cout<<"*************************"<<std::endl;
        std::cout<<x.size()<<std::endl;
        std::cout<<y.size()<<std::endl;
        std::cout<<z.size()<<std::endl;
    }
    
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void ComputeDiagnostics(std::string& outputFolder, TMRSDataTransfer& sim_data, std::set<int>& bcflux, TPZMultiphysicsCompMesh* mixed_operator) {
    std::ofstream out(outputFolder + "fluxintegral.txt");
    std::cout << "\n**************************** Domain Diagnostics **************************** " << std::endl;
    out << "\n**************************** Domain Diagnostics **************************** " << std::endl;
    auto result = computeIntegralOfNormalFlux(bcflux, mixed_operator);
    for(auto it: result)
    {
        out << "Integral for matid " << it.first << " " << it.second << std::endl;
    }
    std::cout << "\n**************************** Fracture Fluxes Diagnostics **************************** " << std::endl;
    out << "\n**************************** Fracture Fluxes Diagnostics **************************** " << std::endl;
    auto &allfrac = sim_data.mTFracProperties.m_fracprops;
    int gluematid = sim_data.mTFracIntersectProperties.m_FractureGlueId;
    int pressureintersect = sim_data.mTGeometry.m_pressureMatId;
    std::set<int> fracmatids;
    for(auto &it : allfrac) fracmatids.insert(it.first);
    for(auto &it : allfrac)
    {
        FractureQuantities frac(mixed_operator,fracmatids, pressureintersect,gluematid);
        frac.ComputeFluxQuantities(it.first);
        frac.Print(out);
        frac.Print(std::cout);
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void VerifyIfNeumannIsExactlyZero(const int matidNeumann, TPZMultiphysicsCompMesh* mixed_operator) {
    TPZCompMesh* fluxmesh = mixed_operator->MeshVector()[0];
    bool passed = true;
    TPZFMatrix<STATE>& solmat = fluxmesh->Solution();
    for(auto cel : fluxmesh->ElementVec()){
        if(!cel) continue;
        
        TPZGeoEl* gel = cel->Reference();
        if(!gel) continue; // SubCompMesh
        const int matid = gel->MaterialId();
        if(matid != matidNeumann) continue; // For now checking only Neumann in domain boundary
        if(gel->Dimension() != 2) DebugStop();
        
        TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
        if(!intel) DebugStop();
        
        const int ncon = intel->NConnects();
        if(ncon != 1) DebugStop();
        
        TPZConnect& con = intel->Connect(0);
        const int64_t seq = con.SequenceNumber();
        const int64_t firsteq = fluxmesh->Block().Position(seq);
        const int64_t blocksize = fluxmesh->Block().Size(seq);
        
        for (int i = 0; i < blocksize; i++) {
            const REAL solzero = solmat(firsteq+i,0);
            if (!IsZero(solzero)) {
                std::cout << "Position " << firsteq+i << " has a dof with value " << solzero << std::endl;
                passed = false;
            }
        }
        
    }
    if (!passed) {
        std::cout << "\n\nERROR! Zero Neumann filter did not work. Check above for diagnostics" << std::endl;
        DebugStop();
    }
    else {
        std::cout << "\n===> Filter Worked! All Neumann zero are zero." << std::endl;
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FilterZeroNeumann(std::string& outputFolder, TMRSDataTransfer& sim_data, TPZAutoPointer<TPZStructMatrix> strmat, TPZCompMesh* cmesh) {

    std::cout << "\n---------------------- Filtering zero neumann equations ----------------------" << std::endl;
    TPZSimpleTimer timer_filter("Timer Filter Equations");
    
   
    std::set<int64_t> matidset;
    
    // First find all the zero neumann in in the 3d domain
    std::cout << "Domain BC matids: ";
    for(auto &chunk : sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
        const int bc_id = chunk.first;
        const std::pair<int,REAL>& typeAndVal = chunk.second;
        const int bc_type = typeAndVal.first;
        const REAL val = typeAndVal.second;
        if(bc_type == 1 && fabs(val) < ZeroTolerance()){
            std::cout << bc_id << " ";
            matidset.insert(bc_id);
        }
    }
    
    // Then all the zero neumann in the fractures
    std::cout << "\nFracture BC matids: ";
    for (auto& chunk : sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
        const int bc_id   = chunk.first;
        const std::pair<int,REAL>& typeAndVal = chunk.second;
        const int bc_type = typeAndVal.first;
        const REAL val = typeAndVal.second;
        if(bc_type == 1 && fabs(val) < ZeroTolerance()){
            std::cout << bc_id << " ";
            matidset.insert(bc_id);
        }
    }
    
    // Set equations to be filted
    std::set<int64_t> eqset;
    cmesh->GetEquationSetByMat(matidset, eqset);
    if (eqset.size()) {
        strmat->EquationFilter().ExcludeEquations(eqset);
    }
    
    int count = 0;
    for(auto cel : cmesh->ElementVec()){
        TPZSubCompMesh* subcmesh = dynamic_cast<TPZSubCompMesh*>(cel);
        if (subcmesh) {
            std::cout << "\n\t------- Submesh " << count << " -------" << std::endl;
//            if(count == 1 || count == 3){
//                count++;
//                std::cout << "===> Skipping filter" << std::endl;
//                continue;
//            }
            count++;
#ifdef PZDEBUG
            {
                std::string filename = outputFolder + "submesh_" + to_string(count) + ".vtk";
                std::ofstream out(filename);
                TPZVTKGeoMesh::PrintCMeshVTK(subcmesh, out);
            }
#endif
            std::set<int64_t> eqsetsub;
            subcmesh->GetEquationSetByMat(matidset, eqsetsub);
            const int64_t ninteq = subcmesh->NumInternalEquations();
            const int64_t neq = subcmesh->TPZCompMesh::NEquations();
            if (eqsetsub.size()) {
                subcmesh->Analysis()->StructMatrix()->EquationFilter().SetNumEq(neq); // Setting again bcz els could have been condensed
                subcmesh->Analysis()->StructMatrix()->EquationFilter().ExcludeEquations(eqsetsub); // will set as active all eqs (internal and external) that are not zero neumann
                auto& activeeq = subcmesh->Analysis()->StructMatrix()->EquationFilter().GetActiveEquations();
                
                std::cout << "size eqsetsub = " << eqsetsub.size() << " | eqsetsub = " << eqsetsub;
                std::cout << "size activeeq = " << activeeq.size() << " | active eq = " << activeeq << std::endl;
                auto& excludedeq = subcmesh->Analysis()->StructMatrix()->EquationFilter().GetExcludedEquations();
                std::cout << "size excludedeq = " << excludedeq.size() << " | excluded eq = " << excludedeq << std::endl;
                const int64_t nexteq = neq - ninteq;
#ifdef PZDEBUG
                // All the external equations must be present in the active equations
                int64_t count = 0;
                for(auto& eq : activeeq){
                    if(eq > ninteq-1) count++;
                }
                if(count < nexteq) DebugStop();
#endif
                // Now we want to set as active only the ones that are internal and not zero Neumann
                // Since activeeq is ordered, we simply pick equations until we reach the first non internal equation
                TPZVec<int64_t> activeinternal(ninteq);
                int64_t i = 0;
                for(auto& eq : activeeq){
                    if(eq > ninteq-1) break;
                    activeinternal[i++] = eq;
                }
                activeinternal.Resize(i);
                subcmesh->Analysis()->StructMatrix()->EquationFilter().Reset();
                subcmesh->Analysis()->StructMatrix()->EquationFilter().SetActiveEquations(activeinternal); // sets new filter
                // Note that, in subcmesh, after we create the data structure of the matred matrices, we need to the
                // set the external equations as active again
            }
        }
    }
    
    std::cout << "\n==> Total Filter time: " << timer_filter.ReturnTimeDouble()/1000. << " seconds" << std::endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
