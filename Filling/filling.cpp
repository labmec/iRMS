// C++ includes
// nothing

// PZ includes
#include <filesystem>

#include "TMRSApproxSpaceGenerator.h"
#include "TPZGenGrid3D.h"
#include "TPZSimpleTimer.h"
#include "TPZVTKGeoMesh.h"
#include "imrs_config.h"
#include "json.hpp"
#include "pzintel.h"
#include "pzlog.h"
#include "pzsmanal.h"

// ----- Namespaces -----
using namespace std;
namespace fs = std::filesystem;

// ----- End of namespaces -----

// ----- Global vars -----
const int glob_n_threads = 0;

// ----- Functions -----
TPZGeoMesh* ReadMeshFromGmsh(TMRSDataTransfer& sim_data);
void FillDataTransfer(std::string filenameBase, TMRSDataTransfer& sim_data);

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
int main(int argc, char* argv[]) {
  string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
  string logpath = basemeshpath + "/../Filling/log4cxx.cfg";
  TPZLogger::InitializePZLOG(logpath);
  if (mainlogger.isDebugEnabled()) {
    std::stringstream sout;
    sout << "\nLogger for Filling problem target\n"
         << endl;
    ;
    LOGPZ_DEBUG(mainlogger, sout.str())
  }
#endif

  // =========> Read the json file and fill the data transfer object
  TMRSDataTransfer sim_data;
  sim_data.mTNumerics.m_four_approx_spaces_Q = true;
  sim_data.mTNumerics.m_mhm_mixed_Q = false;
  sim_data.mTNumerics.m_need_merge_meshes_Q = false;
  sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
  FillDataTransfer(basemeshpath + "/../Filling/module", sim_data);

  // =========> Create GeoMesh
  TPZGeoMesh* gmesh = ReadMeshFromGmsh(sim_data);
  // Print gmesh to vtk
  std::ofstream out("gmesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // =========> Create Multiphysics comp mesh
  TMRSApproxSpaceGenerator aspace;
  aspace.SetDataTransfer(sim_data);
  aspace.SetGeometry(gmesh);
  aspace.BuildMixedMultiPhysicsCompMesh(1);
  TPZMultiphysicsCompMesh* mp_cmesh = aspace.GetMixedOperator();  

  // =========> Create Analysis
  RenumType renumtype = RenumType::EDefault;
  bool UsingPzSparse = true;  // Necessary to use multithread for now...
  bool UsePardiso_Q = true;   // lighting fast!
  cout << "\n---------------------- Creating Analysis (Might optimize bandwidth) ----------------------" << endl;
  TMRSMixedAnalysis* mixAnalisys = new TMRSMixedAnalysis(mp_cmesh, renumtype);
  mixAnalisys->SetDataTransfer(&sim_data);
  UsePardiso_Q = true;
  mixAnalisys->Configure(glob_n_threads, UsePardiso_Q, UsingPzSparse);
  TPZFastCondensedElement::fSkipLoadSolution = false;
  if (sim_data.mTNumerics.m_run_with_transport) {
    aspace.BuildAuxTransportCmesh();
    TPZCompMesh* transport_operator = aspace.GetTransportOperator();
    TMRSSFIAnalysis* sfi_analysis = new TMRSSFIAnalysis(mp_cmesh, transport_operator, renumtype);
    sfi_analysis->SetDataTransferAndBuildAlgDatStruct(&sim_data);
    sfi_analysis->Configure(glob_n_threads, UsePardiso_Q, UsingPzSparse);    
    const int n_steps = sim_data.mTNumerics.m_n_steps;
    const REAL dt = sim_data.mTNumerics.m_dt;

    // Times to report solution
    TPZStack<REAL, 100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    REAL sim_time = 0.0;
    int pos = 0;
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

    TPZFastCondensedElement::fSkipLoadSolution = false;
    const int typeToPPinit = 1;   // 0: both, 1: p/flux, 2: saturation
    const int typeToPPsteps = 2;  // 0: both, 1: p/flux, 2: saturation

    // Looping over time steps
    for (int it = 1; it <= n_steps; it++) {
      sim_time = it * dt;
      sfi_analysis->m_transport_module->SetCurrentTime(dt);
      sfi_analysis->RunTimeStep();
      if (it == 1) {
        sfi_analysis->PostProcessTimeStep(typeToPPinit);
      }
      mp_cmesh->LoadSolution(mp_cmesh->Solution());

      // Only post process based on reporting times
      if (sim_time >= current_report_time) {
        cout << "\n---------------------- SFI Step " << it << " ----------------------" << endl;
        std::cout << "Simulation time:  " << sim_time << std::endl;
        mp_cmesh->UpdatePreviousState(-1.);
        sfi_analysis->PostProcessTimeStep(typeToPPsteps);
        pos++;
        current_report_time = reporting_times[pos];

        REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        std::cout << "Mass report at time : " << sim_time << std::endl;
        std::cout << "Mass integral :  " << mass << std::endl;
      }
      sfi_analysis->m_transport_module->fAlgebraicTransport.VerifyConservation(it);
    }

  }
  else {
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    int dimToPost = 2;
    mixAnalisys->PostProcessTimeStep(dimToPost);
  }


  delete mp_cmesh;
  delete gmesh;
  cout << "-------------------- Simulation Finished --------------------" << endl;
  return 0;
}
// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh* ReadMeshFromGmsh(TMRSDataTransfer& sim_data) {
  // read mesh from gmsh
  TPZGeoMesh* gmesh;
  gmesh = new TPZGeoMesh();
  string file_name = std::string(FRACMESHES) + "/../Filling/" + sim_data.mTGeometry.mGmeshFileName;
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(5);
    stringtoint[2]["dom"] = sim_data.mTGeometry.mDomainNameAndMatId["dom"];

    stringtoint[1]["bcl"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bcl"];
    stringtoint[1]["bcr"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bcr"];
    stringtoint[1]["bct"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bct"];
    stringtoint[1]["bcb"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bcb"];
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransfer(string filenameBase, TMRSDataTransfer& sim_data) {
  using json = nlohmann::json;
  std::string filenamejson = filenameBase + ".json";

  std::ifstream filejson(filenamejson);
  json input = json::parse(filejson, nullptr, true, true);  // to ignore comments in json file

  // ------------------------ Getting number of domains and fractures ------------------------
  if (input.find("Domains") == input.end()) DebugStop();
  if (input.find("Mesh") == input.end()) DebugStop();
  const int ndom = input["Domains"].size();
  std:string mesh = input["Mesh"];
  sim_data.mTGeometry.mGmeshFileName = mesh;
  sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(ndom + 1);
  int countPhi = 0;

  // ------------------------ Reading 3D Domain matids ------------------------
  std::map<int, REAL> idPerm;
  for (auto& domain : input["Domains"]) {
    if (domain.find("matid") == domain.end()) DebugStop();
    if (domain.find("name") == domain.end()) DebugStop();
    if (domain.find("K") == domain.end()) DebugStop();
    if (domain.find("phi") == domain.end()) DebugStop();
    const int matid = domain["matid"];
    const string name = domain["name"];
    const REAL permeability = domain["K"];
    const REAL phi = domain["phi"];
    sim_data.mTGeometry.mDomainNameAndMatId[name] = matid;
    idPerm[matid] = permeability;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(matid, phi, 1.0);
  }

  sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;

  // ------------------------ Reading 3D Domain BC matids ------------------------
  if (input.find("Boundary") == input.end()) DebugStop();
  std::map<int, std::pair<int, REAL>>& BCFlowMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue;
  std::map<int, std::pair<int, REAL>>& BCTransportMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue;
  for (auto& bc : input["Boundary"]) {
    if (bc.find("matid") == bc.end()) DebugStop();
    if (bc.find("type") == bc.end()) DebugStop();
    if (bc.find("value") == bc.end()) DebugStop();
    if (bc.find("name") == bc.end()) DebugStop();
    const int matid = bc["matid"];
    const int type = bc["type"];
    const REAL value = bc["value"];
    const std::string name = bc["name"];
    sim_data.mTBoundaryConditions.mDomainNameAndMatId[name] = matid;
    if (BCFlowMatIdToTypeValue.find(matid) != BCFlowMatIdToTypeValue.end()) DebugStop();
    BCFlowMatIdToTypeValue[matid] = std::make_pair(type, value);
    BCTransportMatIdToTypeValue[matid] = std::make_pair(type, value);
  }

  // Transport properties
  if (input.find("RunWithTransport") != input.end()) {
    sim_data.mTNumerics.m_run_with_transport = input["RunWithTransport"];
    if (sim_data.mTNumerics.m_run_with_transport) {
      if (input.find("DeltaT") == input.end()) DebugStop();
      sim_data.mTNumerics.m_dt = input["DeltaT"];
      if (input.find("NSteps") == input.end()) DebugStop();
      sim_data.mTNumerics.m_n_steps = input["NSteps"];
    }
  } else {
    DebugStop();  // Please set run with tranport in json
                  //        sim_data.mTNumerics.m_n_steps = 100;
                  //        sim_data.mTNumerics.m_dt      = 1.e7;//*day;
  }

  // ------------------------ Setting extra stuff that is still not in JSON ------------------------
  const int D_Type = 0, N_Type = 1, Mixed_Type = 2;
  // sim_data.mTGeometry.mInterface_material_id = 100;
  // sim_data.mTGeometry.mInterface_material_idFracInf = 102;
  // sim_data.mTGeometry.mInterface_material_idFracSup = 101;
  // sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
  // sim_data.mTGeometry.mInterface_material_idFracBound = 104;

  // sim_data.mTGeometry.mSkeletonDiv = 0;
  sim_data.mTNumerics.m_sfi_tol = 0.0001;
  sim_data.mTNumerics.m_res_tol_transport = 0.0001;
  sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
  sim_data.mTNumerics.m_four_approx_spaces_Q = true;
  std::vector<REAL> grav(3, 0.0);
  grav[1] = 0.0;  //-9.8*(1.0e-6); // hor
  sim_data.mTNumerics.m_gravity = grav;
  sim_data.mTNumerics.m_ISLinearKrModelQ = true;
  sim_data.mTNumerics.m_nThreadsMixedProblem = glob_n_threads;
  sim_data.mTNumerics.m_max_iter_sfi = 1;
  sim_data.mTNumerics.m_max_iter_mixed = 1;
  sim_data.mTNumerics.m_max_iter_transport = 1;

  sim_data.mTPostProcess.m_file_name_mixed = "postdarcy.vtk";
  sim_data.mTPostProcess.m_file_name_transport = "posttransport.vtk";
  TPZStack<std::string, 10> scalnames, vecnames, scalnamesTransport;
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
  TPZStack<REAL, 100> reporting_times;
  REAL time = sim_data.mTPostProcess.m_file_time_step;
  int n_reporting_times = (n_steps) / (time / dt) + 1;
  REAL r_time = 0.0;
  for (int i = 1; i <= n_reporting_times; i++) {
    r_time += dt * (time / dt);
    reporting_times.push_back(r_time);
  }
  sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------