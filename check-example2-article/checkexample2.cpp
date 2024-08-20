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
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZVTKGenerator.h"

// ----- Namespaces -----
using namespace std;
namespace fs = std::filesystem;

// ----- End of namespaces -----

// ----- Global vars -----
const int glob_n_threads = 8;
enum EMatId {Enone, Edom1, Edom2, EbcNoFlux, Ebcin, Ebcout, Efrac};

// ----- Functions -----
TPZGeoMesh* ReadMeshFromGmsh(std::string filename);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh);
void SolveSyst(TPZLinearAnalysis &an, TPZCompMesh *cmesh, DecomposeType dtype);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

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

  // =========> Create GeoMesh
  TPZGeoMesh* gmesh = ReadMeshFromGmsh("FractureExample.msh");
  // Print gmesh to vtk
  std::ofstream out("gmesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // =========> Create H1 comp mesh
  TPZCompMesh* cmesh = CreateH1CMesh(gmesh);
    
  // =========> Create Analysis and solve
  TPZLinearAnalysis an(cmesh);
  DecomposeType dtype = ECholesky;
  SolveSyst(an, cmesh, dtype);
  
  // =========> Print results to vtk
  PrintResults(an, cmesh);
  
  delete cmesh;
  delete gmesh;
  cout << "-------------------- Simulation Finished --------------------" << endl;
  return 0;
}
// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh* ReadMeshFromGmsh(std::string filename) {
  // read mesh from gmsh
  TPZGeoMesh* gmesh;
  gmesh = new TPZGeoMesh();
  string file_name = std::string(FRACMESHES) + "/../check-example2-article/" + filename;
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(7);
    stringtoint[2]["Domain1"] = Edom1;
    stringtoint[2]["Domain2"] = Edom2;

    stringtoint[1]["FlowOut"] = Ebcout;
    stringtoint[1]["FlowIn"] = Ebcin;
    stringtoint[1]["NoFlow"] = EbcNoFlux;
    stringtoint[1]["Fracture1"] = Efrac;
    stringtoint[1]["Fracture2"] = Efrac;
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh) {
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
  const int dim = gmesh->Dimension();
  cmesh->SetDefaultOrder(2);
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsContinuous();

  if (dim != 2) DebugStop(); // only testing 2d problem here
  
  // Create domain materials
  const REAL K1 = 1.0, K2 = 1.0e-6;
  TPZDarcyFlow* matdom1 = new TPZDarcyFlow(Edom1, 2);
  matdom1->SetConstantPermeability(K1);
  TPZDarcyFlow* matdom2 = new TPZDarcyFlow(Edom2, 2);
  matdom2->SetConstantPermeability(K2);
  
  cmesh->InsertMaterialObject(matdom1);
  cmesh->InsertMaterialObject(matdom2);
  
  // Create boundary conditions
  TPZFMatrix<STATE> val1(1,1,0.);
  TPZManVector<STATE> val2(1,1.);
  
  const int diri = 0, neu = 1, mixed = 2;
  auto bndIn = matdom1->CreateBC(matdom1, Ebcin, diri, val1, val2);
  cmesh->InsertMaterialObject(bndIn);
  
  val2[0] = 0.;
  auto bndOut = matdom1->CreateBC(matdom1, Ebcout, diri, val1, val2);
  cmesh->InsertMaterialObject(bndOut);
  
  auto bndNoFlux = matdom1->CreateBC(matdom1, EbcNoFlux, neu, val1, val2);
  cmesh->InsertMaterialObject(bndNoFlux);
  
  // Create fracture materials
  const REAL Kfrac = 1.0e3;
  TPZDarcyFlow* matfrac = new TPZDarcyFlow(Efrac, 1);
  matfrac->SetConstantPermeability(Kfrac);
  cmesh->InsertMaterialObject(matfrac);
  
  // Create comp elements
  cmesh->AutoBuild();
  
  return cmesh;
}


// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh) {
 
    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "post";
    constexpr int vtkRes{0};
    TPZVec<std::string> fields = {
        "Pressure",
        "Flux"
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(0);
    vtk.Do();
  
    // Post proc frac
    const std::string plotfile2 = "post_frac";
    auto vtkfrac = TPZVTKGenerator(cmesh, fields, plotfile2, vtkRes, cmesh->Dimension()-1);
    vtkfrac.SetNThreads(0);
    vtkfrac.Do();
  
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void SolveSyst(TPZLinearAnalysis &an, TPZCompMesh *cmesh, DecomposeType dtype) {
    TPZSSpStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;

    step.SetDirect(dtype);
    an.SetSolver(step);
    an.Run();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
