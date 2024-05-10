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

// ----- Global vars -----
const int glob_n_threads = 8;

// ----- End of namespaces -----

// ----- Functions -----
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

enum EMatid {ENone,EDomain,EFarField, EWellbore, EWellboreBC};

// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("imrs");
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
        sout << "\nLogger for WANN problem target\n" << endl;;
        LOGPZ_DEBUG(mainlogger, sout.str())
    }
#endif

    TPZGeoMesh* gmesh = ReadMeshFromGmsh(basemeshpath + "/wann/mesh_rev01.msh"); 
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

    
    return 0;
}
// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        TPZManVector<std::map<std::string,int>,4> stringtoint(13);
        stringtoint[2]["dom"] = EDomain;
        
        stringtoint[1]["far_field"] = EFarField;
        stringtoint[1]["wellbore"] = EWellbore;
        stringtoint[0]["heel"] = EWellboreBC;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);
    }

    return gmesh;
}
