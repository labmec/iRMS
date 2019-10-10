
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "pzgengrid.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pznoderep.h"
#include "pzshapetetra.h"
#include "RSimulatorConfiguration.h"
#include "pzshapepiram.h"


// Geometry generation
#include "../TMRSApproxSpaceGenerator.h"

TPZGeoMesh * ReadGeometry(std::string geometry_file);
void PrintGeometry(TPZGeoMesh * geometry, std::string name);

int main(){
 
    std::string geometry_file = "reservoir.msh";
    std::string name = "reservoir";
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);
    aspace.PrintGeometry(name);
    
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 0;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(2);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(0.1);
    
    sim.c_inlet = 1.0;
    
    int bc_non_flux = -1;
    int bc_inlet  = -2;
    int bc_non_flux2 = -3;
    int bc_outlet = -4;
    
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_non_flux2);
    sim.gamma_dim.push_back(1);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(1);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 1.0;
    REAL p_outlet = 0.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    sim.type.push_back(bc_type_D);
    
    sim.vals.push_back(qn);
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(qn);
    sim.vals.push_back(p_outlet);
    
    
    RSimulatorConfiguration caseSim(sim);
    caseSim.Run();
    
    return 0;
}


