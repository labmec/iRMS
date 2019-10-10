
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
#include "../TMRSDataTransfer.h"

TPZAnalysis * CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,  bool must_opt_band_width_Q, int n_threads, bool UsePardiso_Q);

int main(){
 
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_wbregion_p"] = 2;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_wbregion_i"] = 3;
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(6);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(4,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(5,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(6,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(7,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[4] = std::make_tuple(10,0,20.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[5] = std::make_tuple(11,0,30.0);
    
    std::string geometry_file = "reservoir.msh";
    std::string name = "reservoir";
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);
    aspace.PrintGeometry(name);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    TPZMultiphysicsCompMesh * mixed_operator = aspace.MixedMultiPhysicsCompMesh(order);
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    TPZAnalysis * analysis =  CreateAnalysis(mixed_operator,  must_opt_band_width_Q, n_threads, UsePardiso_Q);
    
    analysis->Assemble();
    analysis->Solve();
    mixed_operator->LoadSolutionFromMultiPhysics();
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    scalnames.Push("p");
    
    int div = 0;
    int dim = aspace.GetGeometry()->Dimension();
    std::string fileresult("flux_and_p.vtk");
    analysis->DefineGraphMesh(dim,scalnames,vecnames,fileresult);
    analysis->PostProcess(div,dim);
    
    return 0;
    
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


TPZAnalysis * CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,  bool must_opt_band_width_Q, int n_threads, bool UsePardiso_Q){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh_mult, must_opt_band_width_Q);
    
    if (UsePardiso_Q) {
        
        TPZSymetricSpStructMatrix matrix(cmesh_mult);
        matrix.SetNumThreads(n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);

    }else{
        
        TPZSkylineStructMatrix matrix(cmesh_mult);
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);

    }
    
    return analysis;
    
}
