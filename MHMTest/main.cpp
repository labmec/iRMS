
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
#include "TPZGenGrid2D.h"
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

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TPZExtendGridDimension.h"
#include "TPZFastCondensedElement.h"
#include "TPZReservoirTools.h"
#include "pzcondensedcompel.h"

#include "TPZAlgebraicDataTransfer.h"

#include <libInterpolate/Interpolate.hpp>
#include <libInterpolate/AnyInterpolator.hpp>

TMRSDataTransfer Setting2D();
TPZFMatrix<STATE> TimeForward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag);

void SimpleTest();
//
int main(){
    InitializePZLOG();
    SimpleTest();
    return 0;
}


void SimpleTest(){
    
    TMRSDataTransfer sim_data  = Setting2D();
    
    std::string geometry_file = "gmsh/simple_2D_coarse.msh";
    std::string name = "simplemesh.vtk";
    
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);
    aspace.CreateUniformMesh(10, 100, 1, 10);
    aspace.GenerateMHMUniformMesh(0);
    aspace.PrintGeometry(name);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    
    {
        std::ofstream out("fluxmesh.txt");
        mixed_operator->MeshVector()[0]->Print(out);
    }
//    TPZAlgebraicTransport transport;
//    TPZAlgebraicDataTransfer transfer;
//    transfer.SetMeshes(*mixed_operator, *transport_operator);
//    transfer.BuildTransportDataStructure(transport);
//    
   
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
  
//    sfi_analysis->RunTimeStep();
    
    
//    exit(0);
//    sfi_analysis->m_mixed_module->RunTimeStep();
  

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    
    
    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    
    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];
    TPZFMatrix<REAL> solution_n;
    solution_n = sfi_analysis->m_transport_module->Solution();
    
//    solution_n *= 0.0;
   
    for (int it = 1; it <= n_steps; it++) {
        sim_time = it*dt;
        sfi_analysis->m_transport_module->SetCurrentTime(sim_time);
        sfi_analysis->RunTimeStep();
       
        if (sim_time >=  current_report_time) {
            std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
            sfi_analysis->PostProcessTimeStep();
            pos++;
            current_report_time =reporting_times[pos];
        }
    }
}


TMRSDataTransfer Setting2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.Interface_material_id = 100;
   
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 1100.0;
    REAL pressure_out = 100.0;
    
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(3);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(2,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(4,D_Type,pressure_in);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(5,D_Type,pressure_out);
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
      sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,D_Type,pressure_out);
 
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-3,bc_outlet,0.0);
     sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_outlet,1.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.2;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 800.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 0.00001;
    sim_data.mTNumerics.m_corr_tol_mixed = 0.000001;
    sim_data.mTNumerics.m_res_tol_transport = 0.00001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00001;
    sim_data.mTNumerics.m_n_steps = 100;
    sim_data.mTNumerics.m_dt      = 0.1;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
//    scalnames.Push("p");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 0.1;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
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


