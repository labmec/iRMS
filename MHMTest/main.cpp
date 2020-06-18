
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
#include "TMRSPropertiesFunctions.h"
#include <libInterpolate/Interpolate.hpp>
#include <libInterpolate/AnyInterpolator.hpp>


TMRSDataTransfer Setting2D();
TMRSDataTransfer Setting3D();
TMRSDataTransfer SettingUNISIM();
TPZGeoMesh * CreateGeoMeshWithTopeAndBase(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ);
void ModifyTopeAndBase(TPZGeoMesh * gmesh, std::string filename);
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers);
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
void SimpleTest();
void SimpleTest3D();
void UNISIMTest();
void PostProcessResProps(TPZMultiphysicsCompMesh *cmesh, TPZAlgebraicTransport *alg);

//
int main(){
    InitializePZLOG();
    SimpleTest();
//    SimpleTest3D();
//    UNISIMTest();
    return 0;
}


void SimpleTest(){
    
    TMRSDataTransfer sim_data  = Setting2D();
    
    std::string geometry_file = "gmsh/simple_2D_coarse.msh";
    std::string name = "simplemesh.vtk";
    
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);

    aspace.CreateUniformMesh(10, 10, 10, 10);
    aspace.GenerateMHMUniformMesh(2);

    aspace.PrintGeometry(name);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh();
    
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
   
    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_kappa(TMRSPropertiesFunctions::EPiecewiseFunction);
    reservoir_properties.set_function_type_phi(TMRSPropertiesFunctions::EConstantFunction);
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::ECircleLevelSetFunction);

    auto kx = reservoir_properties.Create_Kx();
    auto ky = reservoir_properties.Create_Ky();
    auto kz = reservoir_properties.Create_Kz();
    auto phi = reservoir_properties.Create_phi();
    auto s0 = reservoir_properties.Create_s0();
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kx,ky,kz,phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
  
    // Render a graphical map
    PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    
    
    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    
    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];
    
    // Print initial condition
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    sfi_analysis->PostProcessTimeStep();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;
    for (int it = 1; it <= n_steps; it++) {
        sim_time = it*dt;
        sfi_analysis->m_transport_module->SetCurrentTime(dt);
        sfi_analysis->RunTimeStep();
        
       
        if (sim_time >=  current_report_time) {
            std::cout << "Time step number:  " << it << std::endl;
            std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
            sfi_analysis->PostProcessTimeStep();
            pos++;
            current_report_time =reporting_times[pos];
            
            REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
            std::cout << "Mass report at time : " << sim_time << std::endl;
            std::cout << "Mass integral :  " << mass << std::endl;
            
        }
        sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
    }
    
    std::cout  << "Number of transportr equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
}
void SimpleTest3D(){
    
    TMRSDataTransfer sim_data  = Setting3D();
    
    std::string geometry_file = "gmsh/simple_2D_coarse.msh";
    std::string name = "simplemesh.vtk";
    
    TMRSApproxSpaceGenerator aspace;
//    aspace.LoadGeometry(geometry_file);
    
    aspace.CreateUniformMesh(2, 100, 1, 10,1,10.0);
    aspace.GenerateMHMUniformMesh(2);
    
    aspace.PrintGeometry(name);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
      bool must_opt_band_width_Q = true;
      int n_threads = 0;
      bool UsePardiso_Q = true;
      aspace.BuildMixedMultiPhysicsCompMesh(order);
      TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
      TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh();
      
      aspace.BuildTransportMultiPhysicsCompMesh();
      TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
     
      TMRSPropertiesFunctions reservoir_properties;
      reservoir_properties.set_function_type_kappa(TMRSPropertiesFunctions::EPiecewiseFunction);
      reservoir_properties.set_function_type_phi(TMRSPropertiesFunctions::EConstantFunction);
      reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::ECircleLevelSetFunction);

      auto kx = reservoir_properties.Create_Kx();
      auto ky = reservoir_properties.Create_Ky();
      auto kz = reservoir_properties.Create_Kz();
      auto phi = reservoir_properties.Create_phi();
      auto s0 = reservoir_properties.Create_s0();
      
      TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kx,ky,kz,phi,s0);
      sfi_analysis->SetDataTransfer(&sim_data);
      sfi_analysis->Configure(n_threads, UsePardiso_Q);
    
      // Render a graphical map
      PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
      
      int n_steps = sim_data.mTNumerics.m_n_steps;
      REAL dt = sim_data.mTNumerics.m_dt;
      
      
      TPZStack<REAL,100> reporting_times;
      reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
      
      REAL sim_time = 0.0;
      int pos =0;
      REAL current_report_time = reporting_times[pos];
      
      // Print initial condition
      sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
      sfi_analysis->PostProcessTimeStep();
      REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
      std::cout << "Mass report at time : " << 0.0 << std::endl;
      std::cout << "Mass integral :  " << initial_mass << std::endl;
      for (int it = 1; it <= n_steps; it++) {
          sim_time = it*dt;
          sfi_analysis->m_transport_module->SetCurrentTime(dt);
          sfi_analysis->RunTimeStep();
          
         
          if (sim_time >=  current_report_time) {
              std::cout << "Time step number:  " << it << std::endl;
              std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
              sfi_analysis->PostProcessTimeStep();
              pos++;
              current_report_time =reporting_times[pos];
              
              REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
              std::cout << "Mass report at time : " << sim_time << std::endl;
              std::cout << "Mass integral :  " << mass << std::endl;
              
          }
          sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
      }
      
      std::cout  << "Number of transportr equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
}
void UNISIMTest(){
    std::string geometry_file2D ="gmsh/Contorno.msh";
    int nLayers = 2;
    bool is3DQ = true;
    bool print3DMesh = true;
    TPZGeoMesh *gmesh = CreateGeoMeshWithTopeAndBase( geometry_file2D,  nLayers, print3DMesh, is3DQ);
 
    TMRSApproxSpaceGenerator aspace;
    
    TMRSDataTransfer sim_data  = SettingUNISIM();
    
    aspace.SetGeometry(gmesh);
    std::string name="NewMesh";
    std::cout<< gmesh->NElements();
    aspace.GenerateMHMUniformMesh(0);
    aspace.PrintGeometry(name);
    
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh();

    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();

    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_kappa(TMRSPropertiesFunctions::EPiecewiseFunction);
    reservoir_properties.set_function_type_phi(TMRSPropertiesFunctions::EConstantFunction);
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::ECircleLevelSetFunction);

    auto kx = reservoir_properties.Create_Kx();
    auto ky = reservoir_properties.Create_Ky();
    auto kz = reservoir_properties.Create_Kz();
    auto phi = reservoir_properties.Create_phi();
    auto s0 = reservoir_properties.Create_s0();

    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kx,ky,kz,phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);

    // Render a graphical map
    PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;


    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];

    // Print initial condition
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    sfi_analysis->PostProcessTimeStep();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;
    for (int it = 1; it <= n_steps; it++) {
      sim_time = it*dt;
      sfi_analysis->m_transport_module->SetCurrentTime(dt);
      sfi_analysis->RunTimeStep();
      
     
      if (sim_time >=  current_report_time) {
          std::cout << "Time step number:  " << it << std::endl;
          std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
          sfi_analysis->PostProcessTimeStep();
          pos++;
          current_report_time =reporting_times[pos];
          
          REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
          std::cout << "Mass report at time : " << sim_time << std::endl;
          std::cout << "Mass integral :  " << mass << std::endl;
          
      }
      sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
    }

    std::cout  << "Number of transportr equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
    
}

TMRSDataTransfer Setting2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.Interface_material_id = 100;
   
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 10.0;
    REAL pressure_out = 10.0;
    
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(3);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(2,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(4,D_Type,pressure_in);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(5,D_Type,pressure_out);
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
//      sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,D_Type,pressure_out);
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,N_Type,zero_flux);
      sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,D_Type,pressure_out);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,N_Type,zero_flux);
 
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
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.001;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 500.0;

    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 30;
    sim_data.mTNumerics.m_sfi_tol = 0.001;
    sim_data.mTNumerics.m_res_tol_mixed = 0.00001;
    sim_data.mTNumerics.m_corr_tol_mixed = 0.00001;
    sim_data.mTNumerics.m_res_tol_transport = 0.00001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00001;
    sim_data.mTNumerics.m_n_steps = 100;
    REAL day = 86400;
    sim_data.mTNumerics.m_dt      = 0.01*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = -9.81*(1.0e-6);
    sim_data.mTNumerics.m_gravity = grav;
    
    
    
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
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
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
TMRSDataTransfer Setting3D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["RockMatrix"] = 1;
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
    sim_data.mTFluidProperties.mWaterViscosity = 1.0;
    sim_data.mTFluidProperties.mOilViscosity = 1.0;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 800.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 30;
    sim_data.mTNumerics.m_res_tol_mixed = 0.000001;
    sim_data.mTNumerics.m_corr_tol_mixed = 0.000001;
    sim_data.mTNumerics.m_res_tol_transport = 0.000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.000001;
    sim_data.mTNumerics.m_n_steps = 30;
    sim_data.mTNumerics.m_dt      = 10.0;
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
    sim_data.mTPostProcess.m_file_time_step = 10.0;
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
TMRSDataTransfer SettingUNISIM(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["RockMatrix"] = 1;
//    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["RockMatrix"] = 1;
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
    sim_data.mTFluidProperties.mWaterViscosity = 1.0;
    sim_data.mTFluidProperties.mOilViscosity = 1.0;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 800.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 30;
    sim_data.mTNumerics.m_res_tol_mixed = 0.000001;
    sim_data.mTNumerics.m_corr_tol_mixed = 0.000001;
    sim_data.mTNumerics.m_res_tol_transport = 0.000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.000001;
    sim_data.mTNumerics.m_n_steps = 50;
    sim_data.mTNumerics.m_dt      = 1000.0;
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
    sim_data.mTPostProcess.m_file_time_step = 1000.0;
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
TPZGeoMesh * CreateGeoMeshWithTopeAndBase(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ){
    
    //    std::string  geometry_file = "gmsh/Contorno.msh";
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["RockMatrix"] = 1;
    dim_name_and_physical_tag[2]["RockMatrix2"] = 1;
    dim_name_and_physical_tag[1]["Injectors"] = -2;
    dim_name_and_physical_tag[1]["Productors"] = -4;
    dim_name_and_physical_tag[1]["ZeroFlux"] = -1;
    
    TPZGmshReader Geometry;
    TPZGeoMesh * gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);
    double w = 200.0;
    
    std::string name2D("mesh2d.vtk");
    
    
    int topID= -3;
    int baseID = -1;
    TPZGeoMesh * returnedMesh = nullptr;
    if (Is3DQ) {
        TPZExtendGridDimension extend(gmesh2d, w);
        returnedMesh = extend.ExtendedMesh(nLayers,topID,baseID);
        ModifyTopeAndBase2(returnedMesh ,nLayers);
        
    }
    if (!Is3DQ){
        std::string filename1 = "Reservoir/tope_unisim2.txt";
        ModifyTopeAndBase(gmesh2d,filename1 );
        returnedMesh =gmesh2d;
        
    }
    return returnedMesh;
    
}
void ModifyTopeAndBase(TPZGeoMesh * gmesh, std::string filename){
    
    
    std::vector<double> x, y, z;
    ReadData(filename, true, x, y, z);
    
    _2D::ThinPlateSplineInterpolator <double> interp;
    interp.setData(x,y,z);
    
    int nCoordinates = gmesh->NodeVec().NElements();
    double sum=0.0;
    for (auto val:z) {
        sum += val;
    }
    double val_storage= sum / z.size();
    
    for (int icoord=0; icoord< nCoordinates; icoord++) {
        TPZGeoNode node = gmesh->NodeVec()[icoord];
        TPZVec<REAL> co(3);
        node.GetCoordinates(co);
        double val_interp =interp(co[0],co[1]);
        
        if (val_interp==0.0) {
            co[2] =val_storage;
        }
        if (val_interp>1.0) {
            co[2] =val_interp;
        }
        
        gmesh->NodeVec()[icoord].SetCoord(co);
    }
    
    
}
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers){
    std::string filename1 = "Reservoir/tope_unisim2.txt";
    std::string filename2 = "Reservoir/base_unisim2.txt";
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
            }
            if (baseinterpol==0) {
                baseinterpol = val_base;
            }

            if (ilay==1) {
                valinter=topeinterpol;
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
            if (ilay==nlayers+1) {
                valinter = baseinterpol;
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
            if (ilay>1   && ilay < nlayers+1) {
                valinter = topeinterpol + (ilay-1)*(baseinterpol - topeinterpol)/(nlayers);
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
//            if(co[2]==0.0){
//                DebugStop();
//            }
//

        }
    }
}
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z){
    
    bool modpoints = true;
    std::ifstream file;
    file.open(name);
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
                    x.push_back(a - 311214.135108);
                    y.push_back(b - 7478429.28008);
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

void PostProcessResProps(TPZMultiphysicsCompMesh *cmesh, TPZAlgebraicTransport *alTransport){
    
    
     TPZAnalysis *an = new TPZAnalysis(cmesh,false);
     int ncells = alTransport->fCellsData.fVolume.size();
     TPZFMatrix<STATE> phi(ncells,1,0.1);
     TPZFMatrix<STATE> Kx(ncells,1,1.0);
     TPZFMatrix<STATE> Ky(ncells,1,2.0);
     TPZFMatrix<STATE> Kz(ncells,1,3.0);
    
    for (int ivol = 0; ivol<ncells; ivol++) {
        REAL phival = alTransport->fCellsData.fporosity[ivol];
        REAL Kxval = alTransport->fCellsData.fKx[ivol];
        REAL Kyval = alTransport->fCellsData.fKy[ivol];
        REAL Kzval = alTransport->fCellsData.fKz[ivol];
        phi(ivol,0) =phival;
        Kx(ivol,0) =Kxval;
        Ky(ivol,0) =Kyval;
        Kz(ivol,0) =Kzval;
    }
//    cmesh->MeshVector()[0]->Solution().Resize(ncells, 1);
//    cmesh->MeshVector()[1]->Solution().Resize(ncells, 1);
//    cmesh->MeshVector()[2]->Solution().Resize(ncells, 1);
//    cmesh->MeshVector()[3]->Solution().Resize(ncells, 1);
     cmesh->MeshVector()[0]->Solution() = phi;
     cmesh->MeshVector()[1]->Solution() = Kx;
     cmesh->MeshVector()[2]->Solution() = Ky;
     cmesh->MeshVector()[3]->Solution() = Kz;
    
     TPZStack<std::string,10> scalnames, vecnames;
     vecnames.Push("q");
     scalnames.Push("Porosity");
     scalnames.Push("Permeability_x");
     scalnames.Push("Permeability_y");
     scalnames.Push("Permeability_z");
     std::string file("Props.vtk");
     an->DefineGraphMesh(2,scalnames,vecnames,file);
     an->PostProcess(0,2);
    
}
