
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
#include "TRMSpatialPropertiesMap.h"
#include <libInterpolate/Interpolate.hpp>
#include <libInterpolate/AnyInterpolator.hpp>

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

TMRSDataTransfer SettingPlumeGravity2D();
TMRSDataTransfer SettingGravity2D();
TMRSDataTransfer Setting3D();
TMRSDataTransfer SettingUNISIM();

TMRSDataTransfer SettingPaper2D();
TMRSDataTransfer SettingPaper3D();

TPZGeoMesh * CreateGeoMeshWithTopeAndBase(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ);
TPZGeoMesh * CreateGeoMeshMHM3DTest(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ);

void ModifyTopeAndBase(TPZGeoMesh * gmesh, std::string filename);
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers);
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);

void PlumeGravity2D();
void Gravity2D();
void SimpleTest3D();
void UNISIMTest();

void PaperTest2D();
void PaperTest3D();

void PostProcessResProps(TPZMultiphysicsCompMesh *cmesh, TPZAlgebraicTransport *alg);

TMRSDataTransfer SettingSimple2D();
void SimpleTest2D();
//
int main(){
    InitializePZLOG();
    
    // Gravity cases.
//    PlumeGravity2D();
    Gravity2D();
    
    // Refinement effect examples
//    PaperTest2D();
//    PaperTest3D();

    // Box meshes examples
//    SimpleTest2D();
//    SimpleTest3D();
    
    // Complex geometry example
//    UNISIMTest();
    return 0;
}
void SimpleTest2D(){
    TMRSDataTransfer sim_data  = SettingSimple2D();
    
    TMRSApproxSpaceGenerator aspace;
    aspace.CreateUniformMesh(10, 10, 10, 10);
    std::string name = "2D_geo";
    aspace.PrintGeometry(name);
    
    aspace.ApplyUniformRefinement(2);
    std::string name_ref = "2D_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    
    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_kappa(TMRSPropertiesFunctions::EConstantFunction);
    reservoir_properties.set_function_type_phi(TMRSPropertiesFunctions::EConstantFunction);
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::EConstantFunction);
    
    auto kx = reservoir_properties.Create_Kx();
    auto ky = reservoir_properties.Create_Ky();
    auto kz = reservoir_properties.Create_Kz();
    auto phi = reservoir_properties.Create_phi();
    auto s0 = reservoir_properties.Create_s0();
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kx,ky,kz,phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
    
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
    sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
    sfi_analysis->PostProcessTimeStep();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;
    
    for (int it = 1; it <= n_steps; it++) {
        TPZFastCondensedElement::fSkipLoadSolution = false;
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
        TPZFastCondensedElement::fSkipLoadSolution = true;
    }
    
    std::cout  << "Number of transport equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
    
}

void PlumeGravity2D(){
    
  TMRSDataTransfer sim_data  = SettingPlumeGravity2D();
      
  TMRSApproxSpaceGenerator aspace;
  aspace.CreateUniformMesh(1, 1, 1280, 10);
  std::string name = "g_segregation_geo";
  aspace.PrintGeometry(name);
  
  aspace.ApplyUniformRefinement(0);
  std::string name_ref = "g_segregation_ref_geo";
  aspace.PrintGeometry(name_ref);
  aspace.SetDataTransfer(sim_data);
  
  int order = 1;
  bool must_opt_band_width_Q = true;
  int n_threads = 0;
  bool UsePardiso_Q = true;
  aspace.BuildMixedMultiPhysicsCompMesh(order);
  TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
  
  aspace.BuildTransportMultiPhysicsCompMesh();
  TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
 
  TMRSPropertiesFunctions reservoir_properties;
  reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::EPiecewiseFunction);

  auto kappa_phi = reservoir_properties.Create_Kappa_Phi();
  auto s0 = reservoir_properties.Create_s0();
  
  TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kappa_phi,s0);
  sfi_analysis->SetDataTransfer(&sim_data);
  sfi_analysis->Configure(n_threads, UsePardiso_Q);

  TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);

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
  sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
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
  }
  
  std::cout  << "Number of transport equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
    
}

void Gravity2D(){
    
    #ifdef USING_BOOST
        boost::posix_time::ptime toverhead1 = boost::posix_time::microsec_clock::local_time();
    #endif
    
    TMRSDataTransfer sim_data  = SettingGravity2D();
    
    TMRSApproxSpaceGenerator aspace;
    aspace.CreateUniformMesh(10, 10, 10, 10);
    std::string name = "g_segregation_geo";
    aspace.PrintGeometry(name);
    
    aspace.ApplyUniformRefinement(3);
    std::string name_ref = "g_segregation_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
   
    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_kappa(TMRSPropertiesFunctions::EConstantFunction);
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
  
    // Common section
            
    // Render a graphical map
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
    PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
    std::cout << "Spatial properties are transferred." << std::endl;

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;

    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos = 0;
    REAL current_report_time = reporting_times[pos];

    // Mass integral - Injection - Production data
    TPZFNMatrix<200,REAL> time_mass(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_inj(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_prod(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_fluxInlet(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_fluxOutlet(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_sfi_interations(n_steps+1,2,0.0);
    // Print initial condition
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-2);
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-4);
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
    sfi_analysis->PostProcessTimeStep();
    
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
    std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
    std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);
    REAL fluxInlet_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-2);
    REAL fluxOutlet_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-4);
    
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;

    time_mass(0,0) = 0.0;
    time_mass(0,1) = initial_mass;

    time_inj(0,0) = 0.0;
    time_inj(0,1) = inj_data.first;
    time_inj(0,2) = inj_data.second;

    time_prod(0,0) = 0.0;
    time_prod(0,1) = prod_data.first;
    time_prod(0,2) = prod_data.second;
    
    time_fluxInlet(0,0) =0.0;
    time_fluxInlet(0,1) =fluxInlet_data;

    time_fluxOutlet(0,0) =0.0;
    time_fluxOutlet(0,1) =fluxOutlet_data;
    
    time_sfi_interations(0,0) = 0;
    time_sfi_interations(0,1) = 0;
    
#ifdef USING_BOOST
    boost::posix_time::ptime toverhead2 = boost::posix_time::microsec_clock::local_time();
    auto delta_overhead = toverhead2-toverhead1;
    std::cout << "Overhead time " << delta_overhead << std::endl;;
#endif

#ifdef USING_BOOST
    boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
#endif

    for (int it = 1; it <= n_steps; it++) {
        
        sim_time = it*dt;
        if (sim_time >= current_report_time)
        {
            TPZFastCondensedElement::fSkipLoadSolution = false;
        }
        sfi_analysis->m_transport_module->SetCurrentTime(dt);
        
        #ifdef USING_BOOST
            boost::posix_time::ptime tsfi_step1 = boost::posix_time::microsec_clock::local_time();
        #endif
        
        sfi_analysis->RunTimeStep();
    
        #ifdef USING_BOOST
            boost::posix_time::ptime tsfi_step2 = boost::posix_time::microsec_clock::local_time();
            auto delta_sfi_step = tsfi_step2-tsfi_step1;
            std::cout << "SFI step time " << delta_sfi_step << std::endl;;
        #endif
        
        if (sim_time >=  current_report_time) {
          std::cout << "Time step number:  " << it << std::endl;
          std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
          sfi_analysis->m_mixed_module->LoadSolution();
          TPZMultiphysicsCompMesh *mphys = dynamic_cast<TPZMultiphysicsCompMesh *>(sfi_analysis->m_mixed_module->Mesh());
          if(!mphys) DebugStop();
            
          mphys->LoadSolutionFromMultiPhysics();
            
          sfi_analysis->PostProcessTimeStep();
          pos++;
          current_report_time = reporting_times[pos];
        }
        
        REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
        std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
        std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);

        REAL inletflow = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-2);
        REAL outletflow = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-4);
        std::cout << "Mass report at time : " << sim_time << std::endl;
        std::cout << "Mass integral :  " << mass << std::endl;

        time_mass(it,0) = sim_time;
        time_mass(it,1) = mass;

        time_inj(it,0) = sim_time;
        time_inj(it,1) = inj_data.first;
        time_inj(it,2) = inj_data.second;

        time_prod(it,0) = sim_time;
        time_prod(it,1) = prod_data.first;
        time_prod(it,2) = prod_data.second;

        time_fluxInlet(it,0)=sim_time;
        time_fluxInlet(it,1) =inletflow;

        time_fluxOutlet(it,0)=sim_time;
        time_fluxOutlet(it,1) =outletflow;
        
        time_sfi_interations(it,0) = sim_time;
        time_sfi_interations(it,1) = sfi_analysis->m_k_iteration;
    }

#ifdef USING_BOOST
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = tsim2-tsim1;
    std::cout << "Overall SFI simulation time " << deltat << std::endl;;
#endif


    // Writing relevant output

    std::ofstream outdata_mass("imrs_mass.txt");
    time_mass.Print("mass = ",outdata_mass,EMathematicaInput);

    std::ofstream outdata_inj("imrs_inj.txt");
    time_inj.Print("inj = ",outdata_inj,EMathematicaInput);

    std::ofstream outdata_prod("imrs_prod.txt");
    time_prod.Print("prod = ",outdata_prod,EMathematicaInput);
    
    std::ofstream outdata_fluxin("imrs_fluxinlet.txt");
    time_fluxInlet.Print("fluxint = ",outdata_fluxin,EMathematicaInput);
    
    std::ofstream outdata_fluxout("imrs_fluxoutlet.txt");
    time_fluxOutlet.Print("fluxint = ",outdata_fluxout,EMathematicaInput);
    
    std::ofstream outdata_iterations("imrs_iterations.txt");
    time_sfi_interations.Print("iters = ",outdata_iterations,EMathematicaInput);
}

void PaperTest2D(){
    
    #ifdef USING_BOOST
        boost::posix_time::ptime tspatialmap1 = boost::posix_time::microsec_clock::local_time();
    #endif
    
    // total 60-x-220-x-85
    std::vector<size_t> n_blocks = {220,60,1};// first layer
    std::vector<REAL> size_blocks = {1000.0/220.0,100.0/60.0,1.0};
    std::vector<REAL> translation = {-500.0,-50.0,-0.5};
    std::vector<size_t> SAMe_blocks = {5,5,1}; // keep it small
    std::string perm_data = "maps/spe_perm.dat";
    std::string phi_data  = "maps/spe_phi.dat";
    TRMSpatialPropertiesMap properties_map;
    properties_map.SetCartesianMeshData(n_blocks,size_blocks,perm_data,phi_data,SAMe_blocks,translation);
    
    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::EConstantFunction);
    
    auto kappa_phi = reservoir_properties.Create_Kappa_Phi(properties_map);
    auto s0 = reservoir_properties.Create_s0();
   
    #ifdef USING_BOOST
        boost::posix_time::ptime tspatialmap2 = boost::posix_time::microsec_clock::local_time();
        auto delta_spatialmap = tspatialmap2-tspatialmap1;
        std::cout << "Spatial map construction time " << delta_spatialmap << std::endl;;
    #endif

    #ifdef USING_BOOST
        boost::posix_time::ptime toverhead1 = boost::posix_time::microsec_clock::local_time();
    #endif
    
    TMRSDataTransfer sim_data  = SettingPaper2D();
    std::string geometry_file = "gmsh/reservoir_2d_paper.msh";
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["Reservoir"] = 1;
    dim_name_and_physical_tag[2]["Wellbore_p"] = 1;
    dim_name_and_physical_tag[2]["Wellbore_i"] = 1;
    dim_name_and_physical_tag[1]["producers"] = -4;
    dim_name_and_physical_tag[1]["injectors"] = -2;
    dim_name_and_physical_tag[1]["Reservoir_south"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_West"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_north"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_East"] = -1;
    
    TPZGmshReader Geometry;
    TPZGeoMesh * gmesh;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(geometry_file);
    Geometry.PrintPartitionSummary(std::cout);
    
    TMRSApproxSpaceGenerator aspace;
    aspace.SetGeometry(gmesh);
    std::string name = "paper_2d_test_geo";
    aspace.PrintGeometry(name);
    
    aspace.ApplyUniformRefinement(3);
    std::string name_ref = "paper_2d_test_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();

    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();

    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kappa_phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    
    // Common section
        
    // Render a graphical map
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
    PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
    std::cout << "Spatial properties are transferred." << std::endl;
    std::cout << "Memory used by SpatialPropertiesMap is released. " << std::endl;
    properties_map.Clear();

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;

    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos = 0;
    REAL current_report_time = reporting_times[pos];

    // Mass integral - Injection - Production data
    TPZFNMatrix<200,REAL> time_mass(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_inj(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_prod(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_fluxInlet(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_fluxOutlet(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_sfi_interations(n_steps+1,2,0.0);
    // Print initial condition
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-2);
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-4);
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
    sfi_analysis->PostProcessTimeStep();
    
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
    std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
    std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);
    REAL fluxInlet_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-2);
    REAL fluxOutlet_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-4);
    
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;

    time_mass(0,0) = 0.0;
    time_mass(0,1) = initial_mass;

    time_inj(0,0) = 0.0;
    time_inj(0,1) = inj_data.first;
    time_inj(0,2) = inj_data.second;

    time_prod(0,0) = 0.0;
    time_prod(0,1) = prod_data.first;
    time_prod(0,2) = prod_data.second;
    
    time_fluxInlet(0,0) =0.0;
    time_fluxInlet(0,1) =fluxInlet_data;

    time_fluxOutlet(0,0) =0.0;
    time_fluxOutlet(0,1) =fluxOutlet_data;
    
    time_sfi_interations(0,0) = 0;
    time_sfi_interations(0,1) = 0;
    
#ifdef USING_BOOST
    boost::posix_time::ptime toverhead2 = boost::posix_time::microsec_clock::local_time();
    auto delta_overhead = toverhead2-toverhead1;
    std::cout << "Overhead time " << delta_overhead << std::endl;;
#endif

#ifdef USING_BOOST
    boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
#endif

    for (int it = 1; it <= n_steps; it++) {
        
        sim_time = it*dt;
        sfi_analysis->m_transport_module->SetCurrentTime(dt);
        
        #ifdef USING_BOOST
            boost::posix_time::ptime tsfi_step1 = boost::posix_time::microsec_clock::local_time();
        #endif
        
        sfi_analysis->RunTimeStep();
    
        #ifdef USING_BOOST
            boost::posix_time::ptime tsfi_step2 = boost::posix_time::microsec_clock::local_time();
            auto delta_sfi_step = tsfi_step2-tsfi_step1;
            std::cout << "SFI step time " << delta_sfi_step << std::endl;;
        #endif
             
      if (sim_time >=  current_report_time) {
          TPZFastCondensedElement::fSkipLoadSolution = false;

          std::cout << "Time step number:  " << it << std::endl;
          std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
          sfi_analysis->m_mixed_module->LoadSolution();
          TPZMultiphysicsCompMesh *mphys = dynamic_cast<TPZMultiphysicsCompMesh *>(sfi_analysis->m_mixed_module->Mesh());
          if(!mphys) DebugStop();
            
          mphys->LoadSolutionFromMultiPhysics();
            
          sfi_analysis->PostProcessTimeStep();
          pos++;
          current_report_time = reporting_times[pos];          
          TPZFastCondensedElement::fSkipLoadSolution = true;

      }
        
        REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
        std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
        std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);
        REAL inletflow = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-2);
        REAL outletflow = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-4);
        
        std::cout << "Mass report at time : " << sim_time << std::endl;
        std::cout << "Mass integral :  " << mass << std::endl;

        time_mass(it,0) = sim_time;
        time_mass(it,1) = mass;

        time_inj(it,0) = sim_time;
        time_inj(it,1) = inj_data.first;
        time_inj(it,2) = inj_data.second;

        time_prod(it,0) = sim_time;
        time_prod(it,1) = prod_data.first;
        time_prod(it,2) = prod_data.second;

        time_fluxInlet(it,0)=sim_time;
        time_fluxInlet(it,1) =inletflow;

        time_fluxOutlet(it,0)=sim_time;
        time_fluxOutlet(it,1) =outletflow;
        
        time_sfi_interations(it,0) = sim_time;
        time_sfi_interations(it,1) = sfi_analysis->m_k_iteration;

    }

#ifdef USING_BOOST
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = tsim2-tsim1;
    std::cout << "Overall SFI simulation time " << deltat << std::endl;;
#endif


    // Writing relevant output

    std::ofstream outdata_mass("imrs_mass.txt");
    time_mass.Print("mass = ",outdata_mass,EMathematicaInput);

    std::ofstream outdata_inj("imrs_inj.txt");
    time_inj.Print("inj = ",outdata_inj,EMathematicaInput);

    std::ofstream outdata_prod("imrs_prod.txt");
    time_prod.Print("prod = ",outdata_prod,EMathematicaInput);
    
    std::ofstream outdata_fluxin("imrs_fluxinlet.txt");
    time_fluxInlet.Print("fluxint = ",outdata_fluxin,EMathematicaInput);
    
    std::ofstream outdata_fluxout("imrs_fluxoutlet.txt");
    time_fluxOutlet.Print("fluxint = ",outdata_fluxout,EMathematicaInput);
    
    std::ofstream outdata_iterations("imrs_iterations.txt");
    time_sfi_interations.Print("iters = ",outdata_iterations,EMathematicaInput);
}

void PaperTest3D(){
    
    std::string geometry_file = "gmsh/reservoir_2d_paper.msh";
    int n_layers = 4;
    bool is3D_Q = true;
    bool printMesh_Q = true;
    gRefDBase.InitializeAllUniformRefPatterns();
    TPZGeoMesh *gmesh = CreateGeoMeshMHM3DTest( geometry_file,  n_layers, printMesh_Q, is3D_Q);

    TMRSApproxSpaceGenerator aspace;
    TMRSDataTransfer sim_data  = SettingPaper3D();
    sim_data.mTGeometry.mSkeletonDiv =1;

    aspace.SetGeometry(gmesh);
    std::string name = "paper_3d_test_geo";
    aspace.PrintGeometry(name);
    aspace.ApplyUniformRefinement(2);
    std::string name_ref = "paper_3d_test_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);

    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();


    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();

  
    
    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::EConstantFunction);

    auto kappa_phi = reservoir_properties.Create_Kappa_Phi();
    auto s0 = reservoir_properties.Create_s0();

    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kappa_phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);

    // Render a graphical map
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
        PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
    std::cout << "Spatial properties are transferred." << std::endl;
    std::cout << "Memory used by SpatialPropertiesMap is released. " << std::endl;
   

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;


    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];

     // Mass integral - Injection - Production data
    TPZFNMatrix<200,REAL> time_mass(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_inj(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_prod(n_steps+1,3,0.0);

    // Print initial condition
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-2);
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-4);
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
    sfi_analysis->PostProcessTimeStep();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();

    sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
    std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
    std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;

    time_mass(0,0) = 0.0;
    time_mass(0,1) = initial_mass;

    time_inj(0,0) = 0.0;
    time_inj(0,1) = inj_data.first;
    time_inj(0,2) = inj_data.second;

    time_prod(0,0) = 0.0;
    time_prod(0,1) = prod_data.first;
    time_prod(0,2) = prod_data.second;

    for (int it = 1; it <= n_steps; it++) {
    // TPZFastCondensedElement::fSkipLoadSolutionfSkipLoadSolution = false;
     sim_time = it*dt;
     sfi_analysis->m_transport_module->SetCurrentTime(dt);
     sfi_analysis->RunTimeStep();
     

     if (sim_time >=  current_report_time) {
         std::cout << "Time step number:  " << it << std::endl;
         std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
         sfi_analysis->PostProcessTimeStep();
         pos++;
         current_report_time = reporting_times[pos];
         
         REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
         sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
         std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
         std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);
         std::cout << "Mass report at time : " << sim_time << std::endl;
         std::cout << "Mass integral :  " << mass << std::endl;
         
         time_mass(it,0) = sim_time;
         time_mass(it,1) = mass;
         
         time_inj(it,0) = sim_time;
         time_inj(it,1) = inj_data.first;
         time_inj(it,2) = inj_data.second;
         
         time_prod(it,0) = sim_time;
         time_prod(it,1) = prod_data.first;
         time_prod(it,2) = prod_data.second;
         
     }
      // TPZFastCondensedElement::fSkipLoadSolutionfSkipLoadSolution = true;
    }

    std::cout  << "Number of transport equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;

    // Writing relevant output

    std::ofstream outdata_mass("imrs_mass.txt");
    time_mass.Print("mass = ",outdata_mass,EMathematicaInput);

    std::ofstream outdata_inj("imrs_inj.txt");
    time_inj.Print("inj = ",outdata_inj,EMathematicaInput);

    std::ofstream outdata_prod("imrs_prod.txt");
    time_prod.Print("prod = ",outdata_prod,EMathematicaInput);
    
}

void SimpleTest3D(){
    
    TMRSDataTransfer sim_data  = Setting3D();
    
    std::string geometry_file = "gmsh/simple_2D_coarse.msh";
    std::string name = "simplemesh.vtk";
    
    TMRSApproxSpaceGenerator aspace;
//    aspace.LoadGeometry(geometry_file);
    double nx = 2;
    double ny = 2;
    double nz = 1;
    double L = 10;
    double H = 10;
    double W = 10;
    aspace.CreateUniformMesh(nx, L, ny, H,nz,W);
    aspace.ApplyUniformRefinement(3);
    
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

    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_kappa(TMRSPropertiesFunctions::EConstantFunction);
    reservoir_properties.set_function_type_phi(TMRSPropertiesFunctions::EConstantFunction);
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::EConstantFunction);
    
    auto kx = reservoir_properties.Create_Kx();
    auto ky = reservoir_properties.Create_Ky();
    auto kz = reservoir_properties.Create_Kz();
    auto phi = reservoir_properties.Create_phi();
    auto s0 = reservoir_properties.Create_s0();
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kx,ky,kz,phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
    
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
    sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
    sfi_analysis->PostProcessTimeStep();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;
    
    for (int it = 1; it <= n_steps; it++) {
       // TPZFastCondensedElement::fSkipLoadSolutionfSkipLoadSolution = false;
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
       // TPZFastCondensedElement::fSkipLoadSolutionfSkipLoadSolution = false;
    }
    
    std::cout  << "Number of transport equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
}
void UNISIMTest(){
    
#ifdef USING_BOOST
    boost::posix_time::ptime tspatialmap1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // spatial properties
    int64_t n_cells = 38466;
    std::string grid_data = "maps/corner_grid_coordinates.dat";
    std::string props_data = "maps/corner_grid_props.dat";
    TRMSpatialPropertiesMap properties_map;
    std::vector<size_t> SAMe_blocks = {5,5,5};
    properties_map.SetCornerGridMeshData(n_cells, grid_data, props_data, SAMe_blocks);

    TMRSPropertiesFunctions reservoir_properties;
    reservoir_properties.set_function_type_s0(TMRSPropertiesFunctions::EConstantFunction);

    auto kappa_phi = reservoir_properties.Create_Kappa_Phi(properties_map);
    auto s0 = reservoir_properties.Create_s0();
    
#ifdef USING_BOOST
    boost::posix_time::ptime tspatialmap2 = boost::posix_time::microsec_clock::local_time();
    auto delta_spatialmap = tspatialmap2-tspatialmap1;
    std::cout << "Spatial map construction time " << delta_spatialmap << std::endl;;
#endif

#ifdef USING_BOOST
    boost::posix_time::ptime toverhead1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    std::string geometry_file2D ="gmsh/UNISIMT4R8P2p5.msh";
    int nLayers = 5;
    bool is3DQ = true;
    bool print3DMesh = true;
    gRefDBase.InitializeAllUniformRefPatterns();
    TPZGeoMesh *gmesh = CreateGeoMeshWithTopeAndBase( geometry_file2D,  nLayers, print3DMesh, is3DQ);
 
    TMRSApproxSpaceGenerator aspace;
    
    TMRSDataTransfer sim_data  = SettingUNISIM();
    
    aspace.SetGeometry(gmesh);
    std::string name="unisim_geo";
    aspace.PrintGeometry(name);
    aspace.ApplyUniformRefinement(2);
    std::string name_ref = "unisim_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();

    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kappa_phi,s0);
    sfi_analysis->SetDataTransfer(&sim_data);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);


    // Common section
    
    // Render a graphical map
    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
    PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
    std::cout << "Spatial properties are transferred." << std::endl;
    std::cout << "Memory used by SpatialPropertiesMap is released. " << std::endl;
    properties_map.Clear();

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;

    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos = 0;
    REAL current_report_time = reporting_times[pos];

    // Mass integral - Injection - Production data
    TPZFNMatrix<200,REAL> time_mass(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_inj(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_prod(n_steps+1,3,0.0);
    TPZFNMatrix<200,REAL> time_fluxInlet(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_fluxOutlet(n_steps+1,2,0.0);
    TPZFNMatrix<200,REAL> time_sfi_interations(n_steps+1,2,0.0);
    // Print initial condition
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-2);
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-4);
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    sfi_analysis->SetMixedMeshElementSolution(sfi_analysis->m_mixed_module->Mesh());
    
    
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
    std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
    std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);
    REAL fluxInlet_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-2);
    REAL fluxOutlet_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-4);
    
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;

    time_mass(0,0) = 0.0;
    time_mass(0,1) = initial_mass;

    time_inj(0,0) = 0.0;
    time_inj(0,1) = inj_data.first;
    time_inj(0,2) = inj_data.second;

    time_prod(0,0) = 0.0;
    time_prod(0,1) = prod_data.first;
    time_prod(0,2) = prod_data.second;
    
    time_fluxInlet(0,0) =0.0;
    time_fluxInlet(0,1) =fluxInlet_data;

    time_fluxOutlet(0,0) =0.0;
    time_fluxOutlet(0,1) =fluxOutlet_data;
    
    time_sfi_interations(0,0) = 0;
    time_sfi_interations(0,1) = 0;
    
#ifdef USING_BOOST
    boost::posix_time::ptime toverhead2 = boost::posix_time::microsec_clock::local_time();
    auto delta_overhead = toverhead2-toverhead1;
    std::cout << "Overhead time " << delta_overhead << std::endl;;
#endif

#ifdef USING_BOOST
    boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
#endif

    for (int it = 1; it <= n_steps; it++) {
        
        sim_time = it*dt;
        sfi_analysis->m_transport_module->SetCurrentTime(dt);
        
        #ifdef USING_BOOST
            boost::posix_time::ptime tsfi_step1 = boost::posix_time::microsec_clock::local_time();
        #endif
        
        sfi_analysis->RunTimeStep();
    
        #ifdef USING_BOOST
            boost::posix_time::ptime tsfi_step2 = boost::posix_time::microsec_clock::local_time();
            auto delta_sfi_step = tsfi_step2-tsfi_step1;
            std::cout << "SFI step time " << delta_sfi_step << std::endl;;
        #endif
        
        if (sim_time >=  current_report_time) {
          std::cout << "Time step number:  " << it << std::endl;
          std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
          sfi_analysis->m_mixed_module->LoadSolution();
          TPZMultiphysicsCompMesh *mphys = dynamic_cast<TPZMultiphysicsCompMesh *>(sfi_analysis->m_mixed_module->Mesh());
          if(!mphys) DebugStop();
            
          mphys->LoadSolutionFromMultiPhysics();
            
          sfi_analysis->PostProcessTimeStep();
          pos++;
          current_report_time = reporting_times[pos];
        }
        
        REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
        sfi_analysis->m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(sim_data.mTNumerics.m_ISLinearKrModelQ);
        std::pair<REAL, REAL> inj_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-2);
        std::pair<REAL, REAL> prod_data = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxWaterOilIntegralbyID(-4);

        REAL inletflow = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-2);
        REAL outletflow = sfi_analysis->m_transport_module->fAlgebraicTransport.FLuxIntegralbyID(-4);
        std::cout << "Mass report at time : " << sim_time << std::endl;
        std::cout << "Mass integral :  " << mass << std::endl;

        time_mass(it,0) = sim_time;
        time_mass(it,1) = mass;

        time_inj(it,0) = sim_time;
        time_inj(it,1) = inj_data.first;
        time_inj(it,2) = inj_data.second;

        time_prod(it,0) = sim_time;
        time_prod(it,1) = prod_data.first;
        time_prod(it,2) = prod_data.second;

        time_fluxInlet(it,0)=sim_time;
        time_fluxInlet(it,1) =inletflow;

        time_fluxOutlet(it,0)=sim_time;
        time_fluxOutlet(it,1) =outletflow;
        
        time_sfi_interations(it,0) = sim_time;
        time_sfi_interations(it,1) = sfi_analysis->m_k_iteration;
    }

#ifdef USING_BOOST
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = tsim2-tsim1;
    std::cout << "Overall SFI simulation time " << deltat << std::endl;;
#endif


    // Writing relevant output

    std::ofstream outdata_mass("imrs_mass.txt");
    time_mass.Print("mass = ",outdata_mass,EMathematicaInput);

    std::ofstream outdata_inj("imrs_inj.txt");
    time_inj.Print("inj = ",outdata_inj,EMathematicaInput);

    std::ofstream outdata_prod("imrs_prod.txt");
    time_prod.Print("prod = ",outdata_prod,EMathematicaInput);
    
    std::ofstream outdata_fluxin("imrs_fluxinlet.txt");
    time_fluxInlet.Print("fluxint = ",outdata_fluxin,EMathematicaInput);
    
    std::ofstream outdata_fluxout("imrs_fluxoutlet.txt");
    time_fluxOutlet.Print("fluxint = ",outdata_fluxout,EMathematicaInput);
    
    std::ofstream outdata_iterations("imrs_iterations.txt");
    time_sfi_interations.Print("iters = ",outdata_iterations,EMathematicaInput);
    
}
TMRSDataTransfer SettingSimple2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressurein = 25.0;
    REAL pressureout = 10.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressurein);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,D_Type,pressureout);
    
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
    sim_data.mTFluidProperties.mOilDensity = 800.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 50;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.000001;
    sim_data.mTNumerics.m_n_steps = 100;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.01*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
//    grav[1] = -9.81*0.0*(1.0e-6);
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 10;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    //    scalnames.Push("p");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("p");
        scalnames.Push("p_average");
        scalnames.Push("kxx");
        scalnames.Push("kyy");
        scalnames.Push("kzz");
        scalnames.Push("lambda");
        
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

TMRSDataTransfer SettingPlumeGravity2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
   
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_on_top = 0.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,N_Type,zero_flux);
      sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,D_Type,pressure_on_top);
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
    sim_data.mTFluidProperties.mOilDensity = 700.0;

    // Numerical controls
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 50;
    sim_data.mTNumerics.m_sfi_tol = 0.000001;
    sim_data.mTNumerics.m_res_tol_transport = 0.00000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00000001;
    sim_data.mTNumerics.m_n_steps = 10;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.1*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 9.81*(1.0e-6);
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 10;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
        scalnames.Push("kxx");
        scalnames.Push("kyy");
        scalnames.Push("kzz");
        scalnames.Push("lambda");

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

TMRSDataTransfer SettingGravity2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
   
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_on_top = 10.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,N_Type,zero_flux);
      sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,D_Type,pressure_on_top);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,N_Type,zero_flux);
 
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 0.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-3,bc_inlet,sat_in);
     sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_outlet,1.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.001;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 800.0;

    // Numerical controls
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 50;
    sim_data.mTNumerics.m_sfi_tol = 0.000001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_n_steps = 100;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 1.0*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = -9.81*(1.0e-6);
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 10;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
        scalnames.Push("kxx");
        scalnames.Push("kyy");
        scalnames.Push("kzz");
        scalnames.Push("lambda");
        scalnames.Push("p");

    }
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = 10*sim_data.mTPostProcess.m_file_time_step;
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
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 25.0;
    REAL pressure_out = 10.0;
    
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
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.001;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 500.0;
    
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 40;
    
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.00001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_n_steps = 20;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.001*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[2] = -9.8*(1.0e-6);
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
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
    REAL time = 2*sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    return sim_data;
}

TMRSDataTransfer SettingPaper2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["Reservoir"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
   
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 25.0;
    REAL pressure_out = 10.0;
        
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-4,D_Type,pressure_out);
    
 
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-4,bc_outlet,0.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.002;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 800.0;


    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 3;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 50;

    sim_data.mTGeometry.mSkeletonDiv = 3;
    sim_data.mTNumerics.m_sfi_tol = 0.000001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_n_steps = 240;

    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 20.0*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = -9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
//    vecnames.Push("q");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("p");
//        scalnames.Push("g_average");
//        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = 20*sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    return sim_data;
}

TMRSDataTransfer SettingPaper3D(){
    
   TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Reservoir"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 15.0;
    REAL pressure_out = 10.0;
    
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
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_outlet,0.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.001;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 500.0;
    
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 50;
    
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.00001;
    sim_data.mTNumerics.m_res_tol_transport = 0.00000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00000001;
    sim_data.mTNumerics.m_n_steps = 200;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 20.0*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = -9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 10;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
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
    REAL time = 2*sim_data.mTPostProcess.m_file_time_step;
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
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 25.0;
    REAL pressure_out = 10.0;
    
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
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.002;
    sim_data.mTFluidProperties.mWaterDensity = 1000.0;
    sim_data.mTFluidProperties.mOilDensity = 800.0;


    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 30;

    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.00001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00001;
    sim_data.mTNumerics.m_n_steps = 1300;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 10*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[2] = -9.8*(1.0e-6);
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 8;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
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
    REAL time = 10*sim_data.mTPostProcess.m_file_time_step;
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
        extend.SetElType(1);
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

TPZGeoMesh * CreateGeoMeshMHM3DTest(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["Reservoir"] = 1;
    dim_name_and_physical_tag[2]["Wellbore_p"] = 1;
    dim_name_and_physical_tag[2]["Wellbore_i"] = 1;
    dim_name_and_physical_tag[1]["producers"] = -4;
    dim_name_and_physical_tag[1]["injectors"] = -2;
    dim_name_and_physical_tag[1]["Reservoir_south"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_West"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_north"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_East"] = -1;
    
    TPZGmshReader Geometry;
    TPZGeoMesh * gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);
    double w = 25.0;
    
    std::string name2D("mesh2d.vtk");
    
    int topID= -3;
    int baseID = -1;
    TPZGeoMesh * returnedMesh = nullptr;
    TPZExtendGridDimension extend(gmesh2d, w);
    extend.SetElType(1);
    returnedMesh = extend.ExtendedMesh(nLayers,topID,baseID);
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
//    std::string filename2 = "Reservoir/base_unisimMOD.txt";
    
     std::string filename1 = "Reservoir/topeMOD.txt";
    std::string filename2 = "Reservoir/baseMOD.txt";
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

//                if (ipoint==0) {
//                    bool find = 0;
//                    int i =1;
//                    while (find==0) {
//                        TPZGeoNode node = gmesh->NodeVec()[ipoint+i];
//                        TPZVec<REAL> coaux(3);
//                        node.GetCoordinates(coaux);
//                        topeinterpol =interpTope(coaux[0],coaux[1]);
//                        if (topeinterpol!=0) {
//                            find=1;
//                        }
//                        i=i+1;
////                        topeinterpol =coaux[2];
//                    }
//                }
//                if (ipoint != 0) {
//                    TPZGeoNode node = gmesh->NodeVec()[ipoint-1];
//                    TPZVec<REAL> coaux(3);
//                    node.GetCoordinates(coaux);
//                    REAL valz = coaux[2];
//                    topeinterpol =coaux[2];
//
//                }
                
                
            }
            if (baseinterpol==0) {
                
                baseinterpol = val_base;
                if (co[0]>1000.00) {
                   baseinterpol = val_base-80;
                }
//                std::cout<<"{"<<co[0]<<","<<co[1]<<"};"<<std::endl;
//                    if (ipoint==npointsPerLayer) {
//                        bool find = 0;
//                        int i =1;
//                        while (find==0) {
//                            TPZGeoNode node = gmesh->NodeVec()[ipoint+i];
//                            TPZVec<REAL> coaux(3);
//                            node.GetCoordinates(coaux);
//                            baseinterpol =interpBase(coaux[0],coaux[1]);
//                            if (baseinterpol!=0) {
//                                find=1;
//                            }
//                            i=i+1;
//    //                        topeinterpol =coaux[2];
//                        }
//                    }
//                    if (ipoint != npointsPerLayer && ilay!=1) {
//                        TPZGeoNode node = gmesh->NodeVec()[ipoint-1];
//                        TPZVec<REAL> coaux(3);
//                        node.GetCoordinates(coaux);
//                        REAL valz = coaux[2];
//                        baseinterpol =coaux[2];
//
//                    }


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

void PostProcessResProps(TPZMultiphysicsCompMesh *cmesh, TPZAlgebraicTransport *alTransport){
    
    
     TPZAnalysis *an = new TPZAnalysis(cmesh,false);
     int ncells = alTransport->fCellsData.fVolume.size();
     TPZFMatrix<STATE> phi(ncells,1,0.1);
     TPZFMatrix<STATE> Kx(ncells,1,1.0);
     TPZFMatrix<STATE> Ky(ncells,1,2.0);
     TPZFMatrix<STATE> Kz(ncells,1,3.0);
    
    for (int ivol = 0; ivol<ncells; ivol++) {
        int eq = alTransport->fCellsData.fEqNumber[ivol];
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
     int dim = cmesh->Reference()->Dimension();
     TPZStack<std::string,10> scalnames, vecnames;
//     vecnames.Push("q");
     scalnames.Push("Porosity");
     scalnames.Push("Permeability_x");
     scalnames.Push("Permeability_y");
     scalnames.Push("Permeability_z");
     std::string file("Props.vtk");

     an->DefineGraphMesh(dim,scalnames,vecnames,file);
     an->PostProcess(0,dim);

    
}
