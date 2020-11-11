
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif


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


TMRSDataTransfer SettingPaper2D();
void PaperTest2D();

void SimpleTest2D();
//
int main(){
    PaperTest2D();
}

void PaperTest2D(){
    
   
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
    
    aspace.ApplyUniformRefinement(1);
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
    
    auto kappa_phi = reservoir_properties.Create_Kappa_Phi();
    auto s0 = reservoir_properties.Create_s0();

//    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q,kappa_phi,s0);
TMRSSFIAnalysis *sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->SetDataTransfer(&sim_data);
    bool usingpzSparse = false;
    sfi_analysis->Configure(n_threads, UsePardiso_Q, usingpzSparse);
// not use
//    TPZMultiphysicsCompMesh *AuxPosProcessProps = aspace.BuildAuxPosProcessCmesh(sfi_analysis->fAlgebraicDataTransfer);
//    PostProcessResProps(AuxPosProcessProps, &sfi_analysis->m_transport_module->fAlgebraicTransport);
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
    
    // Print initial condition
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-2);
    sfi_analysis->m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-4);
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
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
    
   // TPZFastCondensedElement::fSkipLoadSolutionfSkipLoadSolution = true;

#ifdef USING_BOOST
    boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
#endif

    for (int it = 1; it <= n_steps; it++) {
      // TPZFastCondensedElement::fSkipLoadSolution = false;
      sim_time = it*dt;
      sfi_analysis->m_transport_module->SetCurrentTime(sim_time);
      sfi_analysis->RunTimeStep();
      
     
      if (sim_time >=  current_report_time) {
          std::cout << "Time step number:  " << it << std::endl;
          std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
         
          sfi_analysis->m_mixed_module->LoadSolution();
          TPZMultiphysicsCompMesh *mphys = dynamic_cast<TPZMultiphysicsCompMesh *>(sfi_analysis->m_mixed_module->Mesh());
          if(!mphys) DebugStop();
//          mphys->LoadSolutionFromMultiPhysics();
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

#ifdef USING_BOOST
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = tsim2-tsim1;
    std::cout << "Total timestepping time " << deltat << std::endl;
#endif

    std::cout  << "Number of transport equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
    
    // Writing relevant output
    
    std::ofstream outdata_mass("imrs_mass.txt");
    time_mass.Print("mass = ",outdata_mass,EMathematicaInput);

    std::ofstream outdata_inj("imrs_inj.txt");
    time_inj.Print("inj = ",outdata_inj,EMathematicaInput);
    
    std::ofstream outdata_prod("imrs_prod.txt");
    time_prod.Print("prod = ",outdata_prod,EMathematicaInput);
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
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 800.0;


    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 2;
    sim_data.mTNumerics.m_max_iter_transport = 10;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_sfi_tol = 0.001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 120;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 10.0 ;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.000045;//-9.81;
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
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
