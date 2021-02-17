
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TPZAlgebraicDataTransfer.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


TMRSDataTransfer SettingSimple2D();
void SimpleTest2D();
//
int main(){
    SimpleTest2D();
    return 0;
}
void SimpleTest2D(){

    //sim_data contains the simulation configuration, including contour conditions and muneric controls
    TMRSDataTransfer sim_data  = SettingSimple2D();
   
    //The problem geometry is a rectangle with LxY dimentions and elements number in x y and nxXn
    TMRSApproxSpaceGenerator aspace;
    aspace.CreateUniformMesh(10, 10, 10, 10);
    
    std::string name = "2D_geo";
    aspace.PrintGeometry(name);
    aspace.ApplyUniformRefinement(1);
    std::string name_ref = "2D_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    //This parameter should be always "true"
    bool UsePardiso_Q = true;
    
    //Multiphysic mesh creation of the mixed problem
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->SetDataTransfer(&sim_data);

    //If the parameter "UsingPzSparse" is true, it uses the pz sparse matrix, otherwise it uses eigen sparse matrix
    bool usingpzSparse = false;
    
      //The parallelism is just implemented for the "UsingPzSparse=True" case, with eigen for now is running in serial (the next task to do)
    sfi_analysis->Configure(n_threads, UsePardiso_Q, usingpzSparse);
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    
    
    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    
    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];
    
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;
    
    TPZFastCondensedElement::fSkipLoadSolution = true;
    for (int it = 1; it <= n_steps; it++) {
        sim_time = it*dt;
        if (sim_time >=  current_report_time) {
            TPZFastCondensedElement::fSkipLoadSolution = false;
        }
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
            TPZFastCondensedElement::fSkipLoadSolution = true;
           
            
            
            
        }
    }
    
    std::cout  << "Number of transport equations = " << sfi_analysis->m_transport_module->Solution().Rows() << std::endl;
    
}

//TMRSPropertiesFunctions
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
    sim_data.mTFluidProperties.mReferencePressure =1.013e5;
    sim_data.mTFluidProperties.mWaterViscosity = 0.8;
    sim_data.mTFluidProperties.mOilViscosity = 0.7;
    sim_data.mTPetroPhysics.mWaterViscosity = sim_data.mTFluidProperties.mWaterViscosity;
    sim_data.mTPetroPhysics.mOilViscosity = sim_data.mTFluidProperties.mOilViscosity;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilCompressibility = 1.0e-9;
    sim_data.mTFluidProperties.mWaterCompressibility = 1.0e-9;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 2;
    sim_data.mTNumerics.m_max_iter_transport = 10;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_sfi_tol = 0.001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 5;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.01 ;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
     sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4Space;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-0.0005;//-9.81;//Needs to be multiplied by a constant
    sim_data.mTNumerics.m_gravity = grav;
    
    //Update to m_ISLinearKrModelQ=true if is required a linear relative permeabilities model
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    //    scalnames.Push("p");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
//        scalnames.Push("Pressure");
        scalnames.Push("p_average");
        scalnames.Push("g_average");
        //        scalnames.Push("kyy");
        //        scalnames.Push("kzz");
        //        scalnames.Push("lambda");
        
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
