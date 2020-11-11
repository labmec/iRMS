
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void SimpleTest2DHDiv();
TMRSDataTransfer SettingSimple2DHdiv();
void  ForcingFunction (const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

int main(){
     SimpleTest2DHDiv();
}

void SimpleTest2DHDiv(){
    
    //sim_data contains the simulation configuration, including contour conditions and muneric controls
    TMRSDataTransfer sim_data  = SettingSimple2DHdiv();
    TMRSApproxSpaceGenerator aspace;
    
    //The problem geometry is a rectangle with LxY dimentions and elements number in x y and nxXny
    int nx=25;
    int ny=25;
    REAL L =10.0;
    REAL Y =10.0;
    aspace.CreateUniformMesh(nx, L, ny, Y);
    

    std::string name = "2D_geo";
//    aspace.PrintGeometry(name);
    aspace.ApplyUniformRefinement(0);
    std::string name_ref = "2D_ref_geo";
//    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    
    //This parameter should be always "true"
    bool UsePardiso_Q = true;
    
    //Multiphysic mesh creation of the mixed problem
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    TMRSMixedAnalysis *mixedAnal = new TMRSMixedAnalysis(mixed_operator,must_opt_band_width_Q);
    
    //If the parameter "UsingPzSparse" is true, it uses the pz sparse matrix, otherwise it uses eigen sparse matrix
    bool UsingPzSparse = false;

    //The parallelism is just implemented for the "UsingPzSparse=True" case, with eigen for now is running in serial (the next task to do)
    mixedAnal->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    mixedAnal->SetDataTransfer(&sim_data);
    mixedAnal->Assemble();
    size_t n_dof = mixedAnal->Solver().Matrix()->Rows();

    #ifdef USING_BOOST
        boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
    #endif
    mixedAnal->Solve();
    #ifdef USING_BOOST
        boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
        auto deltat = tsim2-tsim1;
        std::cout << "Overal solve calling time " << deltat << std::endl;
    #endif
    std::cout << "Number of dof = " << n_dof << std::endl;
    mixed_operator->UpdatePreviousState(-1);
    mixedAnal->fsoltransfer.TransferFromMultiphysics();
    mixedAnal->PostProcessTimeStep();
  
}

TMRSDataTransfer SettingSimple2DHdiv(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressurein = 25;
    REAL pressureout = 10;
    
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
    sim_data.mTNumerics.m_res_tol_transport = 0.00001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00001;
    sim_data.mTNumerics.m_n_steps = 1;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.03 ;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q          = false;
    std::vector<REAL> grav(3,0.0);
    grav[1] = -0.0;//-9.81; //Needs to be multiplied by a constant
    sim_data.mTNumerics.m_gravity = grav;
    
    //Update to m_ISLinearKrModelQ=true if is required a linear relative permeabilities model
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
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
void  ForcingFunction (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    double fx= 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    disp[0]=fx;
}
