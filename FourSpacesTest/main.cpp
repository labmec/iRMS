
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
//#include "TPZCondensedElement.h"
#include "pzcondensedcompel.h"
#include <libInterpolate/Interpolate.hpp>   
#include <libInterpolate/AnyInterpolator.hpp>

void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix);
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
    aspace.CreateUniformMesh(2, 100, 1, 10);
    
    aspace.PrintGeometry(name);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    TPZCompMesh *mixed = dynamic_cast<TPZCompMesh*>(mixed_operator);
   
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    std::ofstream file("mixed.vtk");
    
    TPZVTKGeoMesh::PrintCMeshVTK(mixed_operator, file);
    
 
    
//    aspace.LinkMemory(mixed_operator, transport_operator);
//    mixed_operator->ComputeNodElCon();

    std::ofstream file1("mixed_before.txt");
    mixed_operator->Print(file1);
    CreatedCondensedElements(mixed_operator, false, true);
    mixed_operator->CleanUpUnconnectedNodes();
    std::ofstream file2("mixed_after.txt");
    mixed_operator->Print(file2);
//    
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    sfi_analysis->SetDataTransfer(&sim_data);
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

    solution_n *= 0.0;
   
    
    
    for (int it = 1; it <= n_steps; it++) {
        sfi_analysis->RunTimeStepWithOutMemory(solution_n);
        solution_n.Print(std::cout);
        sim_time = it*dt;
        sfi_analysis->m_transport_module->SetCurrentTime(sim_time);
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
    

    //Relative permermeabilities
    TRSLinearInterpolator krw, kro ;
    std::string name_krw("PetroPhysics/krw_linear.txt");
    std::string name_kro("PetroPhysics/krow_linear.txt");
    
    krw.ReadData(name_krw,true);
    kro.ReadData(name_kro,true);
    
    kro.SetLeftExtension(TRSLinearInterpolator::ELinear);
    krw.SetLeftExtension(TRSLinearInterpolator::ELinear);
    kro.SetRightExtension(TRSLinearInterpolator::ELinear);
    krw.SetRightExtension(TRSLinearInterpolator::ELinear);
    
    //    kro.SetInterpolationType(TRSLinearInterpolator::InterpType::THermite);
//    krw.SetInterpolationType(TRSLinearInterpolator::InterpType::THermite);
    
    sim_data.mTPetroPhysics.mLayer_Krw_RelPerModel.resize(1);
    sim_data.mTPetroPhysics.mLayer_Kro_RelPerModel.resize(1);
    sim_data.mTPetroPhysics.mLayer_Krw_RelPerModel[0] = krw;
    sim_data.mTPetroPhysics.mLayer_Kro_RelPerModel[0] = kro;
    
    // Fractional flows composition
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double sw, double p)> fw = [](TRSLinearInterpolator & krw, TRSLinearInterpolator  &kro, double sw, double p) {
        
        double mu_w = 1.0;
        double mu_o = 1.0;
        double Bw = 1.0;
        double Bo = 1.0;
        std::tuple<double, double> krw_t = krw.ValDeriv(sw);
        std::tuple<double, double> kro_t = kro.ValDeriv(sw);
        double krwv = std::get<0>(krw_t);
        double krov = std::get<0>(kro_t);
        double dkrw_dswv = std::get<1>(krw_t);
        double dkro_dswv = std::get<1>(kro_t);
        double lwv  = krwv/(mu_w*Bw);
        double dlw_dswv  = dkrw_dswv/(mu_w*Bw);
        double lov  = krov/(mu_o*Bo);
        double dlo_dswv  = dkro_dswv/(mu_o*Bo);
        double lv   = lwv + lov;
        double dl_dswv   = dlw_dswv + dlo_dswv;
        double fwv  = lwv/lv;
        double dfw_dswv  = (dlw_dswv/lv) - lwv*(dl_dswv/(lv*lv));
        std::tuple<double, double, double> fw_t(fwv, dfw_dswv, 0.0);
        return fw_t;
    };
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double sw, double p)> fo = [](TRSLinearInterpolator & krw, TRSLinearInterpolator & kro, double sw, double p) {
        
        double mu_w = 1.0;
        double mu_o = 1.0;
        double Bw = 1.0;
        double Bo = 1.0;
        
        std::tuple<double, double> krw_t = krw.ValDeriv(sw);
        std::tuple<double, double> kro_t = kro.ValDeriv(sw);
        double krwv = std::get<0>(krw_t);
        double krov = std::get<0>(kro_t);
        double dkrw_dswv = std::get<1>(krw_t);
        double dkro_dswv = std::get<1>(kro_t);
        double lwv  = krwv/(mu_w*Bw);
        double dlw_dswv  = dkrw_dswv/(mu_w*Bw);
        double lov  = krov/(mu_o*Bo);
        double dlo_dswv  = dkro_dswv/(mu_o*Bo);
        double lv   = lwv + lov;
        double dl_dswv   = dlw_dswv + dlo_dswv;
        double fov  = lov/lv;
        double dfo_dswv  = (dlo_dswv/lv) - lov*(dl_dswv/(lv*lv));
        std::tuple<double, double, double> fo_t(fov,dfo_dswv, 0.0);
        return fo_t;
    };
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double sw, double p)> lambda = [](TRSLinearInterpolator & krw, TRSLinearInterpolator & kro, double sw, double p) {
        
        
        double mu_w = 1.0;
        double mu_o = 1.0;
        double Bw = 1.0;
        double Bo = 1.0;
        
        std::tuple<double, double> krw_t = krw.ValDeriv(sw);
        std::tuple<double, double> kro_t = kro.ValDeriv(sw);
        double krwv = std::get<0>(krw_t);
        double krov = std::get<0>(kro_t);
        double dkrw_dswv = std::get<1>(krw_t);
        double dkro_dswv = std::get<1>(kro_t);
        double lwv  = krwv/(mu_w*Bw);
        double dlw_dswv  = dkrw_dswv/(mu_w*Bw);
        double lov  = krov/(mu_o*Bo);
        double dlo_dswv  = dkro_dswv/(mu_o*Bo);
        double lv   = lwv + lov;
        double dl_dswv   = dlw_dswv + dlo_dswv;
        std::tuple<double, double, double> l_t(lv, dl_dswv, 0.0);
        return l_t;
    };
    
    sim_data.mTMultiphaseFunctions.mLayer_fw[0] = fw;
    sim_data.mTMultiphaseFunctions.mLayer_fo[0] = fo;
    sim_data.mTMultiphaseFunctions.mLayer_lambda[0] = lambda;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 5;
    sim_data.mTNumerics.m_max_iter_transport = 5;
    sim_data.mTNumerics.m_max_iter_sfi = 5;
    sim_data.mTNumerics.m_res_tol_mixed = 0.01;
    sim_data.mTNumerics.m_corr_tol_mixed = 0.01;
    sim_data.mTNumerics.m_res_tol_transport = 0.01;
    sim_data.mTNumerics.m_corr_tol_transport = 0.01;
    sim_data.mTNumerics.m_n_steps = 50;
    sim_data.mTNumerics.m_dt      = 0.01;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = false;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    scalnames.Push("p");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 0.01;
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

/// created condensed elements for the elements that have internal nodes
void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix)
{
    //    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        if (KeepOneLagrangian) {
            int count = 0;
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    count++;
                    if(count == 1 && c.NState() == 1)
                    {
                        break;
                    } else if(count == 2 && c.NState() == 2)
                    {
                        break;
                    } else if(count == 3 && c.NState() == 3)
                    {
                        break;
                    }
                    
                }
            }
        }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
        {
            TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
            cond->SetPermeability(4.0);
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}




