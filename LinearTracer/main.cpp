
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

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"


TMRSDataTransfer Setting2D();
TMRSDataTransfer SettingSimple2D();
TMRSDataTransfer Setting3D();
TMRSDataTransfer SettingSimpleMHM2D();
TMRSDataTransfer SettingSimple2DQuads();
void MHMSimpleTest();
void SimpleTest();
//void InsertMaterialObjects(TPZMHMixedMeshControl &control);

void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);
TPZGeoMesh *MHMMesh();
int main(){
    InitializePZLOG();
//     SimpleTest();
   MHMSimpleTest();
    return 0;
}



void SimpleTest(){
    bool is_3D_Q = false;
    bool is_2D_Coarse_Q = true;
    
    TMRSDataTransfer sim_data;
    
    std::string geometry_file, name;
    if (is_3D_Q) {
        geometry_file = "gmsh/TopeMeshPT.msh";
        name = "reservoir_3d";
        sim_data = Setting3D();
    }else{
        if (is_2D_Coarse_Q) {
            // geometry_file = "gmsh/reservoir_2d_coarse.msh";
            geometry_file = "gmsh/simple_2D_coarse.msh";
            geometry_file = "gmsh/reservoir_2d_coarse.msh";
            name = "reservoir_2d";
            sim_data = Setting2D();
        }
        else{
            //            geometry_file = "gmsh/reservoir_2d_coarse.msh";
            geometry_file = "gmsh/reservoir_2d_fine.msh";
        }
        
    }
    
    
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);
//    aspace.CreateUniformMesh(2, 2,2,2);
    aspace.GenerateMHMUniformMesh(1);
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
    std::ofstream file("mixed.vtk");
    
    TPZVTKGeoMesh::PrintCMeshVTK(mixed_operator, file);
    
    aspace.LinkMemory(mixed_operator, transport_operator);
    //    aspace.BuildMixedSCStructures();
    mixed_operator->ComputeNodElCon();
    TPZCompMeshTools::GroupElements(mixed_operator);
    TPZCompMeshTools::CreatedCondensedElements(mixed_operator, true);
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    sfi_analysis->SetDataTransfer(&sim_data);
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    
    // Fill time steps vector
    
    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    
    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];
    
    for (int it = 1; it <= n_steps; it++) {
        sfi_analysis->RunTimeStep();
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

void MHMSimpleTest(){
   
    
    TMRSDataTransfer sim_data;
    sim_data = SettingSimpleMHM2D();
   
    sim_data = SettingSimpleMHM2D();
    TMRSApproxSpaceGenerator aspace;
    aspace.CreateUniformMesh(2, 2,2,2);
    aspace.GenerateMHMUniformMesh(2);
    aspace.SetDataTransfer(sim_data);
    
    
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
//    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    std::ofstream file("mixed.vtk");
    
    TPZCompMesh *fluxcmesh_initial = mixed_operator->MeshVector()[0];
    TPZCompMesh *pressure_initial = mixed_operator->MeshVector()[1];
    std::ofstream fileflux("initial_flux.txt");
    std::ofstream filepressure("initial_pressure.txt");
    fluxcmesh_initial->Print(fileflux);
    pressure_initial->Print(filepressure);
    std::ofstream multimixed("multimixed.txt");
    mixed_operator->Print(multimixed);
    
    TPZVTKGeoMesh::PrintCMeshVTK(transport_operator, file);
    {
        std::ofstream transport("transport.txt");
        transport_operator->Print(transport);
    }
    
//    aspace.LinkMemory(mixed_operator, transport_operator);
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,true);
    sfi_analysis->Configure(0, true);
    sfi_analysis->SetDataTransfer(&sim_data);

    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;

    // Fill time steps vector

    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;

    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];

    for (int it = 1; it <= n_steps; it++) {
        sfi_analysis->RunTimeStep();
        sim_time = it*dt;
        sfi_analysis->m_transport_module->SetCurrentTime(sim_time);
        if (sim_time >=  current_report_time) {
            std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
            sfi_analysis->PostProcessTimeStep();
            pos++;
            current_report_time =reporting_times[pos];

        }
    }
    
//    TPZAnalysis an(mixed_operator,true);
//
//    TPZSymetricSpStructMatrix strmat(mixed_operator);
//    strmat.SetNumThreads(8);
//
//    an.SetStructuralMatrix(strmat);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELDLt);
//    an.SetSolver(step);
//    std::cout << "Assembling\n";
//    an.Assemble();
//    std::ofstream filemate("MatrixCoarse.txt");
//
//
//
//    std::cout << "Solving\n";
//    an.Solve();
//    std::cout << "Finished\n";
//    an.LoadSolution(); // compute internal dofs
//
//
//    TPZStack<std::string> scalar, vectors;
//    TPZManVector<std::string,10> scalnames(1), vecnames(1);
//    vecnames[0]  = "Flux";
//    scalnames[0] = "Pressure";
//
//    std::string name_coarse("Mixedresults.vtk");
//
//    an.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
//    an.PostProcess(0,2);

    
    
//    std::cout<<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"<<std::endl;
//    mixed_operator->Print();
//
//    std::ofstream mult("malhamultph.txt");
//    mixed_operator->Print(mult);
//    std::cout<<"Nel: "<<mixed_operator->NElements()<<std::endl;
//    TPZAnalysis an(mixed_operator,true);
//
//    TPZSymetricSpStructMatrix strmat(mixed_operator);
//    strmat.SetNumThreads(8);
//
//    an.SetStructuralMatrix(strmat);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELDLt);
//    an.SetSolver(step);
//    std::cout << "Assembling\n";
//    an.Assemble();
//    std::ofstream filemate("MatrixCoarse.txt");
//
//
//
//    std::cout << "Solving\n";
//    an.Solve();
//    std::cout << "Finished\n";
////    an.LoadSolution(); // compute internal dofs
//
//
//    TPZStack<std::string> scalar, vectors;
//    TPZManVector<std::string,10> scalnames(1), vecnames(1);
//    vecnames[0]  = "state";
//    scalnames[0] = "state";
//
//    std::string name_coarse("Mixedresults.vtk");
//
//    an.DefineGraphMesh(2, scalnames, vecnames, name_coarse);
//    an.PostProcess(0,2);
    
    
}

TMRSDataTransfer Setting2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_wbregion_p"] = 2;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_wbregion_i"] = 3;
    
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 30.0;
    REAL pressure_out = 20.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(6);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(4,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(5,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(6,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(7,N_Type,zero_flux);
    
    //well pressure
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[4] = std::make_tuple(10,D_Type,pressure_out);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[5] = std::make_tuple(11,D_Type,pressure_in);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(6);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(4,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(5,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(6,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(7,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[4] = std::make_tuple(10,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[5] = std::make_tuple(11,bc_inlet,sat_in);

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
    sim_data.mTNumerics.m_max_iter_mixed = 20;
    sim_data.mTNumerics.m_max_iter_transport = 20;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_n_steps = 2;
    sim_data.mTNumerics.m_dt      = 1.0;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("u_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 1.0;
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
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["matrix"] = 1;
//    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["d_wbregion_p"] = 2;
//    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["d_wbregion_i"] = 3;
    
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 30.0;
    REAL pressure_out = 20.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(5);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(2,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(3,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(4,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(7,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[4] = std::make_tuple(8,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[5] = std::make_tuple(9,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[6] = std::make_tuple(10,N_Type,zero_flux);
//    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[7] = std::make_tuple(11,N_Type,zero_flux);
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(5,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[4] = std::make_tuple(6,D_Type,pressure_out);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(5);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(2,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(3,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(4,bc_inlet,0.0);
//    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(7,bc_inlet,0.0);
//    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[4] = std::make_tuple(8,bc_inlet,0.0);
//    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[5] = std::make_tuple(9,bc_inlet,0.0);
//    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[6] = std::make_tuple(10,bc_inlet,0.0);
//    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[7] = std::make_tuple(11,bc_inlet,0.0);
    
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[4] = std::make_tuple(6,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(5,bc_inlet,sat_in);
    
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
    sim_data.mTNumerics.m_max_iter_mixed = 20;
    sim_data.mTNumerics.m_max_iter_transport = 20;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_n_steps = 10;
    sim_data.mTNumerics.m_dt      = 1.0;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q          = false;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator_3d.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator_3d.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("u_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 1.0;
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


TMRSDataTransfer SettingSimple2D(){
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    REAL zero_flux=0.0;
    REAL pressure_in = 30.0;
    REAL pressure_out = 20.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(2,N_Type,zero_flux);
    
    //domain pressure
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(5,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(4,D_Type,pressure_out);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(2,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(5,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(4,bc_outlet,0.0);
    
    //Relative permermeabilities
    TRSLinearInterpolator krw, kro ;
    std::string name_krw("PetroPhysics/krw_linear.txt");
    std::string name_kro("PetroPhysics/krow_linear.txt");
    
    krw.ReadData(name_krw,true);
    kro.ReadData(name_kro,true);
    

//    kro.SetLeftExtension(TRSLinearInterpolator::Enone,1.0);
//    krw.SetLeftExtension(TRSLinearInterpolator::Enone,1.0);
//    kro.SetRightExtension(TRSLinearInterpolator::Enone,1.0);
//    krw.SetRightExtension(TRSLinearInterpolator::Enone,1.0);
    
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
    sim_data.mTNumerics.m_max_iter_mixed = 10;
    sim_data.mTNumerics.m_max_iter_transport = 10;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_n_steps = 10;
    sim_data.mTNumerics.m_dt      = 1.0;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
   
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("u_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 1.0;
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
TMRSDataTransfer SettingSimpleMHM2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 100.0;
    REAL pressure_out = 1000.0;
    
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
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-3,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_inlet,sat_in);

    
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
    sim_data.mTNumerics.m_max_iter_mixed = 20;
    sim_data.mTNumerics.m_max_iter_transport = 20;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_n_steps = 5;
    sim_data.mTNumerics.m_dt      = 0.0001;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("u_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 0.0001;
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
TMRSDataTransfer SettingSimple2DQuads(){
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    REAL zero_flux=0.0;
    REAL pressure_in = 10.0;
    REAL pressure_out = 30.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    
    //domain pressure
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,D_Type,pressure_out);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-3,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_outlet,0.0);
    
    //Relative permermeabilities
    TRSLinearInterpolator krw, kro ;
    std::string name_krw("PetroPhysics/krw_linear.txt");
    std::string name_kro("PetroPhysics/krow_linear.txt");
    
    krw.ReadData(name_krw,true);
    kro.ReadData(name_kro,true);
    

//    kro.SetLeftExtension(TRSLinearInterpolator::Enone,1.0);
//    krw.SetLeftExtension(TRSLinearInterpolator::Enone,1.0);
//    kro.SetRightExtension(TRSLinearInterpolator::Enone,1.0);
//    krw.SetRightExtension(TRSLinearInterpolator::Enone,1.0);
    
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
    sim_data.mTNumerics.m_max_iter_mixed = 10;
    sim_data.mTNumerics.m_max_iter_transport = 10;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_mixed = 1.0e-7;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_corr_tol_transport = 1.0e-7;
    sim_data.mTNumerics.m_n_steps = 10;
    sim_data.mTNumerics.m_dt      = 1.0;
    sim_data.mTNumerics.m_four_approx_spaces_Q = false;
    sim_data.mTNumerics.m_mhm_mixed_Q          = false;
   
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("u_average");
    }
    sim_data.mTPostProcess.m_file_time_step = 1.0;
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
