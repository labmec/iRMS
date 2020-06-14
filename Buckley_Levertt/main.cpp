
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
#include "pzshapepiram.h"

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"

void BuckleyLeverett(REAL swf, REAL xswf, TPZVec<REAL> &pt,TPZVec<REAL> q, REAL phi,REAL time, TPZVec<REAL> &Saturation, TMRSDataTransfer sim_data);
void BuckleyLeverettCase( TMRSDataTransfer simdata);
void SetSolution(REAL swf, REAL xswf, TPZCompMesh *cmesh, TMRSDataTransfer simdata, REAL sim_time);
REAL S_Newton(REAL x, REAL t, REAL u, REAL Swr, REAL Sor, REAL phi, REAL s_shok, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta, REAL epsilon);
REAL dfdsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta);
REAL df2dsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta);

TMRSDataTransfer SettingSimple2D();



int main(){
    
    
    BuckleyLeverettCase(SettingSimple2D());
    
    
    return 0;
}


TMRSDataTransfer SettingSimple2D(){
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    REAL zero_flux=0.0;
    REAL pressure_in = 1.0;
    REAL pressure_out = 0.0;
    
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
    std::string name_krw("PetroPhysics/krw_quadratic2.txt");
    std::string name_kro("PetroPhysics/krow_quadratic2.txt");
    
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
    
    kro.SetInterpolationType(TRSLinearInterpolator::InterpType::THermite);
    krw.SetInterpolationType(TRSLinearInterpolator::InterpType::THermite);
    
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
        //        double fwv = ((sw * sw)/(1.0 + 2.0*(sw * sw - sw)));
        //        double dfw_dswv =(-2.0*(-1.0 + sw) * sw)/(pow(1.0 + 2.0 * (-1.0 + sw) * sw,2.0));
        
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
        //        double fov=  pow(-1.0 + sw,2.0)/(1.0 + 2.0*(-1.0 + sw)*sw);
        //        double dfo_dswv = (2.0*(-1.0 + sw)*sw)/pow(1.0 + 2.0*(-1.0 + sw)*sw,2.0);
        
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
        //        double lv = pow(1.0 - sw,2.0) + pow(sw,2.0);
        //        double dl_dswv = -2.0 + 4.0*sw;
        
        if (sw>1.0 || sw<0.0) {
            DebugStop();
        }
        
        
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
    sim_data.mTNumerics.m_n_steps = 200;
    sim_data.mTNumerics.m_dt      = 0.0005;
    
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    scalnames.Push("Sw_exact");
    sim_data.mTPostProcess.m_file_time_step = 0.0005;
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


void BuckleyLeverett(REAL swf, REAL x_front, TPZVec<REAL> &pt,TPZVec<REAL> q, REAL phi,REAL time, TPZVec<REAL> &Saturation, TMRSDataTransfer sim_data){
    
    
    
    
    REAL val_x = pt[0];
    if (pt[0] > x_front) {
        Saturation[0] =0.0;
        return;
    }
    if (pt[0] < x_front){
        Saturation[0] = S_Newton(pt[0], time, q[0], 0, 0, 0.1, swf, 1.0, 1.0, 1.0, 1.0, 0.001);
        double swval =Saturation[0];
        if (pt[0]==0.0) {
            Saturation[0] = 1.0;
        }
        
        return;
    }
    
    //
}
void BuckleyLeverettCase(TMRSDataTransfer sim_data){
    bool is_3D_Q = false;
    bool is_2D_Coarse_Q = true;
    
    std::string geometry_file;
    geometry_file = "gmsh/simple_2D_fine.msh";
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);
    
    aspace.CreateUniformMesh(100,1,1,1);
    
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    
    TPZAnalysis *an = new TPZAnalysis(transport_operator);
    TPZCompMesh *cmesh = transport_operator->MeshVector()[2];
    
    std::ofstream filetxt("bl.txt");
    cmesh->Print(filetxt);
    
    const TPZVec<std::string> scalnames(1);
    const TPZVec<std::string> vecnames;
    scalnames[0] = "state";
    std::string plotfile("buckleyLevertt.vtk");
    //
    TRSLinearInterpolator & Krw = sim_data.mTPetroPhysics.mLayer_Krw_RelPerModel[0];
    TRSLinearInterpolator & Kro = sim_data.mTPetroPhysics.mLayer_Kro_RelPerModel[0];
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & lambda = sim_data.mTMultiphaseFunctions.mLayer_lambda[0];
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & fw = sim_data.mTMultiphaseFunctions.mLayer_fw[0];
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & fo = sim_data.mTMultiphaseFunctions.mLayer_fo[0];
    
    // calculo de la saturacion en el frente
    int n_val = 1000;
    double delta_S = 1.0/n_val;
    TPZFMatrix<double> sw_deriv(n_val,4,0.0);
    REAL swf = 0.0;
    REAL dfswswf = 0.0;
    
    REAL swf_secante = 0;
    for (int isat = 1; isat < n_val; isat++) {
        double current_sw =  delta_S*isat;
        std::tuple<double, double, double> fw_l = fw(Krw,Kro,current_sw,0.0);
        double deriv = (std::get<0>(fw_l))/(current_sw);
        if (  deriv > swf_secante) {
            swf_secante = deriv;
            swf = current_sw;
            dfswswf = (std::get<1>(fw_l));
            
        }
    }
    
    // calculo de la posicion en el frente
    TPZVec<REAL> q(1,1.0);
    REAL phi =0.1;
    
    
    //
    TPZVec<REAL> reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    int n_reportingtimes = reporting_times.size();
    for (int i=0; i< n_reportingtimes; i++) {
        
        double x_front = (q[0] * reporting_times[i] /phi) * dfswswf;
        SetSolution(swf, x_front, cmesh, sim_data, reporting_times[i]);
        an->DefineGraphMesh(2, scalnames, vecnames, plotfile);
        an->PostProcess(0, 2);
        
    }
    
    
    
    
    
    
    return 0;
    
}
void SetSolution(REAL swf, REAL xswf,TPZCompMesh *cmesh,TMRSDataTransfer simdata, REAL sim_time){
    int gdl = cmesh->Solution().Rows();
    TPZFMatrix<STATE> sol2;
    sol2.Resize(gdl, 1);
    
    int nel = cmesh->NElements();
    for (int i=0; i<nel; i++) {
        TPZCompEl *cel = cmesh->Element(i);
        int n_conects = cel->NConnects();
        for(int j=0; j<n_conects; j++){
            int sec_number = cel->Connect(j).fSequenceNumber;

            //TEST
            if(sec_number > -1)
            {
                int64_t pos = cmesh->Block().Position(sec_number);
                TPZFMatrix<REAL> cooridnates;
                TPZVec<REAL> cooridnates_vec(3,0);
                cel->Reference()->NodesCoordinates(cooridnates);
                cooridnates_vec[0] = cooridnates(0,j);
                cooridnates_vec[1] = cooridnates(1,j);
                cooridnates_vec[2] = cooridnates(2,j);
                TPZVec<REAL> Saturation(1,0.0);
                TPZVec<REAL> q(1,1.0);
                
                BuckleyLeverett(swf, xswf,cooridnates_vec, q, 0.1, sim_time, Saturation, simdata);
                REAL sw_val =Saturation[0];
                sol2(pos,0)=Saturation[0];
            }
        }
    }
    
    cmesh->LoadSolution(sol2);
    
    
}
REAL S_Newton(REAL x, REAL t, REAL u, REAL Swr, REAL Sor, REAL phi, REAL s_shok, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta, REAL epsilon)
{
    REAL S_trial = ((1.0-Sor) + s_shok)/2.0;
    REAL jac    = 0.0;
    REAL r      = 1.0;
    REAL delta_S;
    REAL ds;
    REAL S_k = S_trial;
    int max = 20;
    int it = 0;
    
    
    
    while ( fabs(r) > epsilon && it < max){
        
        ds = dfdsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
        
        r = x - (u * t * ds)/(phi * rho_alpha);
        
        jac = - (u * t * df2dsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta))/(phi * rho_alpha);
        
        delta_S = -r/jac;
        S_k += delta_S;
        if (S_k < -1.5 || S_k > 1.5) {
            DebugStop();
        }
        ds = dfdsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
        r = x - (u * t * ds)/(phi * rho_alpha);
        it++;
        
    }
    
    return S_k;
}

REAL dfdsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta)
{
    REAL dfwdSwS;
    
    
    dfwdSwS = (2.0*mu_alpha*mu_beta*rho_alpha*rho_beta*(-1.0 + Sor + Sw)*(Sw - Swr)*(-1.0 + Sor + Swr))/std::pow(mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0) + mu_beta*rho_alpha*std::pow(Sw - Swr,2.0),2.0);
    
    return (dfwdSwS);
}
REAL df2dsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta)
{
    REAL dfw2dSwS2;
    
    
    dfw2dSwS2 = (2.0*mu_alpha*mu_beta*rho_alpha*rho_beta*(-1.0 + Sor + Swr)*(-(mu_beta*rho_alpha*std::pow(Sw - Swr,2.0)*(-3.0 + 3.0*Sor + 2.0*Sw + Swr)) + mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0)*(-1.0 + Sor - 2.0*Sw + 3.0*Swr)))/
    std::pow(mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0) + mu_beta*rho_alpha*std::pow(Sw - Swr,2.0),3.0);
    
    return dfw2dSwS2;
}
