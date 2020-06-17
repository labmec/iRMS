//
//  AlgebraicTransport.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "TPZAlgebraicTransport.h"
#include "pzelmat.h"
/// Default constructor
TPZAlgebraicTransport::TPZAlgebraicTransport(){
    
    
}

/// Copy constructor
TPZAlgebraicTransport::TPZAlgebraicTransport(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
    fNFluxCoefficients =other.fNFluxCoefficients;
    inletmatid=other.inletmatid;
    outletmatid=other.outletmatid;
    interfaceid=other.interfaceid;
    fgravity=other.fgravity;
}

/// Assignement constructor
const TPZAlgebraicTransport & TPZAlgebraicTransport::operator=(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
    inletmatid=other.inletmatid;
    outletmatid=other.outletmatid;
    interfaceid=other.interfaceid;
    fgravity=other.fgravity;
    return *this;
}
TPZAlgebraicTransport::~TPZAlgebraicTransport(){
    
}



void TPZAlgebraicTransport::CalcLambdas(){

}
void TPZAlgebraicTransport::CalcDensities(){

}

double TPZAlgebraicTransport::CalcLambda(double sw, double paverage, double densityo, double densityw){
    double kro = (1 - sw)*(1 - sw);
    double krw = sw*sw;
    double muo = 0.05;
    double muw = 0.01;
    double lambda;
    lambda = ((densityo*kro)/muo) + ((densityw*krw)/muw);
    return lambda;
}

double TPZAlgebraicTransport::CalcDensity(double paverage,double compress, double reff, double pref){
    double density = reff*(1+ compress*(paverage-pref));
    return density;
}



void TPZAlgebraicTransport::BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh)
{

}

void TPZAlgebraicTransport::Contribute(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef){
    ek(0,0) = fCellsData.fVolume[index];
    REAL sat =fCellsData.fSaturation[index];
    REAL satLast =fCellsData.fSaturationLastState[index];
    ef(0) = fCellsData.fVolume[index]*(sat-satLast);
}
void TPZAlgebraicTransport::ContributeInterface(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceid].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[interfaceid].fIntegralFlux[index];
   
    REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
    REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];
    REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    REAL dfwSw_R = fCellsData.fDerivativeWfractionalflow[lr_index.second];
    REAL phi_L = fCellsData.fporosity[lr_index.first];
    REAL phi_R = fCellsData.fporosity[lr_index.second];
    REAL beta =0.0;
    //upwind
    if (fluxint>0.0) {
        beta = 1.0;
    }
    
    ef(0) = +1.0*fdt*(beta*fw_L*(1/phi_L) + (1-beta)*fw_R*(1/phi_R))*fluxint;
    ef(1) = -1.0*fdt*(beta*fw_L*(1/phi_L)  + (1-beta)*fw_R*(1/phi_R))*fluxint;
    
    ek(0,0) = +1.0*dfwSw_L *(1/phi_L)* fdt * beta * fluxint;
    ek(0,1) = +1.0*dfwSw_R*(1/phi_R) * fdt *(1-beta) * fluxint;
    ek(1,0) = -1.0*dfwSw_L*(1/phi_L) * fdt* beta * fluxint;
    ek(1,1) = -1.0*dfwSw_R*(1/phi_R)*fdt*(1-beta)*fluxint;
    
    // Gravity fluxes contribution
    ContributeInterfaceIHU(index, ek, ef);
    
}

void TPZAlgebraicTransport::ContributeInterfaceIHU(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceid].fLeftRightVolIndex[index];
    std::tuple<REAL, REAL, REAL> normal = fInterfaceData[interfaceid].fNormalFaceDirection[index];
    
    std::vector<REAL> n(3,0.0);
    n[0] = std::get<0>(normal);
    n[1] = std::get<1>(normal);
    n[2] = std::get<2>(normal);
    
    REAL g_dot_n = n[0]*fgravity[0]+n[1]*fgravity[1]+n[2]*fgravity[2];

    REAL lambdaL = fCellsData.flambda[lr_index.first];
    REAL lambdaR = fCellsData.flambda[lr_index.second];
    
    REAL fwL = fCellsData.fWaterfractionalflow[lr_index.first];
    REAL fwR = fCellsData.fWaterfractionalflow[lr_index.second];
    REAL foL = fCellsData.fOilfractionalflow[lr_index.first];
    REAL foR = fCellsData.fOilfractionalflow[lr_index.second];
    
    REAL dfwL = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    REAL dfwR = fCellsData.fDerivativeWfractionalflow[lr_index.second];
    REAL dfoL = fCellsData.fDerivativeOfractionalflow[lr_index.first];
    REAL dfoR = fCellsData.fDerivativeOfractionalflow[lr_index.second];
    
    REAL lambda_wL = fwL * lambdaL;
    REAL lambda_wR = fwR * lambdaR;
    REAL lambda_oL = foL * lambdaL;
    REAL lambda_oR = foR * lambdaR;
    
    REAL dlambda_wL = fCellsData.fdlambdawdsw[lr_index.first];
    REAL dlambda_wR = fCellsData.fdlambdawdsw[lr_index.second];
    REAL dlambda_oL = fCellsData.fdlambdaodsw[lr_index.first];
    REAL dlambda_oR = fCellsData.fdlambdaodsw[lr_index.second];
    
    // The upwinding logic should be the same for each function
    std::pair<REAL, std::pair<REAL, REAL>> fstarL = f_star(foL, fwR, dfoL, dfwR, g_dot_n);
    std::pair<REAL, std::pair<REAL, REAL>> fstarR = f_star(fwR, foL, dfoR, dfwL, g_dot_n);
    
    REAL rho_wL = fCellsData.fDensityWater[lr_index.first];
    REAL rho_wR = fCellsData.fDensityWater[lr_index.second];
    REAL rho_oL = fCellsData.fDensityOil[lr_index.first];
    REAL rho_oR = fCellsData.fDensityOil[lr_index.second];
    REAL rho_ratio_wL = ((rho_wL - rho_oL)/(rho_wL - rho_oL));
    REAL rho_ratio_wR = ((rho_wR - rho_oR)/(rho_wR - rho_oR));
    REAL rho_ratio_oL = ((rho_oL - rho_oL)/(rho_wL - rho_oL));
    REAL rho_ratio_oR = ((rho_oR - rho_oR)/(rho_wR - rho_oR));
    
    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starL = lambda_star(lambda_wL, lambda_wR, dlambda_wL, dlambda_wR, g_dot_n, rho_ratio_wL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starR = lambda_star(lambda_wR, lambda_wL, dlambda_wR, dlambda_wL, g_dot_n, rho_ratio_wR);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starL = lambda_star(lambda_oL, lambda_oR, dlambda_oL, dlambda_oR, g_dot_n, rho_ratio_oL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starR = lambda_star(lambda_oR, lambda_oL, dlambda_oR, dlambda_oL, g_dot_n, rho_ratio_oR);
    
    // Harmonic permeability mean
    REAL Kx_L =  fCellsData.fKx[lr_index.first];
    REAL Ky_L =  fCellsData.fKy[lr_index.first];
    REAL Kz_L =  fCellsData.fKz[lr_index.first];
    
    REAL Kx_R =  fCellsData.fKx[lr_index.first];
    REAL Ky_R =  fCellsData.fKy[lr_index.first];
    REAL Kz_R =  fCellsData.fKz[lr_index.first];
    
    
    REAL K_x = 2.0*(Kx_L * Kx_R)/(Kx_L + Kx_R);
    REAL K_y = 2.0*(Ky_L * Ky_R)/(Ky_L + Ky_R);
    REAL K_z = 2.0*(Kz_L * Kz_R)/(Kz_L + Kz_R);
    
    
    // Beacuse we assume diagonal abs. perm tensor
    REAL K_times_g_dot_n = K_x*n[0]*fgravity[0]+K_y*n[1]*fgravity[1]+K_z*n[2]*fgravity[2];
    
    ef(0) += fstarL.first * (lamba_w_starL.first + lamba_o_starL.first) * K_times_g_dot_n*(rho_wL - rho_oL);
    ef(1) += fstarR.first * (lamba_w_starR.first + lamba_o_starR.first) * K_times_g_dot_n*(rho_wR - rho_oR);
    
    REAL dGLdSL = fstarL.second.first * (lamba_w_starL.first + lamba_o_starL.first) + fstarL.first * (lamba_w_starL.second.first + lamba_o_starL.second.first);
    REAL dGLdSR = fstarL.second.second * (lamba_w_starL.first + lamba_o_starL.first) + fstarL.first * (lamba_w_starL.second.second + lamba_o_starL.second.second);
    REAL dGRdSL = fstarR.second.first * (lamba_w_starR.first + lamba_o_starR.first) + fstarR.first * (lamba_w_starR.second.first + lamba_o_starR.second.first);
    REAL dGRdSR = fstarR.second.second * (lamba_w_starR.first + lamba_o_starR.first) + fstarR.first * (lamba_w_starR.second.second + lamba_o_starR.second.second);
    
    ek(0,0) += dGLdSL * K_times_g_dot_n*(rho_wL - rho_oL);
    ek(0,1) += dGLdSR * K_times_g_dot_n*(rho_wL - rho_oL);
    ek(1,0) += dGRdSL * K_times_g_dot_n*(rho_wR - rho_oR);
    ek(1,1) += dGRdSR * K_times_g_dot_n*(rho_wR - rho_oR);
    
}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::f_star(REAL fL, REAL fR, REAL dfL, REAL dfR, REAL g_dot_n){
    REAL fstar, dfstardsL, dfstardsR;
    if( g_dot_n < 0.0){
        fstar = fR * fL;
        dfstardsL = fR * dfL;
        dfstardsR = dfR * fL;
    }else{
        fstar = fL * fR;
        dfstardsL = dfL * fR;
        dfstardsR = fL * dfR;
    }
    return std::make_pair(fstar, std::make_pair(dfstardsL, dfstardsR));
}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::lambda_star(REAL lambda_L, REAL lambda_R, REAL dlambda_L, REAL dlambda_R, REAL g_dot_n, REAL rho_ratio){
    REAL lambda_star, dlambda_starL, dlambda_starR;
    if( g_dot_n < 0.0){
        lambda_star = rho_ratio * lambda_R  + (1-rho_ratio) * lambda_L;
        dlambda_starL = (1-rho_ratio) * dlambda_L;
        dlambda_starR = rho_ratio * dlambda_R;
    }else{
        lambda_star = rho_ratio * lambda_L  + (1-rho_ratio) * lambda_R;
        dlambda_starL = rho_ratio * dlambda_L;
        dlambda_starR = (1-rho_ratio) * dlambda_R;
    }
    return std::make_pair(lambda_star, std::make_pair(dlambda_starL, dlambda_starR));
}

void TPZAlgebraicTransport::ContributeBCInletInterface(int index, TPZFMatrix<double> &ef){
   
    int s_inlet =1.0;
    REAL fluxint  = fInterfaceData[inletmatid].fIntegralFlux[index];
    std::pair<int, int> lr_index = fInterfaceData[inletmatid].fLeftRightVolIndex[index];
    REAL phi = fCellsData.fporosity[lr_index.first];
    ef(0,0) = 1.0*fdt*s_inlet*fluxint*(1/phi);
}
void TPZAlgebraicTransport::ContributeBCOutletInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef){
  
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[outletmatid].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[outletmatid].fIntegralFlux[index];
    REAL fw_L= fCellsData.fWaterfractionalflow[lr_index.first];
    REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    REAL phi = fCellsData.fporosity[lr_index.first];
    ef(0,0) = 1.0*fdt*fw_L*fluxint*(1.0/phi);
    ek(0,0) = dfwSw_L*fdt*fluxint*(1.0/phi);
}

void TPZAlgebraicTransport::TCellData::Print(std::ostream &out){
    
    int nels = this->fVolume.size();
    for (int iel =0; iel< nels; iel++) {
        out<<"Material_Id: "<<this->fMatId<<std::endl;
        out<<"Center Cord: ";
        for (int ic=0 ; ic<fCenterCordinate[iel].size(); ic++) {
            out<<fCenterCordinate[iel][ic]<<" ";
        }
        out<<std::endl;
        out<<"Volume : "<<this->fVolume[iel]<<std::endl;
        out<<"Oil density: "<<this->fDensityOil[iel]<<std::endl;
        out<<"Water density: "<<this->fDensityWater[iel]<<std::endl;
        out<<"Pressure: "<<this->fPressure[iel]<<std::endl;
        out<<"Saturation: "<<this->fSaturation[iel]<<std::endl;
        out<<"WaterFracFlow: "<<this->fWaterfractionalflow[iel]<<std::endl;
        out<<"dWaterFracFlowdSw: "<<this->fDerivativeWfractionalflow[iel]<<std::endl;
        out<<"OilFracFlow: "<<this->fOilfractionalflow[iel]<<std::endl;
        out<<"dOilFracFlowdSw: "<<this->fDerivativeOfractionalflow[iel]<<std::endl;
        out<<"Lambda: "<<this->flambda[iel]<<std::endl;
        out<<std::endl;
    }

}

std::pair<std::vector<REAL>,std::vector<REAL>> TPZAlgebraicTransport::fwAndfoVal(REAL sw, REAL muw,REAL muo, bool isLinearQ){
    std::vector<REAL> fwData(2), foData(2);
//    REAL Krw = sw*sw ;
//    REAL Kro = (1-sw)*(1-sw) ;
    
    REAL fw = (muo*sw*sw)/(muw*(sw-1.0)*(sw-1.0) + (muo*sw*sw));
    REAL num = -2.0*(muo*muw*(sw-1.0)*sw);
    REAL dem = ((muw*(sw-1.0)*(sw-1.0))+(muo*sw*sw))*((muw*(sw-1.0)*(sw-1.0))+(muo*sw*sw));
    fwData[0] = fw;
    if(isLinearQ){
        fwData[1]=(muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
    }else{
        fwData[1] = num/dem;
    }

    foData[0] = 1.0;
    foData[1]=1.0;
  
    std::pair<std::vector<REAL>,std::vector<REAL>> fracflows = std::make_pair(fwData, foData);
    return fracflows;
}

std::pair<std::vector<REAL>,std::vector<REAL>> TPZAlgebraicTransport::LinearfwAndfoVal(REAL sw, REAL muw,REAL muo){
    std::vector<REAL> fwData(2), foData(2);
    
    REAL Krw = sw;
    REAL Kro = (1-sw);
    REAL fw = (Krw/muw)/((Krw/muw)+(Kro/muo));
    fwData[0] = fw;
    fwData[1] = (muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
    foData[0] = 1.0;
    foData[1]=1.0;
    
    std::pair<std::vector<REAL>,std::vector<REAL>> fracflows = std::make_pair(fwData, foData);
    return fracflows;
    
}
void TPZAlgebraicTransport::UpdateIntegralFlux(int matid){
    
    if(fInterfaceData.find(matid) == fInterfaceData.end()) DebugStop();
    int nels = fInterfaceData[matid].fCoefficientsFlux.size();
    if(nels == 0) DebugStop();
    fInterfaceData[matid].fIntegralFlux= fInterfaceData[matid].fCoefficientsFlux[0];
    for (int i = 1; i<nels; i++) {
        int np =fInterfaceData[matid].fCoefficientsFlux[i].size();
        for (int index = 0; index<np; index++) {
            REAL val = fInterfaceData[matid].fCoefficientsFlux[i][index];
            fInterfaceData[matid].fIntegralFlux[index] +=val;
        }

    }

    
}
void TPZAlgebraicTransport::TInterfaceDataTransport::Print(std::ostream &out){
    
    int ninterfaces = this->fFluxSign.size();
    out<<"Material_ID: "<<this->fMatid<<std::endl;
    for (int iinter = 0 ; iinter <ninterfaces; iinter++) {
        out<<"Left_Index: "<<this->fLeftRightVolIndex[iinter].first<<std::endl;
        out<<"Right_Index: "<<this->fLeftRightVolIndex[iinter].second<<std::endl;
        int ncflux = fCoefficientsFlux.size();
        out << "fCoefficientsFlux :";
        for (int icoe=0; icoe<fCoefficientsFlux.size(); icoe++) {
            out<<fCoefficientsFlux[icoe][iinter]<<" ";
        }
        out<<std::endl;
        out<<"IntegralFluxFunctions: "<<fIntegralFluxFunctions[iinter]<<std::endl;
        out<<"IntegralFlux: "<<fIntegralFlux[iinter]<<std::endl;
        out<<"fFluxSign: "<<fFluxSign[iinter]<<std::endl;
        out<<"fNormalDirection: ";
        out<<std::get<0>(fNormalFaceDirection[iinter])<<" ";
        out<<std::get<1>(fNormalFaceDirection[iinter])<<" ";
        out<<std::get<2>(fNormalFaceDirection[iinter])<<" ";
        out<<std::endl;
    }
   

}
void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambda(bool isLinearQ){
    int nvols = this->fVolume.size();
    for (int ivol =0 ; ivol< nvols; ivol++) {
        std::pair<std::vector<REAL>,std::vector<REAL>> fWandFo;
        REAL sw = this->fSaturation[ivol];
        REAL Krw, Kro, fw, dfwdsw, fo, dfodsw, dlwdsw,dlodsw;
        REAL muw = fViscosity[0];
        REAL muo = fViscosity[1];
        if (isLinearQ) {
            Krw = sw;
            Kro = (1-sw);
            fw = (Krw/muw)/((Krw/muw)+(Kro/muo));
            fo = (Kro/muo)/((Krw/muw)+(Kro/muo));
            dfwdsw = (muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
            dfodsw = -1.0*dfwdsw;
            dlwdsw = (1.0/muw) ;
            dlodsw = -1.0*(1.0/muo);
        }
        else{
            Krw = sw*sw;
            Kro = (1-sw)*(1-sw);
            fw =(muw*Krw)/((muo*Kro) + (muw*Krw));
            fo =(muo*Kro)/((muo*Kro) + (muw*Krw));
            REAL num = -2.0*(muo*muw*(sw-1.0)*sw);
            REAL dem = ((muw*(sw-1.0)*(sw-1.0))+(muo*sw*sw))*((muw*(sw-1.0)*(sw-1.0))+(muo*sw*sw));
            dfwdsw = num/dem;
            dfodsw = -1.0*dfwdsw;
            
            dlwdsw = (2.0*sw/muw) ;
            dlodsw =(-2*(1.0-sw)/muo);
        }
        this->fWaterfractionalflow[ivol] = fw;
        this->fDerivativeWfractionalflow[ivol] =dfwdsw;
        this->fOilfractionalflow[ivol] = fo;
        this->fDerivativeOfractionalflow[ivol] = dfodsw;
        this->flambda[ivol] = (Krw/muw)+(Kro/muo);
        this->fdlambdawdsw[ivol] = dlwdsw;
        this->fdlambdaodsw[ivol] = dlodsw;
     

    }
}

void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambdaQuasiNewton(){
    int nvols = this->fVolume.size();
    REAL muw= fViscosity[0];
    REAL muo= fViscosity[1];
    for (int ivol =0 ; ivol< nvols; ivol++) {
        std::pair<std::vector<REAL>,std::vector<REAL>> fWandFo;
        REAL sw =this->fSaturation[ivol];
        REAL krw, kro, fw, fo,dfwdsw, dlambdadsw;
        krw = sw*sw;
        kro = (1-sw)*(1-sw);
        fw = (muw*krw)/((muw*kro) + (muo*krw));
        fo = (muo*kro)/((muw*kro) + (muo*krw));
        dfwdsw =(muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
        dlambdadsw = (2.0*sw/muw) + (-2.0*(1.0-sw)/muo);
        this->fWaterfractionalflow[ivol] = fw;
        this->fDerivativeWfractionalflow[ivol] =dfwdsw;
        this->fOilfractionalflow[ivol] =fo;
        this->fDerivativeOfractionalflow[ivol] = -1.0*dfwdsw;
        this->flambda[ivol] = (krw/(fViscosity[0]))+(kro/(fViscosity[1]));
//        this->fdlambdadsw[ivol] = dlambdadsw;
    }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturations(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        fSaturation[icell] = sw(eq_number);
       
    }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturationsTo(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        sw(eq_number) = fSaturation[icell];
       
    }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturationsLastState(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        fSaturationLastState[icell] = sw(eq_number);
    }
}
void TPZAlgebraicTransport::TCellData::UpdateMixedDensity(){
    int ncells = fVolume.size();
    for (int i =0; i< ncells; i++) {
        REAL mixedDen = (fWaterfractionalflow[i]*fDensityWater[i])+fOilfractionalflow[i]*fDensityOil[i];
//        std::cout<<"fw: "<<fWaterfractionalflow[i]<<std::endl;
//        std::cout<<"fo: "<<fOilfractionalflow[i]<<std::endl;
//        std::cout<<"rhow: "<<fDensityWater[i]<<std::endl;
//        std::cout<<"rhoo: "<<fDensityOil[i]<<std::endl;
        fMixedDensity[i] = mixedDen;
    }
}
