//
//  AlgebraicTransport.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "TPZAlgebraicTransport.h"
#include "pzelmat.h"
//#include <Eigen/PardisoSupport>

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
    fHasPropQ = other.fHasPropQ;
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
    fHasPropQ = other.fHasPropQ;
    return *this;
}
TPZAlgebraicTransport::~TPZAlgebraicTransport(){
    
}

void TPZAlgebraicTransport::BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh)
{

}

void TPZAlgebraicTransport::Contribute(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef){

//    REAL sat =fCellsData.fSaturation[index];
//    REAL satLast = fCellsData.fSaturationLastState[index];
//    REAL phi = fCellsData.fporosity[index];
//    ef(0) = fCellsData.fVolume[index]*phi*(sat-satLast)/fdt;
//    ek(0,0) = fCellsData.fVolume[index]*phi/fdt;
    
    REAL sat =fCellsData.fSaturation[index];
    REAL satLast = fCellsData.fSaturationLastState[index];
    REAL phi = fCellsData.fporosity[index];
    ef(0) = fCellsData.fVolume[index]*phi*(sat-satLast);
    ek(0,0) = fCellsData.fVolume[index]*phi;
}

void TPZAlgebraicTransport::ContributeResidual(int index, TPZFMatrix<double> &ef){

//    REAL sat =fCellsData.fSaturation[index];
//    REAL satLast = fCellsData.fSaturationLastState[index];
//    REAL phi = fCellsData.fporosity[index];
//    ef(0) = fCellsData.fVolume[index]*phi*(sat-satLast)/fdt;
    
    REAL sat =fCellsData.fSaturation[index];
    REAL satLast = fCellsData.fSaturationLastState[index];
    REAL phi = fCellsData.fporosity[index];
    ef(0) = fCellsData.fVolume[index]*phi*(sat-satLast);
}

void TPZAlgebraicTransport::ContributeInterface(int index, TPZFMatrix<double> &ek,TPZFMatrix<double> &ef, int interfaceId){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceId].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[interfaceId].fIntegralFlux[index];
    REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
    REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];
    REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    REAL dfwSw_R = fCellsData.fDerivativeWfractionalflow[lr_index.second];
    
    REAL beta =0.0;
    //upwind
//    if (interfaceId==103 || interfaceId==102){
//        REAL Aux =fw_L;
//        fw_L =fw_R;
//        fw_R =Aux;
//    }
    if (fluxint>0.0) {
        beta = 1.0;
    }
    
    ef(0) = +1.0*(beta*fw_L + (1-beta)*fw_R)*fluxint * fdt;
    ef(1) = -1.0*(beta*fw_L  + (1-beta)*fw_R)*fluxint* fdt;
    
    ek(0,0) = +1.0*dfwSw_L  * beta * fluxint* fdt;
    ek(0,1) = +1.0*dfwSw_R * (1-beta) * fluxint* fdt;
    ek(1,0) = -1.0*dfwSw_L* beta * fluxint * fdt;
    ek(1,1) = -1.0*dfwSw_R * (1-beta)*fluxint* fdt;
   
    // Gravity fluxes contribution
//    ContributeInterfaceIHU(index, ek, ef);
    
}

void TPZAlgebraicTransport::ContributeInterfaceResidual(int index, TPZFMatrix<double> &ef, int interfaceID){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceID].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[interfaceID].fIntegralFlux[index];
   
    REAL fw_L = fCellsData.fWaterfractionalflow[lr_index.first];
    REAL fw_R = fCellsData.fWaterfractionalflow[lr_index.second];
    
    REAL beta =0.0;
    //upwind
//    if(interfaceID==102){
//        beta=1;
//    }
  
    if (fluxint>0.0) {
        beta = 1.0;
    }
    
    ef(0) = +1.0*(beta*fw_L + (1-beta)*fw_R)*fluxint* fdt;
    ef(1) = -1.0*(beta*fw_L  + (1-beta)*fw_R)*fluxint* fdt;
    
    // Gravity fluxes contribution
//    ContributeInterfaceIHUResidual(index, ef);
    
#ifdef PZDEBUG
    if(std::isnan(Norm(ef)))
    {
        std::cout << __PRETTY_FUNCTION__ << " nan" << std::endl;
    }
#endif
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
    
    std::pair<REAL, REAL> fwL = {fCellsData.fWaterfractionalflow[lr_index.first], fCellsData.fDerivativeWfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> fwR = {fCellsData.fWaterfractionalflow[lr_index.second],
        fCellsData.fDerivativeWfractionalflow[lr_index.second]};
    std::pair<REAL, REAL> foL = {fCellsData.fOilfractionalflow[lr_index.first],
            fCellsData.fDerivativeOfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> foR = {fCellsData.fOilfractionalflow[lr_index.second],
                fCellsData.fDerivativeOfractionalflow[lr_index.second]};
    
    std::pair<REAL, REAL> lambda_wL = {fwL.first * lambdaL, fCellsData.fdlambdawdsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_wR = {fwR.first * lambdaR, fCellsData.fdlambdawdsw[lr_index.second]};
    std::pair<REAL, REAL> lambda_oL = {foL.first * lambdaL, fCellsData.fdlambdaodsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_oR = {foR.first * lambdaR, fCellsData.fdlambdaodsw[lr_index.second]};
    
    REAL rho_wL = fCellsData.fDensityWater[lr_index.first];
    REAL rho_wR = fCellsData.fDensityWater[lr_index.second];
    REAL rho_oL = fCellsData.fDensityOil[lr_index.first];
    REAL rho_oR = fCellsData.fDensityOil[lr_index.second];
    
    // The upwinding logic should be the same for each function
    std::pair<REAL, std::pair<REAL, REAL>> fstarL = f_star(foL, foR, fwL, fwR, g_dot_n);
    std::pair<REAL, std::pair<REAL, REAL>> fstarR = f_star(foR, foL, fwR, fwL, -g_dot_n);
    
//    REAL rho_ratio_wL = ((rho_wL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_wR = ((rho_wR - rho_oR)/(rho_wR - rho_oR));
//    REAL rho_ratio_oL = ((rho_oL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_oR = ((rho_oR - rho_oR)/(rho_wR - rho_oR));
    REAL rho_ratio_wL = 1.;
    REAL rho_ratio_wR = 1.;
    REAL rho_ratio_oL = 0.;
    REAL rho_ratio_oR = 0.;

    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starL = lambda_w_star(lambda_wL, lambda_wR, g_dot_n, rho_ratio_wL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starR = lambda_w_star(lambda_wR, lambda_wL, -g_dot_n, rho_ratio_wR);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starL = lambda_o_star(lambda_oL, lambda_oR, g_dot_n, rho_ratio_oL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starR = lambda_o_star(lambda_oR, lambda_oL, -g_dot_n, rho_ratio_oR);
    
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
    REAL K_times_g_dot_n = (K_x*n[0]*fgravity[0]+K_y*n[1]*fgravity[1]+K_z*n[2]*fgravity[2]);
    
    REAL res1 = fstarL.first * (lamba_w_starL.first + lamba_o_starL.first) * K_times_g_dot_n * (rho_wL - rho_oL);
    REAL res2 = fstarR.first * (lamba_w_starR.first + lamba_o_starR.first) * K_times_g_dot_n * (rho_wR - rho_oR);
    ef(0) += res1;
    ef(1) -= res2;
    
    REAL dGLdSL = fstarL.second.first * (lamba_w_starL.first + lamba_o_starL.first) + fstarL.first * (lamba_w_starL.second.first + lamba_o_starL.second.first);
    REAL dGLdSR = fstarL.second.second * (lamba_w_starL.first + lamba_o_starL.first) + fstarL.first * (lamba_w_starL.second.second + lamba_o_starL.second.second);
    REAL dGRdSL = fstarR.second.first * (lamba_w_starR.first + lamba_o_starR.first) + fstarR.first * (lamba_w_starR.second.first + lamba_o_starR.second.first);
    REAL dGRdSR = fstarR.second.second * (lamba_w_starR.first + lamba_o_starR.first) + fstarR.first * (lamba_w_starR.second.second + lamba_o_starR.second.second);
    
    ek(0,0) += dGLdSL * K_times_g_dot_n * (rho_wL - rho_oL);
    ek(0,1) += dGLdSR * K_times_g_dot_n * (rho_wL - rho_oL);
    
    ek(1,1) -= dGRdSL * K_times_g_dot_n * (rho_wR - rho_oR);
    ek(1,0) -= dGRdSR * K_times_g_dot_n * (rho_wR - rho_oR);
}

void TPZAlgebraicTransport::ContributeInterfaceIHUResidual(int index, TPZFMatrix<double> &ef){
    
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[interfaceid].fLeftRightVolIndex[index];
    std::tuple<REAL, REAL, REAL> normal = fInterfaceData[interfaceid].fNormalFaceDirection[index];
    
    std::vector<REAL> n(3,0.0);
    n[0] = std::get<0>(normal);
    n[1] = std::get<1>(normal);
    n[2] = std::get<2>(normal);
    
    REAL g_dot_n = n[0]*fgravity[0]+n[1]*fgravity[1]+n[2]*fgravity[2];
    
    REAL lambdaL = fCellsData.flambda[lr_index.first];
    REAL lambdaR = fCellsData.flambda[lr_index.second];
    
    std::pair<REAL, REAL> fwL = {fCellsData.fWaterfractionalflow[lr_index.first], fCellsData.fDerivativeWfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> fwR = {fCellsData.fWaterfractionalflow[lr_index.second],
        fCellsData.fDerivativeWfractionalflow[lr_index.second]};
    std::pair<REAL, REAL> foL = {fCellsData.fOilfractionalflow[lr_index.first],
            fCellsData.fDerivativeOfractionalflow[lr_index.first]};
    std::pair<REAL, REAL> foR = {fCellsData.fOilfractionalflow[lr_index.second],
                fCellsData.fDerivativeOfractionalflow[lr_index.second]};
    
    std::pair<REAL, REAL> lambda_wL = {fwL.first * lambdaL, fCellsData.fdlambdawdsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_wR = {fwR.first * lambdaR, fCellsData.fdlambdawdsw[lr_index.second]};
    std::pair<REAL, REAL> lambda_oL = {foL.first * lambdaL, fCellsData.fdlambdaodsw[lr_index.first]};
    std::pair<REAL, REAL> lambda_oR = {foR.first * lambdaR, fCellsData.fdlambdaodsw[lr_index.second]};
    
    REAL rho_wL = fCellsData.fDensityWater[lr_index.first];
    REAL rho_wR = fCellsData.fDensityWater[lr_index.second];
    REAL rho_oL = fCellsData.fDensityOil[lr_index.first];
    REAL rho_oR = fCellsData.fDensityOil[lr_index.second];
    
    // The upwinding logic should be the same for each function
    std::pair<REAL, std::pair<REAL, REAL>> fstarL = f_star(foL, foR, fwL, fwR, g_dot_n);
    std::pair<REAL, std::pair<REAL, REAL>> fstarR = f_star(foR, foL, fwR, fwL, -g_dot_n);
    
//    REAL rho_ratio_wL = ((rho_wL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_wR = ((rho_wR - rho_oR)/(rho_wR - rho_oR));
//    REAL rho_ratio_oL = ((rho_oL - rho_oL)/(rho_wL - rho_oL));
//    REAL rho_ratio_oR = ((rho_oR - rho_oR)/(rho_wR - rho_oR));
    REAL rho_ratio_wL = 1.;
    REAL rho_ratio_wR = 1.;
    REAL rho_ratio_oL = 0.;
    REAL rho_ratio_oR = 0.;

    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starL = lambda_w_star(lambda_wL, lambda_wR, g_dot_n, rho_ratio_wL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_w_starR = lambda_w_star(lambda_wR, lambda_wL, -g_dot_n, rho_ratio_wR);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starL = lambda_o_star(lambda_oL, lambda_oR, g_dot_n, rho_ratio_oL);
    std::pair<REAL, std::pair<REAL, REAL>> lamba_o_starR = lambda_o_star(lambda_oR, lambda_oL, -g_dot_n, rho_ratio_oR);
    
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
    REAL K_times_g_dot_n = (K_x*n[0]*fgravity[0]+K_y*n[1]*fgravity[1]+K_z*n[2]*fgravity[2]);
    
    REAL res1 = fstarL.first * (lamba_w_starL.first + lamba_o_starL.first) * K_times_g_dot_n * (rho_wL - rho_oL);
    REAL res2 = fstarR.first * (lamba_w_starR.first + lamba_o_starR.first) * K_times_g_dot_n * (rho_wR - rho_oR);
    ef(0) += res1;
    ef(1) -= res2;

}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::f_star(std::pair<REAL, REAL> foL, std::pair<REAL, REAL> foR, std::pair<REAL, REAL> fwL, std::pair<REAL, REAL> fwR, REAL g_dot_n){
    REAL fstar, dfstardsL, dfstardsR;
    if( g_dot_n < 0.0){
        fstar = foL.first * fwR.first;
        dfstardsL = foL.second * fwR.first;
        dfstardsR = foL.first * fwR.second;
    }else{
        fstar = foR.first * fwL.first;
        dfstardsL = foR.first * fwL.second;
        dfstardsR = foR.second * fwL.first;
    }
    return std::make_pair(fstar, std::make_pair(dfstardsL, dfstardsR));
}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::lambda_o_star(std::pair<REAL, REAL> lambda_L, std::pair<REAL, REAL> lambda_R, REAL g_dot_n, REAL rho_ratio){
    REAL lambda_star, dlambda_starL, dlambda_starR;
    if( g_dot_n < 0.0){
        lambda_star = rho_ratio * lambda_L.first  + (1-rho_ratio) * lambda_R.first;
        dlambda_starL = rho_ratio * lambda_L.second;
        dlambda_starR = (1-rho_ratio) * lambda_R.second;
    }else{
        lambda_star = rho_ratio * lambda_R.first  + (1-rho_ratio) * lambda_L.first;
        dlambda_starL = (1-rho_ratio) * lambda_L.second;
        dlambda_starR = rho_ratio * lambda_R.second;
    }
    return std::make_pair(lambda_star, std::make_pair(dlambda_starL, dlambda_starR));
}

std::pair<REAL, std::pair<REAL, REAL>> TPZAlgebraicTransport::lambda_w_star(std::pair<REAL, REAL> lambda_L, std::pair<REAL, REAL> lambda_R, REAL g_dot_n, REAL rho_ratio){
    REAL lambda_star, dlambda_starL, dlambda_starR;
    if( g_dot_n > 0.0){
        lambda_star = rho_ratio * lambda_L.first  + (1-rho_ratio) * lambda_R.first;
        dlambda_starL = rho_ratio * lambda_L.second;
        dlambda_starR = (1-rho_ratio) * lambda_R.second;
    }else{
        lambda_star = rho_ratio * lambda_R.first  + (1-rho_ratio) * lambda_L.first;
        dlambda_starL = (1-rho_ratio) * lambda_L.second;
        dlambda_starR = rho_ratio * lambda_R.second;
    }
    return std::make_pair(lambda_star, std::make_pair(dlambda_starL, dlambda_starR));
}

void TPZAlgebraicTransport::ContributeBCInletInterface(int index, TPZFMatrix<double> &ef, int inId){
   
    REAL s_inlet = 0.01;
    REAL fluxint  = fInterfaceData[inId].fIntegralFlux[index];
    ef(0,0) = 1.0*s_inlet*fluxint* fdt;
}
void TPZAlgebraicTransport::ContributeBCOutletInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef, int outID){
  
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[outID].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[outID].fIntegralFlux[index];
    REAL fw_L= fCellsData.fWaterfractionalflow[lr_index.first];
    REAL dfwSw_L = fCellsData.fDerivativeWfractionalflow[lr_index.first];
    ef(0,0) = 1.0*fw_L*fluxint* fdt;
    ek(0,0) = dfwSw_L*fluxint* fdt;
}

void TPZAlgebraicTransport::ContributeBCOutletInterfaceResidual(int index, TPZFMatrix<double> &ef){
  
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[outletmatid].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[outletmatid].fIntegralFlux[index];
    REAL fw_L= fCellsData.fWaterfractionalflow[lr_index.first];
    ef(0,0) = 1.0*fw_L*fluxint* fdt;
}

void TPZAlgebraicTransport::TCellData::Print(std::ostream &out){
    
    int nels = this->fVolume.size();
    for (int iel =0; iel< nels; iel++) {
        out<<"Material_Id: "<<this->fMatId<<std::endl;
        out<<"Center Cord: ";
        for (int ic=0 ; ic<fCenterCoordinate[iel].size(); ic++) {
            out<<fCenterCoordinate[iel][ic]<<" ";
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

void TPZAlgebraicTransport::UpdateIntegralFlux(int matid){
    
    if(fInterfaceData.find(matid) == fInterfaceData.end()) DebugStop();
    int nels = fInterfaceData[matid].fCoefficientsFlux.size();
    if(nels == 0) return;
    std::vector<REAL> val =fInterfaceData[matid].fCoefficientsFlux[0];
    fInterfaceData[matid].fIntegralFlux=val;
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
void TPZAlgebraicTransport::TCellData::SetDataTransfer(TMRSDataTransfer *simdata){
    fsim_data = simdata;
}
void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambda(bool isLinearQ){
    
    if (!isLinearQ) {
        fsim_data->mTPetroPhysics.CreateQuadraticKrModel();
    }
    else{
        fsim_data->mTPetroPhysics.CreateLinearKrModel();
    }
    
    auto labdaWf = fsim_data->mTPetroPhysics.mLambdaW;
    auto labdaOf = fsim_data->mTPetroPhysics.mLambdaO;
    auto lambdaTotalf = fsim_data->mTPetroPhysics.mLambdaTotal;
    auto fwf = fsim_data->mTPetroPhysics.mFw;
    auto fof = fsim_data->mTPetroPhysics.mFo;
    
    int nvols = this->fVolume.size();
    for (int ivol =0 ; ivol< nvols; ivol++) {

        REAL sw = this->fSaturation[ivol];
        auto fwfvalderiv = fwf(sw);
        auto fovalderiv  = fof(sw);
        auto lambdaWvalderiv = labdaWf(sw);
        auto lambdaOvalderiv = labdaOf(sw);
        auto lambdaTotalvalderiv = lambdaTotalf(sw);
        
        this->fWaterfractionalflow[ivol] = std::get<0>(fwfvalderiv);
        this->fDerivativeWfractionalflow[ivol] = std::get<1>(fwfvalderiv);
        this->fOilfractionalflow[ivol] = std::get<0>(fovalderiv);
        this->fDerivativeOfractionalflow[ivol] = std::get<1>(fovalderiv);
        this->flambda[ivol] = std::get<0>(lambdaTotalvalderiv);
        this->fdlambdawdsw[ivol] =std::get<1>(lambdaWvalderiv);
        this->fdlambdaodsw[ivol] = std::get<1>(lambdaOvalderiv);
    
    }
    
}

void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambdaQuasiNewton(){
    
        fsim_data->mTPetroPhysics.CreateQuadraticKrModel();
        int nvols = this->fVolume.size();
        auto labdaWf = fsim_data->mTPetroPhysics.mLambdaW;
        auto labdaOf = fsim_data->mTPetroPhysics.mLambdaO;
        auto lambdaTotalf = fsim_data->mTPetroPhysics.mLambdaTotal;
        auto fwf = fsim_data->mTPetroPhysics.mFw;
        auto fof = fsim_data->mTPetroPhysics.mFo;
    
        fsim_data->mTPetroPhysics.CreateQuadraticKrModel();
        auto fwLinearf =fsim_data->mTPetroPhysics.mFo;
    
        for (int ivol =0 ; ivol< nvols; ivol++) {
            
            REAL sw = this->fSaturation[ivol];
            auto fwfvalderiv = fwf(sw);
            auto fwflinearvalderiv = fwLinearf(sw);
            auto fovalderiv  = fof(sw);
            auto lambdaWvalderiv = labdaWf(sw);
            auto lambdaOvalderiv = labdaOf(sw);
            auto lambdaTotalvalderiv = lambdaTotalf(sw);
            
            
            this->fWaterfractionalflow[ivol] = std::get<0>(fwfvalderiv);
            this->fOilfractionalflow[ivol] = std::get<0>(fovalderiv);
            this->fDerivativeOfractionalflow[ivol] = std::get<1>(fovalderiv);
            this->flambda[ivol] = std::get<0>(lambdaTotalvalderiv);
            this->fdlambdawdsw[ivol] =std::get<1>(lambdaWvalderiv);
            this->fdlambdaodsw[ivol] = std::get<1>(lambdaOvalderiv);
            
            fsim_data->mTPetroPhysics.CreateLinearKrModel();
            
            this->fDerivativeWfractionalflow[ivol] =std::get<1>(fwflinearvalderiv);
        }
//    int nvols = this->fVolume.size();
//    REAL muw= fViscosity[0];
//    REAL muo= fViscosity[1];
//    for (int ivol =0 ; ivol< nvols; ivol++) {
//        std::pair<std::vector<REAL>,std::vector<REAL>> fWandFo;
//        REAL sw =this->fSaturation[ivol];
//        REAL krw, kro, fw, fo,dfwdsw, dlambdadsw, dlambdaodsw;
//        krw = sw*sw;
//        kro = (1-sw)*(1-sw);
//        fw = (muo*krw)/((muw*kro) + (muo*krw));
//        fo = (muw*kro)/((muw*kro) + (muo*krw));
//        dfwdsw =(muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
//        dlambdadsw = (2.0*sw/muw) ;
//        dlambdaodsw =  (-2.0*(1.0-sw)/muo);
//        this->fWaterfractionalflow[ivol] = fw;
//        this->fDerivativeWfractionalflow[ivol] =dfwdsw;
//        this->fOilfractionalflow[ivol] =fo;
//        this->fDerivativeOfractionalflow[ivol] = -1.0*dfwdsw;
//        this->flambda[ivol] = (krw/muw)+(kro/muo);
//        this->fdlambdawdsw[ivol] = dlambdadsw;
//        this->fdlambdaodsw[ivol] = dlambdaodsw;
//    }
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
        fMixedDensity[i] = mixedDen;
    }
}
void TPZAlgebraicTransport::TCellData::UpdateDensities(){
    UpdateDensitiesLastState();
    int ncells = fVolume.size();
    auto fWaterDensityF = fsim_data->mTFluidProperties.mWaterDensityF;
    auto fOilDensityF = fsim_data->mTFluidProperties.mOilDensityF;
    if(fWaterDensityF){
    for (int icell = 0; icell<ncells; icell++) {
        REAL pressure = fPressure[icell];
        auto densityWvalderiv = fWaterDensityF(pressure);
        auto densityOvalderiv = fOilDensityF(pressure);
        fDensityWater[icell] = std::get<0>(densityWvalderiv);
        fdDensityWaterdp[icell]= std::get<1>(densityWvalderiv);
        fDensityOil[icell] = std::get<0>(densityOvalderiv);
        fdDensityOildp[icell] = std::get<1>(densityOvalderiv);
    }
        
    }
    else{
        for (int icell = 0; icell<ncells; icell++) {
            fDensityWater[icell] = fReferenceDensity[0];
            fDensityOil[icell] = fReferenceDensity[1];
        }
    }
    
}
void TPZAlgebraicTransport::TCellData::UpdateDensitiesLastState(){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        fDensityWaterLastState[icell] = fDensityWater[icell];
        fDensityOilLastState[icell] = fDensityOil[icell];
    }
    
}



REAL TPZAlgebraicTransport::CalculateMass(){
    int ncells = fCellsData.fVolume.size();
    REAL intMass = 0.0;
    for (int icel = 0; icel < ncells; icel++) {
        REAL sat = fCellsData.fSaturation[icel];
        REAL phi = fCellsData.fporosity[icel];
        REAL vol = fCellsData.fVolume[icel];
        intMass += sat*phi*vol;
    }
    return intMass;
}
std::pair<REAL, REAL> TPZAlgebraicTransport::FLuxWaterOilIntegralbyID(int mat_id){
    
    REAL WaterIntegral =0.0;
    REAL OilIntegral = 0.0;
    int ninter = fInterfaceData[mat_id].fIntegralFlux.size();
    for (int iface =0 ; iface < ninter; iface++) {
        int LeftElIndex = fInterfaceData[mat_id].fLeftRightVolIndex[iface].first;
        REAL fracFluxWater = fCellsData.fWaterfractionalflow[LeftElIndex];
        REAL fracFluxOil = fCellsData.fOilfractionalflow[LeftElIndex];
        REAL FluxInttegral = fInterfaceData[mat_id].fIntegralFlux[iface];
        WaterIntegral += FluxInttegral*fracFluxWater;
        OilIntegral += FluxInttegral*fracFluxOil;
        
    }
    return std::make_pair(WaterIntegral, OilIntegral);
}
REAL TPZAlgebraicTransport::FLuxIntegralbyID(int mat_id){
    
    REAL FluxInt =0.0;
    
    int ninter = fInterfaceData[mat_id].fIntegralFlux.size();
    for (int iface =0 ; iface < ninter; iface++) {
        int LeftElIndex = fInterfaceData[mat_id].fLeftRightVolIndex[iface].first;
        REAL FluxInttegral = fInterfaceData[mat_id].fIntegralFlux[iface];
        FluxInt += FluxInttegral;
    }
    return FluxInt;
}
void TPZAlgebraicTransport::VerifyElementFLuxes(){
    
    int nels = fCellsData.fVolume.size();
    std::vector<int> nInterfacesByElement(nels);
    std::vector<int> numZeroFluxByElement(nels);
    std::vector<REAL> SumFluxByElement(nels);
    for( auto interdata: fInterfaceData){
     
        TInterfaceDataTransport &transport  =interdata.second;
        int neles = transport.fLeftRightVolIndex.size();
        for(int iel = 0; iel<neles; iel++){
            int leftIndex = transport.fLeftRightVolIndex[iel].first;
            int rightIndex = transport.fLeftRightVolIndex[iel].second;
            REAL fluxInt =transport.fIntegralFlux[iel];
            if(abs(fluxInt)<1.0e-10){
                numZeroFluxByElement[leftIndex]++;
                numZeroFluxByElement[rightIndex]++;
            }
            
            SumFluxByElement[leftIndex] += fluxInt;
//            std::cout<<transport.fIntegralFlux[iel]<<std::endl;
            nInterfacesByElement[leftIndex] ++;
            if(rightIndex<0){
                continue;
            }
            nInterfacesByElement[rightIndex] ++;
            SumFluxByElement[rightIndex] += -1.0*fluxInt;
        }
    }
    for(int iel = 0; iel<nInterfacesByElement.size(); iel++){
//        std::cout<<"Iel: "<<iel<<" nInterfaces: "<<nInterfacesByElement[iel]<<std::endl;
//        std::cout<<"Integral Flux: "<<SumFluxByElement[iel]<<std::endl;
        if(abs(SumFluxByElement[iel])>1.0e-10){
            std::cout<<"The sum of the flows on each element must be zero. Element:  "<<iel<<" has a value of"<<SumFluxByElement[iel]<<std::endl;
            DebugStop();
        }
        
        if(nInterfacesByElement[iel]==numZeroFluxByElement[iel]){
            std::cout<<"Error: The sum of the flows over the elements is zero because all the flows on the interfaces are zero."<<std::endl;
            DebugStop();
        }
    }
    std::cout<<"The sum of the flows over the elements is zero. ¡This is correct!"<<std::endl;
}
