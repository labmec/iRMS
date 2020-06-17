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
    fgravity.resize(3,0.0);
    fgravity[1] = -10.0;//-9.81;
    
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
    REAL ndotg = (std::get<0>(normal))*fgravity[0]+(std::get<1>(normal))*fgravity[1]+(std::get<2>(normal))*fgravity[2];
    REAL rhowL = fCellsData.fDensityWater[lr_index.first];
    REAL rhowR = fCellsData.fDensityWater[lr_index.second];
    REAL rhooL = fCellsData.fDensityOil[lr_index.first];
    REAL rhooR = fCellsData.fDensityOil[lr_index.second];
    REAL lambdaL = fCellsData.flambda[lr_index.first];
    REAL lambdaR = fCellsData.flambda[lr_index.second];
        
    //    beta = 0.0;
    //    REAL gravfluxL = fCellsData.fKz[lr_index.first]*lambdaL*(rhowL-rhooL);
    //    REAL gravfluxR = fCellsData.fKz[lr_index.second]*lambdaR*(rhowR-rhooR);
    //    if (ndotg>0.0) {
    //        beta = 1.0;
    //    }
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
        REAL Krw, Kro, fw, dfwdsw, fo, dfodsw;
        REAL muw = fViscosity[0];
        REAL muo = fViscosity[1];
        if (isLinearQ) {
            Krw = sw;
            Kro = (1-sw);
            fw = (Krw/muw)/((Krw/muw)+(Kro/muo));
            fo = (Kro/muo)/((Krw/muw)+(Kro/muo));
            dfwdsw = (muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
            dfodsw = -1.0*dfwdsw;
           
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
        }
        this->fWaterfractionalflow[ivol] = fw;
        this->fDerivativeWfractionalflow[ivol] =dfwdsw;
        this->fOilfractionalflow[ivol] = fo;
        this->fDerivativeOfractionalflow[ivol] = dfodsw;
        this->flambda[ivol] = (Krw/muw)+(Kro/muo);

    }
}

void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambdaQuasiNewton(){
    int nvols = this->fVolume.size();
    REAL muw= fViscosity[0];
    REAL muo= fViscosity[1];
    for (int ivol =0 ; ivol< nvols; ivol++) {
        std::pair<std::vector<REAL>,std::vector<REAL>> fWandFo;
        REAL sw =this->fSaturation[ivol];
        REAL krw, kro, fw, fo,dfwdsw;
        krw = sw*sw;
        kro = (1-sw)*(1-sw);
        fw = (muw*krw)/((muw*kro) + (muo*krw));
        fo = (muo*kro)/((muw*kro) + (muo*krw));
        dfwdsw =(muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
       
        this->fWaterfractionalflow[ivol] = fw;
        this->fDerivativeWfractionalflow[ivol] =dfwdsw;
        this->fOilfractionalflow[ivol] =fo;
        this->fDerivativeOfractionalflow[ivol] = -1.0*dfwdsw;
    
        this->flambda[ivol] = (krw/(fViscosity[0]))+(kro/(fViscosity[1]));
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
