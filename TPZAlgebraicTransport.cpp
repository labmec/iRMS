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
    std::vector<REAL> leftinf= fCellsData.fWaterfractionalflow[lr_index.first];
    std::vector<REAL> rightinf= fCellsData.fWaterfractionalflow[lr_index.second];
    REAL fw_L = leftinf[0];
    REAL fw_R = rightinf[0];
    REAL dfwSw_L = leftinf[1];
    REAL dfwSw_R = rightinf[1];
    REAL beta =0.0;
    //upwind
    if (fluxint>0.0) {
        beta = 1.0;
    }
    
    ef(0) = +1.0*fdt*(beta*fw_L + (1-beta)*fw_R)*fluxint;
    ef(1) = -1.0*fdt*(beta*fw_L + (1-beta)*fw_R)*fluxint;
    
    ek(0,0) = +1.0*dfwSw_L * fdt * beta * fluxint;
    ek(0,1) = +1.0*dfwSw_R * fdt *(1-beta) * fluxint;
    ek(1,0) = -1.0*dfwSw_L * fdt* beta * fluxint;
    ek(1,1) = -1.0*dfwSw_R*fdt*(1-beta)*fluxint;


}
void TPZAlgebraicTransport::ContributeBCInletInterface(int index, TPZFMatrix<double> &ef){
   
    int s_inlet =1.0;
    REAL fluxint  = fInterfaceData[inletmatid].fIntegralFlux[index];
    ef(0,0) = 1.0*fdt*s_inlet*fluxint;
}
void TPZAlgebraicTransport::ContributeBCOutletInterface(int index,TPZFMatrix<double> &ek, TPZFMatrix<double> &ef){
  
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[outletmatid].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[outletmatid].fIntegralFlux[index];
    std::vector<REAL> leftinf= fCellsData.fWaterfractionalflow[lr_index.first];
    std::vector<REAL> rightinf= fCellsData.fWaterfractionalflow[lr_index.second];
    REAL fw_L = leftinf[0];
    REAL dfwSw_L = leftinf[1];
    ef(0,0) = 1.0*fdt*fw_L*fluxint;
    ek(0,0) = dfwSw_L*fdt*fluxint;
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
        out<<"WaterFracFlow: "<<this->fWaterfractionalflow[iel][0]<<std::endl;
        out<<"dWaterFracFlowdSw: "<<this->fWaterfractionalflow[iel][1]<<std::endl;
        out<<"OilFracFlow: "<<this->fOilfractionalflow[iel][0]<<std::endl;
        out<<"dOilFracFlowdSw: "<<this->fOilfractionalflow[iel][1]<<std::endl;
        out<<"Lambda: "<<this->flambda[iel]<<std::endl;
        out<<std::endl;
    }

}

std::pair<std::vector<REAL>,std::vector<REAL>> TPZAlgebraicTransport::fwAndfoVal(REAL sw, REAL muw,REAL muo){
    std::vector<REAL> fwData(2), foData(2);
//    REAL Krw = sw*sw ;
//    REAL Kro = (1-sw)*(1-sw) ;
    
    REAL fw = (muo*sw*sw)/(muw*(sw-1.0)*(sw-1.0) + (muo*sw*sw));
    REAL num = -2.0*(muo*muw*(sw-1.0)*sw);
    REAL dem = ((muw*(sw-1.0)*(sw-1.0))+(muo*sw*sw))*((muw*(sw-1.0)*(sw-1.0))+(muo*sw*sw));
    fwData[0] = fw;
//    fwData[1] = num/dem;
    fwData[1]=(muo*muw)/ ((muw + muo*sw - muw*sw)*(muw + muo*sw - muw*sw));
    foData[0] = 1.0;
    foData[1]=1.0;
  
//    REAL Krw = sw;
//    REAL Kro = (1-sw);
//    REAL fw = (Krw)/(Krw+Kro);
//    REAL fo = (Kro)/(Krw+Kro);
//    fwData[0] = fw;
//    fwData[1] = 1.0;
//    foData[0] = fo;
//    foData[1]=-1.0;

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
        REAL sw =this->fSaturation[ivol];
        REAL krw, kro;
        if (isLinearQ) {
            fWandFo=
            LinearfwAndfoVal(fSaturation[ivol], fViscosity[0], fViscosity[1]);
            krw = sw;
            kro = (1-sw);
        }
        else{
        fWandFo=
        fwAndfoVal(fSaturation[ivol], fViscosity[0], fViscosity[1]);
        krw = sw*sw;
        kro = (1-sw)*(1-sw);
            
        }
        this->fWaterfractionalflow[ivol] = fWandFo.first;
        this->fOilfractionalflow[ivol] = fWandFo.second;
       

 

        this->flambda[ivol] = (krw/(fViscosity[0]))+(kro/(fViscosity[1]));
//        this->flambda[ivol] = (krw)+(kro);

    }
}

void TPZAlgebraicTransport::TCellData::UpdateSaturations(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        fSaturation[icell] = sw(eq_number);

    }
}
void TPZAlgebraicTransport::TCellData::UpdateSaturationsLastState(TPZFMatrix<STATE> &sw){
    int ncells = fVolume.size();
    for (int icell = 0; icell<ncells; icell++) {
        int eq_number = fEqNumber[icell];
        fSaturationLastState[icell] = sw(eq_number);
        
    }
}
