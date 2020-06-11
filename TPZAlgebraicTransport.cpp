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
}

/// Assignement constructor
const TPZAlgebraicTransport & TPZAlgebraicTransport::operator=(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
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
void TPZAlgebraicTransport::Assamble(TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef ){

}
void TPZAlgebraicTransport::AssambleResidual(TPZFNMatrix<100, REAL> &ef){
    
}
void TPZAlgebraicTransport::Contribute(int index, TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef){
    ek(0,0) = fCellsData.fVolume[index];
}
void TPZAlgebraicTransport::ContributeInterface(int index, TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef){
    
    int ninternInterrf = fInterfaceData[100].fFluxSign.size();
    if (ninternInterrf<=0) {
        DebugStop();
    }
    std::pair<int64_t, int64_t> lr_index = fInterfaceData[100].fLeftRightVolIndex[index];
    REAL fluxint  = fInterfaceData[100].fIntegralFlux[index];
    REAL dfwSwL = fCellsData.fOilfractionalflow[lr_index.first][1];
    REAL dfwSwR = fCellsData.fOilfractionalflow[lr_index.second][1];
    REAL beta =0;
    //upwind
    if (fluxint>0) {
        beta = 1.0;
    }
    
    ek(0,0) = dfwSwL * beta * fluxint;
    ek(0,1) = dfwSwR * (1-beta) * fluxint;
    ek(1,0) = -1.0 * dfwSwL * beta * fluxint;
    ek(1,1) = -1.0*dfwSwR*(1-beta)*fluxint;
    for (auto inter_bymat : fInterfaceData) {
        int ninter = inter_bymat.second.fIntegralFlux.size();
        for (int inter = 0; inter<ninter; inter++) {
            REAL fluxint = inter_bymat.second.fIntegralFlux[inter];
            
            REAL dfwSwL = fCellsData.fOilfractionalflow[lr_index.first][1];
            REAL dfwSwR = fCellsData.fOilfractionalflow[lr_index.second][1];
           
            REAL beta =0;
            //upwind
            if (fluxint>0) {
                beta = 1.0;
            }
            
            ek(0,0) = dfwSwL * beta * fluxint;
            ek(0,1) = dfwSwR * (1-beta) * fluxint;
            ek(1,0) = -1.0 * dfwSwL * beta * fluxint;
            ek(1,1) = -1.0*dfwSwR*(1-beta)*fluxint;
        }
    }
    ek.Print(std::cout);
}
void TPZAlgebraicTransport::ContributeBCInterface(int index, TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef){
    
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

std::pair<std::vector<REAL>,std::vector<REAL>> TPZAlgebraicTransport::fwAndfoVal(REAL Sw, REAL rhoW,REAL rhoO){
    std::vector<REAL> fwData(2), foData(2);
    REAL Krw = Sw*Sw;
    REAL Kro = (1-Sw)*(1-Sw);
    REAL fw = (Krw)/(Krw+Kro);
    REAL fo = (Kro)/(Krw+Kro);
    fwData[0] = fw;
    fwData[1] = (2.0*(Sw)/(Kro+Krw)) - (Krw*(4*Sw - 2))/((Krw+Kro)*(Krw+Kro));
    foData[0] = fo;
    foData[1]=(2.0*(Sw-1)/(Kro+Krw)) - (Kro*(4*Sw - 2))/((Krw+Kro)*(Krw+Kro));
    std::pair<std::vector<REAL>,std::vector<REAL>> fracflows = std::make_pair(fwData, foData);
    return fracflows;
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
void TPZAlgebraicTransport::TCellData::UpdateFractionalFlowsAndLambda(){
    int nvols = this->fVolume.size();
    for (int ivol =0 ; ivol< nvols; ivol++) {
        std::pair<std::vector<REAL>,std::vector<REAL>> fWandFo=
        fwAndfoVal(fSaturation[ivol], fDensityWater[ivol], fDensityOil[ivol]);
        this->fWaterfractionalflow[ivol] = fWandFo.first;
        this->fOilfractionalflow[ivol] = fWandFo.second;
        REAL sw =this->fSaturation[ivol];
        REAL krw = sw*sw;
        REAL kro = (1-sw)*(1-sw);
        this->flambda[ivol] = (krw/(fViscosity[0]))+(kro/(fViscosity[1]));
    }
}
