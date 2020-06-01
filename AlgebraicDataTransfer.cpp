//
//  AlgebraicDataTransfer.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "AlgebraicDataTransfer.h"


/// Default constructor
AlgebraicDataTransfer::AlgebraicDataTransfer(){
    
}

/// Copy constructor
AlgebraicDataTransfer::AlgebraicDataTransfer(const AlgebraicDataTransfer & other){
    fSaturation= other.fSaturation;
    fPressure = other.fPressure;
    fDensityOil = other.fDensityOil;
    fDensityWater = other.fDensityWater;
    fCompressibility = other.fCompressibility;
    flambdas = other.flambdas;
}

/// Assignement constructor
const AlgebraicDataTransfer & AlgebraicDataTransfer::operator=(const AlgebraicDataTransfer & other){
    fSaturation= other.fSaturation;
    fPressure = other.fPressure;
    fDensityOil = other.fDensityOil;
    fDensityWater = other.fDensityWater;
    fCompressibility = other.fCompressibility;
    flambdas = other.flambdas;
    return *this;
    
}

/// Default desconstructor
 AlgebraicDataTransfer::~AlgebraicDataTransfer(){
    
}

void AlgebraicDataTransfer::CalcLambdas(){
     int nlamdas = fPressure.size();
    flambdas.resize(nlamdas);
    for (int ilambda =0; ilambda<nlamdas; ilambda++) {
       flambdas[ilambda] = CalcLambda(fSaturation[ilambda], fPressure[ilambda], fDensityOil[ilambda], fDensityWater[ilambda] );
    }
}
void AlgebraicDataTransfer::CalcDensities(){
    int ndensities = fPressure.size();
    fDensityOil.resize(ndensities);
    fDensityWater.resize(ndensities);
    double rhoo_ref = 760.00;
    double rhow_ref = 1000.00;
    double p_ref = 1.013E5;
    for (int idensity = 0; idensity< ndensities; idensity++) {
        fDensityWater[idensity] = CalcDensity(fPressure[idensity], fCompressibility[0],rhow_ref, p_ref);
          fDensityOil[idensity] = CalcDensity(fPressure[idensity],fCompressibility[1], rhoo_ref, p_ref);
    }
}
double AlgebraicDataTransfer::CalcLambda(double sw, double paverage, double densityo, double densityw){
    double kro = (1 - sw)*(1 - sw);
    double krw = sw*sw;
    double muo = 0.05;
    double muw = 0.01;
    double lambda;
    lambda = ((densityo*kro)/muo) + ((densityw*krw)/muw);
    
    return lambda;
}
double AlgebraicDataTransfer::CalcDensity(double paverage,double compress, double reff, double pref){
    double density = reff*(1+ compress*(paverage-pref));
    return density;
}

void AlgebraicDataTransfer::SetSaturation(std::vector<double> saturation){
    fSaturation = saturation;
}
void AlgebraicDataTransfer::SetPressure(std::vector<double> pressure){
    fPressure =pressure;
}
void AlgebraicDataTransfer::SetDensities(std::vector<double> densityo, std::vector<double> densityw){
    fDensityOil=densityo;
    fDensityWater=densityw;
}
void AlgebraicDataTransfer::SetCompressibility(std::vector<double> compressibility){
    fCompressibility = compressibility;
}

std::vector<double> AlgebraicDataTransfer::GetSaturation(){
    return fSaturation;
}
std::vector<double> AlgebraicDataTransfer::GetPressure(){
    return fPressure;
}
std::vector<double> AlgebraicDataTransfer::GetWaterDensity(){
    return fDensityWater;
}
std::vector<double> AlgebraicDataTransfer::GetOilDensity(){
    return fDensityOil;
}
std::vector<double> AlgebraicDataTransfer::GetCompressibility(){
    return fCompressibility;
}
