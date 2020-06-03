//
//  AlgebraicTransport.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "AlgebraicTransport.h"

/// Default constructor
AlgebraicTransport::AlgebraicTransport(){
    
}

/// Copy constructor
AlgebraicTransport::AlgebraicTransport(const AlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fCoefficientsFlux = other.fCoefficientsFlux;
    fIntegralFluxFunctions = other.fIntegralFluxFunctions;
    fIntegralFlux = other.fIntegralFlux;
    fFluxSing = other.fFluxSing;
    fNormalFaceDirection = other.fNormalFaceDirection;
    fSaturation= other.fSaturation;
    fPressure = other.fPressure;
    fDensityOil = other.fDensityOil;
    fDensityWater = other.fDensityWater;
    fCompressibility = other.fCompressibility;
    flambdas = other.flambdas;
}

/// Assignement constructor
const AlgebraicTransport & AlgebraicTransport::operator=(const AlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fCoefficientsFlux = other.fCoefficientsFlux;
    fIntegralFluxFunctions = other.fIntegralFluxFunctions;
    fIntegralFlux = other.fIntegralFlux;
    fFluxSing = other.fFluxSing;
    fNormalFaceDirection = other.fNormalFaceDirection;
    fSaturation= other.fSaturation;
    fPressure = other.fPressure;
    fDensityOil = other.fDensityOil;
    fDensityWater = other.fDensityWater;
    fCompressibility = other.fCompressibility;
    flambdas = other.flambdas;
    return *this;
}


void AlgebraicTransport::SetCoefficientsFlux(std::vector< std::vector<REAL> > &fluxCoefficients){
    fCoefficientsFlux =fluxCoefficients;
}
void AlgebraicTransport::SetIntegralFluxFunctions(std::vector< std::vector<REAL> > & fluxIntegralFunctions){
    fIntegralFluxFunctions = fluxIntegralFunctions;
}
void AlgebraicTransport::SetIntegralFlux(std::vector<REAL> &fluxIntegral){
    fIntegralFlux = fluxIntegral;
}
void AlgebraicTransport::SetSingFlux(std::vector<REAL> &fluxSing){
    fFluxSing = fluxSing;
}
void AlgebraicTransport::SetNormalFaceDirection(std::vector<REAL> &NormalFaceDirection){
    fNormalFaceDirection = NormalFaceDirection;
    
}

std::vector< std::vector<REAL> > &AlgebraicTransport::GetCoefficientsFlux(){
    return fCoefficientsFlux;
}
std::vector< std::vector<REAL> > &AlgebraicTransport::GetIntegralFluxFunctions(){
    return fIntegralFluxFunctions;
}
std::vector<REAL> &AlgebraicTransport::GetIntegralFlux(){
    return fIntegralFlux;
}

void AlgebraicTransport::CalcLambdas(){
    int nlamdas = fPressure.size();
    flambdas.resize(nlamdas);
    for (int ilambda =0; ilambda<nlamdas; ilambda++) {
        flambdas[ilambda] = CalcLambda(fSaturation[ilambda], fPressure[ilambda], fDensityOil[ilambda], fDensityWater[ilambda] );
    }
}
void AlgebraicTransport::CalcDensities(){
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

double AlgebraicTransport::CalcLambda(double sw, double paverage, double densityo, double densityw){
    double kro = (1 - sw)*(1 - sw);
    double krw = sw*sw;
    double muo = 0.05;
    double muw = 0.01;
    double lambda;
    lambda = ((densityo*kro)/muo) + ((densityw*krw)/muw);
    
    return lambda;
}

double AlgebraicTransport::CalcDensity(double paverage,double compress, double reff, double pref){
    double density = reff*(1+ compress*(paverage-pref));
    return density;
}

void AlgebraicTransport::SetSaturation(std::vector<double> &saturation){
    fSaturation = saturation;
}
void AlgebraicTransport::SetPressure(std::vector<double> &pressure){
    fPressure =pressure;
}
void AlgebraicTransport::SetDensities(std::vector<double> &densityoil, std::vector<double> &densitywater){
    fDensityOil=densityoil;
    fDensityWater=densitywater;
}
void AlgebraicTransport::SetCompressibility(std::vector<double> &compressibility){
    fCompressibility = compressibility;
}

std::vector<double> &AlgebraicTransport::GetSaturation(){
    return fSaturation;
}

std::vector<double> &AlgebraicTransport::GetPressure(){
    return fPressure;
}

std::vector<double> &AlgebraicTransport::GetWaterDensity(){
    return fDensityWater;
}

std::vector<double> &AlgebraicTransport::GetOilDensity(){
    return fDensityOil;
}

std::vector<double> &AlgebraicTransport::GetCompressibility(){
    return fCompressibility;
}

std::vector<REAL> &AlgebraicTransport::GetSingFlux(){
    return fFluxSing;
}
std::vector<REAL> &AlgebraicTransport::GetNormalFaceDirection(){
    return fNormalFaceDirection;
}

void AlgebraicTransport::BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh)
{
    DebugStop();
}

