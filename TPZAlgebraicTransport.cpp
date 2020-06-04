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
const TPZAlgebraicTransport & TPZAlgebraicTransport::operator=(const TPZAlgebraicTransport & other){
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
TPZAlgebraicTransport::~TPZAlgebraicTransport(){
    
}

void TPZAlgebraicTransport::SetCoefficientsFlux(std::vector< std::vector<REAL> > &fluxCoefficients){
    fCoefficientsFlux =fluxCoefficients;
}
void TPZAlgebraicTransport::SetIntegralFluxFunctions(std::vector< std::vector<REAL> > & fluxIntegralFunctions){
    fIntegralFluxFunctions = fluxIntegralFunctions;
}
void TPZAlgebraicTransport::SetIntegralFlux(std::vector<REAL> &fluxIntegral){
    fIntegralFlux = fluxIntegral;
}
void TPZAlgebraicTransport::SetSingFlux(std::vector<REAL> &fluxSing){
    fFluxSing = fluxSing;
}
void TPZAlgebraicTransport::SetNormalFaceDirection(std::vector<REAL> &NormalFaceDirection){
    fNormalFaceDirection = NormalFaceDirection;
    
}

std::vector< std::vector<REAL> > &TPZAlgebraicTransport::GetCoefficientsFlux(){
    return fCoefficientsFlux;
}
std::vector< std::vector<REAL> > &TPZAlgebraicTransport::GetIntegralFluxFunctions(){
    return fIntegralFluxFunctions;
}
std::vector<REAL> &TPZAlgebraicTransport::GetIntegralFlux(){
    return fIntegralFlux;
}

void TPZAlgebraicTransport::CalcLambdas(){
    int nlamdas = fPressure.size();
    flambdas.resize(nlamdas);
    for (int ilambda =0; ilambda<nlamdas; ilambda++) {
        flambdas[ilambda] = CalcLambda(fSaturation[ilambda], fPressure[ilambda], fDensityOil[ilambda], fDensityWater[ilambda] );
    }
}
void TPZAlgebraicTransport::CalcDensities(){
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

void TPZAlgebraicTransport::SetSaturation(std::vector<double> &saturation){
    fSaturation = saturation;
}
void TPZAlgebraicTransport::SetPressure(std::vector<double> &pressure){
    fPressure =pressure;
}
void TPZAlgebraicTransport::SetDensities(std::vector<double> &densityoil, std::vector<double> &densitywater){
    fDensityOil=densityoil;
    fDensityWater=densitywater;
}
void TPZAlgebraicTransport::SetCompressibility(std::vector<double> &compressibility){
    fCompressibility = compressibility;
}

std::vector<double> &TPZAlgebraicTransport::GetSaturation(){
    return fSaturation;
}

std::vector<double> &TPZAlgebraicTransport::GetPressure(){
    return fPressure;
}

std::vector<double> &TPZAlgebraicTransport::GetWaterDensity(){
    return fDensityWater;
}

std::vector<double> &TPZAlgebraicTransport::GetOilDensity(){
    return fDensityOil;
}

std::vector<double> &TPZAlgebraicTransport::GetCompressibility(){
    return fCompressibility;
}

std::vector<REAL> &TPZAlgebraicTransport::GetSingFlux(){
    return fFluxSing;
}
std::vector<REAL> &TPZAlgebraicTransport::GetNormalFaceDirection(){
    return fNormalFaceDirection;
}

void TPZAlgebraicTransport::BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh)
{
    TPZMultiphysicsCompMesh *transp = &transportmesh;
    int neltransport = transp->NElements();
    for (int iel = 0; iel< neltransport; iel++) {
        TPZCompEl * cel = transp->Element(iel);
        if (!cel) {
            continue;
        }
        if (cel->Dimension() != transp->Dimension()) {
            continue;
        }
        TPZElementMatrix ek;
        TPZElementMatrix ef;
        cel->CalcStiff(ek, ef);
        double vol_val = ek.fMat.Get(0, 0);
        fVolumes.push_back(vol_val);
    }
    int nvols = fVolumes.size();
    
    //these could be default values
    fSaturation.resize(nvols,0.0);
    fDensityOil.resize(nvols,800.00);
    fDensityWater.resize(nvols,1000.00);
    fPressure.resize(nvols, 1.03E5);
    flambdas.resize(nvols,1.0);
    fCompressibility.resize(2);
    fCompressibility[0] = 0.0; //Water compressibility
    fCompressibility[1] = 0.0; //Oil compressibility
}

