//
//  TMRSDataTransfer.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#include "TMRSDataTransfer.h"


TMRSDataTransfer::TMRSDataTransfer() : mTGeometry(), mTPetroPhysics(), mTFluidProperties(),mTReservoirProperties(), mTBoundaryConditions(), mTNumerics(), mTPostProcess(), mTFracProperties(){
    
}


TMRSDataTransfer::TMRSDataTransfer(const TMRSDataTransfer &other){
    mTGeometry = other.mTGeometry;
    mTPetroPhysics = other.mTPetroPhysics;
    mTFluidProperties = other.mTFluidProperties;
    mTReservoirProperties = other.mTReservoirProperties;
    mTBoundaryConditions = other.mTBoundaryConditions;
    mTNumerics = other.mTNumerics;
    mTPostProcess = other.mTPostProcess;
    mSimulationName = other.mSimulationName;
    mTFracProperties =other.mTFracProperties;
    mTFracIntersectProperties =other.mTFracIntersectProperties;
}


TMRSDataTransfer & TMRSDataTransfer::operator=(const TMRSDataTransfer &other){
    
    if (this != & other) // prevent self-assignment
    {
        mTGeometry = other.mTGeometry;
        mTPetroPhysics = other.mTPetroPhysics;
        mTFluidProperties = other.mTFluidProperties;
        mTReservoirProperties = other.mTReservoirProperties;
        mTBoundaryConditions = other.mTBoundaryConditions;
        mTNumerics = other.mTNumerics;
        mTPostProcess = other.mTPostProcess;
        mSimulationName = other.mSimulationName;
        mTFracProperties =other.mTFracProperties;
        mTFracIntersectProperties =other.mTFracIntersectProperties;
    }
    return *this;
}

TMRSDataTransfer::~TMRSDataTransfer(){
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void TMRSDataTransfer::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TMRSDataTransfer::Read(TPZStream &buf, void *context){
    DebugStop();
}

int TMRSDataTransfer::ClassId() const{
    DebugStop();
}
void TMRSDataTransfer::TFluidProperties::CreateLinearDensityFunction(){
    mWaterDensityF = [this](REAL &p){
        REAL dp = p - mReferencePressure;
        REAL rho = mWaterDensityRef*(1+ mWaterCompressibility*dp);
        REAL drho_dp = mWaterCompressibility*mWaterDensityRef;
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
    };
    mOilDensityF = [this](REAL &p){
        REAL dp = p - mReferencePressure;
        REAL rho = mOilDensityRef*(1+ mOilCompressibility*dp);
        REAL drho_dp = mOilCompressibility*mOilDensityRef;
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
    };
}
void TMRSDataTransfer::TFluidProperties::CreateExponentialDensityFunction(){
    mWaterDensityF = [this](REAL &p){
        REAL dp = p - mReferencePressure;
        REAL rho = mWaterDensityRef * std::exp(std::exp( (mWaterCompressibility*(dp))));
        
        REAL drho_dp = mWaterCompressibility*mWaterDensityRef * std::exp(std::exp( (mWaterCompressibility*(dp))));
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
        
    };
    
    mOilDensityF = [this](REAL &p){
        REAL dp = p - mReferencePressure;
        REAL rho = mOilDensityRef * std::exp(((mOilCompressibility*(dp))));
        
        REAL drho_dp = mOilCompressibility* mOilDensityRef * std::exp(((mOilCompressibility*(dp))));
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
    };
}

void TMRSDataTransfer::TPetroPhysics::CreateLinearKrModel(){
    mKrw = [](REAL &sw){
        REAL krw = sw;
        REAL dkrw = 1;
        std::tuple<REAL, REAL> valderiv(krw, dkrw);
        return valderiv;
    };
    mKro = [](REAL &sw){
        REAL krw = (1-sw);
        REAL dkrw = -1;
        std::tuple<REAL, REAL> valderiv(krw, dkrw);
        return valderiv;
    };
    UpdateLambdasAndFracFlows();
}
void TMRSDataTransfer::TPetroPhysics::CreateQuadraticKrModel(){
    
    mKrw = [](REAL &sw){
        REAL krw = sw*sw;
        REAL dkrw = 2*sw;
        std::tuple<REAL, REAL> valderiv(krw, dkrw);
        return valderiv;
    };
    mKro = [](REAL &sw){
        REAL krw = (1-sw)*(1-sw);
        REAL dkrw = -2*(1-sw);
        std::tuple<REAL, REAL> valderiv(krw, dkrw);
        return valderiv;
    };
    UpdateLambdasAndFracFlows();
}

void TMRSDataTransfer::TPetroPhysics::UpdateLambdasAndFracFlows(){
    mLambdaW = [this](REAL &sw){
        std::tuple<REAL, REAL> krwvalderiv = mKrw(sw);
        auto krw = std::get<0>(krwvalderiv);
        REAL dkrw = std::get<1>(krwvalderiv);
        REAL lambdaw = krw/mWaterViscosity;
        REAL dlambdadsw = dkrw/mWaterViscosity;
        std::tuple<REAL, REAL> valderiv(lambdaw, dlambdadsw);
        return valderiv;
    };
    mLambdaO = [this](REAL &sw){
        std::tuple<REAL, REAL> krwvalderiv = mKro(sw);
        auto kro = std::get<0>(krwvalderiv);
        REAL dkro = std::get<1>(krwvalderiv);
        REAL lambdao = kro/mOilViscosity;
        REAL dlambdadso = dkro/mOilViscosity;
        std::tuple<REAL, REAL> valderiv(lambdao, dlambdadso);
        return valderiv;
    };
    
    mLambdaTotal = [this](REAL &sw){
        std::tuple<REAL, REAL> lambdaWvalderiv = mLambdaW(sw);
        std::tuple<REAL, REAL> lambdaOvalderiv = mLambdaO(sw);
        REAL lw = std::get<0>(lambdaWvalderiv);
        REAL dlwdsw = std::get<1>(lambdaWvalderiv);
        REAL lo = std::get<0>(lambdaOvalderiv);
        REAL dlodsw = std::get<1>(lambdaOvalderiv);
        REAL lambdaTotal = lw+lo;
        REAL dlambdaTotaldsw = dlwdsw+dlodsw;
        std::tuple<REAL, REAL> valderiv(lambdaTotal, dlambdaTotaldsw);
        return valderiv;
    };
    mFw = [this](REAL &sw){
        std::tuple<REAL, REAL> lambdaWvalderiv = mLambdaW(sw);
        std::tuple<REAL, REAL> lambdaTotalvalderiv = mLambdaTotal(sw);
        REAL lw = std::get<0>(lambdaWvalderiv);
        REAL dlwdsw = std::get<1>(lambdaWvalderiv);
        REAL ltotal = std::get<0>(lambdaTotalvalderiv);
        REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
        REAL fracflow =lw/ltotal;
        REAL dfracflowdsw = (dlwdsw/ltotal) - ((lw*dltotaldsw)/(ltotal*ltotal));
        std::tuple<REAL, REAL> valderiv(fracflow, dfracflowdsw);
        return valderiv;
    };
    mFo = [this](REAL &sw){
        std::tuple<REAL, REAL> lambdaOvalderiv = mLambdaO(sw);
        std::tuple<REAL, REAL> lambdaTotalvalderiv = mLambdaTotal(sw);
        REAL lo = std::get<0>(lambdaOvalderiv);
        REAL dlodso = std::get<1>(lambdaOvalderiv);
        REAL ltotal = std::get<0>(lambdaTotalvalderiv);
        REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
        REAL fracflow =lo/ltotal;
        REAL dfracflowdsw = (dlodso/ltotal) - ((lo*dltotaldsw)/(ltotal*ltotal));
        std::tuple<REAL, REAL> valderiv(fracflow, dfracflowdsw);
        return valderiv;
    };
    
    
}
