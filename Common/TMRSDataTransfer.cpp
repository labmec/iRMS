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
       
        REAL dp =0.0;// p - mReferencePressure;
        REAL rho = mWaterDensityRef;//*(1+ mWaterCompressibility*dp);
        REAL drho_dp = 0.0;//mWaterCompressibility*mWaterDensityRef;
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
    };
    mOilDensityF = [this](REAL &p){
        REAL dp =0.0;// p - mReferencePressure;
        REAL rho = mOilDensityRef;//*(1+ mOilCompressibility*dp);
        REAL drho_dp = 0.0;//mOilCompressibility*mOilDensityRef;
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
    };
}
void TMRSDataTransfer::TFluidProperties::CreateExponentialDensityFunction(){
    mWaterDensityF = [this](REAL &p){
        REAL dp = 0.0;//p - mReferencePressure;
        REAL rho = mWaterDensityRef;// * std::exp(std::exp( (mWaterCompressibility*(dp))));
        
        REAL drho_dp = 0.0;// mWaterCompressibility*mWaterDensityRef * std::exp(std::exp( (mWaterCompressibility*(dp))));
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
        
    };
    
    mOilDensityF = [this](REAL &p){
        REAL dp = 0.0;// p - mReferencePressure;
        REAL rho = mOilDensityRef;// * std::exp(((mOilCompressibility*(dp))));
        
        REAL drho_dp = 0.0;//mOilCompressibility* mOilDensityRef * std::exp(((mOilCompressibility*(dp))));
        std::tuple<REAL, REAL> valderiv(rho, drho_dp);
        return valderiv;
    };
}

void TMRSDataTransfer::TPetroPhysics::CreateLinearKrModel(){
    mKrw = [](REAL &sw){
        REAL krw = sw;
        REAL dkrw = 1;
        REAL d2krw=0.0;
        std::tuple<REAL, REAL, REAL> valderiv(krw, dkrw,d2krw);
        return valderiv;
    };
    mKro = [](REAL &sw){
        REAL krw = (1-sw);
        REAL dkrw = -1;
        REAL d2krw=0.0;
        std::tuple<REAL, REAL, REAL> valderiv(krw, dkrw,d2krw);
        return valderiv;
    };
    UpdateLambdasAndFracFlows();
}
void TMRSDataTransfer::TPetroPhysics::CreateQuadraticKrModel(){
    
    mKrw = [](REAL &sw){
        REAL krw = sw*sw;
        REAL dkrw = 2*sw;
        REAL d2krw= 2;
        std::tuple<REAL, REAL,REAL> valderiv(krw, dkrw,d2krw);
        return valderiv;
    };
    mKro = [](REAL &sw){
        REAL krw = (1-sw)*(1-sw);
        REAL dkrw = -2*(1-sw);
        REAL d2krw = 2.0;
        std::tuple<REAL, REAL, REAL> valderiv(krw, dkrw,d2krw );
        return valderiv;
    };
    UpdateLambdasAndFracFlows();
}

void TMRSDataTransfer::TPetroPhysics::UpdateLambdasAndFracFlows(){
    mLambdaW = [this](REAL &sw){
        std::tuple<REAL, REAL,REAL> krwvalderiv = mKrw(sw);
        auto krw = std::get<0>(krwvalderiv);
        REAL dkrw = std::get<1>(krwvalderiv);
        REAL d2krw = std::get<2>(krwvalderiv);
        
        REAL lambdaw = krw/mWaterViscosity;
        REAL dlambdadsw = dkrw/mWaterViscosity;
        REAL d2lambdadsw = d2krw/mWaterViscosity;
        
        std::tuple<REAL, REAL, REAL> valderiv(lambdaw, dlambdadsw,d2lambdadsw);
        return valderiv;
    };
    mLambdaO = [this](REAL &sw){
        std::tuple<REAL, REAL, REAL> krwvalderiv = mKro(sw);
        auto kro = std::get<0>(krwvalderiv);
        REAL dkro = std::get<1>(krwvalderiv);
        REAL d2kro = std::get<2>(krwvalderiv);
        
        REAL lambdao = kro/mOilViscosity;
        REAL dlambdadso = dkro/mOilViscosity;
        REAL d2lambdadso = d2kro/mOilViscosity;
        
        std::tuple<REAL, REAL, REAL> valderiv(lambdao, dlambdadso, d2lambdadso);
        return valderiv;
    };
    
    mLambdaTotal = [this](REAL &sw){
        std::tuple<REAL, REAL, REAL> lambdaWvalderiv = mLambdaW(sw);
        std::tuple<REAL, REAL, REAL> lambdaOvalderiv = mLambdaO(sw);
        REAL lw = std::get<0>(lambdaWvalderiv);
        REAL dlwdsw = std::get<1>(lambdaWvalderiv);
        REAL d2lwdsw = std::get<2>(lambdaWvalderiv);
        
        REAL lo = std::get<0>(lambdaOvalderiv);
        REAL dlodsw = std::get<1>(lambdaOvalderiv);
        REAL d2lodsw = std::get<2>(lambdaOvalderiv);
        
        REAL lambdaTotal = lw+lo;
        REAL dlambdaTotaldsw = dlwdsw + dlodsw;
        REAL d2lambdaTotaldsw = d2lwdsw + d2lodsw;
        
        std::tuple<REAL, REAL, REAL> valderiv(lambdaTotal, dlambdaTotaldsw, d2lambdaTotaldsw);
        return valderiv;
    };
    mFw = [this](REAL &sw){
        std::tuple<REAL, REAL, REAL> lambdaWvalderiv = mLambdaW(sw);
        std::tuple<REAL, REAL, REAL> lambdaTotalvalderiv = mLambdaTotal(sw);
        REAL lw = std::get<0>(lambdaWvalderiv);
        REAL dlwdsw = std::get<1>(lambdaWvalderiv);
        REAL d2lwdsw = std::get<2>(lambdaWvalderiv);
        
        REAL ltotal = std::get<0>(lambdaTotalvalderiv);
        REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
        REAL d2ltotaldsw = std::get<2>(lambdaTotalvalderiv);
        
        REAL fracflow =lw/ltotal;
        REAL dfracflowdsw = (dlwdsw/ltotal) - ((lw*dltotaldsw)/(ltotal*ltotal));
    
        REAL termA=(2.0*(lw*dltotaldsw*dltotaldsw)/(ltotal*ltotal*ltotal));
        REAL termB= ((-2.0*dltotaldsw*dlwdsw)/(ltotal*ltotal));
        REAL termC= ((-1.0*lw*d2ltotaldsw)/(ltotal*ltotal));
        REAL termD= d2lwdsw/ltotal;
        
        REAL d2fracflowdsw = termA+termB+termC + termD;
        
        std::tuple<REAL, REAL, REAL> valderiv(fracflow, dfracflowdsw,d2fracflowdsw);
        return valderiv;
    };
    mFo = [this](REAL &sw){
        std::tuple<REAL, REAL, REAL> lambdaOvalderiv = mLambdaO(sw);
        std::tuple<REAL, REAL, REAL> lambdaTotalvalderiv = mLambdaTotal(sw);
        REAL lo = std::get<0>(lambdaOvalderiv);
        REAL dlodso = std::get<1>(lambdaOvalderiv);
        REAL d2lodso = std::get<2>(lambdaOvalderiv);
        
        REAL ltotal = std::get<0>(lambdaTotalvalderiv);
        REAL dltotaldsw = std::get<1>(lambdaTotalvalderiv);
        REAL d2ltotaldsw = std::get<2>(lambdaTotalvalderiv);
        
        REAL fracflow =lo/ltotal;
        REAL dfracflowdsw = (dlodso/ltotal) - ((lo*dltotaldsw)/(ltotal*ltotal));
        
        REAL termA=(2.0*(lo*dltotaldsw*dltotaldsw)/(ltotal*ltotal*ltotal));
        REAL termB= ((-2.0*dltotaldsw*dlodso)/(ltotal*ltotal));
        REAL termC= ((-1.0*lo*d2ltotaldsw)/(ltotal*ltotal));
        REAL termD= d2lodso/ltotal;
        
        REAL d2fracflowdsw = termA+termB+termC + termD;
        
        std::tuple<REAL, REAL, REAL> valderiv(fracflow, dfracflowdsw, d2fracflowdsw);
        return valderiv;
    };
    
    
}
