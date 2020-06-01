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
}

/// Assignement constructor
const AlgebraicTransport & AlgebraicTransport::operator=(const AlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fCoefficientsFlux = other.fCoefficientsFlux;
    fIntegralFluxFunctions = other.fIntegralFluxFunctions;
    fIntegralFlux = other.fIntegralFlux;
    fFluxSing = other.fFluxSing;
    fNormalFaceDirection = other.fNormalFaceDirection;
    return *this;
}


    void AlgebraicTransport::SetCoefficientsFlux(std::vector< std::vector<int> > fluxCoefficients){
        fCoefficientsFlux =fluxCoefficients;
    }
    void AlgebraicTransport::SetIntegralFluxFunctions(std::vector< std::vector<int> > fluxIntegralFunctions){
        fIntegralFluxFunctions = fluxIntegralFunctions;
    }
    void AlgebraicTransport::SetIntegralFlux(std::vector<int > fluxIntegral){
        fIntegralFlux = fluxIntegral;
    }
    void AlgebraicTransport::SetSingFlux(std::vector<int> fluxSing){
        fFluxSing = fluxSing;
    }
    void AlgebraicTransport::SetNormalFaceDirection(std::vector<int> NormalFaceDirection){
        fNormalFaceDirection = NormalFaceDirection;
    
    }
   
    std::vector< std::vector<int> >  AlgebraicTransport::GetCoefficientsFlux(){
        return fCoefficientsFlux;
    }
    std::vector< std::vector<int> > AlgebraicTransport::GetIntegralFluxFunctions(){
        return fIntegralFluxFunctions;
    }
    std::vector< int> AlgebraicTransport::GetIntegralFlux(){
        return fIntegralFlux;
    }
    std::vector< int > AlgebraicTransport::GetSingFlux(){
        return fFluxSing;
    }
    std::vector<int> AlgebraicTransport::GetNormalFaceDirection(){
        return fNormalFaceDirection;
    }

