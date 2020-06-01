//
//  AlgebraicTransport.hpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#ifndef AlgebraicTransport_h
#define AlgebraicTransport_h
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>
class AlgebraicTransport {
    
private:
    int fNFluxCoefficients;
    //vector for each multiplying coefficient of fluxes
    std::vector< std::vector<int> > fCoefficientsFlux;
    //Integral of the flux functions associated with faces
    std::vector< std::vector<int> > fIntegralFluxFunctions;
    //  Integral of the flux
    std::vector<int>  fIntegralFlux;
    //Sign with respect to the natural orientation of the flux (+1 or -1)
    std::vector<int> fFluxSing;
    //Direction of the normal to the face
    std::vector<int> fNormalFaceDirection;
    
public:
     /// Default constructor
    AlgebraicTransport();
    
    /// Copy constructor
    AlgebraicTransport(const AlgebraicTransport & other);
    
    /// Assignement constructor
    const AlgebraicTransport & operator=(const AlgebraicTransport & other);
    
    /// Default desconstructor
    ~AlgebraicTransport();
    
    
    void SetCoefficientsFlux(std::vector< std::vector<int> > fluxCoefficients);
    void SetIntegralFluxFunctions(std::vector< std::vector<int> > fluxIntegralFunctions);
    void SetIntegralFlux(std::vector< int> fluxIntegral);
    void SetSingFlux(std::vector< int > fluxIntegral);
    void SetNormalFaceDirection(std::vector<int> fNormalFaceDirection);
    
    std::vector< std::vector<int> >  GetCoefficientsFlux();
    std::vector< std::vector<int> > GetIntegralFluxFunctions();
    std::vector< int> GetIntegralFlux();
    std::vector< int > GetSingFlux();
    std::vector<int> GetNormalFaceDirection();
};

#endif /* AlgebraicTransport_h */
