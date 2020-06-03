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

#include "TPZMultiphysicsCompMesh.h"

class AlgebraicTransport {
    
public:
    int fNFluxCoefficients;
    //vector for each multiplying coefficient of fluxes
    std::vector< std::vector<REAL> > fCoefficientsFlux;
    //Integral of the flux functions associated with faces
    std::vector< std::vector<REAL> > fIntegralFluxFunctions;
    //  Integral of the flux
    std::vector<REAL>  fIntegralFlux;
    //Sign with respect to the natural orientation of the flux (+1 or -1)
    std::vector<REAL> fFluxSing;
    //Direction of the normal to the face
    std::vector<REAL> fNormalFaceDirection;
    
    // CELL DATA
    // Volume of the cells
    std::vector<REAL> fVolumes;
    
    std::vector<double> fSaturation;
    std::vector<double> fPressure;
    std::vector<double> fDensityOil;
    std::vector<double> fDensityWater;
    std::vector<double> fCompressibility;
    std::vector<double> flambdas;

public:
     /// Default constructor
    AlgebraicTransport();
    
    /// Copy constructor
    AlgebraicTransport(const AlgebraicTransport & other);
    
    /// Assignement constructor
    const AlgebraicTransport & operator=(const AlgebraicTransport & other);
    
    /// Default desconstructor
    ~AlgebraicTransport();
    
    void BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh);
    
    void SetCoefficientsFlux(std::vector< std::vector<REAL> > &fluxCoefficients);
    void SetIntegralFluxFunctions(std::vector< std::vector<REAL> > &fluxIntegralFunctions);
    void SetIntegralFlux(std::vector<REAL> &fluxIntegral);
    void SetSingFlux(std::vector<REAL> &fluxIntegral);
    void SetNormalFaceDirection(std::vector<REAL> &NormalFaceDirection);
    
    std::vector< std::vector<REAL> > &GetCoefficientsFlux();
    std::vector< std::vector<REAL> > &GetIntegralFluxFunctions();
    std::vector<REAL> &GetIntegralFlux();
    std::vector<REAL> &GetSingFlux();
    std::vector<REAL> &GetNormalFaceDirection();
    
    void CalcLambdas();
    void CalcDensities();
    double CalcLambda(double sw, double paverage, double densityo, double densityw);
    double CalcDensity(double paverage,double compress, double reff,  double pref);
    
    void SetSaturation(std::vector<double> &saturation);
    void SetPressure(std::vector<double> &pressure);
    void SetDensities(std::vector<double> &densityoil, std::vector<double> &densitywater);
    void SetCompressibility(std::vector<double> &compressibility);
    
    std::vector<double> &GetSaturation();
    std::vector<double> &GetPressure();
    std::vector<double> &GetWaterDensity();
    std::vector<double> &GetOilDensity();
    std::vector<double> &GetCompressibility();
    

};

#endif /* AlgebraicTransport_h */
