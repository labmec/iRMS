//
//  AlgebraicDataTransfer.hpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#ifndef AlgebraicDataTransfer_h
#define AlgebraicDataTransfer_h

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
class AlgebraicDataTransfer {
    private:
        std::vector<double> fSaturation;
        std::vector<double> fPressure;
        std::vector<double> fDensityOil;
        std::vector<double> fDensityWater;
        std::vector<double> fCompressibility;
        std::vector<double> flambdas;
    public:
    
    /// Default constructor
    AlgebraicDataTransfer();
    
    /// Copy constructor
    AlgebraicDataTransfer(const AlgebraicDataTransfer & other);
    
    /// Assignement constructor
    const AlgebraicDataTransfer & operator=(const AlgebraicDataTransfer & other);
    
    /// Default desconstructor
    ~AlgebraicDataTransfer();
    
    void CalcLambdas();
    void CalcDensities();
    double CalcLambda(double sw, double paverage, double densityo, double densityw);
    double CalcDensity(double paverage,double compress, double reff,  double pref);
    
    void SetSaturation(std::vector<double> saturation);
    void SetPressure(std::vector<double> pressure);
    void SetDensities(std::vector<double> densityo, std::vector<double> densityw);
    void SetCompressibility(std::vector<double> compressibility);
    
    std::vector<double> GetSaturation();
    std::vector<double> GetPressure();
    std::vector<double> GetWaterDensity();
    std::vector<double> GetOilDensity();
    std::vector<double> GetCompressibility();
    
};
#endif /* AlgebraicDataTransfer_h */
