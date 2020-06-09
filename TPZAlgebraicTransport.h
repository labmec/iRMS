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

class TPZAlgebraicTransport {
    
public:
    
    struct TInterfaceDataTransport
    {
        int64_t gelIndex;
        //vector for each multiplying coefficient of fluxes
        std::vector<REAL> fCoefficientsFlux;
        //Integral of the flux functions associated with faces
        std::vector<REAL> fIntegralFluxFunctions;
        //  Integral of the flux
        std::vector<REAL>  fIntegralFlux;
        //Sign with respect to the natural orientation of the flux (+1 or -1)
        REAL fFluxSing;
        //Direction of the normal to the face
        std::vector<REAL> fNormalFaceDirection;
        
        TInterfaceDataTransport(){
            
        }
        void Print(std::ostream &out);
    };
    
    // CELL DATA
    struct TCellData{
        int  index;
        REAL fVolume;
        REAL fSaturation;
        REAL fPressure;
        REAL fDensityOil;
        REAL fDensityWater;
        REAL flambda;
        
        TCellData() : index(-1), fVolume(-1), fSaturation(0.0), fPressure(-1), fDensityOil(800.00),fDensityWater(1000.00), flambda(1)
        {
            
        }
        
        void Print(std::ostream &out);
    };
    // fCompressibility[0] = water, fCompressibility[1]=oil, fCompressibility[2]=gas etc.
    std::vector<REAL> fCompressibility;
    std::vector<REAL> fReferencePressures;
    std::vector<REAL> fReferenceDensity;
    
    int fNFluxCoefficients;
    
    // Cells data structure, one material at a time
    std::map<int, std::vector<TCellData>> fCellsData;
    
    // Interface data structure, by material, element and side
    std::map<int, std::map<int, std::vector<TInterfaceDataTransport>>> fInterfaceData;

    
public:
    
     /// Default constructor
    TPZAlgebraicTransport();
    
    /// Copy constructor
    TPZAlgebraicTransport(const TPZAlgebraicTransport & other);
    
    /// Assignement constructor
    const TPZAlgebraicTransport & operator=(const TPZAlgebraicTransport & other);
    

    /// Default desconstructor
    ~TPZAlgebraicTransport();
    
    void BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh);
    void Assamble(TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef );
    void AssambleResidual(TPZFNMatrix<100, REAL> &ef);
    void Contribute();
    void ContributeInterface();
    void ContributeBCInterface();
   
    
    void CalcLambdas();
    void CalcDensities();
    double CalcLambda(double sw, double paverage, double densityo, double densityw);
    double CalcDensity(double paverage,double compress, double reff,  double pref);
    
    

};

#endif /* AlgebraicTransport_h */
