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
        int64_t fMatid;
        //vector for each multiplying coefficient of fluxes
        std::vector<std::vector<REAL> > fCoefficientsFlux;
        
        // left right volume index in the AlgebraicTransport data structure
        std::vector<std::pair<int64_t, int64_t>> fLeftRightVolIndex;
        
        //Integral of the flux functions associated with faces
        std::vector<REAL> fIntegralFluxFunctions;
        //  Integral of the flux
        std::vector<REAL>  fIntegralFlux;
        //Sign with respect to the natural orientation of the flux (+1 or -1)
        std::vector<REAL> fFluxSign;
        //Direction of the normal to the face
        std::vector<std::tuple<REAL,REAL,REAL>> fNormalFaceDirection;
        
        TInterfaceDataTransport() : fMatid(0), fCoefficientsFlux(0), fIntegralFluxFunctions(0), fIntegralFlux(0),fFluxSign(0), fNormalFaceDirection(0) {
           
        }
        TInterfaceDataTransport(const TInterfaceDataTransport &copy){
            fMatid = copy.fMatid;
            fCoefficientsFlux = copy.fCoefficientsFlux;
            fIntegralFluxFunctions = copy.fIntegralFluxFunctions;
            fIntegralFlux = copy.fIntegralFlux;
            fFluxSign= copy.fFluxSign;
            fNormalFaceDirection = copy.fNormalFaceDirection;
            
        }
        TInterfaceDataTransport &operator=(const TInterfaceDataTransport &copy)
        {
            fMatid = copy.fMatid;
            fCoefficientsFlux = copy.fCoefficientsFlux;
            fIntegralFluxFunctions = copy.fIntegralFluxFunctions;
            fIntegralFlux = copy.fIntegralFlux;
            fFluxSign= copy.fFluxSign;
            fNormalFaceDirection = copy.fNormalFaceDirection;
            
            return *this;
        }
        void Print(std::ostream &out);
    };
    
    // CELL DATA
    struct TCellData{
        int  fMatId;
        std::vector<REAL> fVolume;
        std::vector<REAL> fSaturation;
        std::vector<REAL> fPressure;
        std::vector<REAL> fDensityOil;
        std::vector<REAL> fDensityWater;
        std::vector<REAL> flambda;
        std::vector<REAL> fporosity;
        //Fractional flow and its derivative
        std::vector<std::vector<REAL>> fWaterfractionalflow;
        std::vector<std::vector<REAL>> fOilfractionalflow;
        std::vector<std::vector<REAL>> fCenterCordinate;
        
        // fCompressibility[0] = water, fCompressibility[1]=oil, fCompressibility[2]=gas etc.
        std::vector<REAL> fCompressibility;
        std::vector<REAL> fViscosity;
        std::vector<REAL> fReferencePressures;
        std::vector<REAL> fReferenceDensity;
        
        TCellData() : fMatId(-1), fVolume(0), fSaturation(0), fPressure(0), fDensityOil(0),fDensityWater(0), flambda(0),fporosity(0), fWaterfractionalflow(0),fOilfractionalflow(0), fCenterCordinate(0),
        fCompressibility(0),fViscosity(0),fReferencePressures(0),
        fReferenceDensity(0)
        {
            
        }
        TCellData(const TCellData &copy)
        {
            fVolume = copy.fVolume;
            fSaturation= copy.fSaturation;
            fPressure = copy.fPressure;
            fDensityOil = copy.fDensityOil;
            fDensityWater = copy.fDensityWater;
            flambda = copy.flambda;
            fWaterfractionalflow = copy.fWaterfractionalflow;
            fOilfractionalflow = copy.fOilfractionalflow;
            fCenterCordinate = copy.fCenterCordinate;
            fporosity = copy.fporosity;
            fCompressibility = copy.fCompressibility;
            fViscosity = copy.fViscosity;
            fReferencePressures= copy.fReferencePressures;
            fReferenceDensity = copy.fReferenceDensity;
        }
        TCellData &operator=(const TCellData &copy)
        {
            fVolume = copy.fVolume;
            fSaturation= copy.fSaturation;
            fPressure = copy.fPressure;
            fDensityOil = copy.fDensityOil;
            fDensityWater = copy.fDensityWater;
            flambda = copy.flambda;
            fWaterfractionalflow = copy.fWaterfractionalflow;
            fOilfractionalflow = copy.fOilfractionalflow;
            fCenterCordinate = copy.fCenterCordinate;
            fporosity = copy.fporosity;
            fCompressibility = copy.fCompressibility;
            fViscosity = copy.fViscosity;
            fReferencePressures= copy.fReferencePressures;
            fReferenceDensity = copy.fReferenceDensity;
            return *this;
        }
        void SetNumCells(int64_t ncells)
        {
            fVolume.resize(ncells);
            fSaturation.resize(ncells);
            fporosity.resize(ncells);
            fPressure.resize(ncells);
            fDensityOil.resize(ncells);
            fDensityWater.resize(ncells);
            flambda.resize(ncells);
            fWaterfractionalflow.resize(ncells);
            fOilfractionalflow.resize(ncells);
            fCenterCordinate.resize(ncells);
        }
        void UpdateFractionalFlowsAndLambda();
        void Print(std::ostream &out);
    };
    
    
    REAL fdt;
    int fNFluxCoefficients;
    
    //number of volumetric elements in the transport mesh
    int fNVolumesTransport = 0;
    // Cells data structure, one material at a time
    TCellData fCellsData;
    
    // Interface data structure, by material, element and side
    std::map<int, TInterfaceDataTransport> fInterfaceData;

    
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
    void Contribute(int index, TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef);
    void ContributeInterface(int index, TPZFNMatrix< 100, REAL> &ek,TPZFNMatrix<100, REAL> &ef);
    void ContributeBCInterface(int index, TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef);
    static std::pair<std::vector<REAL>,std::vector<REAL>> fwAndfoVal(REAL Sw, REAL rhoW,REAL rhoO);
    
    void CalcLambdas();
    void CalcDensities();
    double CalcLambda(double sw, double paverage, double densityo, double densityw);
    double CalcDensity(double paverage,double compress, double reff,  double pref);
    
    

};

#endif /* AlgebraicTransport_h */