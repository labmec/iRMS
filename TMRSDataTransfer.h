//
//  TMRSDataTransfer.hpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#ifndef TMRSDataTransfer_hpp
#define TMRSDataTransfer_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include "TMRSSavable.h"
#include "pzmanvector.h"
#include<tuple> // for tuple


/// Object that represents GUI state and store all the required input/output data
class TMRSDataTransfer : public TMRSSavable {
    
public:
    
    /// Default constructor
    TMRSDataTransfer();
    
    /// Copy constructor
    TMRSDataTransfer(const TMRSDataTransfer &other);
    
    // Copy assignment operator
    TMRSDataTransfer &operator=(const TMRSDataTransfer &other);
    
    /// Destructor
    ~TMRSDataTransfer();
    
    /// Write object state
    void Write(TPZStream &buf, int withclassid) const;
    
    /// Read object state
    void Read(TPZStream &buf, void *context);
    
    /// Read object state
    virtual int ClassId() const;
    
    
    class TGeometry : public TMRSSavable {
        
    public:
        
        std::vector<std::map<std::string,int>> mDomainDimNameAndPhysicalTag;
        
        std::vector<std::map<std::string,int>> mDomainFracDimNameAndPhysicalTag;

        
        TGeometry(){
            
            mDomainDimNameAndPhysicalTag.resize(4);
            mDomainFracDimNameAndPhysicalTag.resize(3);
        }
        
        ~TGeometry(){
            
        }
        
    
        TGeometry(const TGeometry &other){
            mDomainDimNameAndPhysicalTag = other.mDomainDimNameAndPhysicalTag;
            mDomainFracDimNameAndPhysicalTag = other.mDomainFracDimNameAndPhysicalTag;
        }
        
        TGeometry &operator=(const TGeometry &other){
            if (this != & other) // prevent self-assignment
            {
                mDomainDimNameAndPhysicalTag = other.mDomainDimNameAndPhysicalTag;
                mDomainFracDimNameAndPhysicalTag = other.mDomainFracDimNameAndPhysicalTag;
            }
            return *this;
        }
        
    };
    
    class TPetroPhysics : public TMRSSavable {
        
    public:
    
        
    };
    
    class TFluidProperties : public TMRSSavable {
        
    public:
        
    };
    
    class TBoundaryConditions : public TMRSSavable {
        
    public:
        
        TPZManVector<std::tuple<int, int, REAL>> mBCMixedPhysicalTagTypeValue;
        
        TPZManVector<std::tuple<int, int, REAL>> mBCMixedFracPhysicalTagTypeValue;
        
        TPZManVector<std::tuple<int, int, REAL>> mBCTransportPhysicalTagTypeValue;
        
        TPZManVector<std::tuple<int, int, REAL>> mBCTransportFracPhysicalTagTypeValue;
        
        TBoundaryConditions(){
            
            mBCMixedPhysicalTagTypeValue.Resize(0);
            
            mBCMixedFracPhysicalTagTypeValue.Resize(0);
            
            mBCTransportPhysicalTagTypeValue.Resize(0);
            
            mBCTransportFracPhysicalTagTypeValue.Resize(0);
        }
        
        ~TBoundaryConditions(){
            
        }
        
        TBoundaryConditions(const TBoundaryConditions &other){
            mBCMixedPhysicalTagTypeValue = other.mBCMixedPhysicalTagTypeValue;
            mBCMixedFracPhysicalTagTypeValue = other.mBCMixedFracPhysicalTagTypeValue;
            mBCTransportPhysicalTagTypeValue = other.mBCTransportPhysicalTagTypeValue;
            mBCTransportFracPhysicalTagTypeValue = other.mBCTransportFracPhysicalTagTypeValue;
        }
        
        TBoundaryConditions &operator=(const TBoundaryConditions &other){
            if (this != & other) // prevent self-assignment
            {
                mBCMixedPhysicalTagTypeValue = other.mBCMixedPhysicalTagTypeValue;
                mBCMixedFracPhysicalTagTypeValue = other.mBCMixedFracPhysicalTagTypeValue;
                mBCTransportPhysicalTagTypeValue = other.mBCTransportPhysicalTagTypeValue;
                mBCTransportFracPhysicalTagTypeValue = other.mBCTransportFracPhysicalTagTypeValue;
            }
            return *this;
        }
        
    };
    
    class TNumerics : public TMRSSavable {
        
    public:
        
    };
    
    class TPostProcess : public TMRSSavable {
        
    public:
        
    };
    
    TGeometry mTGeometry;
    
    TPetroPhysics mTPetroPhysics;
    
    TFluidProperties mTFluidProperties;
    
    TBoundaryConditions mTBoundaryConditions;
    
    TNumerics mTNumerics;
    
    TPostProcess mTPostProcess;
    
};

#endif /* TMRSDataTransfer_h */
