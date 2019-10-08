//
//  TMRSDataTransfer.hpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#ifndef TMRSDataTransfer_hpp
#define TMRSDataTransfer_hpp

#include <stdio.h>
#include "TMRSSavable.h"

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
        
    };
    
    class TPetroPhysics : public TMRSSavable {
        
    public:
        
    };
    
    class TFluidProperties : public TMRSSavable {
        
    public:
        
    };
    
    class TBoundaryConditions : public TMRSSavable {
        
    public:
        
    };
    
    class TNumerics : public TMRSSavable {
        
    public:
        
    };
    
    class TPostProcess : public TMRSSavable {
        
    public:
        
    };
    
};

#endif /* TMRSDataTransfer_h */
