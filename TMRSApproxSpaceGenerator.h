//
//  TMRSApproxSpaceGenerator.h
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#ifndef TMRSApproxSpaceGenerator_h
#define TMRSApproxSpaceGenerator_h

#include <stdio.h>
#include "TMRSSavable.h"
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMixedDarcyFlow.h"
#include "TMRSDataTransfer.h"
#include "TPZTracerFlow.h"
#include "pzl2projection.h"

class TMRSApproxSpaceGenerator : public TMRSSavable {
    
private:
    
    TPZGeoMesh * mGeometry;
    
    TMRSDataTransfer mDataTransfer;
    
    TPZMultiphysicsCompMesh * mMixedOperator;
    
    TPZMultiphysicsCompMesh * mTransportOperator;
    
public:
    
    TMRSApproxSpaceGenerator();
    
    // Copy assignment operator
    TMRSApproxSpaceGenerator &operator=(const TMRSApproxSpaceGenerator &other);
    
    /// Destructor
    ~TMRSApproxSpaceGenerator();
    
    /// Write object state
    void Write(TPZStream &buf, int withclassid) const;
    
    /// Read object state
    void Read(TPZStream &buf, void *context);
    
    /// Read object state
    virtual int ClassId() const;
    
    void SetGeometry(TPZGeoMesh * geometry);
    
    TPZGeoMesh * GetGeometry();
    
    void SetDataTransfer(TMRSDataTransfer & DataTransfer);
    
    TMRSDataTransfer & GetDataTransfer();
    
    void LoadGeometry(std::string geometry_file);
    
    void PrintGeometry(std::string name);
    
    TPZCompMesh * HdivFluxCmesh(int order);
    
    TPZCompMesh * DiscontinuousCmesh(int order = 0);
    
    void BuildMixedMultiPhysicsCompMesh(int order);
    
    void BuildTransportMultiPhysicsCompMesh();
    
    TPZMultiphysicsCompMesh * GetMixedOperator();
    
    TPZMultiphysicsCompMesh * GetTransportOperator();
    
};

#endif /* TMRSApproxSpaceGenerator_h */
