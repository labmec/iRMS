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
#include "TMRSDarcyFlowWithMem_impl.h"
#include "TMRSMultiphaseFlow_impl.h"
#include "TMRSMemory.h"
#include "TMRSDataTransfer.h"
#include "TPZTracerFlow.h"
#include "pzl2projection.h"
#include "TPZCompMeshTools.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TMRSSFIAnalysis.h"
#include "TPZMHMixedMeshControl.h"
//#include "TPZMHMixedMeshWithTransportControl.h"
class TMRSApproxSpaceGenerator : public TMRSSavable {
    
public:
    
    TPZGeoMesh * mGeometry;
    
    TMRSDataTransfer mSimData;
    
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
    
    void LoadGeometry(std::string geometry_file);
    
    void CreateUniformMesh(int nx, REAL L, int ny=0, REAL h=0, int nz=0, REAL w=0);
    void GenerateMHMUniformMesh(int nelref);
    
    void PrintGeometry(std::string name);
    
    TPZGeoMesh * GetGeometry();
    
    void SetDataTransfer(TMRSDataTransfer & DataTransfer);
    
    TMRSDataTransfer & GetDataTransfer();
    
 
    
    TPZCompMesh * HdivFluxCmesh(int order);
    
    TPZCompMesh * DiscontinuousCmesh(int order = 0);
    
    void BuildMixedMultiPhysicsCompMesh(int order);
    
    void BuildMixed2SpacesMultiPhysicsCompMesh(int order);
    
    void BuildMixed4SpacesMultiPhysicsCompMesh(int order);
    
    void BuildMHMMixed2SpacesMultiPhysicsCompMesh();
    
    void BuildMHMMixed4SpacesMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh *CreateMixedOperatorMHM();
    
    void BuildMixedSCStructures();
    
    void BuildTransportMultiPhysicsCompMesh();
    
    void BuildTransport2SpacesMultiPhysicsCompMesh();
    
    void BuildTransport4SpacesMultiPhysicsCompMesh();
    void InsertMaterialObjects(TPZMHMixedMeshControl &control);
    TPZMultiphysicsCompMesh * GetMixedOperator();
    
    TPZMultiphysicsCompMesh * GetTransportOperator();
    
    // Linking the memory between the operators
    void LinkMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void AdjustMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void UnifyMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void FillMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void SetUpdateMaterialMemory(int material_id, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q = true);
    
    static void SetUpdateMemory(int dimension, TMRSDataTransfer & sim_data, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q = true);

    
};

#endif /* TMRSApproxSpaceGenerator_h */
