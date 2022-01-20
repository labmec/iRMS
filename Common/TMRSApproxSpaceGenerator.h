//
//  TMRSApproxSpaceGenerator.h
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#ifndef TMRSApproxSpaceGenerator_h
#define TMRSApproxSpaceGenerator_h

#include "Projection/TPZL2Projection.h"
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
#include "TPZCompMeshTools.h"
#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "TMRSSFIAnalysis.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZHybridizeHDiv.h"

class TPZAlgebraicDataTransfer;

class TMRSApproxSpaceGenerator : public TMRSSavable {
    
public:
    
    TPZGeoMesh * mGeometry;
    
    /// subdomain index of each geometric element
    TPZVec<int64_t> mGelSubdomain;
    
    TMRSDataTransfer mSimData;
    
    TPZMultiphysicsCompMesh * mMixedOperator;
    
    TPZCompMesh * mTransportOperator;
    
    TPZVec<int64_t> mSubdomainIndexGel;
    
    TPZHybridizeHDiv* mHybridizer;
    
    int mMatIDIntersection = -1000000;
 
    ForcingFunctionBCType<STATE> mForcingFunctionBC;
    
public:
    
    TMRSApproxSpaceGenerator();
    
    // Copy assignment operator
    TMRSApproxSpaceGenerator &operator=(const TMRSApproxSpaceGenerator &other);
    
    // Copy constructor
    TMRSApproxSpaceGenerator(const TMRSApproxSpaceGenerator &copy)
    {
        DebugStop();
    }
    
    /// Destructor
    ~TMRSApproxSpaceGenerator();
    
    /// Write object state
    void Write(TPZStream &buf, int withclassid) const;
    
    /// Read object state
    void Read(TPZStream &buf, void *context);
    
    /// Read object state
    virtual int ClassId() const;
    
    void SetGeometry(TPZGeoMesh * geometry);
    
    void SetSubdomainIndexes(TPZVec<int64_t> &subIndexes){
        mSubdomainIndexGel =subIndexes;
    }
    TPZVec<int64_t> GetSubdomainIndexes(){
        return mSubdomainIndexGel;
    }
    
    void SetForcingFunctionBC(ForcingFunctionBCType<STATE> f){
        mForcingFunctionBC = f;
    }
    bool HasForcingFunctionBC() const {
        return (bool)mForcingFunctionBC;
    }
    const ForcingFunctionBCType<STATE> &ForcingFunctionBC() const
    {
        return mForcingFunctionBC;
    }
    ForcingFunctionBCType<STATE> &ForcingFunctionBC()
    {
        return mForcingFunctionBC;
    }
    
    const bool isThereFracIntersection() const;
    void HybridizeIntersections(TPZVec<TPZCompMesh *>& mesh_vec);
    void CreateIntersectionInterfaceElements(TPZVec<TPZCompMesh *>& meshvec_Hybrid);
    void DeleteBCsThatAreOnIntersect(TPZCompMesh* hdivcmesh);
    
    void LoadGeometry(std::string geometry_file);
    
    void LoadGeometry(std::string geometry_file,TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag);
    
    void CreateUniformMesh(int nx, REAL L, int ny=0, REAL h=0, int nz=0, REAL w=0);
    void ApplyUniformRefinement(int nelref);
    void ApplyUniformRefinement();
    
    void PrintGeometry(std::string name,  bool vtkFile=true, bool textfile=false);
    
    TPZGeoMesh * GetGeometry();
    
    void SetDataTransfer(TMRSDataTransfer & DataTransfer);
    
    TMRSDataTransfer & GetDataTransfer();
    
    const TPZHybridizeHDiv* Hybridizer() const {return mHybridizer;}
    TPZHybridizeHDiv* Hybridizer() {return mHybridizer;}
 
    /// create the HDiv computational mesh
    TPZCompMesh * HdivFluxCmesh(int order);
    
    /// create a discontinuous mesh
    TPZCompMesh * DiscontinuousCmesh(int order, char lagrange);
    
    /// create a discontinuous mesh
    TPZCompMesh * TransportCmesh();
    
    void  BuildAuxTransportCmesh();
    
    /// create an HDiv mesh for mortar approximation
    TPZCompMesh *HDivMortarFluxCmesh(char mortarlagrange);
    
    /// create a pressure with mortar elements
    TPZCompMesh *PressureMortarCmesh(char firstlagrangepressure, char lagrangepressure, char lagrangemortar);
    
    /// insert the necessary interface elements
    void InsertInterfaceElements();
    
    // build a discontinuous mesh with the order of the elements according to the
    // algebraic transport mesh. The order is always 0
    TPZCompMesh * DiscontinuousCmesh(TPZAlgebraicDataTransfer &Atransfer);
    
    void BuildMixedMultiPhysicsCompMesh(int order);
    
    void BuildMixed2SpacesMultiPhysicsCompMesh(int order);
    
    void BuildMixed4SpacesMultiPhysicsCompMesh(int order);
    
    /// build a multiphysics mesh corresponding to a zero order mortar approximation
    void BuildMixed4SpacesMortarMesh();
    
    /// insert wrapper elements necessary for creating the (hybridized) mortar spaces
    void InsertGeoWrappersForMortar();
    
    void GeoWrappersForMortarGelSide(TPZGeoElSide &gelside, std::set<int> bcids);
    int  FindNeighSubDomain(TPZGeoElSide &gelside);
    void BuildMHMMixed2SpacesMultiPhysicsCompMesh();
    
    void BuildMHMMixed4SpacesMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh *CreateMixedOperatorMHM();
    
    TPZMultiphysicsCompMesh *BuildAuxPosProcessCmesh(TPZAlgebraicDataTransfer &Atransfer);
    void BuildMixedSCStructures();
    
    void BuildTransportMultiPhysicsCompMesh();
    
    void BuildTransport2SpacesMultiPhysicsCompMesh();
    
    void BuildTransport4SpacesMultiPhysicsCompMesh();
    
    /// return the material ids and boundary condition material ids
    void GetMaterialIds(int dim, std::set<int> &matids, std::set<int> &bcmatids);
    
    void InsertMaterialObjects(TPZMHMixedMeshControl &control);

    TPZMultiphysicsCompMesh * GetMixedOperator();
    
    TPZCompMesh * GetTransportOperator();
    
    // Linking the memory between the operators
    void LinkMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZCompMesh * TransportOperator);
    
    static void AdjustMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void UnifyMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void FillMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator);
    
    static void FillMaterialMemoryDarcy(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZAlgebraicTransport *algebraicTranspor);
    
    static void SetUpdateMaterialMemory(int material_id, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q = true);
    
    static void SetUpdateMemory(int dimension, TMRSDataTransfer & sim_data, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q = true);

    void InitializeFracProperties(TPZMultiphysicsCompMesh * MixedOperator);
   
    //
    void findNeighElementbyMatId(TPZGeoElSide &gelside, std::vector<TPZGeoElSide > &neihside, std::set<int> VolMatIds);
    void CreateElementInterfaces(TPZGeoEl *gel);
    void CreateInterfaces(TPZCompMesh *cmesh);
    void CreateFracInterfaces(TPZGeoEl *gel);
    void CreateInterfaceElements(TPZGeoElSide &gelside, TPZGeoElSide &gelneig, int matid, bool IsAtomic);
    void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);
    void TestSideOrient(TPZCompMesh *MultFlux);
    void TakeElementsbyID(std::map<int, std::vector<TPZGeoEl* >> &, std::vector<int> & );
    void VerifySideOrientsCoarseFine(TPZCompMesh *fluxCmesh);
    void HideTheElements(TPZCompMesh *cmesh);
    void BuildMultiphysicsSpaceWithMemoryByMatId(TPZVec<int> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector);
    
    const int MatIDFracIntesect() const {return mMatIDIntersection;} 
    int& MatIDFracIntesect() {return mMatIDIntersection;}
};

#endif /* TMRSApproxSpaceGenerator_h */
