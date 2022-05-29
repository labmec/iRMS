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
    
private:
        
    void AddMultiphysicsMaterialsToCompMesh(const int order, std::set<int> &MatsWithmem, std::set<int> &MatsWitOuthmem);
    void GetTransportMaterials(std::set<int> &MatsWithmem, std::set<int> &MatsWitOuthmem);
    // associate a lagrange multiplier level to the connects of the meshes
    void SetLagrangeMultiplier4Spaces(TPZVec<TPZCompMesh *>& mesh_vec);
    void AddAtomicMaterials(const int dim, TPZCompMesh* cmesh,
                            std::set<int>& matids,
                            std::set<int>& bcmatids,
                            const bool isInsertBCs = true);    
    void CreateFractureHDivCompMesh(TPZCompMesh* cmesh,
                                    std::set<int>& matids, std::set<int>& bcids,
                                    std::set<int>& matids_dim2, std::set<int>& bcids_dim2);
    
    // order overlapping fracture elements such that their normal follows their position in the original fractures
    void OrderOverlappingFractures();
    
    // order the overlapping fracture elements such that they correspond to the order of the fracture planes
    // create HDivBound glue elements between the fractures
    void OrderFractures(TPZCompMesh *cmesh, TPZVec<TPZGeoElSide> &fracvec);
    
    // create the H(div) spaces of the fracture elements
    void CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh);
    void SplitConnectsAtInterface(TPZCompElSide& compside);
    void AdjustOrientBoundaryEls(TPZCompMesh* cmesh, std::set<int>& buildmatids);
    
    int fInitMatIdForMergeMeshes = -1000;
    
public:
    
    TPZGeoMesh * mGeometry;
    
    TMRSDataTransfer mSimData;
    
    TPZMultiphysicsCompMesh * mMixedOperator;
    
    TPZCompMesh * mTransportOperator;
    
    /// this vector contains the MHM domain index for each geometric element
    TPZVec<int64_t> mSubdomainIndexGel;
    
    TPZHybridizeHDiv* mHybridizer;
 
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
    
    /// Atribute access methods
    const int& InitMatIdForMergeMeshes() const {return fInitMatIdForMergeMeshes;}
    int& InitMatIdForMergeMeshes() {return fInitMatIdForMergeMeshes;}
		
	/// Backwards compatibility attribute in case all fractures have the same matid
//    const int& FractureUniqueMatId() const {return mSimData.mTFracProperties.m_matid;}    
    
    /// For MHM
    /// Sets the geometry based on a fine and a coarse mesh. It creates a list of subdomains based on that
    /// and fills mSubdomainIndexGel vector that will be used to set the macro domains
    void SetGeometry(TPZGeoMesh * gmeshfine, TPZGeoMesh * gmeshcoarse);
    
	/// Checks if there is a skeleton element between volume elements of different domains
	void CheckMeshIntegrity(TPZGeoMesh* gmesh);
	
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
    // split the connects of the fluxmesh, create HDivBound elements and pressure elements
    void HybridizeIntersections(TPZVec<TPZCompMesh *>& mesh_vec);
    // assign a subdomain to the lower level elements
    void IdentifySubdomainForLowdimensionElements(TPZCompMesh *fluxmesh);
    
    // Adjust the neighbouring information such that the boundary of the fracture elements is the first boundary
    // verify if the assigned subdomains are consistent
    void VerifySubdomainIntegrity();
    
    // identify the domain indices of the interface elements
    void SetInterfaceDomains(TPZStack<int64_t> &pressureindices,std::pair<int,int> &interfacematids);

    void CreateIntersectionInterfaceElements(TPZVec<TPZCompMesh *>& meshvec_Hybrid);

    void CreateIntersectionInterfaceElements();
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
    
    /// create hybridized hdiv computational mesh
    TPZCompMesh * HdivFluxMeshHybridized(int order);
    
    /// create a discontinuous mesh
    TPZCompMesh * DiscontinuousCmesh(int order, char lagrange);
    
    /// group the connects of the discontinuous mesh such that all connects of a subdomain are identical
    void GroupConnectsBySubdomain(TPZCompMesh *cmesh);
    
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
    
	
	/// Builds the multiphysics cmesh and atomics cmeshes based on booleans previously set
	/// @param order approximation order
    void BuildMixedMultiPhysicsCompMesh(int order);
    
	
	/// Generates a pressure/flow problem with H(div) space for flow and L2 space for pressure.
	/// May 2022: Only works for problems without fractures and is only tested for 2D domains living in 3D
	/// @param order approximation order
    void BuildMixed2SpacesMultiPhysicsCompMesh(int order);
    
    void BuildMixed4SpacesMultiPhysicsCompMesh(int order);
        
    /// Creates the mesh with 4 spaces and 1 hybridization between the flux elements
    /// First only the 3D elements are hybridized. TODO: hybridize the fracture elements.
    void BuildMixed4SpacesHybridized(int order);
    
    /// build a multiphysics mesh corresponding to a zero order mortar approximation
    void BuildMixed4SpacesMortarMesh();
    
    /// insert wrapper elements necessary for creating the (hybridized) mortar spaces
    void InsertGeoWrappersForMortar();
    
    /// insert wrapper elements necessary for creating the hybridized spaces
    void InsertGeoWrappers();
    
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
    void CreateInterfaceElements(TPZGeoElSide &gelside, TPZGeoElSide &gelneig, int matid);
    void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);
    void TestSideOrient(TPZCompMesh *MultFlux);
    void TakeElementsbyID(std::map<int, std::vector<TPZGeoEl* >> &, std::vector<int> & );
    void VerifySideOrientsCoarseFine(TPZCompMesh *fluxCmesh);
    void HideTheElements(TPZCompMesh *cmesh);
    void BuildMultiphysicsSpaceWithMemoryByMatId(TPZVec<int> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector);
    
    void MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh);
    const bool isFracSim() const {return mSimData.mTGeometry.mDomainFracNameAndMatId.size();}
    
    bool IsFracMatId(int matiD);
    
    bool IsFracBCMatId(int matiD)
    {
        return IsFracMatId(matiD-1);
    }
    bool IsFracIntersectMatId(int matiD)
    {
        return IsFracMatId(matiD-2);
    }
};

#endif /* TMRSApproxSpaceGenerator_h */
