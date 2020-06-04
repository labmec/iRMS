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
#include "TPZMultiphysicsCompMesh.h"
#include "AlgebraicTransport.h"

class AlgebraicDataTransfer {
    
public:
    TPZMultiphysicsCompMesh *fFluxMesh;
    
    TPZMultiphysicsCompMesh *fTransportMesh;
    
    struct TInterfaceWithVolume
    {
        // geometric element index of the interface element
        int64_t fInterface_gelindex;
        // computational element index of the interface element
        int64_t fInterface_celindex;
        // left right geometric element index of the associated volume elements
        std::pair<int64_t, int64_t> fLeftRightGelIndex;
        // left right volume index in the AlgebraicTransport data structure
        std::pair<int64_t, int64_t> fLeftRightVolIndex;
        
        TInterfaceWithVolume() : fInterface_gelindex(-1), fInterface_celindex(-1), fLeftRightGelIndex(-1,-1), fLeftRightVolIndex(-1,-1)
        {
            
        }
    };
    // The indexes of the geometric elements corresponding to internal interfaces
    std::map<int,TPZVec<TInterfaceWithVolume>> fInterfaceGelIndexes;
    
    // The volume elements identified by material id
    std::map<int,TPZVec<int64_t>> fVolumeElements;
    
    
    
    public:
    
    /// Default constructor
    AlgebraicDataTransfer();
    
    /// Copy constructor
    AlgebraicDataTransfer(const AlgebraicDataTransfer & other);
    
    /// Assignement constructor
    const AlgebraicDataTransfer & operator=(const AlgebraicDataTransfer & other);
    
    /// Default desconstructor
    ~AlgebraicDataTransfer();
    
    void SetMeshes(TPZMultiphysicsCompMesh &fluxmesh, TPZMultiphysicsCompMesh &transportmesh)
    {
        fFluxMesh = &fluxmesh;
        fTransportMesh = &transportmesh;
    }
    
    // compute the data transfer data structures between the fluxmesh and transport class
    void BuildTransportDataStructure(AlgebraicTransport &transport);
    
    // Identify the geometric elements corresponding to interface elements. Order them as
    // a function of the number of corner nodes
    void IdentifyInterfaceGeometricElements();
    
    // Identify volume information to the interface data structure (TInterfaceWithVolume)
    void IdentifyVolumeGeometricElements();
    
    // print the datastructure
    void Print(std::ostream &out = std::cout);
};
#endif /* AlgebraicDataTransfer_h */
