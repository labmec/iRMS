//
//  TPZMHMixedMeshControl.hpp
//  PZ
//


#ifndef TPZMHMixedMeshWithTransportControl_hpp
#define TPZMHMixedMeshWithTransportControl_hpp

#include <stdio.h>

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"

/// class for creating TPZMHMM with Mixed Meshes
class TPZMHMixedMeshWithTransportControl : public TPZMHMixedMeshControl
{
    
    TPZAutoPointer<TPZCompMesh> fcmeshPressureAverg;
    TPZAutoPointer<TPZCompMesh> fcmeshFluxAverg ;
    
public:
    
    TPZMHMixedMeshWithTransportControl() : TPZMHMixedMeshControl()
    {
        fcmeshPressureAverg = new TPZCompMesh;
        fcmeshFluxAverg = new TPZCompMesh;
    }
    
//    TPZMHMixedMesh4SpacesControl(int dimension):TPZMHMixedMeshControl(dimension){
//        
//    }
//    
//    
    TPZMHMixedMeshWithTransportControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices):TPZMHMixedMeshControl( gmesh, coarseindices){
        fcmeshPressureAverg = new TPZCompMesh(gmesh);
        fcmeshFluxAverg = new TPZCompMesh(gmesh);
    }
    TPZMHMixedMeshWithTransportControl(const TPZMHMixedMeshWithTransportControl &copy){
        TPZMHMixedMeshControl::operator=(copy);
        this->operator=(copy);
    }
    
    TPZMHMixedMeshWithTransportControl operator=(const TPZMHMixedMeshWithTransportControl &cp) {
        
        fcmeshPressureAverg = cp.fcmeshPressureAverg;
        fcmeshFluxAverg =cp.fcmeshFluxAverg;
        TPZMHMixedMeshControl::operator=(cp);
        return *this;
    }
//    
    TPZMHMixedMeshWithTransportControl(TPZAutoPointer<TPZGeoMesh> gmesh):TPZMHMixedMeshControl(gmesh)
    {
        
    }
//
//    
//    TPZMHMixedMesh4SpacesControl(const TPZMHMixedMesh4SpacesControl &copy) : TPZMHMixedMeshControl(copy)
//    {
//        
//        fFluxMesh = copy.fFluxMesh;
//    }
//    
//    TPZMHMixedMesh4SpacesControl &operator=(const TPZMHMixedMesh4SpacesControl &cp)
//    {
//        fFluxMesh = cp.fFluxMesh;
//        TPZMHMixedMeshControl::operator=(cp);
//        return *this;
//    }
//    
    void BuildComputationalMesh(bool usersubstructure);
    void CreateHDivPressureMHMMesh();
    void CreateAverageFlux();
    void CreateAveragePressure();
    void BuildMultiPhysicsMesh();
//
//    void HideTheElements();
//    
//    int64_t WhichSubdomain(TPZCompEl *cel);
    
};

#endif /* TPZMHMixedMeshChannelControl_hpp */
