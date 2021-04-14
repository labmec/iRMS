//
//  ReservoirTools.hpp
//  Buckley_Levertt
//
//  Created by Jose on 5/26/20.
//
#include "TPZCompMeshTools.h"
#include "pzelchdiv.h"
#include "pzshapepiram.h"
#include "TPZOneShapeRestraint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzelementgroup.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzmultiphysicselement.h"
#include "TPZMeshSolution.h"
#include "TPZFastCondensedElement.h"
#include "TPZMultiphysicsCompMesh.h"

#ifndef ReservoirTools_hpp
#define ReservoirTools_hpp

#include <stdio.h>

class TPZReservoirTools
{
public:
    static void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix);
    static void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix);
    static void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix, std::set<int> matids);
    static void AddDependency(std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons);
    static void TakeFatherSonsCorrespondence(TPZCompMesh *fluxCmesh,  std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons);
    static void TakeElementsbyID(TPZGeoMesh *mGeometry, std::map<int, std::vector<TPZGeoEl* >> &interfaces, std::vector<int> &matIds);
    static void FindCondensed(TPZCompEl *cel, TPZStack<TPZFastCondensedElement *> &condensedelements);
};





#endif /* ReservoirTools_hpp */
