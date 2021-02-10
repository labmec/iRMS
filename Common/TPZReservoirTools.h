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

#ifndef ReservoirTools_hpp
#define ReservoirTools_hpp

#include <stdio.h>

class TPZReservoirTools
{
public:
    static void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix);
    static void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix);
};





#endif /* ReservoirTools_hpp */
