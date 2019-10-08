//
//  TPZMixedDarcyWithFourSpaces.h
//  reservoirlib
//
//  Created by Omar Durán and Jorge Ordoñez on 7/12/19.
//

#ifndef TPZMixedDarcyWithFourSpaces_h
#define TPZMixedDarcyWithFourSpaces_h

#include <stdio.h>
#include "TPZMixedDarcyFlow.h"

class TPZMixedDarcyWithFourSpaces : public TPZMixedDarcyFlow {
    
public:
    
    /** @brief default constructor  */
    TPZMixedDarcyWithFourSpaces();
    
    /** @brief default desconstructor  */
    ~TPZMixedDarcyWithFourSpaces();
    
    /** @brief constructor based on material id and dimension */
    TPZMixedDarcyWithFourSpaces(int mat_id, int dim);
    
    /** @brief constructor based on object copy */
    TPZMixedDarcyWithFourSpaces(const TPZMixedDarcyWithFourSpaces &other);
    
    /** @brief Volumetric contribute with jacobian matrix */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ef) override;
    
    /** @brief Variable index based on variable naming */
    int VariableIndex(const std::string &name) override;
    
    /** @brief size of the current variable (1 -> scalar, 3-> vector, 9 ->  Tensor ) */
    int NSolutionVariables(int var) override;
    
    /** @brief Postprocess required variables multiphysics */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
};


#endif /* TPZMixedDarcyWithFourSpaces_h */
