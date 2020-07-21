//
//  TPZSpStructMatrix_Eeigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//

#ifndef TPZSpStructMatrix_Eeigen_h
#define TPZSpStructMatrix_Eeigen_h

#include <stdio.h>
#include "pzstrmatrix.h"
#include "pzysmp.h"

#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSpStructMatrixEigen : public TPZStructMatrix {
public:
    
    TPZSpStructMatrixEigen(TPZCompMesh *);
    virtual TPZStructMatrix * Clone() override;
    virtual TPZMatrix<STATE> * Create() override;
    
    using TPZStructMatrix::CreateAssemble;
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
//    virtual TPZStructMatrix * Clone() override;
    
public:
    int ClassId() const override;
    
    private :
    TPZSpStructMatrixEigen();
    
    friend TPZPersistenceManager;
};
#endif /* TPZSpStructMatrix_Eeigen_hpp */
