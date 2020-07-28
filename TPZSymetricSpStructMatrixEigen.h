//
//  TPZSymetricSpStructMatrixEigenEigen.h
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#ifndef TPZSymetricSpStructMatrixEigenEigen_h
#define TPZSymetricSpStructMatrixEigenEigen_h

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
class TPZSymetricSpStructMatrixEigen : public TPZStructMatrix {
    
public:
    
    TPZSymetricSpStructMatrixEigen(TPZCompMesh *);
    
    virtual TPZMatrix<STATE> * Create();
    
    virtual TPZMatrix<STATE> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    
    using TPZStructMatrix::CreateAssemble;
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    virtual TPZStructMatrix * Clone();
    
    
    /** Used only for testing */
    static int main();
    
    private :
    
    TPZSymetricSpStructMatrixEigen();
    
    friend TPZPersistenceManager;
};

#endif /* TPZSymetricSpStructMatrixEigenEigen_h */
