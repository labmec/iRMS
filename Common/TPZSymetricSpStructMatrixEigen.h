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
#include "TPZSSpStructMatrix.h"
#include "TPZSSpMatrixEigen.h"
#include "TPZFastCondensedElement.h"
/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSymetricSpStructMatrixEigen : public TPZStructMatrix {
    
    std::vector<Triplet3<REAL> > m_triplets;
    
public:

    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration fAsTotalCalcStifSub = tsim2-tsim2;
    boost::posix_time::time_duration fAsTotalAdkelsSub = tsim2-tsim2;
    boost::posix_time::time_duration fAsTotaAssembleSub = tsim2-tsim2;
    boost::posix_time::time_duration fAsTotaCondensedSub = tsim2-tsim2;
//    std::cout << "Mixed:: Assembly time Global " << deltat << std::endl;

    
    TPZSymetricSpStructMatrixEigen(TPZCompMesh *);
    
    virtual TPZMatrix<STATE> * Create();
    
    virtual TPZMatrix<STATE> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    
    using TPZStructMatrix::CreateAssemble;
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    virtual TPZStructMatrix * Clone();
    void Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    void Serial_AssembleSub(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    void Serial_AssembleGlob(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** Used only for testing */
    static int main();
    
    private :
    
    TPZSymetricSpStructMatrixEigen();
 
    void Swap(int64_t *a, int64_t *b)
    {
        int64_t aux = *a;
        *a = *b;
        *b = aux;
    }
    
    friend TPZPersistenceManager;
};

#endif /* TPZSymetricSpStructMatrixEigenEigen_h */
