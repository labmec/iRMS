//
//  TPZHdivTransfer.cpp
//  reservoirlib
//
//  Created by Jorge Paúl Ordóñez Andrade on 15/07/19.
//

#include "TPZHdivTransfer.h"


template<class TVar>
TPZHdivTransfer<TVar>::TPZHdivTransfer() : TPZMatrix<TVar>() {
}


template<class TVar>
TPZHdivTransfer<TVar>::TPZHdivTransfer(const TPZHdivTransfer<TVar> &cp) : TPZMatrix<TVar>(cp), fIndexes(cp.fIndexes) {
}

template<class TVar>
TPZHdivTransfer<TVar>::TPZHdivTransfer(int64_t rows, int64_t cols, TPZVec<int64_t> &Indexes) : TPZMatrix<TVar>(rows, cols) {
    fIndexes = Indexes;
}

template<class TVar>
TPZHdivTransfer<TVar>::~TPZHdivTransfer(){
    
}

template<class TVar>
void TPZHdivTransfer<TVar>::Print(const char *name, std::ostream &out , const MatrixOutputFormat form) const {
    DebugStop();
}

/**
 * @brief Gather a vector
 * @param y: Compressed vector which stores the information of the full vector
 * @param x: Vector which is going to provide the values for x
 */
template<class TVar>
void TPZHdivTransfer<TVar>::Gather(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x){
    

    
}

/**
 * @brief Scatter a vector
 * @param y: Full vector which stores the information of the compressed vector
 * @param x: Vector which is going to provide the values for y
 */
template<class TVar>
void TPZHdivTransfer<TVar>::Scatter(TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &x){
    DebugStop();
}


template<class TVar>
void TPZHdivTransfer<TVar>::SetIndexes(int64_t rows, int64_t cols, TPZVec<int64_t> &Indexes){
    TPZMatrix<TVar>::Resize(rows, cols);
    fIndexes = Indexes;
    int64_t r = Indexes.size();
    if (r != rows) {
        DebugStop();
    }
    
}


template<class TVar>
TPZVec<int64_t> & TPZHdivTransfer<TVar>::GetIndexes(){
    return fIndexes;
}

/**
 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
 * @param x Is x on the above operation
 * @param y Is y on the above operation
 * @param z Is z on the above operation
 * @param alpha Is alpha on the above operation
 * @param beta Is beta on the above operation
 * @param opt Indicates if is Transpose or not
 */
template<class TVar>
void TPZHdivTransfer<TVar>::MultAdd(const TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z, const TVar alpha, const TVar beta, const int opt) const {
    
    if (opt==0 && (y.Cols() != x.Cols() || y.Rows() != this->Rows() || y.Cols() != z.Cols() || y.Rows() != z.Rows())) {
        DebugStop();
    }
    if (opt && (y.Cols() != x.Cols() || y.Rows() != this->Cols() || y.Cols() != z.Cols() || y.Rows() != z.Rows())) {
        DebugStop();
    }
    int64_t rows = this->Rows();
    int64_t xcols = x.Cols();
    
    TPZMatrix<TVar>::PrepareZ (y,z,beta,opt);
    TVar val = 0.;

    if(opt == 0)
    {
        for (int64_t i = 0; i<rows; i++) {
            for (int64_t j = 0; j<xcols; j++) {
                z(i,j) += alpha*x.GetVal(fIndexes[i],j);
            }
        }
    }
    else
    {
        for (int64_t i = 0; i<rows; i++) {
            for (int64_t j = 0; j<xcols; j++) {
                z(fIndexes[i],j) += alpha*x.GetVal(i,j);
            }
        }
    }
}

template class TPZHdivTransfer<double>;
