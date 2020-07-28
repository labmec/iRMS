//
//  TPZSSpMatrixEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#include "TPZSSpMatrixEigen.h"
/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrixEigen methods.
 */

#include <memory.h>

#include "pzsysmp.h"
#include "pzfmatrix.h"
#include "pzstack.h"

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrixEigen<TVar>::TPZSYsmpMatrixEigen() : TPZRegisterClassId(&TPZSYsmpMatrixEigen::ClassId),
TPZMatrix<TVar>() {

}

template<class TVar>
TPZSYsmpMatrixEigen<TVar>::TPZSYsmpMatrixEigen(const int64_t rows,const int64_t cols ) : TPZRegisterClassId(&TPZSYsmpMatrixEigen::ClassId),
TPZMatrix<TVar>(rows,cols) {
    Eigen::SparseMatrix<REAL> sparse_eigen(rows, cols);
    fsparse_eigen = sparse_eigen;
    fsparse_eigen.setZero();
    

}

template<class TVar>
TPZSYsmpMatrixEigen<TVar>::~TPZSYsmpMatrixEigen() {
    // Deletes everything associated with a TPZSYsmpMatrixEigen
#ifdef DESTRUCTOR
    cerr << "~TPZSYsmpMatrixEigen()\n";
#endif
}

template<class TVar>
TPZSYsmpMatrixEigen<TVar> &TPZSYsmpMatrixEigen<TVar>::operator=(const TPZSYsmpMatrixEigen<TVar> &copy)
{
    TPZMatrix<TVar>::operator=(copy);
    fsparse_eigen = copy.fsparse_eigen;
    fIA =copy.fIA;
    fJA = copy.fJA;
    fA = copy.fA;
    fDiag = copy.fDiag;
#ifdef USING_MKL
//    fPardisoControl = copy.fPardisoControl;
//    fPardisoControl.SetMatrix(this);
#endif
    return *this;
}


template<class TVar>
const TVar &TPZSYsmpMatrixEigen<TVar>::GetVal(const int64_t r,const int64_t c ) const {
    // Get the matrix entry at (row,col) without bound checking
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col ) return fA[ic];
    }
    return this->gZero;
}

/** @brief Put values without bounds checking \n
 *  This method is faster than "Put" if DEBUG is defined.
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & val )
{
    // Get the matrix entry at (row,col) without bound checking
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col )
        {
            fA[ic] = val;
            return 0;
        }
    }
    if (val != (TVar(0.))) {
        DebugStop();
    }
    return 0;
    
}


// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                   TPZFMatrix<TVar> &z,
                                   const TVar alpha,const TVar beta,const int opt) const {
    // computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot share storage
    int64_t  ir, ic;
    int64_t  r = (opt) ? this->Cols() : this->Rows();
    
    // Determine how to initialize z
    if(beta != 0) {
        z = y*beta;
    } else {
        z.Zero();
    }
    
    // Compute alpha * A * x
    int64_t ncols = x.Cols();
    for (int64_t col=0; col<ncols; col++)
    {
        for(int64_t ir=0; ir<this->Rows(); ir++) {
            for(int64_t ic=fIA[ir]; ic<fIA[ir+1]; ic++) {
                int64_t jc = fJA[ic];
                z(ir,col) += alpha * fA[ic] * x.g(jc,col);
                if(jc != ir)
                {
                    z(jc,col) += alpha * fA[ic] * x.g(ir,col);
                }
            }
        }
    }
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::Print(const char *title, std::ostream &out ,const MatrixOutputFormat form) const {
    // Print the matrix along with a identification title
    if(form == EInputFormat) {
        out << "\nTSYsmpMatrix Print: " << title << '\n'
        << "\tNon zero elements    = " << fA.size()  << '\n'
        << "\tRows    = " << this->Rows()  << '\n'
        << "\tColumns = " << this->Cols() << '\n';
        int i;
        out << "\tIA\tJA\tA\n"
        << "\t--\t--\t-\n";
        for(i=0; i<=this->Rows(); i++) {
            out << i      << '\t'
            << fIA[i] << '\t'
            << fJA[i] << '\t'
            << fA[i]  << '\n';
        }
        for(i=this->Rows()+1; i<fIA[this->Rows()]-1; i++) {
            out << i      << "\t\t"
            << fJA[i] << '\t'
            << fA[i]  << '\n';
        }
    } else {
        TPZMatrix<TVar>::Print(title,out,form);
    }
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::ComputeDiagonal() {
    if(!fDiag.size()) fDiag.resize(this->Rows());
    for(int ir=0; ir<this->Rows(); ir++) {
        fDiag[ir] = GetVal(ir,ir);
    }
}
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex){
    int64_t i,j = 0;
    TVar value=0.;
    int64_t ipos,jpos;
    for(i=0;i<sourceindex.NElements();i++){
        for(j=0;j<sourceindex.NElements();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(sourceindex[i],sourceindex[j]);
            if (!isNull(fsparse_eigen, ipos, jpos)) {
                fsparse_eigen.coeffRef(ipos, jpos) += value;
            }
            else{
                std::cout<<"Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;
            }
        }
    }
}
/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
{
    if (!symmetric || nrow != ncol) {
        DebugStop();
    }
    TPZFMatrix<TVar> orig;
    orig.AutoFill(nrow,ncol,symmetric);
    
    TPZVec<int64_t> IA(nrow+1);
    TPZStack<int64_t> JA;
    TPZStack<TVar> A;
    IA[0] = 0;
    TPZVec<std::set<int64_t> > eqs(nrow);
    for (int64_t row=0; row<nrow; row++) {
        eqs[row].insert(row);
        for (int64_t col = 0; col<ncol; col++) {
            REAL test = rand()*1./RAND_MAX;
            if (test > 0.5) {
                eqs[row].insert(col);
                if (symmetric) {
                    eqs[col].insert(row);
                }
            }
        }
    }
    int64_t pos=0;
    for (int64_t row=0; row< nrow; row++) {
        for (std::set<int64_t>::iterator col = eqs[row].begin(); col != eqs[row].end(); col++) {
            if(*col >= row)
            {
                JA.Push(*col);
                A.Push(orig(row,*col));
            }
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}

#ifdef USING_MKL

#include "TPZPardisoControl.h"
/**
 * @name Factorization
 * @brief Those member functions perform the matrix factorization
 * @{
 */


/**
 * @brief Decomposes the current matrix using LDLt. \n
 * The current matrix has to be symmetric.
 * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    Decompose_LDLt();
    return 1;
}
/** @brief Decomposes the current matrix using LDLt. */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_LDLt()
{
    if(this->IsDecomposed() == ELDLt) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
//    fPardisoControl.SetMatrixType(TPZPardisoControl<TVar>::ESymmetric,TPZPardisoControl<TVar>::EIndefinite);
//    fPardisoControl.Decompose();
    fanalysis.factorize(fsparse_eigen);
    this->SetIsDecomposed(ELDLt);
    return 1;
    
}

/** @brief Decomposes the current matrix using Cholesky method. The current matrix has to be symmetric. */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_Cholesky()
{
    if(this->IsDecomposed() == ECholesky) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    
//    fPardisoControl.SetMatrixType(TPZPardisoControl<TVar>::ESymmetric,TPZPardisoControl<TVar>::EPositiveDefinite);
//    fPardisoControl.Decompose();
    
    this->SetIsDecomposed(ECholesky);
    return 1;
}
/**
 * @brief Decomposes the current matrix using Cholesky method.
 * @param singular
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    return Decompose_Cholesky();
}



/** @} */


/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> rhs;
    FromPZtoEigen(x, rhs);
    fanalysis.analyzePattern(fsparse_eigen);
    fanalysis.factorize(fsparse_eigen);
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> ds = fanalysis.solve(rhs);
    std::cout<<ds<<std::endl;
    FromEigentoPZ(x, ds);
//    fPardisoControl.Solve(*b,x);
    *b = x;
    return 1;
}
//template<class TVar>
//void TPZSYsmpMatrixEigen<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
//                                     TPZFMatrix<TVar> &z,
//                                     const TVar alpha,const TVar beta,const int opt) const {
//
//    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> xeigen,yeigen,zeigen;
//    this->FromPZtoEigen(x, xeigen);
//    this->FromPZtoEigen(y, yeigen);
//    zeigen = beta*yeigen + alpha*xeigen;
//    this->FromEigentoPZ(z, zeigen);
//
//
//}
/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is lower triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
//    fPardisoControl.Solve(*b,x);
    *b = x;
    return 1;
}

/**
 * @brief Computes B = Y, where A*Y = B, A is upper triangular.
 * @param b right hand side and result after all
 */
template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}


#endif


template<class TVar>
int TPZSYsmpMatrixEigen<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrixEigen") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template<class TVar>
bool TPZSYsmpMatrixEigen<TVar>::isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col)
{
    for (Eigen::SparseMatrix<REAL>::InnerIterator it(mat, col); it; ++it) {
        if (it.row() == row) return false;
    }
    return false;
}
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
    int nrows = pzmat.Rows();
    int ncols = pzmat.Cols();
    if (nrows<0 || ncols<0) {
        DebugStop();
    }
    eigenmat.resize(nrows, ncols);
    eigenmat.setZero();
    for (int i=0; i< nrows; i++) {
        for (int j=0; j< ncols; j++) {
            eigenmat(i, j) = pzmat.Get(i, j);
        }
    }
}
template<class TVar>
void TPZSYsmpMatrixEigen<TVar>::FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
    int nrows = pzmat.Rows();
    int ncols = pzmat.Cols();
    if (nrows<0 || ncols<0) {
        DebugStop();
    }
    eigenmat.resize(nrows, ncols);
    for (int i=0; i< nrows; i++) {
        for (int j=0; j< ncols; j++) {
            pzmat(i, j)= eigenmat(i, j) ;
        };
    }
}
template class TPZSYsmpMatrixEigen<double>;
template class TPZSYsmpMatrixEigen<float>;
template class TPZSYsmpMatrixEigen<long double>;
