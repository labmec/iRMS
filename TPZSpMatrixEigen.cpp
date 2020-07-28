//
//  TPZSpMatrixEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//

#include "TPZSpMatrixEigen.h"
#include "TPZSSpMatrixEigen.h"
#include "pzysmp.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <memory.h>
#include <string>
#include <map>
#include <pthread.h>

#include "tpzverysparsematrix.h"
#include "pz_pthread.h"
#include "pzstack.h"

using namespace std;

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************
template<class TVar>
void TPZSpMatrixEigen<TVar>::InitializeData(){}



// ****************************************************************************
//
// Constructor
//
// ****************************************************************************

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen(const TPZSpMatrixEigen &cp){
    
}

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen() : TPZRegisterClassId(&TPZSpMatrixEigen::ClassId),
TPZMatrix<TVar>(),  fsparse_eigen(0,0)
{
     
}

template<class TVar>
TPZSpMatrixEigen<TVar> &TPZSpMatrixEigen<TVar>::operator=(const TPZSpMatrixEigen<TVar> &cp) {
    TPZMatrix<TVar>::operator=(cp);
    fsparse_eigen = cp.fsparse_eigen;
    return *this;
}


template<class TVar>
int TPZSpMatrixEigen<TVar>::PutVal(const int64_t row, const int64_t col, const TVar &Value){
//    std::cout<<fsparse_eigen.toDense()<<std::endl;
    if (!isNull(fsparse_eigen, row, col)) {
        fsparse_eigen.coeffRef(row,col)=Value;
    }
    else{
        std::cout<<"Non existing position on sparse matrix: line =" << row << " column =" << col << std::endl;
    }
    
    return 1;
}
template<class TVar>
void TPZSpMatrixEigen<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex){
    
    int64_t i,j;
    TVar value=0.;
    int64_t ipos,jpos;
    for(i=0;i<elmat.Rows();i++){
        for(j=0;j<elmat.Rows();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(i,j);
            if (!isNull(fsparse_eigen, ipos, jpos)) {
                fsparse_eigen.coeffRef(ipos, jpos) += value;
            }
            else{
                std::cout<<"Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;
            }
        }
    }
}

template<class TVar>
void TPZSpMatrixEigen<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex){
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

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen(const int64_t rows,const int64_t cols ) :
TPZRegisterClassId(&TPZSpMatrixEigen::ClassId),TPZMatrix<TVar>(rows,cols) {
  
    Eigen::SparseMatrix<REAL> sparse_eigen(rows, cols);
    fsparse_eigen = sparse_eigen;
    fsparse_eigen.setZero();
    fSymmetric = 0;
   
#ifdef CONSTRUCTOR
    cerr << "TPZSpMatrixEigen(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZSpMatrixEigen<TVar>::~TPZSpMatrixEigen() {
    // Deletes everything associated with a TPZSpMatrixEigen
#ifdef DESTRUCTOR
    cerr << "~TPZSpMatrixEigen()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

template<class TVar>
const TVar & TPZSpMatrixEigen<TVar>::GetVal(const int64_t row,const int64_t col ) const {
    TVar val = fsparse_eigen.coeff(row, col);
    return val;

}

template<class TVar>
void TPZSpMatrixEigen<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                   TPZFMatrix<TVar> &z,
                                   const TVar alpha,const TVar beta,const int opt) const {
    
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> xeigen,yeigen,zeigen;
    this->FromPZtoEigen(x, xeigen);
    this->FromPZtoEigen(y, yeigen);
    zeigen = beta*yeigen + alpha*xeigen;
    this->FromEigentoPZ(z, zeigen);

    
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZSpMatrixEigen<TVar>::Print(const char *title, ostream &out ,const MatrixOutputFormat form) const {
    // Print the matrix along with a identification title
    if(form != EInputFormat) {
        out << "\nSparse Eigen Matrix Print: " << title << '\n';
        out<<fsparse_eigen<<'\n';
    }
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

//template<class TVar>
//void TPZSpMatrixEigen<TVar>::ComputeDiagonal() {
//    if(fDiag.size()) return;
//    int64_t rows = fsparse_eigen.innerSize();
//    fDiag.resize(rows);
//    for(int64_t ir=0; ir<rows; ir++) {
//        fDiag[ir] = fsparse_eigen.coeff(ir,ir);
//    }
//}


template<class TVar>
int TPZSpMatrixEigen<TVar>::Zero()
{
    fsparse_eigen = 0.0*fsparse_eigen;

    return 1;
}



//template<class TVar>
//void *TPZSpMatrixEigen<TVar>::ExecuteMT(void *entrydata)
//{
//    //DebugStop();
//    TPZMThread *data = (TPZMThread *) entrydata;
//    const TPZSpMatrixEigen *mat = data->target;
//    TVar sum;
//    int64_t xcols = data->fX->Cols();
//    int64_t ic,ir,icol;
//    // Compute alpha * A * x
//    if(xcols == 1 && data->fOpt == 0)
//    {
//        for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
//            int64_t icolmin = mat->fIA[ir];
//            int64_t icolmax = mat->fIA[ir+1];
//            const TVar *xptr = &(data->fX->g(0,0));
//            TVar *Aptr = &mat->fA[0];
//            int64_t *JAptr = &mat->fJA[0];
//            for(sum = 0.0, icol=icolmin; icol<icolmax; icol++ ) {
//                sum += Aptr[icol] * xptr[JAptr[icol]];
//            }
//            data->fZ->operator()(ir,0) += data->fAlpha * sum;
//        }
//    }
//    else
//    {
//        for(ic=0; ic<xcols; ic++) {
//            if(data->fOpt == 0) {
//
//                for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
//                    for(sum = 0.0, icol=mat->fIA[ir]; icol<mat->fIA[ir+1]; icol++ ) {
//                        sum += mat->fA[icol] * data->fX->g((mat->fJA[icol]),ic);
//                    }
//                    data->fZ->operator()(ir,ic) += data->fAlpha * sum;
//                }
//            }
//
//            // Compute alpha * A^T * x
//            else
//            {
//                int64_t jc;
//                int64_t icol;
//                for(ir=data->fFirsteq; ir<data->fLasteq; ir++)
//                {
//                    for(icol=mat->fIA[ir]; icol<mat->fIA[ir+1]; icol++ )
//                    {
//                        if(mat->fJA[icol]==-1) break; //Checa a exist�cia de dado ou n�
//                        jc = mat->fJA[icol];
//                        data->fZ->operator()(jc,ic) += data->fAlpha * mat->fA[icol] * data->fX->g(ir,ic);
//                    }
//                }
//
//            }
//        }
//    }
//    return 0;
//}

template<class TVar>
int TPZSpMatrixEigen<TVar>::GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
                                 const int64_t colSize, TPZFMatrix<TVar> & A ) const {
    
    Eigen::SparseMatrix<REAL> sub = fsparse_eigen.block(sRow, sCol, rowSize, colSize);
    A.Resize(sub.innerSize(), sub.outerSize());
    A.Zero();
    for (int col=0; col<sub.outerSize(); col++) {
        for (Eigen::SparseMatrix<REAL>::InnerIterator it(sub, col); it; ++it) {
            A(it.row(),col)=sub.coeffRef(it.row(), col);
        }
    }

    return 0;
}


/*
 * Perform row update of the sparse matrix
 */

/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
//template<class TVar>
//void TPZSpMatrixEigen<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
//{
//    if (symmetric && nrow != ncol) {
//        DebugStop();
//    }
//    TPZFMatrix<TVar> orig;
//    orig.AutoFill(nrow,ncol,symmetric);
//    
//    TPZVec<int64_t> IA(nrow+1);
//    TPZStack<int64_t> JA;
//    TPZStack<TVar> A;
//    IA[0] = 0;
//    TPZVec<std::set<int64_t> > eqs(nrow);
//    for (int64_t row=0; row<nrow; row++) {
//        if(nrow == ncol) eqs[row].insert(row);
//        for (int64_t col = 0; col<ncol; col++) {
//            REAL test = rand()*1./RAND_MAX;
//            if (test > 0.5) {
//                eqs[row].insert(col);
//                if (symmetric) {
//                    eqs[col].insert(row);
//                }
//            }
//        }
//    }
// 
//    for (int64_t row=0; row< nrow; row++) {
//        for (std::set<int64_t>::iterator col = eqs[row].begin(); col != eqs[row].end(); col++) {
//            JA.Push(*col);
//            A.Push(orig(row,*col));
//        }
//        IA[row+1] = JA.size();
//    }
//    TPZMatrix<TVar>::Resize(nrow,ncol);
//    SetData(IA, JA, A);
//}


/**
 * Decomposes the current matrix using LU decomposition.
 */
template<class TVar>
int TPZSpMatrixEigen<TVar>::Decompose_LU(std::list<int64_t> &singular)
{
    return Decompose_LU();
}
template<class TVar>
int TPZSpMatrixEigen<TVar>::Decompose_LU()
{
    if(this->IsDecomposed() == ELU) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    
    m_analysis.analyzePattern(fsparse_eigen);
    m_analysis.factorize(fsparse_eigen);
    
    this->SetIsDecomposed(ELU);
    return 1;
}

template<class TVar>
int TPZSpMatrixEigen<TVar>::Substitution( TPZFMatrix<TVar> *B ) const
{
    TPZFMatrix<TVar> x(*B);
    int nrows = x.Rows();
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> rhs(nrows,1);
    this->FromPZtoEigen(x, rhs);
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> dsol = m_analysis.solve(rhs);
    this->FromEigentoPZ(x, dsol);
    *B = x;
    return 1;
}
template<class TVar>
bool TPZSpMatrixEigen<TVar>::isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col)
{
    for (Eigen::SparseMatrix<REAL>::InnerIterator it(mat, col); it; ++it) {
        if (it.row() == row) return false;
    }
    return true;
}

template<class TVar>
int TPZSpMatrixEigen<TVar>::ClassId() const{
    return Hash("TPZSpMatrixEigen") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template<class TVar>
void TPZSpMatrixEigen<TVar>::FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
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
void TPZSpMatrixEigen<TVar>::FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const{
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
/** @brief Pass the data to the class. */
template<class TVar>
inline void TPZSpMatrixEigen<TVar>::SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A ){
    fsparse_eigen.setZero();
    std::vector<Triplet3<REAL> > triplets(A.size());
    int nrows = IA.size()-1;
    int count =0;
    for (int irow = 0; irow < nrows; irow++) {
        for(int k=IA[irow]; k<IA[irow+1]; k++){
            int row= irow;
            int col = JA[k];
            REAL val = A[k];
            Triplet3<REAL> trip(row, col, val);
            triplets[count] = trip;
            count++;
        }
    }
    fsparse_eigen.setFromTriplets(triplets.begin(), triplets.end());
    //        triplets.clear();
    m_analysis.analyzePattern(fsparse_eigen);
    
    
    if (IA.size() != this->Rows() + 1 ) {
        DebugStop();
    }
    
    if (JA.size() != IA[this->Rows()]) {
        DebugStop();
    }
    
}
;
template class TPZSpMatrixEigen<long double>;
template class TPZSpMatrixEigen<double>;
template class TPZSpMatrixEigen<float>;
