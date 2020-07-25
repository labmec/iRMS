//
//  TPZSpMatrixEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//

#include "TPZSpMatrixEigen.h"
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
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen(const TPZVerySparseMatrix<TVar> &cp) : TPZRegisterClassId(&TPZSpMatrixEigen::ClassId),
TPZMatrix<TVar>
()
{
    
    *this = cp;

}

template<class TVar>
TPZSpMatrixEigen<TVar>::TPZSpMatrixEigen() : TPZRegisterClassId(&TPZSpMatrixEigen::ClassId),
TPZMatrix<TVar>(), fIA(1,0),fJA(),fA(),fDiag(), fsparse_eigen(0,0)
{
     
}

template<class TVar>
TPZSpMatrixEigen<TVar> &TPZSpMatrixEigen<TVar>::operator=(const TPZSpMatrixEigen<TVar> &cp) {
    TPZMatrix<TVar>::operator=(cp);
    fsparse_eigen = cp.fsparse_eigen;
    fIA = cp.fIA;
    fA = cp.fA;
    fJA = cp.fJA;
    fDiag = cp.fDiag;

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
    // Constructs an empty TPZSpMatrixEigen
    //    fSolver = -1;
    Eigen::SparseMatrix<REAL> sparse_eigen(rows, cols);
    fsparse_eigen = sparse_eigen;
    fsparse_eigen.setZero();

    fSymmetric = 0;
    fDiag = 0;
    fA = 0;
    fIA = 0;
    fJA = 0;

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
void TPZSpMatrixEigen<TVar>::MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                     TPZFMatrix<TVar> &z,
                                     const TVar alpha,const TVar beta,const int opt )  {
    // computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot share storage
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || y.Rows() != z.Rows() )
    {
        cout << "\nERROR! in TPZVerySparseMatrix::MultiplyAdd : incompatible dimensions in x, y or z\n";
        return;
    }
    
    int64_t  ir, ic, icol, xcols;
    xcols = x.Cols();
    TVar sum;
    int64_t  r = (opt) ? this->Cols() : this->Rows();
    
    // Determine how to initialize z
    for(ic=0; ic<xcols; ic++) {
        TVar *zp = &(z(0,ic));
        if(beta != 0) {
            const TVar *yp = &(y.g(0,0));
            TVar *zlast = zp+r;
            
            if(&z != &y) {
                memcpy(zp,yp,r*sizeof(TVar));
            }
        } else {
            TVar *zp = &(z(0,0)), *zlast = zp+r;
            while(zp != zlast) {
                *zp = 0.;
                zp ++;
            }
        }
    }
    // Compute alpha * A * x
    if(xcols == 1 && opt == 0)
    {
        if(this->Cols() != x.Rows() || this->Rows() != y.Rows())
        {
            cout << "\nERROR! in TPZSpMatrixEigen::MultiplyAddMT: incompatible dimensions in opt=false\n";
            return;
        }
        for(ir=0; ir<r; ir++) {
            int64_t icolmin = fIA[ir];
            int64_t icolmax = fIA[ir+1];
            const TVar *xptr = &(x.g(0,0));
            TVar *Aptr = &fA[0];
            int64_t *JAptr = &fJA[0];
            for(sum = 0.0, icol=icolmin; icol<icolmax; icol++ ) {
                sum += Aptr[icol] * xptr[JAptr[icol]];
            }
            z(ir,0) += alpha * sum;
        }
    }
    else
    {
        for(ic=0; ic<xcols; ic++) {
            if(opt == 0) {
                
                for(ir=0; ir<this->Rows(); ir++) {
                    for(sum = 0.0, icol=fIA[ir]; icol<fIA[ir+1]; icol++ ) {
                        sum += fA[icol] * x.g((fJA[icol]),ic);
                    }
                    z(ir,ic) += alpha * sum;
                }
            }
            
            // Compute alpha * A^T * x
            else
            {
                if (this->Rows() != x.Rows() || this->Cols() != y.Rows())
                {
                    cout << "\nERROR! in TPZSpMatrixEigen::MultiplyAddMT: incompatible dimensions in opt=true\n";
                    return;
                }
                int64_t jc;
                int64_t icol;
                for(ir=0; ir<this->Rows(); ir++) {
                    for(icol=fIA[ir]; icol<fIA[ir+1]; icol++ ) {
                        if(fJA[icol]==-1) break; //Checa a exist�cia de dado ou n�
                        jc = fJA[icol];
                        TVar aval = fA[icol];
                        //cout << "FA["<<icol<<"] = "<<aval<< " * x["<< ir<<"] ="<< x.Get(ir,ic)<< endl;
                        z(jc,ic) += alpha * aval * x.g(ir,ic);
                    }
                }
            }
        }
    }
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
void TPZSpMatrixEigen<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
                                   TPZFMatrix<TVar> &z,
                                   const TVar alpha,const TVar beta,const int opt) const {
    // computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot share storage
    
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        std::cout << "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" ;
        return;
    }
    if(beta != (double)0. && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        std::cout << "TPZFMatrix::MultAdd matrix y with incompatible dimensions>";
        return;
    }
#endif
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0) {
        z.Zero();
        return;
    }
    
    int64_t  ic, xcols;
    xcols = x.Cols();
    int64_t  r = (opt) ? this->Cols() : this->Rows();
    
    // Determine how to initialize z
    for(ic=0; ic<xcols; ic++) {
        TVar *zp = &(z(0,ic));
        if(beta != 0) {
            const TVar *yp = &(y.g(0,0));
            TVar *zlast = zp+r;
            if(beta != 1.) {
                while(zp < zlast) {
                    *zp = beta * (*yp);
                    zp ++;
                    yp ++;
                }
            }
            else if(&z != &y) {
                memcpy(zp,yp,r*sizeof(REAL));
            }
        } else {
            TVar *zp = &(z(0,ic)), *zlast = zp+r;
            while(zp != zlast) {
                *zp = 0.;
                zp ++;
            }
        }
    }
    /*
     TPZSpMatrixEigen *target;
     int fFirsteq;
     int fLasteq;
     TPZFMatrix<>*fX;
     TPZFMatrix<>*fZ;
     REAL fAlpha;
     int fOpt;
     */
    int numthreads = 2;
    if(opt) numthreads = 1;
    TPZVec<pthread_t> allthreads(numthreads);
    TPZVec<TPZMThread> alldata(numthreads);
    TPZVec<int> res(numthreads);
    int i;
    int eqperthread = r/numthreads;
    int firsteq = 0;
    for(i=0;i<numthreads;i++)
    {
        alldata[i].target = this;
        alldata[i].fFirsteq = firsteq;
        alldata[i].fLasteq = firsteq+eqperthread;
        firsteq += eqperthread;
        if(i==numthreads-1) alldata[i].fLasteq = this->Rows();
        alldata[i].fX = &x;
        alldata[i].fZ = &z;
        alldata[i].fAlpha = alpha;
        alldata[i].fOpt = opt;
        res[i] = PZ_PTHREAD_CREATE(&allthreads[i], NULL,
                                   ExecuteMT, &alldata[i], __FUNCTION__);
    }
    for(i=0;i<numthreads;i++) {
        PZ_PTHREAD_JOIN(allthreads[i], NULL, __FUNCTION__);
    }
    
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

template<class TVar>
void TPZSpMatrixEigen<TVar>::ComputeDiagonal() {
    if(fDiag.size()) return;
    int64_t rows = fsparse_eigen.innerSize();
    fDiag.resize(rows);
    for(int64_t ir=0; ir<rows; ir++) {
        fDiag[ir] = fsparse_eigen.coeff(ir,ir);
    }
}


template<class TVar>
int TPZSpMatrixEigen<TVar>::Zero()
{
    fsparse_eigen = 0.0*fsparse_eigen;
//    fA.Fill(TVar(0.));
//    fDiag.Fill(TVar(0.));
    return 1;
}



template<class TVar>
void *TPZSpMatrixEigen<TVar>::ExecuteMT(void *entrydata)
{
    //DebugStop();
    TPZMThread *data = (TPZMThread *) entrydata;
    const TPZSpMatrixEigen *mat = data->target;
    TVar sum;
    int64_t xcols = data->fX->Cols();
    int64_t ic,ir,icol;
    // Compute alpha * A * x
    if(xcols == 1 && data->fOpt == 0)
    {
        for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
            int64_t icolmin = mat->fIA[ir];
            int64_t icolmax = mat->fIA[ir+1];
            const TVar *xptr = &(data->fX->g(0,0));
            TVar *Aptr = &mat->fA[0];
            int64_t *JAptr = &mat->fJA[0];
            for(sum = 0.0, icol=icolmin; icol<icolmax; icol++ ) {
                sum += Aptr[icol] * xptr[JAptr[icol]];
            }
            data->fZ->operator()(ir,0) += data->fAlpha * sum;
        }
    }
    else
    {
        for(ic=0; ic<xcols; ic++) {
            if(data->fOpt == 0) {
                
                for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
                    for(sum = 0.0, icol=mat->fIA[ir]; icol<mat->fIA[ir+1]; icol++ ) {
                        sum += mat->fA[icol] * data->fX->g((mat->fJA[icol]),ic);
                    }
                    data->fZ->operator()(ir,ic) += data->fAlpha * sum;
                }
            }
            
            // Compute alpha * A^T * x
            else
            {
                int64_t jc;
                int64_t icol;
                for(ir=data->fFirsteq; ir<data->fLasteq; ir++)
                {
                    for(icol=mat->fIA[ir]; icol<mat->fIA[ir+1]; icol++ )
                    {
                        if(mat->fJA[icol]==-1) break; //Checa a exist�cia de dado ou n�
                        jc = mat->fJA[icol];
                        data->fZ->operator()(jc,ic) += data->fAlpha * mat->fA[icol] * data->fX->g(ir,ic);
                    }
                }
                
            }
        }
    }
    return 0;
}
static int  FindCol(int64_t *colf, int64_t *coll, int64_t col)
{
    if(col == *colf) return 0;
    int64_t *begin = colf;
    int64_t *end = coll;
    while (begin != end)
    {
        int64_t dist = (end-begin)/2;
        int64_t *mid = begin+dist;
        if(*mid == col) return (mid-colf);
        else if(*mid > col) end=mid;
        else begin = mid;
    }
    return -1;
}

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


template<class TVar>
void TPZSpMatrixEigen<TVar>::GetSub(const TPZVec<int64_t> &indices,TPZFMatrix<TVar> &block) const
{
    std::map<int64_t,int64_t> indord;
    int64_t i,size = indices.NElements();
    for(i=0; i<size; i++)
    {
        indord[indices[i]] = i;
    }
    std::map<int64_t,int64_t>::iterator itset,jtset;
    for(itset = indord.begin(); itset != indord.end(); itset++)
    {
        int64_t *jfirst = &fJA[0]+fIA[(*itset).first];
        int64_t *jlast = &fJA[0]+fIA[(*itset).first+1]-1;
        //    int row = (*itset).first;
        for(jtset = indord.begin(); jtset != indord.end(); jtset++)
        {
            int64_t col = FindCol(jfirst,jlast,(*jtset).first);
            int64_t dist = jfirst+col-&fJA[0];
            block((*itset).second,(*jtset).second) = fA[dist];
            jfirst += col+1;
        }
    }
}

/*
 * Perform row update of the sparse matrix
 */

/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSpMatrixEigen<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
{
    if (symmetric && nrow != ncol) {
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
        if(nrow == ncol) eqs[row].insert(row);
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
            JA.Push(*col);
            A.Push(orig(row,*col));
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}


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
    Eigen::PardisoLU<Eigen::SparseMatrix<REAL>>  m_analysis;
    m_analysis.analyzePattern(fsparse_eigen);
    m_analysis.factorize(fsparse_eigen);
    
    this->SetIsDecomposed(ELU);
    return 1;
}

template<class TVar>
int TPZSpMatrixEigen<TVar>::Substitution( TPZFMatrix<TVar> *B ) const
{
    
    TPZFMatrix<TVar> x(*B);
//    x.Print(std::cout);
    Eigen::PardisoLU<Eigen::SparseMatrix<REAL>>  m_analysis;
    m_analysis.analyzePattern(fsparse_eigen);
    m_analysis.factorize(fsparse_eigen);
    int nrows = x.Rows();
    Eigen::SparseMatrix<REAL> rhs(nrows,1);
    for (int i=0; i< nrows; i++) {
        rhs.insert(i, 0) = x.Get(i, 0);
    }
//    std::cout<<r_hs.toDense()<<std::endl;
    Eigen::SparseMatrix<REAL> dsol = m_analysis.solve(rhs);
//    std::cout<<dsol.toDense()<<std::endl;
//    fPardisoControl.Solve(*B,x);
    for (int i=0; i< nrows; i++) {
        x(i, 0) = dsol.coeff(i, 0);
    }
    
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
template class TPZSpMatrixEigen<long double>;
template class TPZSpMatrixEigen<double>;
template class TPZSpMatrixEigen<float>;


