//
//  TPZSSpMatrixEigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#ifndef TPZSSpMatrixEigen_hpp
#define TPZSSpMatrixEigen_hpp
#define EIGEN_SUPERLU_SUPPORT

#include <stdio.h>
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZSSpMatrixEigen.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
//#include <Eigen/SparseLDLt>
#include <Eigen/PardisoSupport>

//#include <Eigen/SuperLUSupport>

//#include <Eigen/>
#include "TPZAnalysisAuxEigen.h"
#include "TPZSpMatrixEigen.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include<Eigen/SparseCholesky>

#ifdef USING_MKL

#include "TPZPardisoControl.h"
#endif
/**
 * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
 * @ingroup matrix
 */


template<class TVar>
class TPZSYsmpMatrixEigen : public TPZMatrix<TVar>{
    
#ifdef USING_MKL
    friend class TPZPardisoControl<TVar>;
#endif
    
    public :
    /** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrixEigen();
    /** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrixEigen(const int64_t rows, const int64_t cols );
    /** @brief Copy constructor */
    TPZSYsmpMatrixEigen(const TPZSYsmpMatrixEigen &cp){
 
    }
//    TPZSYsmpMatrixEigen(const TPZSYsmpMatrixEigen<TVar> &cp) :
//    TPZRegisterClassId(&TPZSYsmpMatrixEigen::ClassId),
//    TPZMatrix<TVar>(cp), fIA(cp.fIA), fJA(cp.fJA), fA(cp.fA), fDiag(cp.fDiag)
//    {
//        
//       
//    }
    
    TPZSYsmpMatrixEigen &operator=(const TPZSYsmpMatrixEigen<TVar> &copy);
    
    CLONEDEF(TPZSYsmpMatrixEigen)
    /** @brief Destructor */
    virtual ~TPZSYsmpMatrixEigen();
    
    /** @brief Checks if the current matrix is symmetric */
    virtual int IsSimetric() const  override { return 1; }
    /** @brief Checks if current matrix is square */
    inline int IsSquare() const { return 1;}
    
    /** @brief Zeroes the matrix */
    virtual int Zero() override {
        fsparse_eigen = 0.0*fsparse_eigen;
        
        
        fA.Fill(0.);
        fDiag.Fill(0.);
#ifndef USING_MKL
//        TPZMatrix<TVar>::fDecomposed = ENoDecompose;
#endif
        return 0;
        
    }
    
    /** @brief Zeroes the matrix */
    virtual int Redim(int64_t rows, int64_t cols) override
    {
        if(rows == this->fRow && cols == this->fCol)
        {
            Zero();
        }
        else
        {
            DebugStop();
        }
    }
    
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
    
    /** @brief Get the matrix entry at (row,col) without bound checking */
    virtual const TVar &GetVal(const int64_t row, const int64_t col ) const override;
    
    /** @brief Put values without bounds checking \n
     *  This method is faster than "Put" if DEBUG is defined.
     */
    virtual int PutVal(const int64_t /*row*/,const int64_t /*col*/,const TVar & val ) override;
    
    
    /** @brief Computes z = beta * y + alpha * opt(this)*x */
    /** @note z and x cannot overlap in memory */
    virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                         const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
    
    /** @brief Sets data to the class */
    virtual void SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A );
    
    /// Access function for the coefficients
    TPZVec<TVar> &A()
    {
        return fA;
    }
    
    TPZVec<int64_t> &IA()
    {
        return fIA;
    }
    
    TPZVec<int64_t> &JA()
    {
        return fJA;
    }
    
    /** @brief Print the matrix along with a identification title */
    virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const override;
    
#ifdef USING_MKL
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
    virtual int Decompose_LDLt(std::list<int64_t> &singular) override;
    /** @brief Decomposes the current matrix using LDLt. */
    virtual int Decompose_LDLt() override;
    
    /** @brief Decomposes the current matrix using Cholesky method. The current matrix has to be symmetric. */
    virtual int Decompose_Cholesky() override;
    /**
     * @brief Decomposes the current matrix using Cholesky method.
     * @param singular
     */
    virtual int Decompose_Cholesky(std::list<int64_t> &singular)  override;
    
    
    /** @} */
    
    /**
     * @name Substitutions
     * @brief Substitutions forward and backward
     * @{
     */
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LForward( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LBackward( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
     * @param b right hand side and result after all
     */
    virtual int Subst_Diag( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Forward( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Backward( TPZFMatrix<TVar>* b ) const override;
    
    
    /** @} */
    
    
#endif
public:
    int ClassId() const override;
    
    void ComputeDiagonal();
    
    void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex) override;
    bool isNull( Eigen::SparseMatrix<REAL>& mat, int row, int col);
    Eigen::SparseMatrix<REAL> fsparse_eigen;
    mutable Eigen::SparseLU<Eigen::SparseMatrix<REAL>> fanalysis;
//    mutable Eigen::SuperLU<Eigen::SparseMatrix<REAL>> fanalysis;
//    Eigen::SuperLU<Eigen::SparseMatrix<double> > slu;
//    slu.compute(A);
//    x = slu.solve(b);
    void FromPZtoEigen(const TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat) const;
    
    void FromEigentoPZ( TPZFMatrix<TVar> &pzmat, Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> &eigenmat)const;
    virtual int Decompose_LU() override;
    virtual int Substitution( TPZFMatrix<TVar> * B ) const override;
   
private:
    
    
    TPZVec<int64_t>  fIA;
    TPZVec<int64_t>  fJA;
    TPZVec<TVar> fA;
    
#ifdef USING_MKL
//    TPZPardisoControl<TVar> fPardisoControl;
#endif
    
    TPZVec<TVar> fDiag;
};

template<class TVar>
inline void TPZSYsmpMatrixEigen<TVar>::SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A )
{
    //
    fsparse_eigen.setZero();
    
    int nrows = IA.size()-1;
    std::vector<Triplet3<REAL> > triplets(2*A.size() -nrows);
    int count =0;
    for (int irow = 0; irow < nrows; irow++) {
        for(int k=IA[irow]; k<IA[irow+1]; k++){
            int row= irow;
            int col = JA[k];
            REAL val = A[k];
            Triplet3<REAL> trip(row, col, val);
            if (row!=col) {
                Triplet3<REAL> trip(col, row, val);
                triplets[count] = trip;
                count++;
            }
            triplets[count] = trip;
            count++;
        }
    }
    fsparse_eigen.setFromTriplets(triplets.begin(), triplets.end());
    triplets.clear();
    fanalysis.analyzePattern(fsparse_eigen);
//    std::cout<<fsparse_eigen<<std::endl;
    //
    // Pass the data to the class.
    fIA = IA;
    fJA = JA;
    fA  =  A;
    ComputeDiagonal();
}

#endif /* TPZSSpMatrixEigen_hpp */
