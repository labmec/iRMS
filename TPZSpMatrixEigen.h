//
//  TPZSpMatrixEigen.hpp
//  ALL_BUILD
//
//  Created by Jose on 7/20/20.
//

#ifndef TPZSpMatrixEigen_hpp
#define TPZSpMatrixEigen_hpp

#include <stdio.h>
#include "pzanalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TMRSDataTransfer.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZAlgebraicTransport.h"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/PardisoSupport>
#include "TPZAnalysisAuxEigen.h"
#include "TPZSpMatrixEigen.h"

template<typename StorageIndex = typename Eigen::SparseMatrix<REAL>::StorageIndex >
class Triplet3
{
public:
    
    Triplet3() : m_row(0), m_col(0), m_value(0) {}
    
    Triplet3(const StorageIndex& i, const StorageIndex& j, const REAL& v = REAL(0))
    : m_row(i), m_col(j), m_value(v)
    {}
    
    const StorageIndex& row() const { return m_row; }
    const StorageIndex& col() const { return m_col; }
    const REAL & value() const { return m_value; }
    
protected:
    
    StorageIndex m_row, m_col;
    
    REAL m_value;
};

//
    template<class TVar>
    class TPZSpMatrixEigen : public TPZMatrix<TVar> {
        
        public :
        
        /** @brief An auxiliary structure to hold the data of the subset \n of equations used to multiply in a multi-threaded environment */
        /**
         In future versions this structure should be defined in a derived class
         */
        struct TPZMThread {
            const TPZSpMatrixEigen<TVar> *target;
            int64_t fFirsteq;
            int64_t fLasteq;
            const TPZFMatrix<TVar> *fX;
            TPZFMatrix<TVar> *fZ;
            TVar fAlpha;
            int fOpt;
        };
        
    private:
        
        static void * ExecuteMT(void *entrydata);
        
    public:
        Eigen::PardisoLU<Eigen::SparseMatrix<REAL>>  *f_solver;
        TPZSpMatrixEigen();
        int jose=1;
        TPZSpMatrixEigen(const int64_t rows,const int64_t cols );
        
        TPZSpMatrixEigen(const TPZVerySparseMatrix<TVar> &cp);
        
        TPZSpMatrixEigen &operator=(const TPZSpMatrixEigen<TVar> &copy);
//        TPZSpMatrixEigen &operator=(const TPZVerySparseMatrix<TVar> &cp);
        
        CLONEDEF(TPZSpMatrixEigen)
        
        virtual ~TPZSpMatrixEigen();
        
        /** @brief Fill matrix storage with randomic values */
        /** This method use GetVal and PutVal which are implemented by each type matrices */
        void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
        
        
        
        /** @brief Get the matrix entry at (row,col) without bound checking */
        virtual const TVar &GetVal(const int64_t row,const int64_t col ) const override;
        
        
        bool isNull(const Eigen::SparseMatrix<TVar>& mat, int row, int col);

        int64_t NumTerms()
        {
            return fIA[this->Rows()];
        }
        
        int PutVal(const int64_t row, const int64_t col, const TVar &Value) override;
        
        virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                             const TVar alpha=1.,const TVar beta = 0., const int opt = 0) const override;
        
        virtual void MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                               const TVar alpha=1.,const TVar beta = 0., const int opt = 0);
        
        virtual int GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
                           const int64_t colSize, TPZFMatrix<TVar> & A ) const override;
        
        void GetSub(const TPZVec<int64_t> &indices,TPZFMatrix<TVar> &block) const override;
        
        /** @brief Pass the data to the class. */
        virtual void SetData( int64_t *IA, int64_t *JA, TVar *A );
        
        /** @brief Pass the data to the class. */
        virtual void SetData( TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A );
        
        /** @brief Print the matrix along with a identification title */
        virtual void Print(const char *title, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
        
        /**
         * @name Solvers
         * @brief Linear system solvers. \n
         */
        /** For symmetric decompositions lower triangular matrix is used. \n
         * Solves a system A*X = B returning X in B
         */
        //@{
        /**
         * @brief Solves the linear system using Jacobi method. \n
         * @param numiterations The number of interations for the process.
         * @param F The right hand side of the system.
         * @param result The solution.
         * @param residual Returns F - A*U which is the solution residual.
         * @param scratch Available manipulation area on memory.
         * @param tol The tolerance value.
         * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
         */
        virtual void SolveJacobi(int64_t & numiterations, const TPZFMatrix<TVar> & F, TPZFMatrix<TVar> & result,
                                 TPZFMatrix<TVar> * residual, TPZFMatrix<TVar> & scratch, REAL & tol, const int FromCurrent = 0)  override;
        
        void SolveSOR(int64_t &numiterations, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
                      TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,
                      const REAL overrelax, REAL &tol,
                      const int FromCurrent = 0,const int direction = 1 )  override;
        // @}
        
        /**
         * @brief Add a contribution of a stiffness matrix
         * putting it on destination indexes position
         */
//        virtual void AddKelOld(
//                               TPZFMatrix<TVar> & elmat //! Member stiffness matrix beeing added
//                               , TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
//        );
        
        virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex) override;
        
        virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex) override;
        
        void MultiplyDummy(TPZSpMatrixEigen<TVar> & B, TPZSpMatrixEigen<TVar> & Res);
        
        virtual int Zero() override;
        
        /**
         * @name Factorization
         * @brief Those member functions are related to matrices factorization
         */
        //@{
        /**
         * @brief Decomposes the current matrix using LU decomposition.
         */
        virtual int Decompose_LU(std::list<int64_t> &singular) override;
        virtual int Decompose_LU() override;
        
        //@}
        
        /**
         * @name Substitutions
         * @brief Substitutions forward and backward
         */
        //@{
        /**
         * @brief Computes Forward and Backward substitution for a "LU" decomposed matrix.
         * @param B right hand side and result after all
         */
        virtual int Substitution( TPZFMatrix<TVar> * B ) const override;
        
        //@}
        
    public:
        int ClassId() const override;
        
    private:
        
        void ComputeDiagonal();
        
        /*
         * @brief Perform row update of the sparse matrix
         */
        void RowLUUpdate(int64_t sourcerow, int64_t destrow);
        
 
    public:
        Eigen::SparseMatrix<REAL> fsparse_eigen;
       
        Eigen::SparseMatrix<REAL> rhs;
        protected:
        TPZVec<int64_t>  fIA;
        TPZVec<int64_t>  fJA;
        TPZVec<TVar> fA;
        
        TPZVec<TVar> fDiag;
        
        int   fSymmetric;
        
#ifdef USING_MKL
//        TPZPardisoControl<TVar> fPardisoControl;
#endif
    protected:
        
        /**
         * @brief Implements a initialization method for the sparse structure. It sets the initial value for the fIA and fJA.
         */
        /**
         * -fIA will contain the initial positions for all the equations
         * -fJA will contain (-1) on all its positions
         * -fA will contain 0 on all its value
         */
        void InitializeData();
    };
    
    
    template<class TVar>
    inline void TPZSpMatrixEigen<TVar>::SetData( int64_t *IA, int64_t *JA, TVar *A ) {
        // Pass the data to the class.
        int nel = this->Rows()+1;
        fIA.resize(nel);
        memccpy(&fIA[0], IA, nel, sizeof(int64_t));
        int64_t nval = fIA[nel-1];
        fJA.resize(nval);
        memccpy(&fJA[0], JA, nval, sizeof(int64_t));
        fA.resize(nval);
        memccpy(&fA[0], A, nval, sizeof(TVar));
        ComputeDiagonal();
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
        Eigen::PardisoLU<Eigen::SparseMatrix<REAL>>  *solver;
        f_solver->analyzePattern(fsparse_eigen);
//        std::cout<<fsparse_eigen<<std::endl;
////        fsparse_eigen.insert(0,3)=4.0;
//        std::cout<<fsparse_eigen<<std::endl;
        
       
        std::cout<<"IA VECTOR: "<<std::endl;
        for (int i=0; i< IA.size(); i++) {
            std::cout<<IA[i]<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"JA VECTOR: "<<std::endl;
        for (int i=0; i< JA.size(); i++) {
            std::cout<<JA[i]<<std::endl;
        }
        std::cout<<std::endl;
        std::cout<<"A VECTOR: "<<std::endl;
        for (int i=0; i< A.size(); i++) {
            std::cout<<A[i]<<std::endl;
        }
        std::cout<<std::endl;
        
        std::cout<<"ni= "<<IA.size()<<" nj= "<<JA.size()<<" a= "<<A.size()<<std::endl;
        for(int ival =0; ival < IA.size(); ival++){
            for(int jval =0; jval < JA.size(); jval++){
                
//                tripletsinitial[count] = trip;
                std::cout<<"Triplet: i="<<ival<<" j= "<<jval<<" val= "<<ival+jval<<std::endl;
                count++;
            }
           
        }
//
//       
//        fsparse_eigen.setFromTriplets(tripletsinitial.begin(), tripletsinitial.end());
//        std::cout<<std::endl;
//        std::cout<<":o "<<fsparse_eigen.toDense();
//        sparse_eigen.AnalysePattern();
        if (IA.size() != this->Rows() + 1 ) {
            DebugStop();
        }
        
        if (JA.size() != IA[this->Rows()]) {
            DebugStop();
        }
        
        fIA = IA;
        fJA = JA;
        fA = A;
    }
;

    
#endif


