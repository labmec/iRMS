//
//  TPZFastCondensedElement.cpp
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//


#include "TPZFastestCondensedElement.h"
#include "pzcmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMixedDarcyWithFourSpaces.h"
#include "TPZDarcyFlowWithMem.h"
#include "TPZSSpMatrixEigen.h"

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZFastestCondensedElement::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    
    if(this->fMatrixComputed == false)
    {
     
        TPZDarcyFlowWithMem *matDarcy = dynamic_cast<TPZDarcyFlowWithMem *>(Material());
        if (!matDarcy) {
            DebugStop();
        }
        
        TPZCondensedCompEl::CalcStiff(ek, ef);
//        ShrinkElementMatrix(ek, fEK);
//        ShrinkElementMatrix(ef, fEF);
        this->fMatrixComputed = true;
    }
    
    ek = fEK;
    ef = fEF;
    
    int nrows = ek.fMat.Rows();
    int ncols = ek.fMat.Rows();
    

    REAL Glambda = fMixedDensity;

    ek.fMat *= (1./fLambda);
    for (int icol=0; icol<ncols; icol++) {
        ek.fMat(nrows-1,icol) *= fLambda;
    }
    for (int irow=0; irow<nrows; irow++) {
        ek.fMat(irow,ncols-1) *= fLambda;
    }
    ek.fMat(nrows-1,ncols-1) *=fLambda;
//    ek.fMat(nrows-1,ncols-1) *=fCompressibilityMatrixTerm;
    
    TPZFMatrix<STATE> solvec(fEK.fMat.Rows(),1,0.);
    GetSolutionVector(solvec);
    

    ef.fMat *= -1.0*Glambda;
//    ef.fMat(nrows-1) = fCompressibiilityRhsTerm;
    
    
    
    /** @brief Computes z = alpha * opt(this)*x + beta * y */
    /** @note z and x cannot overlap in memory */
    //    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
    //                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
    STATE alpha = 1.;
    ek.fMat.MultAdd(solvec, ef.fMat, ef.fMat, alpha, 1);
 
}

// extract the solution vector of the condensed element
void TPZFastestCondensedElement::GetSolutionVector(TPZFMatrix<STATE> &solvec)
{
    int nc = fEK.fConnect.size();
    TPZCompMesh *cmesh = Mesh();
    int64_t vecsize = fEK.fMat.Rows();
    int count = 0;
    for(int ic=0; ic<nc; ic++)
    {
        int64_t cindex = fEK.fConnect[ic];
        TPZConnect &c = cmesh->ConnectVec()[cindex];
        int64_t seqnum = c.SequenceNumber();
        int blsize = c.NShape()*c.NState();
        for(int dof=0; dof<blsize; dof++)
        {
            solvec(count+dof,0) = cmesh->Block()(seqnum,0,dof,0);
        }
        count += blsize;
    }
    if(count != vecsize) DebugStop();
}

/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZFastestCondensedElement::CalcResidual(TPZElementMatrix &ef)
{
      TPZElementMatrix ek;
      CalcStiff(ek, ef);

}


