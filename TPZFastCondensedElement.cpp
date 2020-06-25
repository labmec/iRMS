//
//  TPZFastCondensedElement.cpp
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//


#include "TPZFastCondensedElement.h"
#include "pzcmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMixedDarcyWithFourSpaces.h"

bool TPZFastCondensedElement::fSkipLoadSolution = true;

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZFastCondensedElement::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    
    if(this->fMatrixComputed == false)
    {
     
        TPZMixedDarcyWithFourSpaces *matDarcy = dynamic_cast<TPZMixedDarcyWithFourSpaces *>(Material());
        if (!matDarcy) {
            DebugStop();
        }
        
        matDarcy->SetPermeability(fPermeabilityTensor);
        TPZCondensedCompEl::CalcStiff(ek, ef);
        
        ShrinkElementMatrix(ek, fEK);
        ShrinkElementMatrix(ef, fEF);
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
    
    TPZFMatrix<STATE> solvec(fEK.fMat.Rows(),1,0.);
    GetSolutionVector(solvec);
    ef.fMat *= -1.0*Glambda;
    /** @brief Computes z = alpha * opt(this)*x + beta * y */
    /** @note z and x cannot overlap in memory */
//    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
//                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
    STATE alpha = -1.;
    ek.fMat.MultAdd(solvec, ef.fMat, ef.fMat, alpha, 1.);
}

// extract the solution vector of the condensed element
void TPZFastCondensedElement::GetSolutionVector(TPZFMatrix<STATE> &solvec)
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
void TPZFastCondensedElement::CalcResidual(TPZElementMatrix &ef)
{
    //JOSE: ¿Como vai ser calculado o residuo agora?
    TPZCondensedCompEl::CalcResidual(ef);
    
}

void TPZFastCondensedElement::ShrinkElementMatrix(TPZElementMatrix &input, TPZElementMatrix &output)
{
    output.fType = input.fType;
    output.fMesh = input.fMesh;
    int nindep = 0;
    int64_t condense_size = 0;
    int ncon = input.fConnect.size();
    output.fConnect.Resize(ncon);
    int firstnon = -1;
    for (int ic = 0; ic<ncon; ic++) {
        int64_t cindex = input.fConnect[ic];
        TPZConnect &c = input.fMesh->ConnectVec()[cindex];
        if(c.IsCondensed() && firstnon != -1)
        {
            DebugStop();
        }
        if(c.IsCondensed()) continue;
        if(firstnon == -1) firstnon = ic;
        output.fConnect[nindep] = cindex;
        condense_size += c.NShape()*c.NState();
        nindep++;
    }
    output.fConnect.Resize(nindep);
    if(input.fType == TPZElementMatrix::EK)
    {
        output.fMat.Redim(condense_size,condense_size);
    } else
    {
        output.fMat.Redim(condense_size,input.fMat.Cols());
    }
    int row_orig = input.fMat.Rows();
    if(input.fType == TPZElementMatrix::EK)
    {
        input.fMat.GetSub(row_orig-condense_size, row_orig-condense_size, condense_size, condense_size, output.fMat);
    } else
    {
        input.fMat.GetSub(row_orig-condense_size, 0, condense_size,
                          input.fMat.Cols(), output.fMat);
    }
    output.fBlock.SetNBlocks(nindep);
    
    for (int ic = ncon-nindep; ic < ncon; ic++) {
        output.fBlock.Set(ic-ncon+nindep, input.fBlock.Size(ic));
    }
    output.fBlock.Resequence();
    
    
    
}

void TPZFastCondensedElement::SetLambda(REAL lambda){
    fLambda = lambda;
}
REAL TPZFastCondensedElement::GetLambda(){
    return fLambda;
}
void TPZFastCondensedElement::SetMixedDensity(REAL mdensity){
    fMixedDensity = mdensity;
}
REAL TPZFastCondensedElement::GetMixedDensity(){
    return fMixedDensity;
}
void TPZFastCondensedElement::SetPermTensorAndInv(TPZFNMatrix<9, REAL> &PermeabilityTensor, TPZFNMatrix<9, REAL> &InvPerm){
    fPermeabilityTensor = PermeabilityTensor;
    fInvPerm = InvPerm;
}
TPZFMatrix<REAL> &TPZFastCondensedElement::GetPermTensor(){
    return  fPermeabilityTensor;
}

/**
 * @brief Calculates the solution - sol - for the variable var
 * at point qsi, where qsi is expressed in terms of the
 * master element coordinates
 * @param qsi master element coordinate
 * @param var variable name
 * @param sol vetor for the solution
 */
void TPZFastCondensedElement::Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    switch (var) {
        case 7:
            sol[0] = fPermeabilityTensor(0,0);
            break;
        case 8:
            sol[0] = fPermeabilityTensor(1,1);
            break;
        case 9:
            sol[0] = fPermeabilityTensor(2,2);
            break;
        case 10:
            sol[0] = fLambda;
            break;
        default:
            TPZCondensedCompEl::Solution(qsi, var, sol);
            break;
    }
}

/** @brief Loads the solution within the internal data structure of the element */
/**
 * Is used to initialize the solution of connect objects with dependency \n
 * Is also used to load the solution within SuperElements
 */
void TPZFastCondensedElement::LoadSolution()
{
    if(fSkipLoadSolution)
    {
        TPZCompEl::LoadSolution();
    }
    else
    {
        TPZCondensedCompEl::LoadSolution();
    }
}

