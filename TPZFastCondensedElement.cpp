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
    

//    REAL val = fMixedDensity;
    ek.fMat *= (1./fLambda);
    for (int icol=0; icol<ncols; icol++) {
        ek.fMat(nrows-1,icol) *= fLambda;
    }
    for (int irow=0; irow<nrows; irow++) {
        ek.fMat(irow,ncols-1) *= fLambda;
    }
    ek.fMat(nrows-1,ncols-1) *=fLambda;
    
    ef.fMat *= fSource;
    
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
TPZFNMatrix<9, REAL> TPZFastCondensedElement::GetPermTensor(){
    return  fPermeabilityTensor;
}
