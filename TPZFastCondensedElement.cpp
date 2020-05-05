//
//  TPZFastCondensedElement.cpp
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//


#include "TPZFastCondensedElement.h"
#include "pzcmesh.h"
/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZFastCondensedElement::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    if(this->fMatrixComputed == false)
    {
        TPZCondensedCompEl::CalcStiff(ek, ef);
        ShrinkElementMatrix(ek, fEK);
        ShrinkElementMatrix(ef, fEF);
        this->fMatrixComputed = true;
    }
    ek = fEK;
    ef = fEF;
    // THIS IS WRONG - ONLY THE FLUX EQUATIONS WILL BE MULTIPLIED
    ek.fMat *= (1./fPermeability);
    // THIS IS WRONG - Anxiously waiting for the theoretical development
    DebugStop();
    ef.fMat *= fSource;
    
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZFastCondensedElement::CalcResidual(TPZElementMatrix &ef)
{
    // dont know how to implement this. Lets see the theoretical derivation
    DebugStop();
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
    for (int ic = nindep; ic < ncon; ic++) {
        output.fBlock.Set(ic-nindep, input.fBlock.Size(ic));
    }
    output.fBlock.Resequence();
    
    // @TODO check if the rest of the datastructure is empty after calcstiff
    /*
    // @brief Vector of pointers to TPZConnect objects
    TPZStack<int64_t> fConnect;
    // @brief Pointer to a blocked matrix object
    TPZFNMatrix<1000, STATE> fMat;
    // @brief Block structure associated with fMat
    TPZBlock<STATE> fBlock;
    // @brief Vector of all nodes connected to the element
    TPZStack<int64_t> fConstrConnect;
    // @brief Pointer to the constrained matrix object
    TPZFNMatrix<1000, STATE> fConstrMat;
    // @brief Block structure associated with fConstrMat
    TPZBlock<STATE> fConstrBlock;
    
    TPZManVector<int64_t> fDestinationIndex, fSourceIndex;
    
    /// list of one degree of freedom restraints
    std::list<TPZOneShapeRestraint> fOneRestraints;
    */

}

