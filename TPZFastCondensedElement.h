//
//  TPZFastCondensedElement.h
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//

#ifndef TPZFASTCONDENSEDELEMENT_h
#define TPZFASTCONDENSEDELEMENT_h

#include "pzcondensedcompel.h"
#include "pzelmat.h"

class TPZFastCondensedElement : public TPZCondensedCompEl
{
    // this will be the multiplying factor for the condensed stiffness matrix K11
    REAL fPermeability = 1.;
     
    // this constant contains the source term
    REAL fSource = 0.;
    // this flag indicates whether the matrix of the father element has been computed
    bool fMatrixComputed = false;
    
    // reference stifness matrix and rhs
    TPZElementMatrix fEK, fEF;
    
    void ShrinkElementMatrix(TPZElementMatrix &ekinput, TPZElementMatrix &output);
    
public:
    
    TPZFastCondensedElement(TPZCompEl *ref, bool keepmatrix = true) :
        TPZCondensedCompEl(ref,keepmatrix)
    {
        
        
    }
    
    /** @brief create a copy of the condensed computational element in the other mesh */
    TPZFastCondensedElement(const TPZFastCondensedElement &copy, TPZCompMesh &mesh) :
        TPZCondensedCompEl(copy, mesh)
    {
        fEK = copy.fEK;
        fEF = copy.fEF;
        fMatrixComputed = copy.fMatrixComputed;
    }
    
    
    virtual ~TPZFastCondensedElement()
    {
        
    }

    /**
     * @brief Computes the element stifness matrix and right hand side
     * @param ek element stiffness matrix
     * @param ef element load vector
     */
    virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) override;
    
    
    /**
     * @brief Computes the element right hand side
     * @param ef element load vector(s)
     */
    virtual void CalcResidual(TPZElementMatrix &ef) override;
    
    void SetPermeability(REAL perm);
    REAL GetPermeability();

    
};
#endif /* TPZFASTCONDENSEDELEMENT_h */
