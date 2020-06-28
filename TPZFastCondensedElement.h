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
    REAL fLambda = 1.0;
    
    TPZFNMatrix<9, REAL> fPermeabilityTensor;
    TPZFNMatrix<9, REAL> fInvPerm;
    
    //mixture density = rhow*fw + rhoo*fo
    REAL fMixedDensity = 1.;
     
    // this constant contains the source term
    REAL fSource = 0.;
    // this flag indicates whether the matrix of the father element has been computed
    bool fMatrixComputed = false;
    
    // reference stifness matrix and rhs
    TPZElementMatrix fEK, fEF;
    
public:
    static bool fSkipLoadSolution;
    
private:
    
    void ShrinkElementMatrix(TPZElementMatrix &ekinput, TPZElementMatrix &output);
    
    // extract the solution vector of the condensed element
    void GetSolutionVector(TPZFMatrix<STATE> &solvec);
    
    // global indices of the pressure equations
    void PressureEquations(TPZVec<int64_t> &eqs);
    
    // global index of the average pressure equation
    int64_t AveragePressureEquation();
    
    // global indices of the boundary flux equations
    void BoundaryFluxEquations(TPZVec<int64_t> &eqs);
    
    // adjust the multiplying coeficients of the pressure equations
    void AdjustPressureCoefficients();

public:
    
    TPZFastCondensedElement(TPZCompEl *ref, bool keepmatrix = true) :
        TPZCondensedCompEl(ref,keepmatrix)
    {
       
        
    }
    
    
    /// Assignement constructor
    const TPZFastCondensedElement & operator=(const TPZFastCondensedElement & other){
        fLambda = other.fLambda;
        return *this;
    }
    
    /** @brief create a copy of the condensed computational element in the other mesh */
    TPZFastCondensedElement(const TPZFastCondensedElement &copy, TPZCompMesh &mesh) :
        TPZCondensedCompEl(copy, mesh)
    {
        fEK = copy.fEK;
        fEF = copy.fEF;
        fLambda = copy.fLambda;
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
    
    /** @brief Loads the solution within the internal data structure of the element */
    /**
     * Is used to initialize the solution of connect objects with dependency \n
     * Is also used to load the solution within SuperElements
     */
    virtual void LoadSolution() override;
    

    void SetLambda(REAL lambda);
    REAL GetLambda();
    
    void SetMixedDensity(REAL density);
    REAL GetMixedDensity();

    /**
     * @brief Calculates the solution - sol - for the variable var
     * at point qsi, where qsi is expressed in terms of the
     * master element coordinates
     * @param qsi master element coordinate
     * @param var variable name
     * @param sol vetor for the solution
     */
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
    
    void SetPermTensorAndInv(TPZFNMatrix<9, REAL> &PermeabilityTensor, TPZFNMatrix<9, REAL> &InvPerm);
    TPZFMatrix<REAL> &GetPermTensor();
    
};
#endif /* TPZFASTCONDENSEDELEMENT_h */
