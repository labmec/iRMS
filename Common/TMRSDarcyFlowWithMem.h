//
//  TMRSDarcyFlowWithMem.h
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSDarcyFlowWithMem_h
#define TMRSDarcyFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"

template <class TMEM>
class TMRSDarcyFlowWithMem : public TPZMatWithMem<TMEM> {
    
    /// Dimension
    int m_dimension;
    
    /// Scale factor for pressure
    STATE m_scale_pressure = 1.0;//1.e-6;

    /// Scale factor for flux variable
    STATE m_scale_flux = 1.0;
    
    /// Directive that stands for the use of four approximations spaces (iterative method)
    bool m_is_four_spaces_Q;
    
    TMRSDataTransfer mSimData;
    
public:
    
    /// Default constructor
    TMRSDarcyFlowWithMem();
    
    /// Constructor based on a material id
    TMRSDarcyFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFlowWithMem(const TMRSDarcyFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFlowWithMem &operator=(const TMRSDarcyFlowWithMem & other);
    
    /// Default destructor
    ~TMRSDarcyFlowWithMem();
    
    /// Set the required data at each integration point
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;
    
    /// Returns the name of the material
    std::string Name() override {
        return "TMRSDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() override
    {
        return new TMRSDarcyFlowWithMem(*this);
    }
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) override;
    
    /// returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var) override;
    
    /// Returns the solution associated with the var index based on a finite element approximation (Used for TPZPostProcAnalysis)
    void Solution(TPZMaterialData &datavec, int var, TPZVec<REAL> &Solout) override {
        DebugStop();
    }
    
    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    // Contribute Methods being used
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeFourSpaces(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
};

#endif /* TMRSDarcyFlowWithMem_h */
