//
//  TPZDarcyFlowWithMem.h
//
//  Created by José Villegas on 07/07/20.
//

#ifndef TPZDarcyFlowWithMem_h
#define TPZDarcyFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"
#include "TPZAlgebraicTransport.h"
#include "TPZDarcyMemory.h"


class TPZDarcyFlowWithMem : public TPZMatWithMem<TPZDarcyMemory> {
    
    std::vector<REAL> m_gravity;
    /// Dimension
    int m_dimension;
    
    /// Scale factor for pressure
    STATE m_scale_pressure = 1.0;//1.e-6;

    /// Scale factor for flux variable
    STATE m_scale_flux = 1.0;
    
    /// Directive that stands for the use of four approximations spaces (iterative method)
    bool m_is_four_spaces_Q;
    
    TPZAlgebraicTransport *fAlgebraicTransport;
    
public:
    TMRSDataTransfer mSimData;
    /// Default constructor
    TPZDarcyFlowWithMem();
    
    /// Constructor based on a material id
    TPZDarcyFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TPZDarcyFlowWithMem(const TPZDarcyFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TPZDarcyFlowWithMem &operator=(const TPZDarcyFlowWithMem & other);
    
    /// Default destructor
    ~TPZDarcyFlowWithMem();
    
    void SetAlgebraicTransport(TPZAlgebraicTransport *algebraicTransport){
        fAlgebraicTransport=algebraicTransport;
    }
    TPZAlgebraicTransport *GetAlgebraicTransport(){
        return fAlgebraicTransport;
    }
    
    void SetGravity(std::vector<REAL> gravity){
        m_gravity =gravity;
    }
    std::vector<REAL> GetGravity( ){
        return m_gravity;
    }
    /// Set the required data at each integration point
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;
    
    /// Returns the name of the material
    std::string Name() override {
        return "TPZDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() override
    {
        return new TPZDarcyFlowWithMem(*this);
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

#endif /* TPZDarcyFlowWithMem_h */
