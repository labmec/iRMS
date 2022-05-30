//
//  RSDarcyGlueFractureFlowWithMem.h
//
//  Created by José Villegas on 30/05/22.
//

#ifndef TMRSDarcyFractureGlueFlowWithMem_h
#define TMRSDarcyFractureGlueFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"
#include "TMRSDarcyFlowWithMem.h"
#include "TPZMaterialDataT.h"

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "TPZMatWithMem.h"
template <class TMEM>
class TMRSDarcyFractureGlueFlowWithMem : public TMRSDarcyFlowWithMem<TMEM> {
    using TBase =TMRSDarcyFlowWithMem<TMEM>;
    
public:
    
    /// Default constructor
    TMRSDarcyFractureGlueFlowWithMem();
    
    /// Constructor based on a material id
    TMRSDarcyFractureGlueFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureGlueFlowWithMem(const TMRSDarcyFractureGlueFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureGlueFlowWithMem &operator=(const TMRSDarcyFractureGlueFlowWithMem & other);
    
    /// Default destructor
    ~TMRSDarcyFractureGlueFlowWithMem();
    
    /// Set the required data at each integration point
    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec)const override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Returns the name of the material
    std::string Name() const override {
        return "TMRSDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return this->m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TMRSDarcyFractureFlowWithMem<TMEM>(*this);
    }
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) const override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) const override;
    
    /// returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var) const override;
    
  
    
    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    
    
    
    // Contribute Methods being used
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    void ContributeFourSpaces( const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)override;
    
};

#endif /* TMRSDarcyFractureGlueFlowWithMem_h */
