//
//  TMRSMultiphaseFlow.h
//
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSMultiphaseFlow_h
#define TMRSMultiphaseFlow_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzdiscgal.h"
#include "TMRSDataTransfer.h"
#include <tuple>
#include <functional>

template <class TMEM>
class TMRSMultiphaseFlow : public TPZMatWithMem<TMEM,TPZDiscontinuousGalerkin> {
  
private:
    
    /// Dimension
    int m_dimension;
    
    TMRSDataTransfer mSimData;
    
public:
    
    /// Default constructor
    TMRSMultiphaseFlow();
    
    /// Constructor based on a material id
    TMRSMultiphaseFlow(int matid, int dimension);
    
    /// Constructor based on a TMRSMultiphaseFlow object
    TMRSMultiphaseFlow(const TMRSMultiphaseFlow &other);
    
    /// Assignment operator
    TMRSMultiphaseFlow &operator=(const TMRSMultiphaseFlow &other);
    
    /// Default destructor
    ~TMRSMultiphaseFlow();
    
    /// Set the required data at each integration point
    virtual void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override ;
    
    /// Set the required data at each integration point
    virtual void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data) override;
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right) ;
    
    /// Returns the name of the material
    virtual std::string Name() override{
        return "TMRSMultiphaseFlow";
    }
    
    /// Returns the integrable dimension of the material
    int Dimension() const override {return m_dimension;}
    
    /// Sets material dimension
    void SetDimension(int dim) { m_dimension = dim; }
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;} // Deprecated, must to be removed
    
    /// Returns material copied form this object
    virtual TPZMaterial *NewMaterial() override
    {
        return new TMRSMultiphaseFlow(*this);
    }
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) override;
    
    /// Returns the number of variables associated with varindex
    int NSolutionVariables(int var) override;
    
    /// Not used contribute methods
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)  override { DebugStop();}
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)  override{DebugStop();}
    void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)  {DebugStop();}
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) {DebugStop();}
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)  {DebugStop();}
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {DebugStop();}
    
    
    // Contribute Methods being used
    
    /// Returns the solution associated with the var index
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) ;
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ef) ;
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) ;
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /// Unique identifier for serialization purposes
    int ClassId() const override;
    
     /// Save the element data to a stream
    void Write(TPZStream &buf, int withclassid);
    
    /// Read the element data from a stream
    void Read(TPZStream &buf, void *context) override;
    
    void BuckleyLeverett(TPZVec<REAL> pt,TPZVec<REAL> q, REAL phi,REAL time, TPZVec<REAL> &Saturation, TMRSDataTransfer sim_data);
    
};

#endif /* TMRSMultiphaseFlow_h */
