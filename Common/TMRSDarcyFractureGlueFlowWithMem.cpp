//
//  TMRSDarcyFractureGlueFlowWithMem.cpp
//
//  Created by José Villegas on 30/05/22.

#include "TMRSMemory.h"

#include "TMRSDarcyFractureGlueFlowWithMem.h"

template <class TMEM>
TMRSDarcyFractureGlueFlowWithMem<TMEM>::TMRSDarcyFractureGlueFlowWithMem() : TBase() {
}

template <class TMEM>
TMRSDarcyFractureGlueFlowWithMem<TMEM>::TMRSDarcyFractureGlueFlowWithMem(int mat_id, int dimension) : TBase(mat_id, dimension) {
    
    
}

template <class TMEM>
TMRSDarcyFractureGlueFlowWithMem<TMEM>::TMRSDarcyFractureGlueFlowWithMem(const TMRSDarcyFractureGlueFlowWithMem & other) : TMRSDarcyFlowWithMem<TMEM>(other){
}

template <class TMEM>
TMRSDarcyFractureGlueFlowWithMem<TMEM> & TMRSDarcyFractureGlueFlowWithMem<TMEM>::operator=(const TMRSDarcyFractureGlueFlowWithMem & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    return *this;
}

template <class TMEM>
TMRSDarcyFractureGlueFlowWithMem<TMEM>::~TMRSDarcyFractureGlueFlowWithMem(){
    
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fDeformedDirections = true;
    }
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::SetDataTransfer(TMRSDataTransfer & SimData){
    this->mSimData = SimData;
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::Print(std::ostream &out) const{
    TPZMaterial::Print(out);
}

template <class TMEM>
int TMRSDarcyFractureGlueFlowWithMem<TMEM>::VariableIndex(const std::string &name) const{
    if(!strcmp("Flux",name.c_str()))            return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    if(!strcmp("div_q",name.c_str()))           return  3;
    if(!strcmp("kappa",name.c_str()))           return  4;
    if(!strcmp("g_average",name.c_str()))        return  5;
    if(!strcmp("p_average",name.c_str()))        return  6;
    return TPZMaterial::VariableIndex(name);
}

template <class TMEM>
int TMRSDarcyFractureGlueFlowWithMem<TMEM>::NSolutionVariables(int var) const{
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    return TBase::NSolutionVariables(var);
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    TMRSDarcyFlowWithMem<TMEM>::Solution(datavec,var,Solout);
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int qb = 0;
    int pb = 1;
 
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
    // Get the data at integrations points
    long gp_index = datavec[qb].intGlobPtIndex;
    if(gp_index < 0) DebugStop();
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    REAL kappaNormal = memory.m_kappa_normal;
    std::pair<int, int> fractureindexes = memory.m_fracindexes;
    
    
    
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    

    
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){

   
    
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
   
}

template <class TMEM>
void TMRSDarcyFractureGlueFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    
   
    
}

template class TMRSDarcyFractureGlueFlowWithMem<TMRSMemory>;
