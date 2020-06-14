//
//  TPZMixedDarcyWithFourSpaces.cpp
//  reservoirlib
//
//  Created by Omar Durán on 7/12/19.
//

#include "TPZMixedDarcyWithFourSpaces.h"


TPZMixedDarcyWithFourSpaces::TPZMixedDarcyWithFourSpaces() : TPZMixedDarcyFlow() {
    
}

TPZMixedDarcyWithFourSpaces::~TPZMixedDarcyWithFourSpaces(){
    
}

TPZMixedDarcyWithFourSpaces::TPZMixedDarcyWithFourSpaces(int mat_id, int dim) : TPZMixedDarcyFlow(mat_id,dim){
    
}

TPZMixedDarcyWithFourSpaces::TPZMixedDarcyWithFourSpaces(const TPZMixedDarcyWithFourSpaces &other) : TPZMixedDarcyFlow(other){
    
}

void TPZMixedDarcyWithFourSpaces::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    TPZMixedDarcyFlow::Contribute(datavec, weight, ek, ef);
    
    int qb = 0;
    int pb = 1;
    int g_avgb = 2;
    int p_avgb = 3;
    
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    
    int nphi_gb = datavec[g_avgb].phi.Rows();
    int nphi_ub = datavec[p_avgb].phi.Rows();
    if(nphi_q+nphi_p+nphi_gb+nphi_ub != ek.Rows())
    {
        DebugStop();
    }
    
    TPZManVector<STATE> force(1);
    if(fForcingFunction) {
        fForcingFunction->Execute(datavec[qb].x,force);
    }
    
    for(int ip=0; ip<nphi_p; ip++)
    {
        ef(nphi_q+ip,0) += weight * force[0]*phi_ps(ip,0);
        ek(nphi_q+ip,nphi_q+nphi_p) += phi_ps(ip,0)*weight;
        ek(nphi_q+nphi_p,nphi_q+ip) += phi_ps(ip,0)*weight;
    }
    ek(nphi_q+nphi_p+1,nphi_q+nphi_p) += -weight;
    ek(nphi_q+nphi_p,nphi_q+nphi_p+1) += -weight;
    
}

void TPZMixedDarcyWithFourSpaces::Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ekfake, ef);
}


int TPZMixedDarcyWithFourSpaces::VariableIndex(const std::string &name) {
    if(!strcmp("g_average",name.c_str()))        return  5;
    if(!strcmp("u_average",name.c_str()))        return  6;
    return TPZMixedDarcyFlow::VariableIndex(name);
}

int TPZMixedDarcyWithFourSpaces::NSolutionVariables(int var) {
    if(var == 5) return 1;
    if(var == 6) return 1;
    return TPZMixedDarcyFlow::NSolutionVariables(var);
}

void TPZMixedDarcyWithFourSpaces::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){

    int g_avgb = 2;
    int p_avgb = 3;
    
    if(var ==5)
    {
        Solout[0] = datavec[g_avgb].sol[0][0];
        return;
    }
    if(var ==6)
    {
        Solout[0] = datavec[p_avgb].sol[0][0];
        return;
    }
    
    TPZMixedDarcyFlow::Solution(datavec, var, Solout);
}
