//
//  TMRSDarcyFlowWithMem_impl.h
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSDarcyFlowWithMem_impl_h
#define TMRSDarcyFlowWithMem_impl_h

#include "TMRSDarcyFlowWithMem.h"

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::TMRSDarcyFlowWithMem() : TPZMatWithMem<TMEM>(), mSimData() {
    m_dimension = 0;
    m_is_four_spaces_Q = false;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::TMRSDarcyFlowWithMem(int mat_id, int dimension) : TPZMatWithMem<TMEM>(mat_id), mSimData(){
    m_dimension = dimension;
    m_is_four_spaces_Q = false;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::TMRSDarcyFlowWithMem(const TMRSDarcyFlowWithMem & other) : TPZMatWithMem<TMEM>(other){
    m_dimension         = other.m_dimension;
    m_scale_pressure    = other.m_scale_pressure;
    m_scale_flux        = other.m_scale_flux;
    m_is_four_spaces_Q  = other.m_is_four_spaces_Q;
    mSimData  = other.mSimData;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM> & TMRSDarcyFlowWithMem<TMEM>::operator=(const TMRSDarcyFlowWithMem & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    m_dimension         = other.m_dimension;
    m_scale_pressure    = other.m_scale_pressure;
    m_scale_flux        = other.m_scale_flux;
    m_is_four_spaces_Q  = other.m_is_four_spaces_Q;
    mSimData  = other.mSimData;
    return *this;
}

template <class TMEM>
TMRSDarcyFlowWithMem<TMEM>::~TMRSDarcyFlowWithMem(){
    
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::FillDataRequirements(TPZVec<TPZMaterialData> &datavec) {
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) {
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fDeformedDirections = true;
    }
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::SetDataTransfer(TMRSDataTransfer & SimData){
    mSimData = SimData;
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Print(std::ostream &out) {
    out << m_dimension << std::endl;
    out << m_scale_pressure << std::endl;
    out << m_scale_flux << std::endl;
    out << m_is_four_spaces_Q << std::endl;
}

template <class TMEM>
int TMRSDarcyFlowWithMem<TMEM>::VariableIndex(const std::string &name) {
    if(!strcmp("Flux",name.c_str()))            return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    if(!strcmp("div_q",name.c_str()))           return  3;
    if(!strcmp("kappa",name.c_str()))           return  4;
    if(!strcmp("g_average",name.c_str()))        return  5;
    if(!strcmp("u_average",name.c_str()))        return  6;
    return TPZMatWithMem<TMEM>::VariableIndex(name);
}

template <class TMEM>
int TMRSDarcyFlowWithMem<TMEM>::NSolutionVariables(int var) {
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    return TPZMatWithMem<TMEM>::NSolutionVariables(var);
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int qb = 0;
    int pb = 1;
    Solout.Resize( this->NSolutionVariables(var));
    TPZManVector<STATE,3> p, q;
    
    q = datavec[qb].sol[0];
    p = datavec[pb].sol[0];
    REAL div_q = datavec[qb].divsol[0][0];
    
    if(var == 1){
        for (int i=0; i < 3; i++)
        {
            Solout[i] = q[i];
        }
        return;
    }
    
    if(var == 2){
        Solout[0] = p[0];
        return;
    }
    
    if(var == 3){
        Solout[0] = div_q;
        return;
    }
    
    if(var == 4){
     
            TPZManVector<double, 3> point;
        
            point = datavec[qb].XCenter;
        
            int val = rand() % 100;
            
            REAL kappa = 1.0;
            if (val<75) {
                kappa =  10000;
                
            }
            else{
                kappa =  10;
            }
        
            Solout[0] = kappa;
        
        return;
    }
    
    if (mSimData.mTNumerics.m_four_approx_spaces_Q) {
        
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
    }
    
    DebugStop();
    
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int qb = 0;
    int pb = 1;
    int sb = 2;
    
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    
    TPZFNMatrix<40, REAL> div_phi = datavec[qb].divphi;
    REAL div_q = datavec[qb].divsol[0][0];
    
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
    STATE sw                 = datavec[sb].sol[0][0];
    
    // Get the data at integrations points
    long gp_index = datavec[qb].intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    TPZFNMatrix<3,STATE> phi_q_i(3,1,0.0), kappa_inv_phi_q_j(3,1,0.0), kappa_inv_q(3,1,0.0);
    
    TRSLinearInterpolator & Krw = mSimData.mTPetroPhysics.mLayer_Krw_RelPerModel[0];
    TRSLinearInterpolator & Kro = mSimData.mTPetroPhysics.mLayer_Kro_RelPerModel[0];
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & lambda = mSimData.mTMultiphaseFunctions.mLayer_lambda[0];
    
    
    // Total mobility
    std::tuple<double, double, double> lambda_t = lambda(Krw,Kro,memory.sw_n(),p);
   REAL lambda_v = std::get<0>(lambda_t);
 
    
    int s_i, s_j;
    int v_i, v_j;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            kappa_inv_q(i,0) += memory.m_kappa_inv(i,j)*(1.0/lambda_v)*q[j];
        }
    }
    
    for (int iq = 0; iq < nphi_q; iq++)
    {
        
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        
        STATE kappa_inv_q_dot_phi_q_i = 0.0;
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }
        
        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
            
            kappa_inv_phi_q_j.Zero();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    kappa_inv_phi_q_j(i,0) += memory.m_kappa_inv(i,j) * phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
                }
            }
            
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
            
            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;
        }
        
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
        }
        
    }
    
    for (int ip = 0; ip < nphi_p; ip++)
    {
        
        ef(ip + first_p) += -1.0 * weight * (div_q) * phi_ps(ip,0);
        
        for (int jq = 0; jq < nphi_q; jq++)
        {
            ek(ip + first_p, jq + first_q) += -1.0 * weight * div_phi(jq,0) * phi_ps(ip,0);
        }
        
    }
    
    if(mSimData.mTNumerics.m_four_approx_spaces_Q){
        ContributeFourSpaces(datavec,weight,ek,ef);
    }
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::ContributeFourSpaces(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    int qb = 0;
    int pb = 1;
    int g_avgb = 2;
    int p_avgb = 3;
    
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    
    int nphi_gb = datavec[g_avgb].phi.Rows();
    int nphi_pb = datavec[p_avgb].phi.Rows();
    if(nphi_q+nphi_p+nphi_gb+nphi_pb != ek.Rows())
    {
        DebugStop();
    }
    
    STATE p     = datavec[pb].sol[0][0];
    STATE g_avg = datavec[g_avgb].sol[0][0];
    STATE p_avg = datavec[p_avgb].sol[0][0];
    
    for(int ip=0; ip<nphi_p; ip++)
    {
        ef(nphi_q+ip,0) += weight * g_avg * phi_ps(ip,0);
        ek(nphi_q+ip,nphi_q+nphi_p) += weight * phi_ps(ip,0);
        
        ek(nphi_q+nphi_p,nphi_q+ip) += weight * phi_ps(ip,0);
    }
    
    ef(nphi_q+nphi_p+1,0) += -weight * g_avg;
    ek(nphi_q+nphi_p+1,nphi_q+nphi_p) += -weight;
    
    ef(nphi_q+nphi_p,0) += weight * (p - p_avg);
    ek(nphi_q+nphi_p,nphi_q+nphi_p+1) += -weight;
    
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->Contribute(datavec, weight, ekfake, ef);
    
    if(TMRSDarcyFlowWithMem<TMEM>::fUpdateMem){
        int qb = 0;
        int pb = 1;
        long gp_index = datavec[pb].intGlobPtIndex;
        TMEM & memory = this->GetMemory().get()->operator[](gp_index);
        TPZVec<REAL> q_n = datavec[qb].sol[0][0];
        memory.m_flux = q_n;
        
        REAL p_n = datavec[pb].sol[0][0];
        memory.m_p = p_n;
    }
    
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}

template <class TMEM>
void TMRSDarcyFlowWithMem<TMEM>::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    REAL gBigNumber = TPZMaterial::gBigNumber;
    int qb = 0;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()(0,0);
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += weight * p_D * phi_qs(iq,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                REAL qn_N = bc_data[0];
                REAL qn = 0.0;
                qn = q[0];

                ef(iq + first_q) += weight * gBigNumber * (qn - qn_N) * phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    ek(iq + first_q,jq + first_q) += weight * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

#endif /* TMRSDarcyFlowWithMem_impl_h */
