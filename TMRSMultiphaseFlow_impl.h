//
//  TMRSMultiphaseFlow_impl.h
//
//
//  Created by Omar Durán on 10/10/19.
//

#include "TMRSMultiphaseFlow.h"

template <class TMEM>
TMRSMultiphaseFlow<TMEM>::TMRSMultiphaseFlow() : TPZMatWithMem<TMEM,TPZDiscontinuousGalerkin>(), mSimData(){
    m_dimension = 0;
}

template <class TMEM>
TMRSMultiphaseFlow<TMEM>::TMRSMultiphaseFlow(int matid, int dimension) : TPZMatWithMem<TMEM,TPZDiscontinuousGalerkin>(matid), mSimData(){
    m_dimension = dimension;
}

template <class TMEM>
TMRSMultiphaseFlow<TMEM>::TMRSMultiphaseFlow(const TMRSMultiphaseFlow &other){
    m_dimension = other.m_dimension;
    mSimData = other.mSimData;
}

template <class TMEM>
TMRSMultiphaseFlow<TMEM> & TMRSMultiphaseFlow<TMEM>::operator=(const TMRSMultiphaseFlow &other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    m_dimension = other.m_dimension;
    mSimData = other.mSimData;
    return *this;
}

template <class TMEM>
TMRSMultiphaseFlow<TMEM>::~TMRSMultiphaseFlow(){

}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::FillDataRequirements(TPZVec<TPZMaterialData> &datavec) {
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) {
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::FillDataRequirementsInterface(TPZMaterialData &data) {
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    if(TPZMatWithMem<TMEM,TPZDiscontinuousGalerkin>::fLinearContext == false){
        data.fNeedsNeighborSol = true;
    }
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right) {
   
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
    data.fNeedsNormal = true;
    int nref_left = datavec_left.size();
    for(int iref = 0; iref<nref_left; iref++){
        datavec_left[iref].SetAllRequirements(false);
        datavec_left[iref].fNeedsSol = true;
        datavec_left[iref].fNeedsNormal = true;
    }
    int nref_right = datavec_right.size();
    for(int iref = 0; iref<nref_right; iref++){
        datavec_right[iref].SetAllRequirements(false);
        datavec_right[iref].fNeedsSol = true;
    }
    
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::Print(std::ostream &out){
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::SetDataTransfer(TMRSDataTransfer & SimData){
    mSimData = SimData;
}

template <class TMEM>
int TMRSMultiphaseFlow<TMEM>::VariableIndex(const std::string &name) {
    if (!strcmp("Sw", name.c_str())) return 0;
    if (!strcmp("So", name.c_str())) return 1;
    if (!strcmp("Sw_exact", name.c_str())) return 2;
    
    return TPZMatWithMem<TMEM,TPZDiscontinuousGalerkin>::VariableIndex(name);
}

template <class TMEM>
int TMRSMultiphaseFlow<TMEM>::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Scalar
        case 2:
            return 1; // Scalar
            
    }
    return TPZMatWithMem<TMEM,TPZDiscontinuousGalerkin>::NSolutionVariables(var);
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int s_b    = 2;
    REAL sw = datavec[s_b].sol[0][0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = sw;
        }
            break;
        case 1:
        {
            Solout[0] = 1.0-sw;
        }
            break;
        case 2:
        {
            TPZVec<STATE> sw(1);
            BuckleyLeveret(datavec[0].x, sw);
            Solout[0] = sw[0];
        }
            break;
    }
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 3 ) {
        std::cout << " Erro. The size of the datavec is different from 3 \n";
        DebugStop();
    }
#endif
    
    int s_b = 2;
    long gp_index = datavec[s_b].intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    TPZFMatrix<REAL>  &phiS =  datavec[s_b].phi;

    
    REAL phi = memory.m_phi;
    REAL sw = memory.m_sw;
    REAL sw_n = datavec[s_b].sol[0][0];


    int n_phi_s = phiS.Rows();
    int firsts_s    = 0;
    
    for (int is = 0; is < n_phi_s; is++)
    {
        ef(is + firsts_s) += 1.0 * weight * phi * (sw_n - sw)  * phiS(is,0);
        
        for (int js = 0; js < n_phi_s; js++)
        {
            ek(is + firsts_s, js + firsts_s) += weight * phi * (phiS(js,0) )* phiS(is,0);
        }
    }
    
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->Contribute(datavec, weight, ek_fake, ef);
    
    
    if(TMRSMultiphaseFlow<TMEM>::fUpdateMem){
        int s_b = 2;
        long gp_index = datavec[s_b].intGlobPtIndex;
        TMEM & memory = this->GetMemory().get()->operator[](gp_index);
        REAL sw_n = datavec[s_b].sol[0][0];
        memory.m_sw = sw_n;
        memory.m_so = 1.0-sw_n;
    }
    
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    int q_b = 0;
//    int p_b = 1;
    int s_b = 2;
    REAL dt = mSimData.mTNumerics.m_dt;
    
    TRSLinearInterpolator & Krw = mSimData.mTPetroPhysics.mLayer_Krw_RelPerModel[0];
    TRSLinearInterpolator & Kro = mSimData.mTPetroPhysics.mLayer_Kro_RelPerModel[0];
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & fw = mSimData.mTMultiphaseFunctions.mLayer_fw[0];
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & fo = mSimData.mTMultiphaseFunctions.mLayer_fo[0];
    

    
    // Getting phis and solution for left material data
    TPZFMatrix<REAL>  &phiS_l =  datavecleft[s_b].phi;
    int n_phi_s_l = phiS_l.Rows();
    REAL p_l = 0.0;//datavecleft[p_b].sol[0][0];
    REAL s_l = datavecleft[s_b].sol[0][0];
    int firsts_s_l    = 0;
    
    // Getting phis and solution for right material data
    TPZFMatrix<REAL>  &phiS_r =  datavecright[s_b].phi;
    int n_phi_s_r = phiS_r.Rows();
    REAL p_r = 0.0;//datavecright[p_b].sol[0][0];
    REAL s_r = datavecright[s_b].sol[0][0];
    int firsts_s_r    = n_phi_s_l;
    
    TPZManVector<REAL,3> n = data.normal;
    TPZManVector<REAL,3> q_l =  datavecleft[q_b].sol[0];
    REAL qn = 0.0;
    for (int i = 0; i < 3; i++) {
        qn += q_l[i]*n[i];
    }
    
    REAL beta = 0.0;
    // upwinding
    if (qn > 0.0) {
        beta = 1.0;
    }

    std::tuple<double, double, double> fw_l = fw(Krw,Kro,s_l,p_l);
    std::tuple<double, double, double> fw_r = fw(Krw,Kro,s_r,p_r);

    REAL fw_lv = std::get<0>(fw_l);
    REAL dfw_dsw_lv = std::get<1>(fw_l);
    
    REAL fw_rv = std::get<0>(fw_r);
    REAL dfw_dsw_rv = std::get<1>(fw_r);
    
    
    for (int is = 0; is < n_phi_s_l; is++) {
        
        ef(is + firsts_s_l) += +1.0 * dt * weight * (beta*fw_lv + (1.0-beta)*fw_rv)*phiS_l(is,0)*qn;
        
        for (int js = 0; js < n_phi_s_l; js++) {
            ek(is + firsts_s_l, js + firsts_s_l) += +1.0* dt * weight * beta * dfw_dsw_lv * phiS_l(js,0) * phiS_l(is,0)*qn;
        }
        
        for (int js = 0; js < n_phi_s_r; js++) {
            ek(is + firsts_s_l, js + firsts_s_r) += +1.0* dt * weight * (1.0-beta) * dfw_dsw_rv * phiS_r(js,0) * phiS_l(is,0)*qn;
        }
        
    }
    
    for (int is = 0; is < n_phi_s_r; is++) {
        
        ef(is + firsts_s_r) += -1.0* dt * weight * (beta*fw_lv + (1.0-beta)*fw_rv)*phiS_r(is,0)*qn;
        
        for (int js = 0; js < n_phi_s_l; js++) {
            ek(is + firsts_s_r, js + firsts_s_l) += -1.0* dt * weight * beta * dfw_dsw_lv * phiS_l(js,0) * phiS_r(is,0)*qn;
        }
        
        for (int js = 0; js < n_phi_s_r; js++) {
            ek(is + firsts_s_r, js + firsts_s_r) += -1.0* dt * weight * (1.0-beta) * dfw_dsw_rv * phiS_r(js,0) * phiS_r(is,0)*qn;
        }
        
    }
    
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ef) {
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->ContributeInterface(data,datavecleft,datavecright, weight, ek_fake, ef);
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
    
    REAL tol = 1.0e-10;
 //   REAL tol = 0.01;
    int q_b = 0;
    int s_b = 2;
    
    
    
    REAL dt = mSimData.mTNumerics.m_dt;
    
    // Getting phis and solution for left material data
    TPZFMatrix<REAL>  &phiS_l =  datavecleft[s_b].phi;
    int n_phi_s_l = phiS_l.Rows();
    REAL s_l = datavecleft[s_b].sol[0][0];
    int firsts_s_l    = 0;
    
    
    TPZManVector<REAL,3> n = data.normal;
    TPZManVector<REAL,3> q_l =  datavecleft[q_b].sol[0];
    REAL qn = 0.0;

    
    for (int i = 0; i < 3; i++) {
        qn += q_l[i]*n[i];
    }

    //
    
    TRSLinearInterpolator & Krw = mSimData.mTPetroPhysics.mLayer_Krw_RelPerModel[0];
    TRSLinearInterpolator & Kro = mSimData.mTPetroPhysics.mLayer_Kro_RelPerModel[0];
    
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & fw = mSimData.mTMultiphaseFunctions.mLayer_fw[0];
    std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> & fo = mSimData.mTMultiphaseFunctions.mLayer_fo[0];
    
    std::tuple<double, double, double> fw_l = fw(Krw,Kro,s_l,0.0);
    

    REAL fw_lv = std::get<0>(fw_l);
    REAL dfw_dsw_lv = std::get<1>(fw_l);
    
    
    //
    switch (bc.Type()) {
            
        case 0 :    // BC inlet
        {
            
            REAL s_inlet = bc.Val2()(0,0);
            if (qn < 0.0 || fabs(qn) < tol) {
                for (int is = 0; is < n_phi_s_l; is++) {
                    ef(is + firsts_s_l) += +1.0* dt * weight * s_inlet * phiS_l(is,0)*qn;
                }
            }else{

                std::cout << "TPZTracerFlow:: Outlet flux in inlet boundary condition qn = " << qn << std::endl;
            }
 
            
        }
            break;
            
        case 1 :    // BC outlet
        {
            
            if (qn > 0.0 || fabs(qn) < tol) {
                for (int is = 0; is < n_phi_s_l; is++) {
                    
                    ef(is + firsts_s_l) += 1.0* dt * weight * fw_lv *phiS_l(is,0)*qn;
                    
                    for (int js = 0; js < n_phi_s_l; js++) {
                        ek(is + firsts_s_l, js + firsts_s_l) +=  1.0* dt * weight * phiS_l(js,0) * phiS_l(is,0)*qn*dfw_dsw_lv;
                    }
                }
            }else{
                std::cout << "TPZTracerFlow:: Inlet flux on outlet boundary condition qn = " << qn << std::endl;
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

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    TPZFMatrix<STATE> ek_fake(ef.Rows(),ef.Rows());
    this->ContributeBCInterface(data,datavecleft, weight, ek_fake, ef, bc);
    
}

template <class TMEM>
int TMRSMultiphaseFlow<TMEM>::ClassId() const{
    DebugStop();
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::Write(TPZStream &buf, int withclassid){
    DebugStop();
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::Read(TPZStream &buf, void *context) {
    DebugStop();
}

template <class TMEM>
void TMRSMultiphaseFlow<TMEM>::BuckleyLeveret(TPZVec<REAL> &pt, TPZVec< STATE > &Saturation ) {
    REAL Porosity = 0.2;
    REAL x = pt[0];
    
    REAL x_shock        = 0.0;
    REAL S_shock        = 0.0;
    REAL epsilon        = 1.0e-8;
    REAL u              = 0.1;
    REAL mu_alpha       = 1.0;
    REAL mu_beta        = 1.0;
    REAL rho_alpha      = 1.0;
    REAL rho_beta       = 1.0;
    REAL Sor            = 0.0;
    REAL Swr            = 0.0;
    double time           = 1.5;
    
    S_shock = Swr + sqrt(mu_alpha*rho_beta*(mu_beta*rho_alpha + mu_alpha*rho_beta)*std::pow(-1.0 + Sor + Swr,2.0))/(mu_beta*rho_alpha + mu_alpha*rho_beta);
    x_shock  =  (u*time)/(Porosity*rho_alpha)*dfdsw(S_shock, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
    
    if(x < x_shock)
    {
        REAL Sw = S_Newton(x, time, u, Swr, Sor, Porosity, S_shock, mu_alpha, mu_beta, rho_alpha, rho_beta, epsilon);
        Saturation[0] = Sw;
        
    }
    else
    {
        Saturation[0] = Swr;
    }
    
    return;
}

 template <class TMEM>
 REAL TMRSMultiphaseFlow<TMEM>::S_Newton(REAL x, REAL t, REAL u, REAL Swr, REAL Sor, REAL phi, REAL s_shok, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta, REAL epsilon)
{
    REAL S_trial = ((1.0-Sor) + s_shok)/2.0;
    REAL jac    = 0.0;
    REAL r      = 1.0;
    REAL delta_S;
    REAL ds;
    REAL S_k = S_trial;
    int max = 20;
    int it = 0;
    
    
    
    while ( fabs(r) > epsilon && it < max){
        
        ds = dfdsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
        r = x - (u * t * ds)/(phi * rho_alpha);
        
        jac = - (u * t * df2dsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta))/(phi * rho_alpha);
        
        delta_S = -r/jac;
        S_k += delta_S;
        
        ds = dfdsw(S_k, Swr, Sor, mu_alpha, mu_beta, rho_alpha, rho_beta);
        r = x - (u * t * ds)/(phi * rho_alpha);
        it++;
        
    }
    
    return S_k;
}

template <class TMEM>
REAL TMRSMultiphaseFlow<TMEM>::dfdsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta)
{
    REAL dfwdSwS;
    
    
    dfwdSwS = (2.0*mu_alpha*mu_beta*rho_alpha*rho_beta*(-1.0 + Sor + Sw)*(Sw - Swr)*(-1.0 + Sor + Swr))/std::pow(mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0) + mu_beta*rho_alpha*std::pow(Sw - Swr,2.0),2.0);
    
    return (dfwdSwS);
}
template <class TMEM>
REAL TMRSMultiphaseFlow<TMEM>::df2dsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta)
{
    REAL dfw2dSwS2;
    
    
    dfw2dSwS2 = (2.0*mu_alpha*mu_beta*rho_alpha*rho_beta*(-1.0 + Sor + Swr)*(-(mu_beta*rho_alpha*std::pow(Sw - Swr,2.0)*(-3.0 + 3.0*Sor + 2.0*Sw + Swr)) + mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0)*(-1.0 + Sor - 2.0*Sw + 3.0*Swr)))/
    std::pow(mu_alpha*rho_beta*std::pow(-1.0 + Sor + Sw,2.0) + mu_beta*rho_alpha*std::pow(Sw - Swr,2.0),3.0);
    
    return dfw2dSwS2;
}
