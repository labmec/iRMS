//
//  TMRSMixedAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSMixedAnalysis.h"

TMRSMixedAnalysis::TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::~TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::TMRSMixedAnalysis(TPZMultiphysicsCompMesh * cmesh_mult, bool must_opt_band_width_Q) : TPZAnalysis(cmesh_mult, must_opt_band_width_Q){
    m_soltransportTransfer.BuildTransferData(cmesh_mult);
}

void TMRSMixedAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
}

TMRSDataTransfer * TMRSMixedAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSMixedAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSMixedAnalysis::Configure(int n_threads, bool UsePardiso_Q){
    
    if (UsePardiso_Q) {
        TPZSymetricSpStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        SetSolver(step);
        
    }else{
        TPZSkylineStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
    Assemble();
}

void TMRSMixedAnalysis::RunTimeStep(){
    

    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    
    int n = m_sim_data->mTNumerics.m_max_iter_mixed;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_mixed;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_mixed;

    AssembleResidual();
    res_norm = Norm(Rhs());
    if (res_norm < res_tol && corr_norm<corr_tol) {
        std::cout << "Already converged solution with res_norm = " << res_norm << std::endl;
        return;
    }
    
    TPZFMatrix<STATE> dx,x(Solution());
  
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        NewtonIteration();
        dx = Solution();
        corr_norm = Norm(dx);
        cmesh->UpdatePreviousState(-1);
        cmesh->LoadSolutionFromMultiPhysics();
        
        
        
      
//
//        m_soltransportTransfer.TransferFromMultiphysics();
        
        AssembleResidual();
        res_norm = Norm(Rhs());
        
        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_Q && stop_criterion_corr_Q) {
//        if (stop_criterion_Q) {
      
            std::cout << "Mixed operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
//            x.Print("x = ",std::cout,EMathematicaInput);
//            Rhs().Print("r = ",std::cout,EMathematicaInput);
            break;
        }
        if (m_k_iteration >= n) {
            std::cout << "Mixed operator not converge " << std::endl;
        }
        
    }
}


void TMRSMixedAnalysis::NewtonIteration(){
    
    Assemble();
//    Rhs() *= -1.0; 
    Solve();
}

void TMRSMixedAnalysis::PostProcessTimeStep(){
    TPZStack<std::string,10> scalnames, vecnames;
    // @TODO:: Locate these variables in mTPostProcess
    
    scalnames = m_sim_data->mTPostProcess.m_scalnames;
    vecnames = m_sim_data->mTPostProcess.m_vecnames;

    int div = 0;
    int dim = Mesh()->Reference()->Dimension();
    std::string file = m_sim_data->mTPostProcess.m_file_name_mixed;
    DefineGraphMesh(dim,scalnames,vecnames,file);
    PostProcess(div,dim);
}
