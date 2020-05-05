//
//  TMRSTransportAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSTransportAnalysis.h"
#include "pzfunction.h"
#include "TPZTracerFlow.h"

TMRSTransportAnalysis::TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::~TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::TMRSTransportAnalysis(TPZMultiphysicsCompMesh * cmesh_mult, bool must_opt_band_width_Q) : TPZAnalysis(cmesh_mult, must_opt_band_width_Q){
    
    m_soltransportTransfer.BuildTransferData(cmesh_mult);
    
}

void TMRSTransportAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
   
}

TMRSDataTransfer * TMRSTransportAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSTransportAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSTransportAnalysis::Configure(int n_threads, bool UsePardiso_Q){
    
    if (UsePardiso_Q) {
        TPZSpStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        SetSolver(step);
        
    }else{
        TPZSkylineNSymStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
}

void TMRSTransportAnalysis::RunTimeStep(){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    
    int n = m_sim_data->mTNumerics.m_max_iter_transport;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_transport;
    AssembleResidual();
    
    res_norm = Norm(Rhs());
    
    
    if (res_norm < res_tol) {
        std::cout << "Already converged solution with res_norm = " << res_norm << std::endl;
        return;
    }
    
    TPZFMatrix<STATE> dx,x(Solution());
    
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        NewtonIteration();
        dx = Solution();
        corr_norm = Norm(dx);
        cmesh->UpdatePreviousState(1);
        Rhs() *=-1.0;
//        m_soltransportTransfer.TransferFromMultiphysics();
        cmesh->LoadSolutionFromMultiPhysics();
        AssembleResidual();
        res_norm = Norm(Rhs());
        
        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_Q && stop_criterion_corr_Q) {
//        if (stop_criterion_Q ) {
            std::cout << "Transport operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            break;
        }

    }
    
}
void TMRSTransportAnalysis::RunTimeStepWithoutMemory(TPZFMatrix<REAL> s_n){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    
    int n = m_sim_data->mTNumerics.m_max_iter_transport;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_transport;
    AssembleResidual();
    
    res_norm = Norm(Rhs());
    
    
    if (res_norm < res_tol) {
        std::cout << "Already converged solution with res_norm = " << res_norm << std::endl;
        return;
    }
    
    TPZFMatrix<STATE> dx,x(Solution());
    
    /// Compute mass matrix M.
    TPZAutoPointer<TPZMatrix<STATE> > M;
    TPZFMatrix<REAL> F_inlet;
    TPZFMatrix<STATE>  M_diag;
    {
        
        bool mass_matrix_Q = true;
        std::set<int> volumetric_mat_ids = {1,2};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing Mass Matrix." << std::endl;
        Assemble();
        std::cout << "Mass Matrix is computed." << std::endl;
        M = this->Solver().Matrix()->Clone();
    }
    int n_rows = M->Rows();
    M_diag.Resize(n_rows,1);
    for (int64_t i = 0; i < n_rows; i++) {
        M_diag(i,0) = M->Get(i, i);
    }
    
    {
        bool mass_matrix_Q = false;
        std::set<int> volumetric_mat_ids = {1,2};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetTimeStep(100);
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing transport operator K = M + T, and F_inlet " << std::endl;
        this->Assemble();
        F_inlet = this->Rhs();
    }
    
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        int64_t n_eq = this->Mesh()->NEquations();
        
        
        {
//            TPZFMatrix<REAL> s_n(n_eq,1,0.0);
//            if (m_k_iteration!=1) {
//                TPZFMatrix<REAL> s_n=cmesh->SolutionN();
//            }
            
            TPZFMatrix<REAL> last_state_mass(n_eq,1,0.0);
            TPZFMatrix<REAL> s_np1;
            
        
                for (int64_t i = 0; i < n_eq; i++) {
                    last_state_mass(i,0) = M_diag(i,0)*s_n(i,0);
                }
                
                this->Rhs() = F_inlet - last_state_mass;
                this->Rhs() *= -1.0;
                
                this->Solve(); /// (LU decomposition)
                s_np1 = this->Solution();
                this->LoadSolution(s_np1);
                
            
        }
        
        
        
//        //
//        NewtonIteration();
        dx = Solution();
        corr_norm = Norm(dx);
//        cmesh->UpdatePreviousState(-1);
//        Rhs() *=-1.0;
//        //        m_soltransportTransfer.TransferFromMultiphysics();
//        cmesh->LoadSolutionFromMultiPhysics();
        AssembleResidual();
        res_norm = Norm(Rhs());
        
        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
//        if (stop_criterion_Q && stop_criterion_corr_Q) {
        if (stop_criterion_Q ) {
            std::cout << "Transport operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            break;
        }
        
    }
    
}

void TMRSTransportAnalysis::NewtonIteration(){
    
    Assemble();
//    Rhs() *= -1.0;
    Solve();
}

void TMRSTransportAnalysis::PostProcessTimeStep(){
    
    TPZStack<std::string,10> scalnames, vecnames;
    // @TODO:: Locate these variables in mTPostProcess
    scalnames.Push("Sw");
    scalnames.Push("So");

    int div = 0;
    int dim = Mesh()->Reference()->Dimension();
    std::string file = m_sim_data->mTPostProcess.m_file_name_transport;
    
    DefineGraphMesh(dim,scalnames,vecnames,file);
    PostProcess(div,dim);
    
}


