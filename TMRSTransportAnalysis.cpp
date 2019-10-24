//
//  TMRSTransportAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSTransportAnalysis.h"
#include "pzfunction.h"

TMRSTransportAnalysis::TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::~TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::TMRSTransportAnalysis(TPZMultiphysicsCompMesh * cmesh_mult, bool must_opt_band_width_Q) : TPZAnalysis(cmesh_mult, must_opt_band_width_Q){
    
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
        
        x+=dx;
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        
        AssembleResidual();
        res_norm = Norm(Rhs());
        
        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_Q && stop_criterion_corr_Q) {
            std::cout << "Transport operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
//            x.Print("x = ",std::cout,EMathematicaInput);
//            Rhs().Print("r = ",std::cout,EMathematicaInput);
            break;
        }
        if (m_k_iteration >= n) {
            std::cout << "Transport operator not converge " << std::endl;
        }
    }
    
}


void TMRSTransportAnalysis::NewtonIteration(){
    Assemble();
    Rhs() *= -1.0;
    Solve();
}

void TMRSTransportAnalysis::PostProcessTimeStep(){
    
    TPZStack<std::string,10> scalnames, vecnames;
    // @TODO:: Locate these variables in mTPostProcess
    scalnames.Push("Sw");
    scalnames.Push("So");
    scalnames.Push("Sw_exact");

    int div = 0;
    int dim = Mesh()->Reference()->Dimension();
    std::string file = m_sim_data->mTPostProcess.m_file_name_transport;
    
    DefineGraphMesh(dim,scalnames,vecnames,file);
    PostProcess(div,dim);
    
}


