//
//  TMRSMixedAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.

#include "TMRSMixedAnalysis.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif
#include "TPZSpStructMatrix_Eigen.h"
#include "TPZSpMatrixEigen.h"

TMRSMixedAnalysis::TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::~TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::TMRSMixedAnalysis(TPZMultiphysicsCompMesh * cmesh_mult, bool must_opt_band_width_Q) : TPZAnalysis(cmesh_mult, must_opt_band_width_Q){
    fsoltransfer.BuildTransferData(cmesh_mult);
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
//        TPZSpStructMatrixEigen matrix(Mesh());
//        matrix.SetNumThreads(n_threads);
//        SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELU);
//        SetSolver(step);

    }else{
        TPZSkylineStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
    
   
    //    Assemble();
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
    
    TPZFMatrix<STATE> dx,x(Solution());
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        NewtonIteration();
        dx = Solution();
        corr_norm = Norm(dx);
        x +=dx;
        cmesh->UpdatePreviousState(-1);
        cmesh->Solution().Print(std::cout);
        fsoltransfer.TransferFromMultiphysics();
        Assemble();
        res_norm = Norm(Rhs());

#ifdef PZDEBUG
        {
 
            if(std::isnan(corr_norm) || std::isnan(res_norm))
            {
                DebugStop();
            }
        }
#endif

        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_Q) {
            std::cout << "Mixed operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
  
            break;
        }
                if (m_k_iteration >= n) {
                    std::cout << "Mixed operator not converge " << std::endl;
                }
        
    }
}


void TMRSMixedAnalysis::NewtonIteration(){
    
#ifdef USING_BOOST
    boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
#endif

    if(mIsFirstAssembleQ == true)
    {
        fStructMatrix->SetNumThreads(m_sim_data->mTNumerics.m_nThreadsMixedProblem);
        int64_t nel = fCompMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCompMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
                sub->Analysis()->StructMatrix()->SetNumThreads(m_sim_data->mTNumerics.m_nThreadsMixedProblem);
            }
        }
        mIsFirstAssembleQ=false;
    }

    
   Assemble();

#ifdef USING_BOOST
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = tsim2-tsim1;
    std::cout << "Mixed:: Assembly time " << deltat << std::endl;
#endif
   
    Solve();
    
#ifdef USING_BOOST
    boost::posix_time::ptime tsim3 = boost::posix_time::microsec_clock::local_time();
    auto deltat2 = tsim3-tsim1;
    std::cout << "Mixed:: Solve time " << deltat2 << std::endl;
#endif
//    TPZMatrix<STATE>*mat = Solver().Matrix().operator->();
//    
//    TPZSpMatrixEigen<STATE> *mateig = dynamic_cast<TPZSpMatrixEigen<STATE> *>(mat);
//    mateig->fsparse_eigen*0.0;
}

void TMRSMixedAnalysis::PostProcessTimeStep(){
    TPZStack<std::string,10> scalnames, vecnames;
    
    scalnames = m_sim_data->mTPostProcess.m_scalnames;
    vecnames = m_sim_data->mTPostProcess.m_vecnames;
    
    int div = 0;
    int dim = Mesh()->Reference()->Dimension();
    std::string file = m_sim_data->mTPostProcess.m_file_name_mixed;
    DefineGraphMesh(dim,scalnames,vecnames,file);
    PostProcess(div,dim);
}
