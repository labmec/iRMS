//
//  TMRSTransportAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSTransportAnalysis.h"
#include "pzfunction.h"
#include "TPZTracerFlow.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif
#include "TPZSpStructMatrix_Eigen.h"
#include "TPZSpMatrixEigen.h"
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
//        TPZSpStructMatrix matrix(Mesh());
//        matrix.SetNumThreads(n_threads);
//        SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELU);
//        SetSolver(step);
        TPZSpStructMatrixEigen matrix(Mesh());
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

void TMRSTransportAnalysis::Assemble(){
    if (m_parallel_execution_Q) {
        fTransportSpMatrix->Assemble();
//        Assemble_parallel();
    }else{
        Assemble_serial();
    }
}

void TMRSTransportAnalysis::Assemble_serial(){
//    TPZAnalysis::Assemble();
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
    if(!this->fSolver){
        DebugStop();
    }
    TPZMatrix<STATE> *mat = 0;
    if(!fSolver->Matrix())
    {
        mat = fStructMatrix->Create();
    
        fSolver->SetMatrix(mat);
    }
    else
    {
        mat = fSolver->Matrix().operator->();
    }
//    mat->Redim(ncells, ncells);
//    fRhs.Redim(ncells,1);
    mat->Zero();
    fRhs.Zero();
    //Volumetric Elements
    for (int ivol = 0; ivol<ncells; ivol++) {
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1,1);
        fAlgebraicTransport.Contribute(ivol, elmat, ef);
        mat->AddSub(eqindex, eqindex, elmat);
        fRhs.AddSub(eqindex, 0,ef);
    }
   
    //Interface Elements
    int interID = m_sim_data->mTGeometry.mInterface_material_id;
    int ninterfaces = fAlgebraicTransport.fInterfaceData[interID].fFluxSign.size();
    if (ninterfaces<1) {
        DebugStop();
    }
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[interID].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];
        TPZVec<int64_t> destinationindex(2);
        destinationindex[0]=lefteq;
        destinationindex[1]=righteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterface(interf,elmat, ef);
        mat->AddKel(elmat, destinationindex);
        fRhs.AddFel(ef, destinationindex);
    }
    
    int inlet_mat_id = -2;
    //INLET
    ninterfaces = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fFluxSign.size();
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_mat_id].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
//        int right = lrindex.second;
        TPZVec<int64_t> destinationindex(1);
        destinationindex[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(interf,ef);
        fRhs.AddFel(ef, destinationindex);
    }
    
    int outlet_mat_id = -4;
    //outlet
    ninterfaces = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fFluxSign.size();
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[outlet_mat_id].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        TPZVec<int64_t> destinationindex(1);
        destinationindex[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterface(interf,elmat,ef);
        mat->AddKel(elmat, destinationindex);
        fRhs.AddFel(ef, destinationindex);
    }
}

void TMRSTransportAnalysis::AssembleResidual(){
    if (m_parallel_execution_Q) {
        AssembleResidual_Eigen();
    }else{
        AssembleResidual_serial();
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

    
    TPZFMatrix<STATE> dx(Solution()),x(Solution());
    TPZFMatrix<STATE> correction(Solution());
    correction.Zero();
    
    ComputeInitialGuess(x); // from the linear problem (tangent and residue)
    bool QN_converge_Q = QuasiNewtonSteps(x,20); // assuming linear operator (tangent)
    if(QN_converge_Q){
        return;
    }
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
       
        NewtonIteration();
        dx = Solution();

        x += dx;
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        fAlgebraicTransport.fCellsData.UpdateSaturations(x);
        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTNumerics.m_ISLinearKrModelQ);
    
        AssembleResidual();
        corr_norm = Norm(dx);
        res_norm = Norm(Rhs());
        
#ifdef PZDEBUG
        {
            if(std::isnan(corr_norm) || std::isnan(res_norm))
            {
                DebugStop();
            }
        }
#endif
        std::cout << "res_norm " << res_norm << " corr_norm " << corr_norm << std::endl;
        stop_criterion_Q = (res_norm < res_tol);
        stop_criterion_corr_Q = (corr_norm < corr_tol);
        if (stop_criterion_Q || stop_criterion_corr_Q) {
            std::cout << "Transport operator: Converged" << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            std::cout << "residue norm = " << res_norm << std::endl;
            break;
        }

    }
    
}

void TMRSTransportAnalysis::ComputeInitialGuess(TPZFMatrix<STATE> &x){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(true);
    
    LoadSolution(x);
    cmesh->LoadSolutionFromMultiPhysics();
    
    NewtonIteration();
    x += Solution();
    LoadSolution(x);
    cmesh->LoadSolutionFromMultiPhysics();
//    PostProcessTimeStep();
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(true);
//    PostProcessTimeStep();
    AssembleResidual();
    REAL res_norm = Norm(Rhs());
    std::cout << "Initial guess residue norm : " <<  res_norm << std::endl;
}

bool TMRSTransportAnalysis::QuasiNewtonSteps(TPZFMatrix<STATE> &x, int n){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    std::cout << "Quasi-Newton process : " <<  std::endl;
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        NewtonIteration();
        
        x += Solution();
        
#ifdef PZDEBUG
        {
            REAL norm = Norm(x);
            if(std::isnan(norm))
            {
                DebugStop();
            }
        }
#endif
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        
        
        if (m_sim_data->mTNumerics.m_ISLinearKrModelQ) {
            fAlgebraicTransport.fCellsData.UpdateSaturations(x);
            fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(true);
        }else{
            fAlgebraicTransport.fCellsData.UpdateSaturations(x);
            fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambdaQuasiNewton();
        }
        
        AssembleResidual();
        REAL res_norm = Norm(Rhs());
        std::cout << " Residue norm : " <<  res_norm << std::endl;
        
        res_norm = Norm(Rhs());
        
#ifdef PZDEBUG
        if(std::isnan(res_norm))
        {
            DebugStop();
        }
#endif
        
        bool stop_criterion_Q = (res_norm < res_tol);
        if (stop_criterion_Q) {
            std::cout << "Transport operator: Converged" << std::endl;
            std::cout << "Quasi-Newton iterations = " << m_k_iteration << std::endl;
            std::cout << "Residue norm = " << res_norm << std::endl;
            return true;
        }
    }
    
    return false;
}

void TMRSTransportAnalysis::NewtonIteration(){
    
    if (m_parallel_execution_Q) {
#ifdef USING_BOOST
    boost::posix_time::ptime tnewton1 = boost::posix_time::microsec_clock::local_time();
#endif

        NewtonIteration_Eigen();
#ifdef USING_BOOST
        boost::posix_time::ptime tnewton2 = boost::posix_time::microsec_clock::local_time();
        auto deltat = tnewton2-tnewton1;
        std::cout << "Transport:: Parallel Newton step time: " << deltat << std::endl;
#endif
    }else{
#ifdef USING_BOOST
        boost::posix_time::ptime tnewtons1 = boost::posix_time::microsec_clock::local_time();
#endif

        NewtonIteration_serial();
        #ifdef USING_BOOST
            boost::posix_time::ptime tnewtons2 = boost::posix_time::microsec_clock::local_time();
            auto deltats = tnewtons2-tnewtons1;
            std::cout << "Transport:: Serial Newton step time: " << deltats << std::endl;
        #endif
    }
}

void TMRSTransportAnalysis::NewtonIteration_serial(){
    
    Assemble();
    Rhs() *= -1.0;
//    Solver().Matrix()->Print("j = ",std::cout,EMathematicaInput);
//    Rhs().Print("r = ",std::cout,EMathematicaInput);
    Solve();

}

void TMRSTransportAnalysis::AnalyzePattern(){
//    Assemble_mass_parallel();
//    Assemble_parallel();
//    m_transmissibility += m_mass;
//    m_analysis.analyzePattern(m_transmissibility);
    fTransportSpMatrix = new TPZAnalysisAuxEigen(&fAlgebraicTransport);
//    fTransportSpMatrix->SetAlgebraicTransport(&fAlgebraicTransport);
    fTransportSpMatrix->AnalyzePattern();
    // Because in some parts this objects are needed.
    Solution().Resize(fTransportSpMatrix->NRhsRows(), 1);
    Rhs().Resize(fTransportSpMatrix->NRhsRows(), 1);

}

void TMRSTransportAnalysis::NewtonIteration_Eigen(){
    
//    Assemble_parallel();
    fTransportSpMatrix->Assemble();
#ifdef PZDEBUG
    {
        auto norm = fTransportSpMatrix->RhsNorm();
        if(std::isnan(norm))
        {
            DebugStop();
        }
    }
#endif

    fTransportSpMatrix->Solve();
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> ds = fTransportSpMatrix->Solution();
    assert(Solution().Rows() == ds.rows());
    
#ifdef USING_TBB
    tbb::parallel_for(size_t(0), size_t(ds.rows()), size_t(1),
        [this,&ds] (size_t & i){
        Solution()(i,0) = ds(i,0);
        }
    );
#else
    for (int i = 0; i < ds.rows(); i++) {
        Solution()(i,0) = ds(i,0);
    }
#endif
#ifdef PZDEBUG
    STATE norm = Norm(Solution());
    if(std::isnan(norm)) DebugStop();
#endif
}


void TMRSTransportAnalysis::AssembleResidual_serial(){
    
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
     if(!this->fSolver){
         DebugStop();
     }
     fRhs.Resize(ncells,1);
     fRhs.Zero();
     
     //Volumetric Elements
     for (int ivol = 0; ivol<ncells; ivol++) {
         int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
         TPZFMatrix<double> ef;
         ef.Resize(1,1);
         fAlgebraicTransport.ContributeResidual(ivol, ef);
         fRhs.AddSub(eqindex, 0,ef);
     }
    
     //Interface Elements
     int interID = m_sim_data->mTGeometry.mInterface_material_id;
     int ninterfaces = fAlgebraicTransport.fInterfaceData[interID].fFluxSign.size();
     if (ninterfaces<1) {
         DebugStop();
     }
     for (int interf = 0; interf<ninterfaces; interf++) {
         std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[interID].fLeftRightVolIndex[interf];
         int left = lrindex.first;
         int right = lrindex.second;
         int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
         int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];
         TPZVec<int64_t> destinationindex(2);
         destinationindex[0]=lefteq;
         destinationindex[1]=righteq;
         TPZFMatrix<double> ef;
         ef.Resize(2, 1);
         fAlgebraicTransport.ContributeInterfaceResidual(interf, ef);
         fRhs.AddFel(ef, destinationindex);
     }
     
     int inlet_mat_id = -2;
     //INLET
     ninterfaces = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fFluxSign.size();
     for (int interf = 0; interf<ninterfaces; interf++) {
         std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_mat_id].fLeftRightVolIndex[interf];
         int left = lrindex.first;
         int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
         TPZVec<int64_t> destinationindex(1);
         destinationindex[0]=lefteq;
         TPZFMatrix<double> ef;
         ef.Resize(1, 1);
         fAlgebraicTransport.ContributeBCInletInterface(interf,ef);
         fRhs.AddFel(ef, destinationindex);
     }
     
     int outlet_mat_id = -4;
     //outlet
     ninterfaces = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fFluxSign.size();
     for (int interf = 0; interf<ninterfaces; interf++) {
         std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[outlet_mat_id].fLeftRightVolIndex[interf];
         int left = lrindex.first;
         int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
         TPZVec<int64_t> destinationindex(1);
         destinationindex[0]=lefteq;
         TPZFMatrix<double> ef;
         ef.Resize(1, 1);
         fAlgebraicTransport.ContributeBCOutletInterfaceResidual(interf,ef);
         fRhs.AddFel(ef, destinationindex);
     }
}

void TMRSTransportAnalysis::AssembleResidual_Eigen(){
 
    fTransportSpMatrix->AssembleResidual();
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> r = fTransportSpMatrix->Rhs().toDense();
    assert(Rhs().Rows() == r.rows());
    #ifdef USING_TBB
        tbb::parallel_for(size_t(0), size_t(r.rows()), size_t(1),
            [this,&r] (size_t & i){
             Rhs()(i,0) = r(i,0);
            }
        );
    #else
        for (int i = 0; i < r.rows(); i++) {
            Rhs()(i,0) = r(i,0);
        }
    #endif
    
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

void TMRSTransportAnalysis::UpdateInitialSolutionFromCellsData(){
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    fAlgebraicTransport.fCellsData.UpdateSaturationsTo(Solution());
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTNumerics.m_ISLinearKrModelQ);
    LoadSolution();
    cmesh->LoadSolutionFromMultiPhysics();
}
