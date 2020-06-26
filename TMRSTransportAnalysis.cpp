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

void TMRSTransportAnalysis::Assemble(){
    if (m_parallel_execution_Q) {
        Assemble_parallel();
    }else{
        Assemble_serial();
    }
}

void TMRSTransportAnalysis::Assemble_serial(){
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
    if(!this->fSolver){
        DebugStop();
    }
    TPZMatrix<STATE> *mat = fStructMatrix->Create();
    fSolver->SetMatrix(mat);
    mat->Resize(ncells, ncells);
    fRhs.Resize(ncells,1);
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
        AssembleResidual_parallel();
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
        
        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
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
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(true);
    
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
        bool stop_criterion_Q = res_norm < res_tol;
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

        NewtonIteration_parallel();
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
    Assemble_mass_parallel();
    Assemble_parallel();
    m_transmissibility += m_mass;
    m_analysis.analyzePattern(m_transmissibility);
    
    // Because in some parts this objects are needed.
    Solution().Resize(m_rhs.rows(), 1);
    Rhs().Resize(m_rhs.rows(), 1);
}

void TMRSTransportAnalysis::NewtonIteration_parallel(){
    
    Assemble_parallel();

    m_transmissibility += m_mass;
    m_rhs *= -1.0;
    m_analysis.factorize(m_transmissibility);
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> ds = m_analysis.solve(m_rhs);
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
}

void TMRSTransportAnalysis::Assemble_mass_parallel(){
    
    int n_cells = fAlgebraicTransport.fCellsData.fVolume.size();
    m_mass =  Eigen::SparseMatrix<REAL>( n_cells, n_cells );
    
    m_mass_triplets.resize(n_cells);
    //Volumetric Elements
    
#ifdef USING_TBB
    tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
        [this] (size_t & ivol){
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1,1);
        fAlgebraicTransport.Contribute(ivol, elmat, ef);
        m_mass_triplets[ivol] = Triplet<REAL>(eqindex,eqindex, elmat(0,0));
        }
    );
#else
    for (int ivol = 0; ivol < n_cells; ivol++) {
      int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
      TPZFMatrix<double> elmat, ef;
      elmat.Resize(1, 1);
      ef.Resize(1,1);
      fAlgebraicTransport.Contribute(ivol, elmat, ef);
      m_mass_triplets[ivol] = Triplet<REAL>(eqindex,eqindex, elmat(0,0));
    }
#endif
    
    m_mass.setFromTriplets( m_mass_triplets.begin(), m_mass_triplets.end() );
    m_mass_triplets.clear();
}

void TMRSTransportAnalysis::Assemble_parallel(){
    
    int n_cells = fAlgebraicTransport.fCellsData.fVolume.size();
    m_transmissibility =  Eigen::SparseMatrix<REAL>( n_cells, n_cells );
    m_rhs =  Eigen::SparseMatrix<REAL>( n_cells, 1 );
    if(!this->fSolver){
        DebugStop();
    }
    
    m_rhs.setZero();
    m_transmissibility.setZero();
    
    int internal_faces_id = m_sim_data->mTGeometry.mInterface_material_id;
    int inlet_faces_id = -2;
    int outlet_faces_id = -4;

    int n_internal_faces = fAlgebraicTransport.fInterfaceData[internal_faces_id].fFluxSign.size();
    int n_inlet_faces = fAlgebraicTransport.fInterfaceData[inlet_faces_id].fFluxSign.size();
    int n_outlet_faces = fAlgebraicTransport.fInterfaceData[outlet_faces_id].fFluxSign.size();

    size_t n_nzeros_trans = n_internal_faces * 4 + n_outlet_faces;
    size_t n_nzeros_res = n_cells + n_internal_faces * 2 + n_inlet_faces + n_outlet_faces;
    m_trans_triplets.resize(n_nzeros_trans);
    m_rhs_triplets.resize(n_nzeros_res);
    
#ifdef USING_TBB
    
    tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
        [this] (size_t & ivol){
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> ef;
        ef.Resize(1,1);
        fAlgebraicTransport.ContributeResidual(ivol, ef);
        m_rhs_triplets[ivol] = Triplet<REAL>(eqindex,0, ef(0,0));
        }
    );
    
    tbb::parallel_for(size_t(0), size_t(n_internal_faces), size_t(1),
        [this,&internal_faces_id,&n_cells] (size_t & iface){
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];

        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterface(iface,elmat, ef);
        size_t i_begin = 2*2*(iface);
        m_trans_triplets[i_begin] = (Triplet<REAL>(indexes[0],indexes[0], elmat(0,0)));
        m_trans_triplets[i_begin+1] = (Triplet<REAL>(indexes[0],indexes[1], elmat(0,1)));
        m_trans_triplets[i_begin+2] = (Triplet<REAL>(indexes[1],indexes[0], elmat(1,0)));
        m_trans_triplets[i_begin+3] = (Triplet<REAL>(indexes[1],indexes[1], elmat(1,1)));
        
        size_t i_rhs_begin = 2*(iface) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet<REAL>(indexes[1],0, ef(1,0));
        }
    );
    
    tbb::parallel_for(size_t(0), size_t(n_inlet_faces), size_t(1),
        [this,&inlet_faces_id,&n_internal_faces,&n_cells] (size_t & iface){
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(iface,ef);
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        }
    );
    
    tbb::parallel_for(size_t(0), size_t(n_outlet_faces), size_t(1),
        [this,&outlet_faces_id,&n_internal_faces,&n_inlet_faces,&n_cells] (size_t & iface){
        std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];

        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterface(iface,elmat, ef);
        size_t i_begin = iface + n_internal_faces * 4;
        m_trans_triplets[i_begin] = (Triplet<REAL>(indexes[0],indexes[0], elmat(0,0)));
        
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces + n_inlet_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        }
    );
    
#else
    
    for (int ivol = 0; ivol < n_cells; ivol++) {
      int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
      TPZFMatrix<double> ef;
      ef.Resize(1,1);
      fAlgebraicTransport.ContributeResidual(ivol, ef);
      m_rhs_triplets[ivol] = Triplet<REAL>(eqindex,0, ef(0,0));
    }
    
    for (int iface = 0; iface < n_internal_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];

        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterface(iface,elmat, ef);
        size_t i_begin = 2*2*(iface);
        m_trans_triplets[i_begin] = (Triplet<REAL>(indexes[0],indexes[0], elmat(0,0)));
        m_trans_triplets[i_begin+1] = (Triplet<REAL>(indexes[0],indexes[1], elmat(0,1)));
        m_trans_triplets[i_begin+2] = (Triplet<REAL>(indexes[1],indexes[0], elmat(1,0)));
        m_trans_triplets[i_begin+3] = (Triplet<REAL>(indexes[1],indexes[1], elmat(1,1)));
        
        size_t i_rhs_begin = 2*(iface) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet<REAL>(indexes[1],0, ef(1,0));
    }
    
    for (int iface = 0; iface < n_inlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(iface,ef);
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
    }
    
    for (int iface = 0; iface < n_outlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];

        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterface(iface,elmat, ef);
        size_t i_begin = iface + n_internal_faces * 4;
        m_trans_triplets[i_begin] = (Triplet<REAL>(indexes[0],indexes[0], elmat(0,0)));
        
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces + n_inlet_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
    }
    
#endif
    
    m_rhs.setFromTriplets( m_rhs_triplets.begin(), m_rhs_triplets.end() );
    m_rhs_triplets.clear();
    
    m_transmissibility.setFromTriplets( m_trans_triplets.begin(), m_trans_triplets.end() );
    m_trans_triplets.clear();
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

void TMRSTransportAnalysis::AssembleResidual_parallel(){
 
    int n_cells = fAlgebraicTransport.fCellsData.fVolume.size();
    m_rhs =  Eigen::SparseMatrix<REAL>( n_cells, 1 );
    
    int internal_faces_id = m_sim_data->mTGeometry.mInterface_material_id;
    int inlet_faces_id = -2;
    int outlet_faces_id = -4;

    int n_internal_faces = fAlgebraicTransport.fInterfaceData[internal_faces_id].fFluxSign.size();
    int n_inlet_faces = fAlgebraicTransport.fInterfaceData[inlet_faces_id].fFluxSign.size();
    int n_outlet_faces = fAlgebraicTransport.fInterfaceData[outlet_faces_id].fFluxSign.size();

    size_t n_nzeros_res = n_cells + n_internal_faces * 2 + n_inlet_faces + n_outlet_faces;
    m_rhs_triplets.resize(n_nzeros_res);
    m_rhs.setZero();
    
#ifdef USING_TBB
    
    tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
        [this] (size_t & ivol){
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> ef;
        ef.Resize(1,1);
        fAlgebraicTransport.ContributeResidual(ivol, ef);
        m_rhs_triplets[ivol] = Triplet<REAL>(eqindex,0, ef(0,0));
        }
    );
  
    
    tbb::parallel_for(size_t(0), size_t(n_internal_faces), size_t(1),
        [this,&internal_faces_id,&n_cells] (size_t & iface){
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];

        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> ef;
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterfaceResidual(iface, ef);
        
        size_t i_rhs_begin = 2*(iface) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet<REAL>(indexes[1],0, ef(1,0));
        }
    );
    
    tbb::parallel_for(size_t(0), size_t(n_inlet_faces), size_t(1),
        [this,&inlet_faces_id,&n_internal_faces,&n_cells] (size_t & iface){
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(iface,ef);
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        }
    );
    
    tbb::parallel_for(size_t(0), size_t(n_outlet_faces), size_t(1),
        [this,&outlet_faces_id,&n_internal_faces,&n_inlet_faces,&n_cells] (size_t & iface){
        std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];

        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double>  ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterfaceResidual(iface, ef);
        
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces + n_inlet_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        }
    );
    
#else
    
    for (int ivol = 0; ivol < n_cells; ivol++) {
      int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
      TPZFMatrix<double> ef;
      ef.Resize(1,1);
      fAlgebraicTransport.ContributeResidual(ivol, ef);
      m_rhs_triplets[ivol] = Triplet<REAL>(eqindex,0, ef(0,0));
    }
    
    for (int iface = 0; iface < n_internal_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];

        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> ef;
      
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterfaceResidual(iface, ef);
        
        size_t i_rhs_begin = 2*(iface) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet<REAL>(indexes[1],0, ef(1,0));
    }
    
    for (int iface = 0; iface < n_inlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(iface,ef);
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
    }
    
    for (int iface = 0; iface < n_outlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];

        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterfaceResidual(iface, ef);
        
        size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces + n_inlet_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet<REAL>(indexes[0],0, ef(0,0));
    }
    
#endif
    
    m_rhs.setFromTriplets( m_rhs_triplets.begin(), m_rhs_triplets.end() );
    m_rhs_triplets.clear();
    
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> r = m_rhs.toDense();
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
