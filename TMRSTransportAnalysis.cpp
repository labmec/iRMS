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

void TMRSTransportAnalysis::Assemble(){
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
    int interID = m_sim_data->mTGeometry.Interface_material_id;
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
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
    if(!this->fSolver){
        DebugStop();
    }
    fRhs.Resize(ncells,1);
    fRhs.Zero();
    
    //Volumetric Elements
    for (int ivol = 0; ivol<ncells; ivol++) {
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1,1);
        fAlgebraicTransport.Contribute(ivol, elmat, ef); // @jose -> res at element level
        fRhs.AddSub(eqindex, 0,ef);
    }
   
    //Interface Elements
    int interID = m_sim_data->mTGeometry.Interface_material_id;
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
        fRhs.AddFel(ef, destinationindex);
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
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda();
    
    QuasiNewtonSteps(x,10); // assuming linear operator
    
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda();
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
       
        NewtonIteration();
        dx = Solution();

        x += dx;
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        fAlgebraicTransport.fCellsData.UpdateSaturations(x);
        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda();
        corr_norm = Norm(dx);

        cmesh->LoadSolutionFromMultiPhysics();
        
        AssembleResidual();
        res_norm = Norm(Rhs());
        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_corr_Q) {
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

    NewtonIteration();
    x += Solution();
//    x.Print("mat=",std::cout,EMathematicaInput);
    LoadSolution(x);
    cmesh->LoadSolutionFromMultiPhysics();
//    PostProcessTimeStep();
}

void TMRSTransportAnalysis::QuasiNewtonSteps(TPZFMatrix<STATE> &x, int n){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    
    std::cout << "Quasi-Newton process : " <<  std::endl;
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        NewtonIteration();

        x += Solution();
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        fAlgebraicTransport.fCellsData.UpdateSaturations(x);
        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambdaQuasiNewton();

        cmesh->LoadSolutionFromMultiPhysics();
//        PostProcessTimeStep();
        
        AssembleResidual();
        REAL res_norm = Norm(Rhs());
        std::cout << " Residue norm : " <<  res_norm << std::endl;
    }
    
}

void TMRSTransportAnalysis::NewtonIteration(){
    
    Assemble();
    Rhs() *= -1.0;
    Solve();
//    this->PostProcessTimeStep();
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
    LoadSolution();
    cmesh->LoadSolutionFromMultiPhysics();
}
