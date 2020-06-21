//
//  TMRSSFIAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSSFIAnalysis.h"
#include "TPZMFSolutionTransfer.h"


TMRSSFIAnalysis::TMRSSFIAnalysis(){
    m_sim_data = nullptr;
    m_k_iteration = 0;
    m_mixed_module = nullptr;
    m_transport_module = nullptr;
    m_x_mixed.Resize(0, 0);
    m_x_transport.Resize(0, 0);
}

TMRSSFIAnalysis::~TMRSSFIAnalysis(){
    
}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
   
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    
     fAlgebraicDataTransfer.TransferPermeabiliyTensor();
    
//    fAlgebraicDataTransfer.TransferPermeabiliyTensor();
   
//    int n_mixed_dof = m_mixed_module->Solution().Rows();
//    int n_transport_dof = m_transport_module->Solution().Rows();
//    m_x_mixed.Resize(n_mixed_dof, 1);
//    m_x_transport.Resize(n_transport_dof, 1);
    
}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<REAL(const TPZVec<REAL> & )> & kx, std::function<REAL(const TPZVec<REAL> & )> & ky, std::function<REAL(const TPZVec<REAL> & )> & kz, std::function<REAL(const TPZVec<REAL> & )> & phi, std::function<REAL(const TPZVec<REAL> & )> & s0){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
   
    fAlgebraicDataTransfer.fkx = kx;
    fAlgebraicDataTransfer.fky = ky;
    fAlgebraicDataTransfer.fkz = kz;
    fAlgebraicDataTransfer.fphi = phi;
    fAlgebraicDataTransfer.fs0 = s0;
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
//    fAlgebraicDataTransfer.CheckDataTransferTransportToMixed();
    fAlgebraicDataTransfer.TransferPermeabiliyTensor();
    
//    fAlgebraicDataTransfer.CheckDataTransferTransportToMixed();
//    fAlgebraicDataTransfer.TransferPermeabiliyTensor();
   
//    int n_mixed_dof = m_mixed_module->Solution().Rows();
//    int n_transport_dof = m_transport_module->Solution().Rows();
//    m_x_mixed.Resize(n_mixed_dof, 1);
//    m_x_transport.Resize(n_transport_dof, 1);
    
}

void TMRSSFIAnalysis::Configure(int n_threads, bool UsePardiso_Q){
    m_mixed_module->Configure(n_threads, UsePardiso_Q);
    m_transport_module->Configure(n_threads, UsePardiso_Q);
   
}

void TMRSSFIAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    m_mixed_module->SetDataTransfer(sim_data);
    m_transport_module->SetDataTransfer(sim_data);
    
    m_transport_module->fAlgebraicTransport.interfaceid = 100;
    m_transport_module->fAlgebraicTransport.inletmatid = -2;
    m_transport_module->fAlgebraicTransport.outletmatid = -4;
   
    m_transport_module->fAlgebraicTransport.fgravity = m_sim_data->mTNumerics.m_gravity;
    
    //Set initial properties
    std::cout<<m_sim_data->mTGeometry.Interface_material_id;
    std::cout<<m_sim_data->mTFluidProperties.mWaterViscosity;
    m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[0] = m_sim_data->mTFluidProperties.mWaterViscosity;
    m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[1] = m_sim_data->mTFluidProperties.mOilViscosity;
    int ncells = m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil.size();
    REAL rhow = m_sim_data->mTFluidProperties.mWaterDensity;
    REAL rhoo = m_sim_data->mTFluidProperties.mOilDensity;
    for (int icell =0; icell<ncells; icell++) {
        m_transport_module->fAlgebraicTransport.fCellsData.fDensityWater[icell]= rhow; m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil[icell]= rhoo;
    }
}

TMRSDataTransfer * TMRSSFIAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSSFIAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSSFIAnalysis::RunTimeStep(){
    
    m_x_mixed = m_mixed_module->Solution();
    m_x_transport = m_transport_module->Solution();
    
    
    int n_iterations = m_sim_data->mTNumerics.m_max_iter_sfi;
    REAL eps_tol = m_sim_data->mTNumerics.m_sfi_tol;
    bool stop_criterion_Q = false;
    REAL error_rel_mixed = 1.0;
    REAL error_rel_transport = 1.0;

    for (int i = 1; i <= n_iterations; i++) {
        
        SFIIteration();
        error_rel_mixed = Norm(m_x_mixed - m_mixed_module->Solution())/Norm(m_mixed_module->Solution());
        
        error_rel_transport = Norm(m_x_transport - m_transport_module->Solution())/Norm(m_transport_module->Solution());
        
        stop_criterion_Q = error_rel_mixed < eps_tol && error_rel_transport < eps_tol;
        if (stop_criterion_Q && i > 1) {
            std::cout << "SFI converged " << std::endl;
            std::cout << "Number of iterations = " << i << std::endl;
            std::cout << "Mixed problem variations = " << error_rel_mixed << std::endl;
            std::cout << "Transport problem variations = " << error_rel_transport << std::endl;
            m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
            break;
        }
        
        m_x_mixed = m_mixed_module->Solution();
        m_x_transport = m_transport_module->Solution();
 
//        m_mixed_module->PostProcessTimeStep();
    }
    
}

void TMRSSFIAnalysis::PostProcessTimeStep(int val){
    if (val == 0) {
        m_mixed_module->PostProcessTimeStep();
        m_transport_module->PostProcessTimeStep();
    }
    if (val == 1) {
        m_mixed_module->PostProcessTimeStep();
    }
    if (val == 2) {
        m_transport_module->PostProcessTimeStep();
    }
   
}

void TMRSSFIAnalysis::SFIIteration(){
    
    {
        
        m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTNumerics.m_ISLinearKrModelQ);
        m_transport_module->fAlgebraicTransport.fCellsData.UpdateMixedDensity();
        fAlgebraicDataTransfer.TransferLambdaCoefficients();
       

        m_mixed_module->RunTimeStep();
        
#ifdef USING_BOOST2
        boost::posix_time::ptime mixed_process_t1 = boost::posix_time::microsec_clock::local_time();
#endif
        
#ifdef USING_BOOST2
        boost::posix_time::ptime mixed_process_t2 = boost::posix_time::microsec_clock::local_time();
        REAL mixed_process_time = boost::numeric_cast<double>((mixed_process_t2-mixed_process_t1).total_milliseconds());
        std::cout << "Mixed approximation performed in :" << setw(10) <<  mixed_process_time/1000.0 << setw(5)   << " seconds." << std::endl;
#endif
    }
    
    //   m_mixed_module->PostProcessTimeStep();
    
  
    fAlgebraicDataTransfer.TransferMixedMeshMultiplyingCoefficients();
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(100);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-2);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(-4);
    m_transport_module->fAlgebraicTransport.fdt = m_transport_module->GetCurrentTime();
    m_transport_module->RunTimeStep();
    
//    m_transport_module->PostProcessTimeStep();
    //    m_transport_module->PostProcessTimeStep();
    //    m_transport_module->Solution() = m_transport_module->Solution() + solution_n;
    //    TransferToMixedModule();        // Transfer to mixed
    
 }


void TMRSSFIAnalysis::TransferToTransportModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
   
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }
    
//     m_mixed_module->m_soltransportTransfer.TransferFromMultiphysics();
    mixed_cmesh->LoadSolutionFromMultiPhysics();

    // flux and pressure are transferred to transport module
    int q_b = 0;
    int p_b = 1;
    TPZFMatrix<STATE> & q_dof = mixed_cmesh->MeshVector()[q_b]->Solution();
    TPZFMatrix<STATE> & p_dof = mixed_cmesh->MeshVector()[p_b]->Solution();
    transport_cmesh->MeshVector()[q_b]->LoadSolution(q_dof);
    transport_cmesh->MeshVector()[p_b]->LoadSolution(p_dof);
    
    if (m_sim_data->mTNumerics.m_four_approx_spaces_Q) {
        // average flux and pressure are transferred to transport module
        int qavg_b = 2;
        int pavg_b = 3;
        TPZFMatrix<STATE> & q_dof = mixed_cmesh->MeshVector()[qavg_b]->Solution();
        TPZFMatrix<STATE> & p_dof = mixed_cmesh->MeshVector()[pavg_b]->Solution();
        transport_cmesh->MeshVector()[qavg_b]->LoadSolution(q_dof);
        transport_cmesh->MeshVector()[pavg_b]->LoadSolution(p_dof);
    }

//    m_transport_module->m_soltransportTransfer.TransferToMultiphysics();
    transport_cmesh->LoadSolutionFromMeshes();

}

void TMRSSFIAnalysis::TransferToMixedModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }

//    m_transport_module->m_soltransportTransfer.TransferFromMultiphysics();
    transport_cmesh->LoadSolutionFromMultiPhysics();

    // Saturations are transferred to mixed module
    int s_b = 2;
    TPZFMatrix<STATE> & s_dof = transport_cmesh->MeshVector()[s_b]->Solution();
    mixed_cmesh->MeshVector()[s_b]->LoadSolution(s_dof);
    

//     m_mixed_module->m_soltransportTransfer.TransferToMultiphysics();
    mixed_cmesh->LoadSolutionFromMeshes();
    
}

void TMRSSFIAnalysis::UpdateMemoryMixedModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    
    if (!mixed_cmesh) {
        DebugStop();
    }
    
    // Updating memory
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, mixed_cmesh, true);
    m_mixed_module->AssembleResidual();
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, mixed_cmesh, false);
    
}

void TMRSSFIAnalysis::UpdateMemoryTransportModule(){
 
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    if (!transport_cmesh) {
        DebugStop();
    }
    
    // Updating memory
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, transport_cmesh, true);
    m_transport_module->AssembleResidual();
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, *m_sim_data, transport_cmesh, false);
}

void TMRSSFIAnalysis::UpdateMemoryInModules(){
//    UpdateMemoryMixedModule();
    UpdateMemoryTransportModule();
}

#include "TPZFastCondensedElement.h"
// transfer the permeability and lambda to the element solution for post processing
void TMRSSFIAnalysis::SetMixedMeshElementSolution(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    cmesh->ElementSolution().Redim(nel, 4);
    int64_t count = 0;
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(submesh)
        {
            SetMixedMeshElementSolution(submesh);
            continue;
        }
        TPZFastCondensedElement *fast = dynamic_cast<TPZFastCondensedElement *>(cel);
        if(!fast) continue;
        count++;
        cmesh->ElementSolution()(el,0) = fast->GetPermTensor()(0,0);
        cmesh->ElementSolution()(el,1) = fast->GetPermTensor()(1,1);
        cmesh->ElementSolution()(el,2) = fast->GetPermTensor()(2,2);
        cmesh->ElementSolution()(el,3) = fast->GetLambda();
    }
}

