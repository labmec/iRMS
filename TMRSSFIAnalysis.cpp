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
    int n_mixed_dof = m_mixed_module->Solution().Rows();
    int n_transport_dof = m_transport_module->Solution().Rows();
    m_x_mixed.Resize(n_mixed_dof, 1);
    m_x_transport.Resize(n_transport_dof, 1);
}

void TMRSSFIAnalysis::Configure(int n_threads, bool UsePardiso_Q){
    m_mixed_module->Configure(n_threads, UsePardiso_Q);
    m_transport_module->Configure(n_threads, UsePardiso_Q);
}

void TMRSSFIAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    m_mixed_module->SetDataTransfer(sim_data);
    m_transport_module->SetDataTransfer(sim_data);
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
    bool stop_criterion_Q = false;
    REAL error_rel_mixed = 1.0;
    REAL error_rel_transport = 1.0;
    REAL eps_tol = 0.01;

    for (int i = 1; i <= n_iterations; i++) {
        
        SFIIteration();
        error_rel_mixed = Norm(m_x_mixed - m_mixed_module->Solution())/Norm(m_mixed_module->Solution());
        error_rel_transport = Norm(m_x_transport - m_transport_module->Solution())/Norm(m_transport_module->Solution());
        
        stop_criterion_Q = error_rel_mixed <= eps_tol && error_rel_transport <= eps_tol;
        if (stop_criterion_Q) {
            std::cout << "SFI converged " << std::endl;
            std::cout << "Number of iterations = " << i << std::endl;
            UpdateMemoryInModules();
            break;
        }
        
        m_x_mixed = m_mixed_module->Solution();
        m_x_transport = m_transport_module->Solution();
    }
    
}

void TMRSSFIAnalysis::PostProcessTimeStep(){
    m_mixed_module->PostProcessTimeStep();
    m_transport_module->PostProcessTimeStep();
}

void TMRSSFIAnalysis::SFIIteration(){
    
    {
#ifdef USING_BOOST
        boost::posix_time::ptime mixed_process_t1 = boost::posix_time::microsec_clock::local_time();
#endif
        m_mixed_module->RunTimeStep();
        
        {
            TPZStack<std::string> scalnames, vecnames;
            std::string file("mixed.vtk");
            scalnames.push_back("Pressure");
            vecnames.push_back("Flux");
            m_mixed_module->DefineGraphMesh(2, scalnames, vecnames, file);
            m_mixed_module->PostProcess(0, 2);
        }
        
      
#ifdef USING_BOOST
        boost::posix_time::ptime mixed_process_t2 = boost::posix_time::microsec_clock::local_time();
        REAL mixed_process_time = boost::numeric_cast<double>((mixed_process_t2-mixed_process_t1).total_milliseconds());
        std::cout << "Mixed approximation performed in :" << setw(10) <<  mixed_process_time/1000.0 << setw(5)   << " seconds." << std::endl;
#endif
    }

    
    TransferToTransportModule();    // Transfer to transport
    m_transport_module->RunTimeStep();
    

    TransferToMixedModule();        // Transfer to mixed
}

void TMRSSFIAnalysis::TransferToTransportModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }
    
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
    
    transport_cmesh->LoadSolutionFromMeshes();
    
}

void TMRSSFIAnalysis::TransferToMixedModule(){
    
    TPZMultiphysicsCompMesh * mixed_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh());
    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }
    
    transport_cmesh->LoadSolutionFromMultiPhysics();
    // Saturations are transferred to mixed module
    int s_b = 2;
    TPZFMatrix<STATE> & s_dof = transport_cmesh->MeshVector()[s_b]->Solution();
    mixed_cmesh->MeshVector()[s_b]->LoadSolution(s_dof);
    //
    
TPZMFSolutionTransfer trans;
   trans.BuildTransferData(m_mixed_module->Mesh());
    
    int count =0;
    int ncorrespond = trans.fmeshTransfers.size();
    for (int icorres=0; icorres<ncorrespond; icorres++) {
        int nmach = trans.fmeshTransfers[icorres].fconnecttransfer.size();
        for (int imach=0; imach<nmach; imach++) {
            std::cout<<"Correspond: "<<icorres<<" match number: "<<imach<<std::endl;
            std::cout<<"inicial: "<<trans.fmeshTransfers[icorres].fconnecttransfer[imach].fblocknumber;
            std::pair<TPZCompMesh*, int64_t> match =trans.fmeshTransfers[icorres].fconnecttransfer[imach].fblockTarget;
            std::cout<<"target: "<<std::get<1>(match)<<std::endl;
            std::cout<<"count: "<<count<<std::endl;
            count++;
        }
    }
    
//
    
    std::ofstream filemixed("mixedoperator.txt");
    mixed_cmesh->Print(filemixed);
    std::ofstream filePressure("pressureperator.txt");
    mixed_cmesh->MeshVector()[1]->Print(filePressure);
    std::ofstream fileflux("fluxperator.txt");
    mixed_cmesh->MeshVector()[0]->Print(fileflux);
    
   
    TPZMFSolutionTransfer soltrans;
    soltrans.BuildTransferData(mixed_cmesh ) ;
    
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
    UpdateMemoryMixedModule();
    UpdateMemoryTransportModule();
}
