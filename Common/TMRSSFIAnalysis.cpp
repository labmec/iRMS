//
//  TMRSSFIAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSSFIAnalysis.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZDarcyMemory.h"
#include "TPZFastCondensedElement.h"
#include "TPZDarcyFlowWithMem.h"

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

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZCompMesh * cmesh_transport, bool must_opt_band_width_Q){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    
    
}

void TMRSSFIAnalysis::BuildAlgebraicDataStructure(){
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    FillProperties();
}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<REAL(const TPZVec<REAL> & )> & kx, std::function<REAL(const TPZVec<REAL> & )> & ky, std::function<REAL(const TPZVec<REAL> & )> & kz, std::function<REAL(const TPZVec<REAL> & )> & phi, std::function<REAL(const TPZVec<REAL> & )> & s0){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    bool propsfromPre = m_sim_data->mTReservoirProperties.fPropsFromPreProcess;
    if (propsfromPre==false) {
            fAlgebraicDataTransfer.fkx = kx;
            fAlgebraicDataTransfer.fky = ky;
            fAlgebraicDataTransfer.fkz = kz;
            fAlgebraicDataTransfer.fphi = phi;
            fAlgebraicDataTransfer.fs0 = s0;
    }

    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    FillMaterialMemoryDarcy(1);
    std::string fileprops("Props.txt");

    std::ofstream file(fileprops);
    if(file && !propsfromPre){
        FillProperties(fileprops, &m_transport_module->fAlgebraicTransport);
    }

}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<std::vector<REAL>(const TPZVec<REAL> & )> & kappa_phi, std::function<REAL(const TPZVec<REAL> & )> & s0){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    
    fAlgebraicDataTransfer.fkappa_phi = kappa_phi;
    fAlgebraicDataTransfer.fs0 = s0;
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    
    FillMaterialMemoryDarcy(1);
    std::string fileprops("Props.txt");
    FillProperties(fileprops, &m_transport_module->fAlgebraicTransport);
}

void TMRSSFIAnalysis::Configure(int n_threads, bool UsePardiso_Q, bool usepz){
    m_mixed_module->Configure(n_threads, UsePardiso_Q, usepz);
    m_transport_module->Configure(n_threads, UsePardiso_Q);
    
}
void TMRSSFIAnalysis::FillMaterialMemoryDarcy(int material_id){
  
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }
    // Initialize integration points memory
    
    TPZCompMesh * cmesh = m_mixed_module->Mesh();
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    TPZDarcyFlowWithMem *mat = dynamic_cast<TPZDarcyFlowWithMem *>(material);

    if (!material)
        DebugStop();

    TPZMatWithMem<TPZDarcyMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPZDarcyMemory> * >(material);
    if (!mat_with_memory)
        DebugStop();

    std::shared_ptr<TPZAdmChunkVector<TPZDarcyMemory>> & memory_vector = mat_with_memory->GetMemory();

//    for(auto &meshit : fAlgebraicDataTransfer.fTransportMixedCorrespondence) {
//        int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
//#ifdef PZDEBUG
//        if(meshit.fTransport == 0)
//            DebugStop();
//
//#endif
//        for (int icell = 0; icell < ncells; icell++) {
//            int64_t algbindex = meshit.fAlgebraicTransportCellIndex[icell];
//            TPZFastCondensedElement * fastCond = meshit.fMixedCell[icell];
//            if(fastCond->Reference()->MaterialId()!=material_id){
//                continue;
//            }
//            TPZCompEl *celcomp = fastCond->ReferenceCompEl();
//
//
//            TPZManVector<int64_t> indices;
//            celcomp->GetMemoryIndices(indices);
//            indices.Print();
//            for (int index = 0; index<indices.size(); index++) {
//                int valIndex = indices[index];
//                TPZDarcyMemory &mem = memory_vector.get()->operator [](valIndex);
//                mem.fTransportCellIndex = algbindex;
//                m_transport_module->fAlgebraicTransport.fCellsData.fKx[algbindex] = 1.0;
//                m_transport_module->fAlgebraicTransport.fCellsData.fKy[algbindex] = 1.0;
//                m_transport_module->fAlgebraicTransport.fCellsData.fKz[algbindex] = 1.0;
//                m_transport_module->fAlgebraicTransport.fCellsData.fporosity[algbindex] = 0.1;
//            }
//        }
//    }
    
}
void TMRSSFIAnalysis::FillProperties(){
    
    bool propsfromPre = m_sim_data->mTReservoirProperties.fPropsFromPreProcess;
    if (propsfromPre==false) {
        if (m_sim_data->mTReservoirProperties.kappa_phi) {
            fAlgebraicDataTransfer.fkappa_phi = m_sim_data->mTReservoirProperties.kappa_phi;
            fAlgebraicDataTransfer.fs0=m_sim_data->mTReservoirProperties.s0;
            FillProperties(&m_transport_module->fAlgebraicTransport);
        }
        else{
            std::vector<REAL> kappa_phi(4,1.0);
//            kappa_phi[3]=1.0;
                FillProperties(&m_transport_module->fAlgebraicTransport, kappa_phi);
        }
        
    }
    else{
        if (m_sim_data->mTReservoirProperties.mPropsFileName=="") {
            std::string propsname ="PreProcess/props/"+ m_sim_data->mSimulationName + "Props_nL_"+ std::to_string(m_sim_data->mTGeometry.mnLayers)  +"_nRef_"+std::to_string(m_sim_data->mTGeometry.mnref)+".txt" ;
            m_sim_data->mTReservoirProperties.mPropsFileName = propsname;
        }
        
        std::ifstream file(m_sim_data->mTReservoirProperties.mPropsFileName);
        
        if(file){
            FillProperties(m_sim_data->mTReservoirProperties.mPropsFileName, &m_transport_module->fAlgebraicTransport);
        }
        else{
            
            std::cout<<"The properties of the reservoir have not been loaded. Please enter the properties in a text file or set the functions in the simulationdata object"<<std::endl;
            std::vector<REAL> kappa_phi(4,1.0e-7);
            kappa_phi[3]=0.1;
            FillProperties(&m_transport_module->fAlgebraicTransport, kappa_phi);
        }
        
    }
}
void TMRSSFIAnalysis::FillProperties(std::string fileprops, TPZAlgebraicTransport *algebraicTransport){
    
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }
    std::vector<REAL> Kx, Ky, Kz, Phi;
    ReadProperties(fileprops, false, Kx, Ky, Kz, Phi);
    
    TPZCompMesh * cmesh = m_mixed_module->Mesh();
    
//    for(auto &meshit : fAlgebraicDataTransfer.fTransportMixedCorrespondence)
//    {
//        int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
//#ifdef PZDEBUG
//        if(meshit.fTransport == 0)
//        {
//            DebugStop();
//        }
//#endif
//        for (int icell = 0; icell < ncells; icell++)
//        {
//            int64_t algbindex = meshit.fAlgebraicTransportCellIndex[icell];
//            TPZFastCondensedElement * fastCond = meshit.fMixedCell[icell];
//            TPZCompEl *celcomp = fastCond->ReferenceCompEl();
//            int geoIndexMixed = celcomp->Reference()->Index();
//            if (Kx[geoIndexMixed] <0) {
//                DebugStop();
//            }
//            algebraicTransport->fCellsData.fKx[algbindex] = Kx[geoIndexMixed];
//            algebraicTransport->fCellsData.fKy[algbindex] = Ky[geoIndexMixed];
//            algebraicTransport->fCellsData.fKz[algbindex] = Kz[geoIndexMixed];
//            algebraicTransport->fCellsData.fporosity[algbindex] = Phi[geoIndexMixed];
//        }
//    }
    m_transport_module->fAlgebraicTransport.fHasPropQ=true;
    
}
void TMRSSFIAnalysis::FillProperties(TPZAlgebraicTransport *algebraicTransport){
    
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }

    
    TPZCompMesh * cmesh = m_mixed_module->Mesh();
    
//    for(auto &meshit : fAlgebraicDataTransfer.fTransportMixedCorrespondence)
//    {
//        int64_t ncells = meshit.fAlgebraicTransportCellIndex.size();
//#ifdef PZDEBUG
//        if(meshit.fTransport == 0)
//        {
//            DebugStop();
//        }
//#endif
//        for (int icell = 0; icell < ncells; icell++)
//        {
//            int64_t algbindex = meshit.fAlgebraicTransportCellIndex[icell];
//            TPZFastCondensedElement * fastCond = meshit.fMixedCell[icell];
//            TPZCompEl *celcomp = fastCond->ReferenceCompEl();
//            TPZGeoEl *gel = celcomp->Reference();
//            int dim= gel->Dimension();
//            TPZVec<REAL> ximasscent(dim);
//            gel->CenterPoint(gel->NSides()-1, ximasscent);
//            std::vector<REAL> center(3,0.0);
//            TPZManVector<REAL,3> coord(3,0.0);
//            gel->X(ximasscent, coord);
//            
//            std::vector<REAL> kappa_phi = m_sim_data->mTReservoirProperties.kappa_phi(coord);
//                algebraicTransport->fCellsData.fKx[algbindex]  = kappa_phi[0];
//                algebraicTransport->fCellsData.fKx[algbindex]  = kappa_phi[1];
//                algebraicTransport->fCellsData.fKx[algbindex]  = kappa_phi[2];
//                algebraicTransport->fCellsData.fKx[algbindex] = kappa_phi[3] + 0.01;
//        }
//    }
    m_transport_module->fAlgebraicTransport.fHasPropQ=true;
}

void TMRSSFIAnalysis::FillProperties(TPZAlgebraicTransport *algebraicTransport, std::vector<REAL> kappa_phi){
    if (!m_mixed_module|| !m_transport_module) {
        DebugStop();
    }
    
    int ncells = algebraicTransport->fCellsData.fVolume.size();
    for (int icell =0; icell<ncells; icell++) {
        algebraicTransport->fCellsData.fKx[icell]  = kappa_phi[0];
        algebraicTransport->fCellsData.fKy[icell]  = kappa_phi[1];
        algebraicTransport->fCellsData.fKz[icell]  = kappa_phi[2];
//        algebraicTransport->fCellsData.fporosity[icell] = kappa_phi[3] ;
    }
    
   
    m_transport_module->fAlgebraicTransport.fHasPropQ=true;
    
}

// TODOJOSE: Modify name
void TMRSSFIAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    m_mixed_module->SetDataTransfer(sim_data);
    m_transport_module->SetDataTransfer(sim_data);
    m_transport_module->fAlgebraicTransport.fCellsData.SetDataTransfer(sim_data);
    
    // Creates the interface data structure which holds data such as integrated fluxes at
    // interfaces and information of neighborhood
    // Also creats the TCellData which holds information such as the saturation, density,
    // porosity, etc. in each element
    BuildAlgebraicDataStructure();
    
    m_transport_module->fAlgebraicTransport.interfaceid = m_sim_data->mTGeometry.mInterface_material_id;
    m_transport_module->fAlgebraicTransport.inletmatid = 2;
    m_transport_module->fAlgebraicTransport.outletmatid = 3;
    m_transport_module->fAlgebraicTransport.fgravity = m_sim_data->mTNumerics.m_gravity;
    
    //Set initial properties
    m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[0] = m_sim_data->mTFluidProperties.mWaterViscosity;
    m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[1] = m_sim_data->mTFluidProperties.mOilViscosity;
    int ncells = m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil.size();
    REAL rhow = m_sim_data->mTFluidProperties.mWaterDensityRef;
    REAL rhoo = m_sim_data->mTFluidProperties.mOilDensityRef;
    for (int icell =0; icell<ncells; icell++) {
        m_transport_module->fAlgebraicTransport.fCellsData.fDensityWater[icell]= rhow; m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil[icell]= rhoo;
    }
    
    m_transport_module->fAlgebraicTransport.fdt = sim_data->mTNumerics.m_dt;
    for (auto cond:sim_data->mTBoundaryConditions.mBCTransportPhysicalTagTypeValue){
        REAL idVal = std::get<0>(cond);
        REAL idType = std::get<1>(cond);
        REAL idValue = std::get<2>(cond);
        std::pair<int, REAL> bccond = std::make_pair(idType, idValue);
        m_transport_module->fAlgebraicTransport.fboundaryCMatVal[idVal] =bccond;
        
    }
    
    m_transport_module->AnalyzePattern();
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
    
    for (m_k_iteration = 1; m_k_iteration <= n_iterations; m_k_iteration++) {
        
        SFIIteration();
        error_rel_mixed = Norm(m_x_mixed - m_mixed_module->Solution())/Norm(m_mixed_module->Solution());
        m_x_transport = m_transport_module->Solution();
        error_rel_transport = Norm(m_x_transport - m_transport_module->Solution())/Norm(m_transport_module->Solution());
        stop_criterion_Q = error_rel_transport < eps_tol; // Stop by saturation variation
        if (stop_criterion_Q && m_k_iteration >= 1) {
            std::cout << "SFI converged " << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
//            std::cout << "Mixed problem variation = " << error_rel_mixed << std::endl;
//            std::cout << "Transport problem variation = " << error_rel_transport << std::endl;
            m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
            break;
        }
     
       
        
        m_x_mixed = m_mixed_module->Solution();
        m_x_transport = m_transport_module->Solution();
        break;
    }
    
    if (!stop_criterion_Q) {
        std::cout << "SFI fail to converge " << std::endl;
        std::cout << "Number of iterations = " << n_iterations << std::endl;
        std::cout << "SFI will continue with : " << std::endl;
        std::cout << "Mixed problem variation = " << error_rel_mixed << std::endl;
        std::cout << "Transport problem variation = " << error_rel_transport << std::endl;
        m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
        return;
    }
    
    
}

void TMRSSFIAnalysis::PostProcessTimeStep(int val){
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
    
#endif
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
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = t2-t1;
    std::cout << "PostProcess time : " << deltat << std::endl;
#endif
    
}

void TMRSSFIAnalysis::SFIIteration(){
    
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTNumerics.m_ISLinearKrModelQ);
    
//    m_transport_module->fAlgebraicTransport.fCellsData.UpdateMixedDensity();
//    fAlgebraicDataTransfer.TransferLambdaCoefficients();
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    auto deltat = t2-t1;
    std::cout << "Transfer from transport to mixed time: " << deltat << std::endl;
#endif
    
    if(isLinear){
        m_mixed_module->RunTimeStep(); // Newton iterations for mixed problem are done here till convergence
        isLinear = false; // so it leaves after this iteration

        
        // Now, we transfer the mixed problem dofs to the transport problem. These will be used
        // to compute the transport of saturations in an algebraic way
        fAlgebraicDataTransfer.TransferMixedMeshMultiplyingCoefficients();
        //TODOJOSE: wrap all these in a function
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(100); //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracInf);
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracSup);
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(103);  //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(104); //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(2); //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(3); //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(4); //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(5); //TODOJOSE: PUT MATID FROM DATATRANSFER
        m_transport_module->fAlgebraicTransport.VerifyElementFLuxes();
    }
    
   
  
//    fAlgebraicDataTransfer.TransferPressures();
//    m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensities();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t4 = boost::posix_time::microsec_clock::local_time();
    deltat = t4-t3;
    std::cout << "Transfer mixed to transport time: " << deltat << std::endl;
#endif
    
    // Solves the transport problem
    m_transport_module->RunTimeStep();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t5 = boost::posix_time::microsec_clock::local_time();
    deltat = t5-t4;
    std::cout << "Total Transport time: " << deltat << std::endl;
#endif
//    TransferToMixedModule(); // Transfer to mixed
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
    TPZFMatrix<STATE> &elsol = cmesh->ElementSolution();
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
        elsol(el,0) = fast->GetPermTensor()(0,0);
        elsol(el,1) = fast->GetPermTensor()(1,1);
        elsol(el,2) = fast->GetPermTensor()(2,2);
        elsol(el,3) = fast->GetLambda();
    }
}
void TMRSSFIAnalysis::ReadProperties(std::string name, bool print_table_Q, std::vector<REAL> &Kx, std::vector<REAL> &Ky, std::vector<REAL> &Kz, std::vector<REAL> &Phi){
    
    
    std::ifstream file;
    file.open(name);
    int i=1;
    
    
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        if(l != '/'){
            i=i+1;
//            int val = i%15;
//            if(val ==0){
                double a, b, c, d;
                if(iss >> a >> b >> c >> d) ;
                Kx.push_back(a);
                Ky.push_back(b);
                Kz.push_back(c);
                Phi.push_back(d);
//            };
        };
    };
    
    if(Kx.size() == 0){
        std::cout<<"No data read."<<std::endl;
        
        DebugStop();
    }
    if(print_table_Q){
        std::cout<<"*************************"<<std::endl;
        std::cout<<"Reading file... ok!"<<std::endl;
        std::cout<<"*************************"<<std::endl;
    }
    file.close();
}
