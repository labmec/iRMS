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
#include "TPZCompElHDivCollapsed.h"
#include <pzshapequad.h>
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "TPZSimpleTimer.h"
#include "pznonlinanalysis.h"


TMRSSFIAnalysis::TMRSSFIAnalysis(){
    m_sim_data = nullptr;
    m_k_iteration = 0;
    m_mixed_module = nullptr;
    m_transport_module = nullptr;
    m_x_mixed.Resize(0, 0);
    m_x_transport.Resize(0, 0);
    freport_data = new std::ofstream("Report_SFI.txt");
    fcurrentError=1.0;
}

TMRSSFIAnalysis::~TMRSSFIAnalysis(){
    
}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZCompMesh * cmesh_transport, bool must_opt_band_width_Q){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,true);
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    freport_data = new std::ofstream("Report_SFI.txt");
    m_k_iteration=0;
    fcurrentError=1.0;
}

void TMRSSFIAnalysis::BuildAlgebraicDataStructure(){
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    FillProperties();
    freport_data = new std::ofstream("Report_SFI.txt");

}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<REAL(const TPZVec<REAL> & )> & kx, std::function<REAL(const TPZVec<REAL> & )> & ky, std::function<REAL(const TPZVec<REAL> & )> & kz, std::function<REAL(const TPZVec<REAL> & )> & phi, std::function<REAL(const TPZVec<REAL> & )> & s0){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    freport_data = new std::ofstream("Report_SFI.txt");
    bool propsfromPre = m_sim_data->mTReservoirProperties.fPropsFromPreProcess;
    if (propsfromPre==false) {
            fAlgebraicDataTransfer.fkx = kx;
            fAlgebraicDataTransfer.fky = ky;
            fAlgebraicDataTransfer.fkz = kz;
            fAlgebraicDataTransfer.fphi = phi;
            fAlgebraicDataTransfer.fs0 = s0;
    }
    freport_data = new std::ofstream("Report_SFI.txt");
    m_k_iteration=0;
    fcurrentError=1.0;
    
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
//    FillMaterialMemoryDarcy(1);
//    std::string fileprops("Props.txt");
//    std::ofstream file(fileprops);
//    if(file && !propsfromPre){
//        FillProperties(fileprops, &m_transport_module->fAlgebraicTransport);
//    }

}

TMRSSFIAnalysis::TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<std::vector<REAL>(const TPZVec<REAL> & )> & kappa_phi, std::function<REAL(const TPZVec<REAL> & )> & s0){
    m_mixed_module = new TMRSMixedAnalysis(cmesh_mixed,must_opt_band_width_Q);
    m_transport_module = new TMRSTransportAnalysis(cmesh_transport,must_opt_band_width_Q);
    freport_data = new std::ofstream("Report_SFI.txt");
    m_k_iteration=0;
    fcurrentError=1.0;
    
    fAlgebraicDataTransfer.SetMeshes(*cmesh_mixed, *cmesh_transport);
    
    fAlgebraicDataTransfer.fkappa_phi = kappa_phi;
    fAlgebraicDataTransfer.fs0 = s0;
    fAlgebraicDataTransfer.BuildTransportDataStructure(m_transport_module->fAlgebraicTransport);
    m_transport_module->PostProcessProps(0);//Permeabilities
    m_transport_module->PostProcessProps(1);//Porosities
   // FillProperties();
    
//    FillMaterialMemoryDarcy(1);
//    std::string fileprops("Props.txt");
//    FillProperties(fileprops, &m_transport_module->fAlgebraicTransport);
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
    
    
    m_transport_module->fAlgebraicTransport.outletmatid =4;
    m_transport_module->fAlgebraicTransport.inletmatid =3;
    m_transport_module->fAlgebraicTransport.fdt = m_sim_data->mTNumerics.m_dt;
    m_transport_module->fAlgebraicTransport.fgravity = m_sim_data->mTNumerics.m_gravity;
    
    return;
    bool propsfromPre = m_sim_data->mTReservoirProperties.fPropsFromPreProcess;
    if (propsfromPre==false) {
        if (m_sim_data->mTReservoirProperties.kappa_phi) {
            fAlgebraicDataTransfer.fkappa_phi = m_sim_data->mTReservoirProperties.kappa_phi;
            fAlgebraicDataTransfer.fs0=m_sim_data->mTReservoirProperties.s0;
            FillProperties(&m_transport_module->fAlgebraicTransport);
        }
        else{
            std::vector<REAL> kappa_phi(4,1.0);
            
            m_transport_module->fAlgebraicTransport.interfaceid = m_sim_data->mTGeometry.mInterface_material_id;
            m_transport_module->fAlgebraicTransport.fgravity = m_sim_data->mTNumerics.m_gravity;
            
            //Set initial properties
            m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[0] = m_sim_data->mTFluidProperties.mWaterViscosity;
            m_transport_module->fAlgebraicTransport.fCellsData.fViscosity[1] = m_sim_data->mTFluidProperties.mOilViscosity;
            int ncells = m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil.size();
            REAL rhow = m_sim_data->mTFluidProperties.mWaterDensityRef;
            REAL rhoo = m_sim_data->mTFluidProperties.mOilDensityRef;
          
            for (int icell =0; icell<ncells; icell++) {
                m_transport_module->fAlgebraicTransport.fCellsData.fDensityWater[icell]= rhow; m_transport_module->fAlgebraicTransport.fCellsData.fDensityOil[icell]= rhoo;
               
                int matid = m_transport_module->fAlgebraicTransport.fCellsData.fMatId[icell];
                bool fountmat = false;
                for (auto i:m_sim_data->mTReservoirProperties.mPorosityAndVolumeScale) {
                    int id = std::get<0>(i);
                    REAL porosity =std::get<1>(i);
                    REAL volfactor =std::get<2>(i);
                  
             
                    if(id == matid){
                        m_transport_module->fAlgebraicTransport.fCellsData.fporosity[icell] = porosity;
                        m_transport_module->fAlgebraicTransport.fCellsData.fVolumefactor[icell] = volfactor;
                        m_transport_module->fAlgebraicTransport.fCellsData.fVolume[icell] *= volfactor;
                        
                        //volumetric
                        
                        //fracture
                        bool isfrac = m_sim_data->mTFracProperties.m_fracprops.find(matid) != m_sim_data->mTFracProperties.m_fracprops.end();
                        REAL valperm=1.0;
                        if(isfrac){
                            auto fracprop = m_sim_data->mTFracProperties.m_fracprops[matid];
                            valperm = fracprop.m_perm;
                            
                        }
                        else{
                            auto kappa = m_sim_data->mTReservoirProperties.m_permeabilitiesbyId;
                            valperm = kappa[matid];
                        }
                        m_transport_module->fAlgebraicTransport.fCellsData.fKx[icell] = valperm;
                        m_transport_module->fAlgebraicTransport.fCellsData.fKy[icell] = valperm;
//                        m_transport_module->fAlgebraicTransport.fCellsData.fKz[icell] = valperm;
                        std::cout<<"Warning: fkz"<<std::endl;
                        m_transport_module->fAlgebraicTransport.fCellsData.fKz[icell] = 1.0e-6;
                    
                        fountmat =true;
                        break;
                    }
                }
                if (!fountmat){
                    DebugStop();
                }
//
                
            }
            //here
            m_transport_module->fAlgebraicTransport.outletmatid =4;
            m_transport_module->fAlgebraicTransport.fdt = m_sim_data->mTNumerics.m_dt;
            //Type 0.- InletCondition
            //Type 1.- OutletCondition
            // This part of the code tries to set the inlet and outlet matids automatically based on where the flux is
            // entering and leaving.
            bool foundinlet = false;
			for (auto& chunk : m_sim_data->mTBoundaryConditions.mBCTransportMatIdToTypeValue) {
				const int idVal   = chunk.first;
				std::pair<int,REAL>& typeAndVal = chunk.second;
				const int idType = typeAndVal.first;
				const REAL idValue   = typeAndVal.second;
				std::pair<int, REAL> bccond = std::make_pair(idType, idValue);
				m_transport_module->fAlgebraicTransport.fboundaryCMatVal[idVal] =bccond;
                //PressureImposed
				if(idType==0){
//                    if(idValue!=0){
                        if (!foundinlet) {
                            m_transport_module->fAlgebraicTransport.inletmatid =idVal;
                            foundinlet = true;
                        }
                        else{
                            int inletant = m_transport_module->fAlgebraicTransport.inletmatid;
                            REAL val = m_transport_module->fAlgebraicTransport.fboundaryCMatVal[inletant].second;
                            int type = m_sim_data->mTBoundaryConditions.mBCTransportMatIdToTypeValue[inletant].first;
                            if(type==1){
                                m_transport_module->fAlgebraicTransport.outletmatid = idVal;
                            }
                            if((idType==0) && (val>idValue) && (type!=1)){
                                m_transport_module->fAlgebraicTransport.inletmatid = inletant;
                                m_transport_module->fAlgebraicTransport.outletmatid = idVal;
                            }
//                            else{
//                                m_transport_module->fAlgebraicTransport.inletmatid = idVal;
//                                m_transport_module->fAlgebraicTransport.outletmatid = inletant;
//                            }
//                        }
                    }
					
				}
				if(idType==1 && idValue!=0){
                    if(idValue<0){
                        if(!foundinlet){
                            m_transport_module->fAlgebraicTransport.inletmatid=idVal;
                            foundinlet = true;
                        }
                        else{
                            int antval =m_transport_module->fAlgebraicTransport.inletmatid;
                            m_transport_module->fAlgebraicTransport.inletmatid=idVal;
                            m_transport_module->fAlgebraicTransport.outletmatid=antval;
                        }
                        
                    }
                    else{
                        
                        m_transport_module->fAlgebraicTransport.outletmatid =idVal;
                        
                    }
					
				}
			}
            
//            kappa_phi[3]=1.0;
//                FillProperties(&m_transport_module->fAlgebraicTransport, kappa_phi);
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
    m_transport_module->fAlgebraicTransport.outletmatid =4;
    m_transport_module->fAlgebraicTransport.inletmatid =3;

    
  
    
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


void TMRSSFIAnalysis::SetDataTransferAndBuildAlgDatStruct(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    
    m_mixed_module->SetDataTransfer(sim_data);
    m_transport_module->SetDataTransfer(sim_data);
    m_transport_module->fAlgebraicTransport.fCellsData.SetDataTransfer(sim_data);
    
    // Creates the interface data structure which holds data such as integrated fluxes at
    // interfaces and information of neighborhood
    // Also creats the TCellData which holds information such as the saturation, density,
    // porosity, etc. in each element
    
    //BuildAlgebraicDataStructure();
    FillProperties();
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
    bool stop_crit1=false;
    bool stop_crit2=false;
    REAL error_rel_mixed = 1.0;
    
    for (m_k_iteration = 1; m_k_iteration <= n_iterations; m_k_iteration++) {
        
        *(m_mixed_module->fmixed_report_data)<<"   "<<m_k_iteration<<"   ";
        *(m_transport_module->ftransport_report_data)<<"   "<<m_k_iteration<<"   ";
        SFIIteration();
        m_mixed_module->VerifyElementFluxes();
        error_rel_mixed = Norm(m_x_mixed - m_mixed_module->Solution())/Norm(m_mixed_module->Solution());
        
//        m_x_transport = m_transport_module->Solution();
        if(Norm(m_transport_module->Solution())==0){
            fcurrentError =Norm(m_x_transport - m_transport_module->Solution());
        }else{
            fcurrentError = Norm(m_x_transport - m_transport_module->Solution())/Norm(m_transport_module->Solution());
        }
        stop_criterion_Q = fcurrentError < eps_tol; // Stop by saturation variation
        *freport_data<<"    "<< m_k_iteration <<"          "<<fcurrentError<<std::endl;
        
        if (stop_criterion_Q && m_k_iteration > 1) {
            std::cout << "SFI converged " << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            std::cout << "Mixed problem variation = " << error_rel_mixed << std::endl;
            std::cout << "Transport problem variation = " << fcurrentError << std::endl;
            m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
//            m_mixed_module->PostProcessTimeStep();
            break;
        }
        
//        else if (m_k_iteration == n_iterations){
//            std::cout << "SFI not converged " << std::endl;
//            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
//            std::cout << "Mixed problem variation = " << error_rel_mixed << std::endl;
//            std::cout << "Transport problem variation = " << fcurrentError << std::endl;
//            DebugStop();
//        }
        
        m_x_mixed = m_mixed_module->Solution();
        m_x_transport = m_transport_module->Solution();
        
//        m_mixed_module->AllZero(m_mixed_module->Mesh());
//        m_mixed_module->Mesh()->UpdatePreviousState(1.0);
        m_mixed_module->fsoltransfer.TransferFromMultiphysics();
//        m_mixed_module->PostProcessTimeStep();
//        break;
    }
    
    if (!stop_criterion_Q) {
        std::cout << "SFI fail to converge " << std::endl;
        std::cout << "Number of iterations = " << n_iterations << std::endl;
        std::cout << "SFI will continue with : " << std::endl;
        std::cout << "Mixed problem variation = " << error_rel_mixed << std::endl;
        std::cout << "Transport problem variation = " << fcurrentError << std::endl;
        m_transport_module->fAlgebraicTransport.fCellsData.fSaturationLastState = m_transport_module->fAlgebraicTransport.fCellsData.fSaturation;
        return;
    }
    std::cout<<"SFI it: " <<m_k_iteration<<std::endl;
    
}

void TMRSSFIAnalysis::PostProcessTimeStep(const int type, const int dim){

    std::cout << "\n---------------------- TMRSSFIAnalysis Post Process ----------------------" << std::endl;
    TPZSimpleTimer timer_pp("Timer SFIAnalysis Post Process");
    if (type == 0) {
        m_mixed_module->PostProcessTimeStep(dim);
        m_transport_module->PostProcessTimeStep();
    }
    if (type == 1) {
        m_mixed_module->PostProcessTimeStep(dim);
//        m_mixed_module->PostProcessTimeStep(dim-1);
    }
    if (type == 2) {
        m_transport_module->PostProcessTimeStep();
    }

    std::cout << "TMRSSFIAnalysis Post Process total time : " << timer_pp.ReturnTimeDouble()/1000. << " seconds" << std::endl;    
}

void TMRSSFIAnalysis::SFIIteration(){
    

    TPZSimpleTimer timer_sfi("Timer SFI Iteration");
    m_transport_module->fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTNumerics.m_ISLinearKrModelQ);
    
    m_transport_module->fAlgebraicTransport.fCellsData.UpdateMixedDensity();
    fAlgebraicDataTransfer.TransferLambdaCoefficients();
    
    if(isLinearTracer){
        
//        if(m_k_iteration<11){
            m_mixed_module->RunTimeStep(); // Newton iterations for mixed problem are done here till convergence
    //      m_mixed_module->PostProcessTimeStep();
            
            UpdateAllFluxInterfaces();
        
            isLinearTracer = false; // so it leaves after this iteration
//        }
    }
    
//    fAlgebraicDataTransfer.TransferPressures();
//    m_transport_module->fAlgebraicTransport.fCellsData.UpdateDensities();
    
    std::cout << "Running transport problem now..." << std::endl;
    // Solves the transport problem
    m_transport_module->RunTimeStep();
    
    std::cout << "\n ==> Total SFIIteration time: " << timer_sfi.ReturnTimeDouble()/1000 << " seconds" << std::endl;
    
    TransferToMixedModule(); // Transfer to mixed
}

void TMRSSFIAnalysis::UpdateAllFluxInterfaces(){
    // to compute the transport of saturations in an algebraic way
    fAlgebraicDataTransfer.TransferMixedMeshMultiplyingCoefficients();
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_id);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracInf);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracSup);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracFrac);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_sim_data->mTGeometry.mInterface_material_idFracBound);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_transport_module->fAlgebraicTransport.inletmatid);
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(m_transport_module->fAlgebraicTransport.outletmatid);

    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(4); //Mat With No Flux
    m_transport_module->fAlgebraicTransport.UpdateIntegralFlux(5); //Mat With No Flux
    m_transport_module->fAlgebraicTransport.VerifyElementFLuxes();
}

void TMRSSFIAnalysis::VerifyElementFluxes(){
    const REAL tol = 1.e-7;
    TPZMultiphysicsCompMesh *mixedmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_mixed_module->Mesh()) ;
    TPZCompMesh *cmesh =mixedmesh->MeshVector()[0];
//    std::ofstream file("fuxmesh.txt");
    const TPZFMatrix<STATE> &meshSol = cmesh->Solution();
//    mixedmesh->MeshVector()[0]->Print(file);
    int nels =mixedmesh->MeshVector()[0]->NElements();
    for (int iel =0; iel<nels-1; iel++) {
        TPZCompEl *cel = mixedmesh->MeshVector()[0]->Element(iel);
//            if(cel->Dimension() != cmesh->Dimension()) continue;
        TPZCompElHDiv<pzshape::TPZShapeCube> *hdivel = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(cel);
        TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivelq = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(cel);
        TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollaps = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivelq);
        if (!hdivel && !hdivCollaps) {
            continue;
        }
        if (!cel) continue;
        int ncon = cel->NConnects();
        int nCorners = cel->Reference()->NCornerNodes();
        int nsides1 = cel->Reference()->NSides(1);
        int nsides = cel->Reference()->NSides();
        REAL sumel=0.0;
        for (int icon=0; icon<ncon-1; icon++) {
            TPZConnect &con = cel->Connect(icon);
            int sideOrient =0;
            if(!hdivCollaps){
                 sideOrient = hdivel->GetSideOrient(nCorners+nsides1+icon);
            }
            else{
                if (hdivCollaps && icon < 5) {
                    sideOrient = hdivCollaps->GetSideOrient(nCorners+icon);
                }
                 
            }
            if (hdivCollaps && icon == 4) {
                ncon++;
                continue;
            }
            if (hdivCollaps && icon == 5) {
                sideOrient = hdivCollaps->GetSideOrient(8);
            }
            if (hdivCollaps && icon == 6) {
                sideOrient = hdivCollaps->GetSideOrient(9);
            }
            
           
            int sequence =con.fSequenceNumber;
            int64_t pos = cmesh->Block().Position(sequence);
            if(sequence==-1) continue;
            int blocksize = cmesh->Block().Size(sequence);
            for(int ieq=0; ieq< blocksize; ieq++)
            {
            sumel += sideOrient*meshSol.GetVal(cmesh->Block().Index(sequence,ieq),0);
            }
        }
        if(std::abs(sumel)> tol ){
            std::cout << "\n\nERROR! Conservation of element index " << cel->Reference()->Index() << " is " << sumel << std::endl;
            DebugStop();
        }
    }
    std::cout << "\n\n===> Nice! All flux elements satisfy conservation up to tolerance " << tol << std::endl;
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
//    TPZMultiphysicsCompMesh * transport_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(m_transport_module->Mesh());
    TPZCompMesh * transport_cmesh = m_transport_module->Mesh();
    
    if (!mixed_cmesh || !transport_cmesh) {
        DebugStop();
    }
    
   // //    m_transport_module->m_soltransportTransfer.TransferFromMultiphysics();
   // transport_cmesh->LoadSolutionFromMultiPhysics(); //Load Solution???
  
    
    // Saturations are transferred to mixed module
    int s_b = 2;
//    TPZFMatrix<STATE> & s_dof = transport_cmesh->MeshVector()[s_b]->Solution();
//    TPZFMatrix<STATE> & s_dof = transport_cmesh->Solution();
    
//    mixed_cmesh->MeshVector()[s_b]->LoadSolution(s_dof);
    
    
    //     m_mixed_module->m_soltransportTransfer.TransferToMultiphysics();
//    mixed_cmesh->LoadSolutionFromMeshes();
    
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
