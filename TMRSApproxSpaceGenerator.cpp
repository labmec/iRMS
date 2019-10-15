//
//  TMRSApproxSpaceGenerator.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#include "TMRSApproxSpaceGenerator.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

using namespace std;

TMRSApproxSpaceGenerator::TMRSApproxSpaceGenerator(){
    mGeometry = nullptr;
    mMixedOperator = nullptr;
    mTransportOperator = nullptr;
}

TMRSApproxSpaceGenerator & TMRSApproxSpaceGenerator::operator=(const TMRSApproxSpaceGenerator &other){
    DebugStop();
}

TMRSApproxSpaceGenerator::~TMRSApproxSpaceGenerator(){
    
}

void TMRSApproxSpaceGenerator::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TMRSApproxSpaceGenerator::Read(TPZStream &buf, void *context){
    DebugStop();
}

int TMRSApproxSpaceGenerator::ClassId() const{
    DebugStop();
}

void TMRSApproxSpaceGenerator::SetGeometry(TPZGeoMesh * geometry){
    mGeometry = geometry;
}

TPZGeoMesh * TMRSApproxSpaceGenerator::GetGeometry(){
    return mGeometry;
}


void TMRSApproxSpaceGenerator::LoadGeometry(std::string geometry_file){
    
    TPZGmshReader Geometry;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    mGeometry = Geometry.GeometricGmshMesh(geometry_file);
    Geometry.PrintPartitionSummary(std::cout);
#ifdef PZDEBUG
    if (!mGeometry)
    {
        std::cout << "The geometrical mesh was not generated." << std::endl;
        DebugStop();
    }
#endif
}

void TMRSApproxSpaceGenerator::PrintGeometry(std::string name)
{
    if (!mGeometry) {
        DebugStop();
    }
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name  << name << "_geometry" << ".txt";
    vtk_name   << name << "_geometry"  << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    mGeometry->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, vtkfile, true);
    
}


TPZCompMesh * TMRSApproxSpaceGenerator::HdivFluxCmesh(int order){
    
    if (!mGeometry) {
        DebugStop();
    }

    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZMixedDarcyFlow * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(order);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mDataTransfer.mTGeometry.mDomainDimNameAndPhysicalTag;
    REAL kappa = 1.0;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TPZMixedDarcyFlow(materia_id,d);
            volume->SetPermeability(kappa);
            cmesh->InsertMaterialObject(volume);
        }
    }

    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mDataTransfer.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        cmesh->InsertMaterialObject(face);
    }
    
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "q_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
}

TPZCompMesh * TMRSApproxSpaceGenerator::DiscontinuousCmesh(int order){
 
    if (!mGeometry) {
        DebugStop();
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZL2Projection * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(order);
    
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mDataTransfer.mTGeometry.mDomainDimNameAndPhysicalTag;
    int nstate = 1;
    TPZVec<STATE> sol;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TPZL2Projection(materia_id,d,nstate, sol);
            cmesh->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    if (order == 0) {
        TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
        TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mDataTransfer.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
        for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
            int bc_id   = get<0>(chunk);
            int bc_type = get<1>(chunk);
            val2(0,0)   = get<2>(chunk);
            TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
            cmesh->InsertMaterialObject(face);
        }
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }else{
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    

    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    if (order == 0) {
        file_name << "s_cmesh" << ".txt";
    }
    else{
        file_name << "p_cmesh" << ".txt";
    }
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
}


void TMRSApproxSpaceGenerator::BuildMixedMultiPhysicsCompMesh(int order){
    
    int dimension = mGeometry->Dimension();
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
    mMixedOperator->SetDefaultOrder(order);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mDataTransfer.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,d);
            mMixedOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mDataTransfer.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mMixedOperator->InsertMaterialObject(face);
    }
    
    
    TPZManVector<TPZCompMesh * ,2> mesh_vec(2);

    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order);
    TPZManVector<int,5> active_approx_spaces(2);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    mMixedOperator->SetDimModel(dimension);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec);
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "mixed_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    mMixedOperator->Print(sout);
#endif
    
}

void TMRSApproxSpaceGenerator::BuildTransportMultiPhysicsCompMesh(){
    
    if (!mMixedOperator || !mGeometry) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,2> mixed_meshvec = mMixedOperator->MeshVector();
    TPZManVector<TPZCompMesh *,3> transport_meshvec(3);
    
    transport_meshvec[0] = mixed_meshvec[0];
    transport_meshvec[1] = mixed_meshvec[1];
    transport_meshvec[2] = DiscontinuousCmesh();


    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
    mTransportOperator->SetDefaultOrder(0);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mDataTransfer.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSMultiphaseFlow<TMRSMemory>(materia_id,d);
            mTransportOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mDataTransfer.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mTransportOperator->InsertMaterialObject(face);
    }
    
    int transport_matid = 100;
    {
        REAL phi = 0.1;
        
        TPZTracerFlow * interface = new TPZTracerFlow(transport_matid,dimension-1);
        interface->SetPorosity(phi);
        mTransportOperator->InsertMaterialObject(interface);
    }
    
    mTransportOperator->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 1;
    mTransportOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,transport_meshvec);
    
#ifdef PZDEBUG
    std::ofstream transport_a("transport_cmesh_after.txt");
    mTransportOperator->Print(transport_a);
#endif
    
    {
        mTransportOperator->Reference()->ResetReference();
        mTransportOperator->LoadReferences();
        
        TPZManVector<std::vector<int64_t>,4> cel_indexes(4);
        
        TPZManVector<int64_t,3> left_mesh_indexes(2,0);
        left_mesh_indexes[0] = 0;
        left_mesh_indexes[1] = 2;
        TPZManVector<int64_t,3> right_mesh_indexes(1,0);
        right_mesh_indexes[0] = 2;
        
        int64_t nel = mTransportOperator->NElements();
        for (int64_t el = 0; el < nel; el++) {
            
            TPZCompEl *cel = mTransportOperator->Element(el);
            if(!cel) DebugStop();
            TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!celmp) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            
            int gel_dim = gel->Dimension();
            cel_indexes[gel_dim].push_back(el);
            
        }
        
        for (auto cel_index: cel_indexes[dimension]) { // Higher dimension case
            TPZCompEl *cel = mTransportOperator->Element(cel_index);
            TPZMultiphysicsElement * celmult = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!celmult) {
                DebugStop();
            }
            
            if (!cel){continue;};
            TPZGeoEl *gel = cel->Reference();
            if (!gel){continue;};
            int nsides = gel->NSides();
            
            for (int iside = gel->NNodes(); iside < nsides; iside++) {
                
                TPZGeoElSide gelside(gel,iside);
                TPZCompElSide celside_l(cel,iside);
                TPZGeoElSide neig = gelside.Neighbour();
                TPZGeoEl *neihel = neig.Element();
                TPZCompElSide celside_r = neig.Reference();
                if ((neihel->Dimension() == gel->Dimension()) && (gel->Id() < neihel->Id()) ) {
                    TPZGeoElBC gbc(gelside,transport_matid);
                    
                    int64_t index;
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                }
                if ((neihel->Dimension() == dimension - 1)) { // BC cases
                    
                    TPZGeoElBC gbc(gelside,neihel->MaterialId());
                    
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                    
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    
                }
                
            }
            
        }
    }
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    mTransportOperator->Print(transport);
#endif
}

void TMRSApproxSpaceGenerator::SetDataTransfer(TMRSDataTransfer & DataTransfer){
    mDataTransfer = DataTransfer;
}

TMRSDataTransfer & TMRSApproxSpaceGenerator::GetDataTransfer(){
    return mDataTransfer;
}

TPZMultiphysicsCompMesh * TMRSApproxSpaceGenerator::GetMixedOperator(){
    return mMixedOperator;
}

TPZMultiphysicsCompMesh * TMRSApproxSpaceGenerator::GetTransportOperator(){
    return mTransportOperator;
}


void TMRSApproxSpaceGenerator::AdjustMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator){
    
    if (!MixedOperator || !TransportOperator) {
        DebugStop();
    }
    
    /// Adjust integration rule
    /// o Stands for reservoir
    /// d Stands for pressure
    
    TPZCompMesh * cmesh_res = MixedOperator;
    TPZCompMesh * cmesh_tra = TransportOperator;
    
    cmesh_tra->LoadReferences();
    int nel_res = cmesh_res->NElements();
    int gmesh_dim = cmesh_tra->Reference()->Dimension();
    
    // Scanning structure
    std::vector<std::pair<int64_t, int64_t>> cel_pairs;
    for (long el = 0; el < nel_res; el++) {
        TPZCompEl *cel_res = cmesh_res->Element(el);
        if (!cel_res) {
            continue;
        }
        
        TPZGeoEl * gel = cel_res->Reference();
        if (!gel) {
            continue;
        }
        
        if (gel->Dimension() != gmesh_dim || gel->HasSubElement()) {
            continue;
        }
        
        /// Finding the other computational element
        TPZCompEl * cel_tra = gel->Reference();
        if (!cel_tra) {
            continue;
        }
        
        int64_t cel_res_index = cel_res->Index();
        int64_t cel_tra_index = cel_tra->Index();
        cel_pairs.push_back(std::make_pair(cel_res_index, cel_tra_index));
        cel_tra->SetFreeIntPtIndices();  // operation involving resize. It is not thread safe.
    }
    
    int nel = cel_pairs.size();
#ifdef USING_TBB
    tbb::parallel_for(size_t(0), size_t(nel), size_t(1), [&cel_pairs,&cmesh_res,&cmesh_tra] (size_t & i)
      {
          int64_t cel_res_index = cel_pairs[i].first;
          int64_t cel_geo_index = cel_pairs[i].second;
          TPZCompEl *cel_res = cmesh_res->Element(cel_res_index);
          TPZCompEl * cel_tra = cmesh_tra->Element(cel_geo_index);
          
          const TPZIntPoints & rule = cel_res->GetIntegrationRule();
          TPZIntPoints * cloned_rule = rule.Clone();
          TPZManVector<int64_t,20> indices;
          cel_res->GetMemoryIndices(indices);
          cel_tra->SetMemoryIndices(indices);
          cel_tra->SetIntegrationRule(cloned_rule);
      }
);
    
#else
    for (long i = 0; i < nel; i++) {
        
        int64_t cel_res_index = cel_pairs[i].first;
        int64_t cel_tra_index = cel_pairs[i].second;
        TPZCompEl *cel_res = cmesh_res->Element(cel_res_index);
        TPZCompEl * cel_tra = cmesh_tra->Element(cel_tra_index);
        
        const TPZIntPoints & rule = cel_res->GetIntegrationRule();
        TPZIntPoints * cloned_rule = rule.Clone();
        TPZManVector<int64_t,20> indices;
        cel_res->GetMemoryIndices(indices);
        cel_tra->SetMemoryIndices(indices);
        cel_tra->SetIntegrationRule(cloned_rule);
    }
#endif
    
#ifdef PZDEBUG
    std::ofstream out_res("Cmesh_res_adjusted.txt");
    cmesh_res->Print(out_res);
#endif
    
#ifdef PZDEBUG
    std::ofstream out_geo("Cmesh_tra_adjusted.txt");
    cmesh_tra->Print(out_geo);
#endif
    
}

void TMRSApproxSpaceGenerator::UnifyMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator) {
    
    if (!MixedOperator || !TransportOperator) {
        DebugStop();
    }
    
    TPZCompMesh * cmesh_o = MixedOperator;
    TPZCompMesh * cmesh_d = TransportOperator;
    
    /// Apply memory link
    TPZMaterial * material_o = cmesh_o->FindMaterial(material_id);
    TPZMaterial * material_d = cmesh_d->FindMaterial(material_id);
    if (!material_o || !material_d) {
        DebugStop();
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory_o = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material_o);
    TPZMatWithMem<TMRSMemory,TPZDiscontinuousGalerkin> * mat_with_memory_d = dynamic_cast<TPZMatWithMem<TMRSMemory,TPZDiscontinuousGalerkin> * >(material_d);
    if (!mat_with_memory_o || !mat_with_memory_d) {
        DebugStop();
    }
    
    mat_with_memory_d->SetMemory(mat_with_memory_o->GetMemory());
}

void TMRSApproxSpaceGenerator::FillMaterialMemory(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator){
    
    if (!MixedOperator || !TransportOperator) {
        DebugStop();
    }
    
    TPZCompMesh * cmesh = MixedOperator;
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    if (!material) {
        DebugStop();
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
    if (!mat_with_memory) {
        DebugStop();
    }
    
    
    // Set initial porosity, permeability, saturations, etc ...
    {
        std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        int ndata = memory_vector->NElements();
        
#ifdef USING_TBB
        tbb::parallel_for(size_t(0), size_t(ndata), size_t(1), [&memory_vector] (size_t & i) {
            TMRSMemory &mem = memory_vector.get()->operator [](i);
            
            mem.m_sw = 0.0;
            mem.m_phi = 0.1;
            
            REAL kappa = 1.0;
            mem.m_kappa.Resize(3, 3);
            mem.m_kappa.Zero();
            mem.m_kappa_inv.Resize(3, 3);
            mem.m_kappa_inv.Zero();
            for (int i = 0; i < 3; i++) {
                mem.m_kappa(i,i) = kappa;
                mem.m_kappa_inv(i,i) = 1.0/kappa;
            }
            
        }
);
        
#else
        for (int i = 0; i < ndata; i++) {
            TMRSMemory &mem = memory_vector.get()->operator [](i);

            
            mem.m_sw = 0.0;
            mem.m_phi = 0.1;
            
            REAL kappa = 1.0;
            mem.m_kappa.Resize(3, 3);
            mem.m_kappa.Zero();
            mem.m_kappa_inv.Resize(3, 3);
            mem.m_kappa_inv.Zero();
            for (int i = 0; i < 3; i++) {
                mem.m_kappa(i,i) = kappa;
                mem.m_kappa_inv(i,i) = 1.0/kappa;
            }

        }
#endif
        
        
    }
    
    
}

void TMRSApproxSpaceGenerator::SetUpdateMaterialMemory(int material_id, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q){
    
    if (!cmesh) {
        DebugStop();
    }
    
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    if (!material) {
        DebugStop();
    }
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
    if (mat_with_memory) {
        mat_with_memory->SetUpdateMem(update_memory_Q);
        return;
    }
    
    TPZMatWithMem<TMRSMemory,TPZDiscontinuousGalerkin> * mat_with_memory_trans = dynamic_cast<TPZMatWithMem<TMRSMemory,TPZDiscontinuousGalerkin> * >(material);
    if (mat_with_memory_trans) {
        mat_with_memory_trans->SetUpdateMem(update_memory_Q);
        return;
    }

}

void TMRSApproxSpaceGenerator::SetUpdateMemory(int dimension, TMRSDataTransfer & sim_data, TPZMultiphysicsCompMesh * cmesh, bool update_memory_Q){
    for (auto item : sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[dimension]) {
        int material_id = item.second;
        SetUpdateMaterialMemory(material_id, cmesh, update_memory_Q);
    }
}


