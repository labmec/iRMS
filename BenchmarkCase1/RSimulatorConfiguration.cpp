//
//  RSimulatorConfiguration.cpp
//  MonophasicTest
//
//  Created by Jose on 27/7/19.
//

#include "RSimulatorConfiguration.h"
#include "MMeshType.h"
#include "TPZGenGrid2D.h"

void Ladoderecho_2D (const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

RSimulatorConfiguration::RSimulatorConfiguration(){
    
}

RSimulatorConfiguration::RSimulatorConfiguration(SimulationCase sim_case){
    fsim_case = sim_case;
}
RSimulatorConfiguration::RSimulatorConfiguration(TPZGeoMesh *gmesh){
    fGmesh = gmesh;
}
void RSimulatorConfiguration::CreateGeomesh(int nx, int ny, double l, double h, MMeshType eltype){
   
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;
    nels[1]=ny;
    
    
    TPZVec<REAL> x0(3,0.0);
    x0[0]=0.1;
    TPZVec<REAL> x1(3,l+0.1);
    x1[1]=h;
    x1[2]=0;
    
    //Setting boundary conditions (negative numbers to recognize them)
    TPZGenGrid2D gen(nels,x0,x1);
    switch(eltype)
    {
        
        case MMeshType::ETriangular:
            gen.SetElementType(MMeshType::ETriangular);
            break;
        case MMeshType::EQuadrilateral:
            gen.SetElementType(MMeshType::EQuadrilateral);
            break;
        default:
            DebugStop();
    }
    
    
    gmesh->SetDimension(2);
  
    gen.Read(gmesh);
    gen.SetRefpatternElements(true);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -2);
    gen.SetBC(gmesh, 6, -3);
    gen.SetBC(gmesh, 7, -4);
    
    gmesh->BuildConnectivity();
    
    fGmesh = gmesh;
    
}

TPZGeoMesh *  RSimulatorConfiguration::GetGeomesh(){
    if (!fGmesh) {
        std::cout<<"a malha geometrica nao esta creada"<<std::endl;
        DebugStop();
    }
    return fGmesh;
}

TPZCompMesh * RSimulatorConfiguration::CreateFluxCmesh(TPZGeoMesh * gmesh, int order){
    
    int dimension = gmesh->Dimension();
    int nvols = fsim_case.omega_ids.size();
    int nbound = fsim_case.gamma_ids.size();
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    cmesh->SetDefaultOrder(order);
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],fsim_case.omega_dim[ivol]);
        volume->SetPermeability(fsim_case.permeabilities[ivol]);
        
  //  TPZMatPoisson3d * volume = new TPZMatPoisson3d(1,2);
        TPZAutoPointer<TPZFunction<STATE> > sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_2D, 10);
        volume->SetForcingFunction(sourceterm);
        
        
        cmesh->InsertMaterialObject(volume);
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {

                val2(0,0)=fsim_case.vals[ibound];
                int condType=fsim_case.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,fsim_case.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
   
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
//    for (int dim=dimension; dim>0; dim--) {
//        cmesh->Reference()->ResetReference();
//        cmesh->SetDimModel(dim);
//        cmesh->SetDefaultOrder(order);
//        std::set<int> matids;
//        for (int ivol=0; ivol<nvols; ivol++) {
//            if(fsim_case.omega_dim[ivol]==dim) matids.insert(fsim_case.omega_ids[ivol]);
//        }
//        for (int ibound=0; ibound<nbound; ibound++) {
//            if(fsim_case.gamma_dim[ibound]==dim) matids.insert(fsim_case.gamma_ids[ibound]);
//        }
//
//        cmesh->SetAllCreateFunctionsHDiv();
//        cmesh->AutoBuild(matids);
//    }

    cmesh->InitializeBlock();
        
    return cmesh;
    
}
TPZCompMesh * RSimulatorConfiguration::CreatePressureCmesh(TPZGeoMesh * gmesh, int order){
    
    int dimension = gmesh->Dimension();
    int nvols = fsim_case.omega_ids.size();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    std::set<int> matids;
    for (int ivol=0; ivol < nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],dimension);
        volume->SetPermeability(fsim_case.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        if (fsim_case.omega_dim[ivol] == dimension) {
            matids.insert(fsim_case.omega_ids[ivol]);
        }
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh->AutoBuild(matids);
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    
    return cmesh;
    
}
TPZCompMesh *RSimulatorConfiguration::CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh, int ref)
{
    TPZCompMesh *q_cmesh = cmesh->MeshVector()[0];
    TPZGeoMesh * geometry = q_cmesh->Reference();
    if (!geometry) {
        DebugStop();
    }
    // UniformRefinement(geometry, ref);
    
    TPZCompMesh *s_cmesh = new TPZCompMesh(geometry);
    
    geometry->ResetReference();
    q_cmesh->LoadReferences();
    
    std::set<int> allmat_ids = {1};
    
//    for(auto frac_set: fractures)
//    {
//        for(auto matid:frac_set.m_id)
//        {
//            allmat_ids.insert(matid);
//        }
//    }
    
    std::set<int> bcmat_ids = {-1,-2,-3,-4};
    
    int nstate = 1;
    TPZVec<STATE> sol(1,0.0);
    
    /// Inserting the materials
    for (auto mat_id: allmat_ids) {
        TPZL2Projection * volume = new TPZL2Projection(mat_id,2,nstate, sol);
        s_cmesh->InsertMaterialObject(volume);
    }
    
    TPZMaterial * material = s_cmesh->FindMaterial(1);
    if (!material) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0);
    TPZFMatrix<STATE> val2(1,1,1);
    
    int typ = 0; // inlet
    /// Inserting the materials
    for (auto mat_id: bcmat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ, val1, val2);
        s_cmesh->InsertMaterialObject(bc);
    }
    geometry->ResetReference();
    
    s_cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int s_order = 0;
    s_cmesh->SetDefaultOrder(s_order);
    TPZStack<TPZGeoEl *> Flux_Elements;
    
    for (auto cel : q_cmesh->ElementVec()) {
        
        if (!cel) {
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        Flux_Elements.push_back(gel);
        //TPZVec<TPZGeoEl *> sons;
        
        //gel->Divide(sons);
        //for (auto geld: sons){
        int gel_index = gel->Index();
        int mat_id = gel->MaterialId();
        if (allmat_ids.find(mat_id) != allmat_ids.end()) {
            TPZGeoEl * geld = geometry->Element(gel_index);
            //CreateTransportElement(s_order,s_cmesh, geld, false);
        }
        if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
            TPZGeoEl * geld = geometry->Element(gel_index);
            // CreateTransportElement(s_order,s_cmesh, geld, true);
        }
        
        
        //}
    }
    
    
    //ref
    TPZManVector<TPZGeoEl*> sons;
    int h_level =0;
    
    for (int level=0; level<h_level; level++) {
        int nel = Flux_Elements.size();
        for (int iel = 0; iel<nel; iel++) {
            TPZGeoEl *gel = Flux_Elements[iel];
            if(!gel || gel->HasSubElement()){continue;}
            gel->Divide(sons);
            int nsons = sons.size();
            for (int i=0; i<nsons; i++) {
                Flux_Elements.push_back(sons[i]);
            }
            
        }
        
    }
    
    int nel = Flux_Elements.NElements();
    for (int iel =0; iel<nel; iel++) {
        TPZGeoEl *gel = Flux_Elements[iel];
        if (gel->HasSubElement()) {
            continue;
        }
        int gel_index = gel->Index();
        int mat_id = gel->MaterialId();
        if (allmat_ids.find(mat_id) != allmat_ids.end()) {
            TPZGeoEl * geld = geometry->Element(gel_index);
            CreateTransportElement(s_order,s_cmesh, geld, false);
        }
        if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
            TPZGeoEl * geld = geometry->Element(gel_index);
            CreateTransportElement(s_order,s_cmesh, geld, true);
        }
        
    }
    
    
    //
    //    for(int i=0; i < h_level; i++)
    //    {
    //        int64_t nels = s_cmesh->NElements();
    //        //TPZGeoMesh *geo_c = s_cmesh->Reference();
    //        for(int64_t elem = 0; elem < nels; elem++)
    //        {
    //            TPZCompEl * cel = s_cmesh->Element(elem);
    //            TPZGeoEl * gel = cel->Reference();
    //            if(!gel || gel->HasSubElement()){
    //                continue;
    //
    //            }
    //            gel->Divide(sons);
    //            for (auto geld: sons){
    //            int gel_index = geld->Index();
    //            int mat_id = gel->MaterialId();
    //            if (allmat_ids.find(mat_id) != allmat_ids.end()) {
    //                TPZGeoEl * geld2 = geometry->Element(gel_index);
    //                CreateTransportElement(s_order,s_cmesh, geld2, false);
    //            }
    //            if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
    //                TPZGeoEl * geld2 = geometry->Element(gel_index);
    //                CreateTransportElement(s_order,s_cmesh, geld2, true);
    //            }
    //        }
    //        }
    //    }
    
    geometry->ResetConnectivities();
    geometry->BuildConnectivity();
    //ref
    
    
    s_cmesh->SetDimModel(0);
    s_cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int ind_0 = -1;
//    for (int i = 0; i<fractures.size(); i++) {
//        if (fractures[i].m_dim == 0) {
//            ind_0 = i;
//            break;
//        }
//    }
////    if(ind_0 >= 0)
////    {
////        /// Create point element
////        for (auto cel : cmesh->ElementVec()) {
////            if (!cel) {
////                continue;
////            }
////
////            TPZGeoEl * gel = cel->Reference();
////            if (!gel) {
////                DebugStop();
////            }
////
////            if (gel->Dimension() != 0) {
////                continue;
////            }
////            int matid = gel->MaterialId();
////            if(bcmat_ids.find(matid) != bcmat_ids.end())
////            {
////                continue;
////            }
////            if(allmat_ids.find(matid) == allmat_ids.end())
////            {
////                continue;
////            }
////
////            TPZMultiphysicsInterfaceElement * int_face = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
////            if (int_face) {
////                continue;
////            }
////
////            int n_connect = cel->NConnects();
////            bool matid_found = fractures[ind_0].m_id.find(gel->MaterialId())  != fractures[ind_0].m_id.end();
////            if (n_connect !=1 && matid_found) {
////                DebugStop();
////            }
////            else if(n_connect == 1 && !matid_found)
////            {
////                std::cout << "I dont understand matid should be included in the fracture data structure\n";
////            }
////            else if (n_connect == 0)
////            {
////                DebugStop();
////            }
////            TPZConnect & c = cel->Connect(0);
////
////            if (c.NElConnected() > 2 && matid_found) {
////                CreateTransportElement(s_order,s_cmesh, gel, false);
////            }
////
////
////        }
//    }
    geometry->ResetReference();
    s_cmesh->InitializeBlock();
    
    std::ofstream fileprint("malla_test.txt");
    s_cmesh->Print(fileprint);
    
    return s_cmesh;
}
void RSimulatorConfiguration::CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC){
    int64_t cel_index;
    int dimension = gel->Dimension();
    cmesh->SetDimModel(dimension);
    if (is_BC) {
        cmesh->SetDimModel(dimension+1);
    }
    TPZCompEl * cel = cmesh->ApproxSpace().CreateCompEl(gel, *cmesh, cel_index);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
    TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
    if (intel){
        intel->PRefine(p_order);
    } else if (intelDisc) {
        intelDisc->SetDegree(p_order);
        intelDisc->SetTrueUseQsiEta();
    } else {
        DebugStop();
    }
    gel->ResetReference();
}

TPZMultiphysicsCompMesh * RSimulatorConfiguration::MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZManVector<TPZCompMesh *> meshvec){
    
    TPZGeoMesh *geometry = mixed->Reference();
    int dimension = geometry->Dimension();
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);
    
    /// Inserting matrix materials

    int n_vols = fsim_case.omega_ids.size();
    for (int i = 0; i < n_vols; i++) {
        int mat_id = fsim_case.omega_ids[i];
        REAL phi = fsim_case.porosities[i];
        TPZTracerFlow * volume = new TPZTracerFlow(mat_id,dimension);
        volume->SetPorosity(phi);
        cmesh->InsertMaterialObject(volume);
    }
    
    int transport_matid = 10;
    {
        REAL phi = 0.1;

        TPZTracerFlow * interface = new TPZTracerFlow(transport_matid,dimension-1);
        interface->SetPorosity(phi);
        cmesh->InsertMaterialObject(interface);
    }
    /// Inserting fracture materials
//    int n_fracs = fracture_data.size();
//    for (int i = 0; i < n_fracs; i++) {
//        for(auto mat_id :fracture_data[i].m_id)
//        {
//            REAL phi = fracture_data[i].m_porosity;
//            REAL d_opening = fracture_data[i].m_d_opening;
//            TPZTracerFlow * volume = new TPZTracerFlow(mat_id,0);
//            volume->SetPorosity(phi);
//            volume->SetFractureCrossLength(d_opening);
//            cmesh->InsertMaterialObject(volume);
//        }
//    }
    
    std::set<int> bc_inlet_mat_ids = {-1,-2,-3};
    std::set<int> bc_outlet_mat_ids = {-4};
    
    TPZMaterial * material = cmesh->FindMaterial(1);
    if (!material) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
    
    int typ_inlet = 0; // inlet
    /// Inserting the materials
      val2(0,0) = fsim_case.c_inlet;
    for (auto mat_id: bc_inlet_mat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ_inlet, val1, val2);
         cmesh->InsertMaterialObject(bc);
     
    }
 
    int typ_outlet = 1; // outlet
    /// Inserting the materials
    val2(0,0) = 0.0;
    for (auto mat_id: bc_outlet_mat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ_outlet, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->SetDimModel(dimension);
    
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,meshvec);
    
    {
        cmesh->Reference()->ResetReference();
        cmesh->LoadReferences();

        TPZManVector<std::vector<int64_t>,4> cel_indexes(4);
        
        TPZManVector<int64_t,3> left_mesh_indexes(2,0);
        left_mesh_indexes[0] = 0;
        left_mesh_indexes[1] = 2;
        TPZManVector<int64_t,3> right_mesh_indexes(1,0);
        right_mesh_indexes[0] = 2;
        
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            
            TPZCompEl *cel = cmesh->Element(el);
            if(!cel) DebugStop();
            TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!celmp) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
            
            int gel_dim = gel->Dimension();
            cel_indexes[gel_dim].push_back(el);
            
        }
        
        for (auto cel_index: cel_indexes[2]) { // 2D case
            TPZCompEl *cel = cmesh->Element(cel_index);
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
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside_l,celside_r);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                }
                if ((neihel->Dimension() == dimension - 1)) { // BC cases
                    
                    TPZGeoElBC gbc(gelside,neihel->MaterialId());
                    
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside_l,celside_r);
                    
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    
                }
                
            }
            
        }
    }
    
    std::cout << "Created multi-physics transport mesh\n";
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    cmesh->Print(transport);
#endif
    
    return cmesh;
}

TPZMultiphysicsCompMesh *RSimulatorConfiguration::CreateMultiPhysicsCompMesh(TPZGeoMesh *geometry){
    
    int dimension = geometry->Dimension();
    int nvols = fsim_case.omega_ids.size();
    int nbound= fsim_case.gamma_ids.size();
    if (nvols<1) {
        std::cout<<"Error: Omega is not defined."<<std::endl;
        DebugStop();
    }
    if (nbound<1) {
        std::cout<<"Error: Gamma is not defined."<<std::endl;
        DebugStop();
    }
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);
    TPZFNMatrix<9,STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(fsim_case.omega_ids[ivol],dimension);
        volume->SetPermeability(fsim_case.permeabilities[ivol]);
 //       TPZMixedPoisson * volume = new TPZMixedPoisson(1,2);
        
        TPZAutoPointer<TPZFunction<STATE> > sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_2D, 10);
     volume->SetForcingFunction(sourceterm);
        
        cmesh->InsertMaterialObject(volume);
        
   
        
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=fsim_case.vals[ibound];
                int condType=fsim_case.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,fsim_case.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    TPZManVector<TPZCompMesh * ,2> mesh_vec(2);
    
    int order_q= fsim_case.order_q;
    int order_p = fsim_case.order_p;
    mesh_vec[0] = CreateFluxCmesh(geometry,order_q);
    mesh_vec[1] = CreatePressureCmesh(geometry, order_p);
    TPZManVector<int,5> active_approx_spaces(2); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    
    
    
    if (fsim_case.IsMHMQ) {
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    //    else{
    //        TPZCompMeshTools::GroupElements(cmesh);
    //        std::cout << "Created grouped elements\n";
    //        bool keepmatrix = false;
    //        bool keeponelagrangian = true;
    //        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    //        std::cout << "Created condensed elements\n";
    //        cmesh->CleanUpUnconnectedNodes();
    //        cmesh->ExpandSolution();
    //    }
    
    std::cout << "Created multi-physics DFN mesh\n";
    
    //meshvec = mesh_vec;
    return cmesh;
    
}

//TPZAnalysis * RSimulatorConfiguration::CreateAnalysis(){
//
//}

TPZAnalysis * RSimulatorConfiguration::CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,  bool must_opt_band_width_Q, int number_threads, bool UsePardiso_Q){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh_mult, true);
    
    if (fsim_case.UsePardisoQ) {
        
        TPZSymetricSpStructMatrix matrix(cmesh_mult);
      
        //        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(fsim_case.n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    
//    if (fsim_case.UseFrontalQ) {
//        
//        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh_mult);
//        matrix.SetDecomposeType(ELDLt);
//        matrix.SetNumThreads(fsim_case.n_threads);
//     
//        analysis->SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELDLt);
//        analysis->SetSolver(step);
//        
//        return analysis;
//    }
    else{
        
        TPZSkylineStructMatrix matrix(cmesh_mult);
        matrix.SetNumThreads(fsim_case.n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    
    return analysis;
    
}


void RSimulatorConfiguration::PrintCmesh(int mesh_index, std::ofstream &file_name){
   
}
void RSimulatorConfiguration::PrintGmesh( std::ofstream &file_name){
    TPZVTKGeoMesh::PrintGMeshVTK(fGmesh, file_name);
}

void RSimulatorConfiguration::PosProcess(){
    
}

//#define Verbose_Q

void RSimulatorConfiguration::Run(){
    CreateGeomesh(10, 1, 2, 1, MMeshType::EQuadrilateral);
    fsim_case.order_p=1;
    fsim_case.order_q=1;
    TPZMultiphysicsCompMesh *c_mult = CreateMultiPhysicsCompMesh(fGmesh);
    
    TPZCompMesh *s_cmesh = CreateTransportMesh(c_mult, 0);
    TPZManVector<TPZCompMesh *> meshvect = c_mult->MeshVector();
    TPZManVector<TPZCompMesh *> meshvec(3);
    meshvec[0]=meshvect[0];
    meshvec[1]=meshvect[1];
    meshvec[2]=s_cmesh;
    TPZMultiphysicsCompMesh *mul_sat = MPTransportMesh(c_mult, meshvec);
    
#ifdef Verbose_Q
    {
        TPZGeoMesh * gmsh = mul_sat->Reference();
        std::string vtk_name = "geometry_sat.vtk";
        std::ofstream vtkfile(vtk_name.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmsh, vtkfile, true);
    }
    
    std::ofstream fileprint("test_sat.txt");
    mul_sat->Print(fileprint);
    std::ofstream file("sat.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(mul_sat, file);
#endif
    
    TPZAnalysis *an = CreateAnalysis(c_mult, false, 1, true);
    
#ifdef Verbose_Q
    std::ofstream fileg("geometry.txt");
    fGmesh->Print(fileg);
    std::ofstream filep("pressure.txt");
    c_mult->MeshVector()[1]->Print(filep);
    std::ofstream flux("flux.txt");
    c_mult->MeshVector()[0]->Print(flux);
    std::ofstream multi("multi.txt");
    c_mult->Print(multi);
#endif
    
    an->Assemble();
        an->Solve();
    c_mult->LoadSolutionFromMultiPhysics();
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    scalnames.Push("p");
    
    int div = 0;
    int dim=fGmesh->Dimension();
    std::string fileresult("flux_and_p.vtk");
    an->DefineGraphMesh(dim,scalnames,vecnames,fileresult);
    an->PostProcess(div,dim);
  
    TPZAnalysis *tracer_an = CreateTransportAnalysis(mul_sat, fsim_case);
    int n_steps = 20;
    REAL dt     = 0.1;
    TPZFMatrix<STATE> M_diag;
        TPZFMatrix<STATE> saturations = TimeForward(tracer_an, n_steps, dt, M_diag);
   
}
void RSimulatorConfiguration::InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh)
{
    
    int transport_matid = 1;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int64_t nel = cmesh->NElements();
    TPZManVector<std::vector<int>,4> gel_vol_index_to_cel(4);
    
    for (int64_t el = 0; el< nel; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) DebugStop();
        TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!celmp) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        
        int gel_dim = gel->Dimension();
        gel_vol_index_to_cel[gel_dim].push_back(el);
        
    }
    
  //  InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[3]);
    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[2]);
    std::ofstream file("mult_antes2D.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file);
    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[1]);
    std::ofstream file2("mult_antes2DP1D.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file2);
    
    cmesh->ComputeNodElCon();
}

//void RSimulatorConfiguration::InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes){
//
//    TPZGeoMesh * geometry = cmesh->Reference();
//    if (!geometry) {
//        DebugStop();
//    }
//
//    int mesh_dim = geometry->Dimension();
//
//    std::set<int> bcmat_ids = {-1,-2,-3,-4};
//    bool needs_all_boundaries_Q = true;
//    TPZManVector<int64_t,3> left_mesh_indexes(2,0);
//    left_mesh_indexes[0] = 0;
//    left_mesh_indexes[1] = 2;
//    TPZManVector<int64_t,3> right_mesh_indexes(1,0);
//    right_mesh_indexes[0] = 2;
//    for (auto cel_index: cel_indexes) {
//        TPZCompEl *cel = cmesh->Element(cel_index);
//        TPZGeoEl *gel = cel->Reference();
//
//        int nsides = gel->NSides();
//        int geldim = gel->Dimension();
//        for (int side = 0; side<nsides; side++) {
//            int sidedim = gel->SideDimension(side);
//            if(sidedim < geldim-1) continue;
//            int mat_id = gel->MaterialId();
//            if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
//                continue;
//            }
//            TPZStack<TPZCompElSide> celstack;
//            TPZCompElSide celside(cel,side);
//            TPZGeoElSide gelside(gel,side);
//            gelside.EqualLevelCompElementList(celstack, 0, 0);
//            if(celstack.size() == 0 && sidedim < mesh_dim && needs_all_boundaries_Q)
//            {
//                DebugStop();
//            }
//            if(celstack.size() == 0){
//                continue;
//            }
//
//            if(sidedim == geldim-1)
//            {
//
//                int count_dim_m_1 = 0;
//                int stack_index_dim_m_1 = 0;
//                int count = 0;
//                int stack_index = 0;
//                bool same_dimension_Q = false;
//                for (int icel = 0; icel < celstack.size(); icel++) {
//                    TPZCompEl * cel_neigh = celstack[icel].Element();
//                    TPZMultiphysicsInterfaceElement * mp_interface_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel_neigh);
//                    int gel_neigh_dim = cel_neigh->Dimension();
//                    if (!mp_interface_cel) {
//                        if (gel_neigh_dim == geldim-1) {
//                            count_dim_m_1++;
//                            stack_index_dim_m_1 = icel;
//                        }else if (gel_neigh_dim == geldim){
//                            count++;
//                            stack_index = icel;
//                            same_dimension_Q = true;
//                        }
//                    }
//                }
//
//                if (same_dimension_Q && count_dim_m_1 == 0) {
//
//                    /// There must be a lower dimensional transport element to equate the fluxes
//                    if(count != 1){
//                        DebugStop();
//                    }
//
//                    TPZGeoEl *neighgel = celstack[stack_index].Element()->Reference();
//                    if(geldim < neighgel->Dimension())
//                    {
//                        // we only create interfaces from lower to higher dimensional elements
//                        // neighgel must be a boundary element
//                        continue;
//                    }
//                    if(geldim == neighgel->Dimension() && gel->Id() > neighgel->Id())
//                    {
//                        // the interface is created from the lower id to the higher id
//                        continue;
//                    }
//
//                    // we need to create the interface
//                    TPZGeoElBC gbc(gelside,transport_matid);
//                    int64_t index;
//                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstack[stack_index]);
//                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
//                    continue;
//                }
//
//                /// There must be a lower dimensional transport element to equate the fluxes
//                if(count_dim_m_1 != 1){
//                    DebugStop();
//                }
//
//                TPZGeoEl *neighgel = celstack[stack_index_dim_m_1].Element()->Reference();
//                if(geldim < neighgel->Dimension())
//                {
//                    // we only create interfaces from lower to higher dimensional elements
//                    // neighgel must be a boundary element
//                    continue;
//                }
//                if(geldim == neighgel->Dimension() && gel->Id() > neighgel->Id())
//                {
//                    // the interface is created from the lower id to the higher id
//                    continue;
//                }
//                {
//                    int mat_interface_id;
//                    int neigh_mat_id = neighgel->MaterialId();
//                    if (bcmat_ids.find(neigh_mat_id) != bcmat_ids.end()) {
//                        mat_interface_id = neigh_mat_id;
//                    }else{
//                        mat_interface_id = transport_matid;
//                    }
//
//                    // we need to create the interface
//                    TPZGeoElBC gbc(gelside,mat_interface_id);
//                    int64_t index;
//                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstack[stack_index_dim_m_1]);
//                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
//                }
//            }
//        }
//    }
//
//}

void RSimulatorConfiguration::InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes){

    TPZGeoMesh * geometry = cmesh->Reference();
    if (!geometry) {
        DebugStop();
    }
    int mesh_dim = geometry->Dimension();

    std::set<int> bcmat_ids = {-1,-2,-3,-4};
    bool needs_all_boundaries_Q = false;
    TPZManVector<int64_t,3> left_mesh_indexes(2,0);
    left_mesh_indexes[0] = 0;
    left_mesh_indexes[1] = 2;
    TPZManVector<int64_t,3> right_mesh_indexes(1,0);
    right_mesh_indexes[0] = 2;
    //Create interface 2D- 2D
    for (auto cel_index: cel_indexes) {
        TPZCompEl *cel = cmesh->Element(cel_index);
        TPZMultiphysicsElement * celmult = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!celmult) {
            DebugStop();
        }
        celmult->CreateInterfaces();
    }
  

}

TPZAnalysis * RSimulatorConfiguration::CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    if (sim_data.UsePardisoQ) {
        
        TPZSpStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        
        return analysis;
    }else{
        
        TPZSkylineNSymStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
        
    }
    
}
TPZFMatrix<STATE> RSimulatorConfiguration::TimeForward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag){
    
    TPZMultiphysicsCompMesh * cmesh_transport = dynamic_cast<TPZMultiphysicsCompMesh *>(tracer_analysis->Mesh());
    
    if (!cmesh_transport) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,3> meshtrvec = cmesh_transport->MeshVector();
    
    /// Compute mass matrix M.
    TPZAutoPointer<TPZMatrix<STATE> > M;
    TPZFMatrix<REAL> F_inlet;
    {
        
        bool mass_matrix_Q = true;
        std::set<int> volumetric_mat_ids = {1,10};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing Mass Matrix." << std::endl;
        tracer_analysis->Assemble();
        std::cout << "Mass Matrix is computed." << std::endl;
        M = tracer_analysis->Solver().Matrix()->Clone();
    }
    
//    M->Print("M = ",std::cout,EMathematicaInput);
    
    int n_rows = M->Rows();
    M_diag.Resize(n_rows,M->Cols());
    
    for (int64_t i = 0; i < n_rows; i++) {
        M_diag(i,0) = M->Get(i, i);
    }
    


   
    int64_t n_eq = tracer_analysis->Mesh()->NEquations();
    TPZFMatrix<STATE> saturations(n_eq,n_steps);
    
    {
        bool mass_matrix_Q = false;
        std::set<int> volumetric_mat_ids = {1,10};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetTimeStep(dt);
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing transport operator K = M + T, and F_inlet " << std::endl;
        tracer_analysis->Assemble();
//        tracer_analysis->Rhs().Print("f = ",std::cout,EMathematicaInput);
 
        F_inlet = tracer_analysis->Rhs();
   
    }
    
    /// Time evolution
    std::string file_reservoir("transport.vtk");
    
    {
        int div = 0;
        TPZFMatrix<REAL> s_n(n_eq,1,0.0);
        TPZFMatrix<REAL> last_state_mass(n_eq,1,0.0);
        TPZFMatrix<REAL> s_np1;
        
        for (int it = 0; it < n_steps; it++) {
            
            for (int64_t i = 0; i < n_eq; i++) {
                last_state_mass(i,0) = M_diag(i,0)*s_n(i,0);
            }
            
            tracer_analysis->Rhs() = F_inlet - last_state_mass;
            tracer_analysis->Rhs() *= -1.0;
            
            tracer_analysis->Solve(); /// (LU decomposition)
            
            s_np1 = tracer_analysis->Solution();
            tracer_analysis->LoadSolution(s_np1);
    
            /// postprocess ...
            TPZStack<std::string,10> scalnames, vecnames;
            scalnames.Push("Sw");
            scalnames.Push("So");
            
            std::map<int,int> volumetric_ids;
            volumetric_ids.insert(std::make_pair(1, 2));
            
            std::map<int,int> fracture_ids;
            fracture_ids.insert(std::make_pair(6, 2));
            
            std::map<int,int> fracture_intersections_ids;
            fracture_intersections_ids.insert(std::make_pair(7, 1));
            
            for (auto data: volumetric_ids) {
                TPZMaterial * mat = cmesh_transport->FindMaterial(data.first);
                TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
                if (!volume) {
                    DebugStop();
                }
                volume->SetDimension(data.second);
            }
            
            int div = 0;
            std::set<int> mat_id_3D;
            mat_id_3D.insert(1);
            
            std::string file_reservoir("saturation.vtk");
            tracer_analysis->DefineGraphMesh(2,mat_id_3D,scalnames,vecnames,file_reservoir);
            tracer_analysis->PostProcess(div,2);
            
            //                std::set<int> mat_id_2D;
            //                mat_id_2D.insert(6);
            //                std::string file_frac("fracture_s.vtk");
            //                tracer_analysis->DefineGraphMesh(1,mat_id_2D,scalnames,vecnames,file_frac);
            //                tracer_analysis->PostProcess(div,1);
            
            
            
            // configuring next time step
            s_n = s_np1;
            for (int64_t i = 0; i < n_eq; i++) {
                saturations(i,it) = cmesh_transport->Solution()(i,0);
            }
            
        }
        
        
    }
    return saturations;
}
void Ladoderecho_2D (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    STATE rw=0.0;
 //   STATE r = Norm(x);
    double fx= (79.5775/((x + rw)*(x - rw)));
    
 //   double fx= 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);        //Force function definition
    
    //    double fx =-4144.653167389283*pow(10,
    //                                      -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x) + 4771.70829943056*
    //    pow(10,-pow(-2*M_PI + 15.*x,2) -
    //        pow(-2*M_PI + 15.*y,2))*pow(-2*M_PI + 15.*x,3) +
    //    4771.70829943056*pow(10,
    //                         -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
    //
    //    double fx = -4144.653167389282*pow(2,2 - pow(-2*M_PI + 15.*x,2) -
    //                                       pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x) + 4771.708299430558*
    //    pow(2,2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(-2*M_PI + 15.*x,3) +
    //    4771.708299430558*pow(2,
    //                          2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
    
    disp[0]=0.0;
    
    
}
