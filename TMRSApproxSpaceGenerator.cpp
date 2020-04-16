//
//  TMRSApproxSpaceGenerator.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#include "TMRSApproxSpaceGenerator.h"
#include "TPZMHMixedMeshWithTransportControl.h"
#include "TPZCompMeshTools.h"
#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

using namespace std;
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);
//void InsertMaterialObjects(TPZMHMixedMeshControl &control);

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
void TMRSApproxSpaceGenerator::CreateUniformMesh(int nx, REAL L, int ny, REAL h, int nz, REAL w){
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y
    
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,0.0);
    x1[0]=L;
    
    if (ny!=0) {
        x1[1]=h;
    }
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGenGrid2D gen(nels,x0,x1);
    gen.SetRefpatternElements(true);
    
    //    if (ny!=0) {
    //        <#statements#>
    //    }
    //    gen.SetElementType(EQuadrilateral);
    gen.Read(gmesh);
    
    if (nz!=0 ) {
        double var = w/nz;
        TPZExtendGridDimension extend(gmesh,var);
        extend.SetElType(1);
        gmesh = extend.ExtendedMesh(nz);
    }
    
    if (nz!=0) {
        for (auto gel:gmesh->ElementVec()) {
            TPZFMatrix<REAL> coordinates;
            gel->NodesCoordinates(coordinates);
            if(coordinates(2,0)==0){
                gel->CreateBCGeoEl(20, -1);
            }
            if(coordinates(2,4)==w){
                gel->CreateBCGeoEl(25, -2);
            }
            
            if(coordinates(0,0)==0.0 ){
                gel->CreateBCGeoEl(24, -3);
            }
            if(coordinates(1,0)==0.0 ){
                gel->CreateBCGeoEl(21, -4);
            }
            
            if(coordinates(0,1)== L ){
                gel->CreateBCGeoEl(22, -5);
            }
            if(coordinates(1,3)==h){
                gel->CreateBCGeoEl(23, -6);
            }
        };
        gmesh->SetDimension(3);
    }
    
    if (ny!=0 && nz==0) {
                gen.SetBC(gmesh, 4, -1);
                gen.SetBC(gmesh, 5, -2);
                gen.SetBC(gmesh, 6, -3);
                gen.SetBC(gmesh, 7, -4);
        gmesh->SetDimension(2);
    }
    
    if (ny==0 && nz==0) {
        double dh = L/nx;
        int Domain_Mat_Id = 1;
        int Inlet_bc_Id = -1;
        int Outletbc_Id = -2;
        TPZVec<REAL> xp(3,0.0);   //Creates vector nodes
        
        gmesh->NodeVec().Resize(nx+1);
        for (int64_t i=0; i<nx+1; i++) {
            xp[0] =(i)*dh;
            gmesh->NodeVec()[i]= TPZGeoNode(i, xp, *gmesh);
        }
        
        TPZVec<int64_t> cornerindexes(2);   //Creates elements
        for (int64_t iel=0; iel<nx; iel++) {
            cornerindexes[0]=iel;
            cornerindexes[1]=iel+1;
            gmesh->CreateGeoElement(EOned, cornerindexes, Domain_Mat_Id, iel);
        }
        gmesh->Element(0)->CreateBCGeoEl(0, Inlet_bc_Id);     //Sets BC
        gmesh->Element(nx-1)->CreateBCGeoEl(1, Outletbc_Id);
        gmesh->SetDimension(1);
        gmesh->BuildConnectivity();
    }
    
    gmesh->BuildConnectivity();
    
    mGeometry = gmesh;
    
   
#ifdef PZDEBUG
    if (!mGeometry)
    {
        std::cout << "The geometrical mesh was not generated." << std::endl;
        DebugStop();
    }
#endif
}

void TMRSApproxSpaceGenerator::GenerateMHMUniformMesh(int nelref){
    
    
    for (int iref=0; iref<nelref; iref++) {
        int nel = mGeometry->NElements();
        for (int iel = 0; iel <nel; iel++) {
            TPZGeoEl *gel = mGeometry->Element(iel);
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            TPZVec<TPZGeoEl *> sons;
            gel->Divide(sons);
        }
    }
   
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
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
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
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
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
    
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
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
        TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
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
    
    if (!mSimData.mTNumerics.m_four_approx_spaces_Q) {
        int ncon = cmesh->NConnects();
        //Set Lagrange multiplier
        for(int i=0; i<ncon; i++){
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
    }

    
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
    
    if (mSimData.mTNumerics.m_four_approx_spaces_Q) {
        BuildMixed4SpacesMultiPhysicsCompMesh(order);
    }else{
        if (mSimData.mTNumerics.m_mhm_mixed_Q) {
            BuildMHMMixed2SpacesMultiPhysicsCompMesh();
        }
        else{
        BuildMixed2SpacesMultiPhysicsCompMesh(order);
        }
    }
}

void TMRSApproxSpaceGenerator::BuildMixed2SpacesMultiPhysicsCompMesh(int order){
    
    int dimension = mGeometry->Dimension();
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
//    TPZMixedDarcyFlow *volume = nullptr;
    mMixedOperator->SetDefaultOrder(order);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,d);
//            volume = new TPZMixedDarcyFlow(materia_id, d);
//             volume->SetPermeability(1.0);
            volume->SetDataTransfer(mSimData);
            mMixedOperator->InsertMaterialObject(volume);
           
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mMixedOperator->InsertMaterialObject(face);
    }
    
    
    TPZManVector<TPZCompMesh *, 3> mesh_vec(3);
    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order);
    mesh_vec[2] = DiscontinuousCmesh();
    TPZManVector<int,5> active_approx_spaces(3);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 0;
    mMixedOperator->SetDimModel(dimension);
//    mMixedOperator->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec);

    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "mixed_cmesh_two_spaces" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    mMixedOperator->Print(sout);
#endif
    
}

void TMRSApproxSpaceGenerator::BuildMixed4SpacesMultiPhysicsCompMesh(int order){
    
    int dimension = mGeometry->Dimension();
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
    mMixedOperator->SetDefaultOrder(order);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,d);
            volume->SetDataTransfer(mSimData);
            mMixedOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mMixedOperator->InsertMaterialObject(face);
    }
    
    
    TPZManVector<TPZCompMesh *, 5> mesh_vec(5);
    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order);
    mesh_vec[2] = DiscontinuousCmesh();
    mesh_vec[3] = DiscontinuousCmesh();
    mesh_vec[4] = DiscontinuousCmesh();
    TPZManVector<int,5> active_approx_spaces(5);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    active_approx_spaces[3] = 1;
    active_approx_spaces[4] = 0;
    mMixedOperator->SetDimModel(dimension);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec);
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "mixed_cmesh_four_spaces" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    mMixedOperator->Print(sout);
#endif
    
}

void TMRSApproxSpaceGenerator::BuildMHMMixed2SpacesMultiPhysicsCompMesh(){
    
    TPZGeoMesh *gmeshcoarse = GetGeometry();
    
        TPZGeoMesh * gmeshauto = gmeshcoarse; //Autopointer2
       
        TPZMHMixedMeshWithTransportControl *mhm = new TPZMHMixedMeshWithTransportControl(gmeshauto);
        TPZVec<int64_t> coarseindices;
        ComputeCoarseIndices(gmeshauto, coarseindices); //operator->()
        
//        gmeshauto->AddInterfaceMaterial(1, 2, interface_mat_id);
//        gmeshauto->AddInterfaceMaterial(2, 1, interface_mat_id);
        
    
        mhm->DefinePartitionbyCoarseIndices(coarseindices);
        std::set<int> matids;
    
        //Insert Material Objects
        {
            int dimension = mGeometry->Dimension();
            TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
            std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
            for (int d = 0; d <= dimension; d++) {
                for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
                    std::string material_name = chunk.first;
                    std::cout << "physical name = " << material_name << std::endl;
                    int materia_id = chunk.second;
                    matids.insert(materia_id);
                    volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,d);
                    volume->SetDataTransfer(mSimData);
                    mhm->CMesh()->InsertMaterialObject(volume);
                }
            }
            mhm->fMaterialIds = matids;
            matids.clear();
            if (!volume) {
                DebugStop();
            }
            
            TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
            TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
            for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
                int bc_id   = get<0>(chunk);
                int bc_type = get<1>(chunk);
                val2(0,0)   = get<2>(chunk);
                matids.insert(bc_id);
                TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
                mhm->CMesh()->InsertMaterialObject(face);
            }
            mhm->fMaterialBCIds = matids;
        }

        
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
    
        mhm->SetInternalPOrder(1);
        mhm->SetSkeletonPOrder(1);
        
        mhm->DivideSkeletonElements(0);
        mhm->DivideBoundarySkeletonElements();
    
        bool substructure = true;
    
        mhm->SetApproxSpaceGenerator(this);
        mhm->BuildComputationalMesh(substructure);
    
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream out("MixedMeshControlHDiv.txt");
            meshcontrol.Print(out);
        }
#endif
        
        std::cout << "MHM Hdiv Computational meshes created\n";
#ifdef PZDEBUG2
        if(1)
        {
            std::ofstream gfile("geometryMHMHdiv.txt");
            gmeshauto->Print(gfile);
            std::ofstream out_mhm("MHM_hdiv.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        
        std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;
    
    TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
    
    TPZMultiphysicsCompMesh *cmeshtest = dynamic_cast<TPZMultiphysicsCompMesh*>(MixedMesh);
    
    mMixedOperator = cmeshtest;

}

void TMRSApproxSpaceGenerator::BuildMHMMixed4SpacesMultiPhysicsCompMesh(){
    DebugStop();
}

void TMRSApproxSpaceGenerator::BuildMixedSCStructures(){
    
    std::cout << "DoF:: Before SC = " << mMixedOperator->NEquations() << std::endl;
    
    mMixedOperator->ComputeNodElCon();
    if (mSimData.mTNumerics.m_four_approx_spaces_Q) {
        int dim = mMixedOperator->Dimension();
        int64_t nel = mMixedOperator->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = mMixedOperator->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
    }
    
    // Created condensed elements for the elements that have internal nodes
    bool KeepOneLagrangianQ = true;
    bool KeepMatrixQ = true;
    TPZCompMeshTools::CreatedCondensedElements(mMixedOperator, KeepOneLagrangianQ, KeepMatrixQ);
    mMixedOperator->ComputeNodElCon();
    mMixedOperator->InitializeBlock();
    std::cout << "DoF:: After SC = " << mMixedOperator->NEquations() << std::endl;
}


void TMRSApproxSpaceGenerator::BuildTransportMultiPhysicsCompMesh(){
    
    if (mSimData.mTNumerics.m_four_approx_spaces_Q) {
        BuildTransport4SpacesMultiPhysicsCompMesh();
    }else{
        BuildTransport2SpacesMultiPhysicsCompMesh();
    }
    
}

void TMRSApproxSpaceGenerator::BuildTransport2SpacesMultiPhysicsCompMesh(){
    
    if (!mMixedOperator || !mGeometry) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,3> mixed_meshvec = mMixedOperator->MeshVector();
    TPZManVector<TPZCompMesh *,3> transport_meshvec(3);
    
    transport_meshvec[0] = mixed_meshvec[0];
    transport_meshvec[1] = mixed_meshvec[1];
    transport_meshvec[2] = DiscontinuousCmesh();
    
    
    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
    mTransportOperator->SetDefaultOrder(0);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSMultiphaseFlow<TMRSMemory>(materia_id,d);
            volume->SetDataTransfer(mSimData);
            mTransportOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mTransportOperator->InsertMaterialObject(face);
        
    }
    
    int transport_matid = 100;
    {
        TMRSMultiphaseFlow<TMRSMemory> * interface = new TMRSMultiphaseFlow<TMRSMemory>(transport_matid,dimension-1);
        interface->SetDataTransfer(mSimData);
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
                bool condition =false;
                while (neig !=gelside && condition != true) {

                    for (int d = 0; d <= dimension; d++) {
                        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
                            int material_id = chunk.second;
                            if (neig.Element()->MaterialId() == material_id) {
                                condition = true;
                            }

                        }
                    }
                    if (condition == false) {
                        for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
                            int bc_id   = get<0>(chunk);
                            if (neig.Element()->MaterialId() == bc_id) {
                                condition = true;
                            }
                        }
                    }


                    if (condition == false) {
                        neig=neig.Neighbour();
                    }
                }

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

void TMRSApproxSpaceGenerator::BuildTransport4SpacesMultiPhysicsCompMesh(){
    
    if (!mMixedOperator || !mGeometry) {
        DebugStop();
    }
    
    TPZManVector<TPZCompMesh *,5> mixed_meshvec = mMixedOperator->MeshVector();
    TPZManVector<TPZCompMesh *,5> transport_meshvec(5);
    
    transport_meshvec[0] = mixed_meshvec[0];
    transport_meshvec[1] = mixed_meshvec[1];
    transport_meshvec[2] = mixed_meshvec[2];
    transport_meshvec[3] = mixed_meshvec[3];
    transport_meshvec[4] = DiscontinuousCmesh();
    
    
    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
    mTransportOperator->SetDefaultOrder(0);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSMultiphaseFlow<TMRSMemory>(materia_id,d);
            volume->SetDataTransfer(mSimData);
            mTransportOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2(0,0)   = get<2>(chunk);
        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mTransportOperator->InsertMaterialObject(face);
    }
    
    int transport_matid = 100;
    {
        TMRSMultiphaseFlow<TMRSMemory> * interface = new TMRSMultiphaseFlow<TMRSMemory>(transport_matid,dimension-1);
        interface->SetDataTransfer(mSimData);
        mTransportOperator->InsertMaterialObject(interface);
    }
    
    mTransportOperator->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(5); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 0;
    active_approx_spaces[3] = 0;
    active_approx_spaces[4] = 1;
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
        left_mesh_indexes[1] = 4;
        TPZManVector<int64_t,3> right_mesh_indexes(1,0);
        right_mesh_indexes[0] = 4;
        
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

void TMRSApproxSpaceGenerator::SetDataTransfer(TMRSDataTransfer & SimData){
    mSimData = SimData;
}

TMRSDataTransfer & TMRSApproxSpaceGenerator::GetDataTransfer(){
    return mSimData;
}

TPZMultiphysicsCompMesh * TMRSApproxSpaceGenerator::GetMixedOperator(){
    return mMixedOperator;
}

TPZMultiphysicsCompMesh * TMRSApproxSpaceGenerator::GetTransportOperator(){
    return mTransportOperator;
}

void TMRSApproxSpaceGenerator::LinkMemory(TPZMultiphysicsCompMesh * MixedOperator, TPZMultiphysicsCompMesh * TransportOperator){
 
    AdjustMemory(MixedOperator, TransportOperator);
    for (auto item : mSimData.mTGeometry.mDomainDimNameAndPhysicalTag[2]) {
        int material_id = item.second;
        UnifyMaterialMemory(material_id, MixedOperator, TransportOperator);
        FillMaterialMemory(material_id, MixedOperator, TransportOperator);
    }
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
//            kappa *= rand() % 40 +1;
            for (int j = 0; j < 3; j++) {
                mem.m_kappa(j,j) = kappa;
                mem.m_kappa_inv(j,j) = 1.0/kappa;
            }
        }
        
        
//        int nels = cmesh->NElements();
//        for (int iel = 0; iel< nels; iel++) {
//            TPZCompEl *cel = cmesh->Element(iel);
//            if(!cel){
//                continue;
//
//            }
//            TPZGeoEl *gel = cel->Reference();
//            if (!gel || gel->HasSubElement()) {
//                continue;
//            }
//
//            if (cel->Dimension()!= cmesh->Dimension()) {
//                continue;
//            }
//            if (!MixedOperator->Element(iel)) {
//                continue;
//            }
//            TPZVec<int64_t> indices;
//            cel->GetMemoryIndices(indices);
//            TPZVec<REAL> qsi(3,0.25);
//            qsi[2]=0.0;
//            TPZVec<REAL> point(3,0.0);
//            gel->X(qsi, point);
////            if (gel->MaterialId()!=2){
////                continue;
////            }
//
////            int val = rand() % 100;
//
//
//            REAL kappa =  1000*(sin(point[0])*sin(point[1]) + 2);
//
//
//
//
////            REAL kappa = 100000.0 + 1*(sin(point[0])*sin(point[1])+2);
//
//
//            for (auto memIndex: indices) {
//                if (memIndex<=0) {
//                    continue;
//                }
//
//
//                TMRSMemory &mem = memory_vector.get()->operator [](memIndex);
//                mem.m_sw = 0.0;
//                mem.m_phi = 0.1;
//                mem.m_kappa.Resize(3, 3);
//                mem.m_kappa.Zero();
//                mem.m_kappa_inv.Resize(3, 3);
//                mem.m_kappa_inv.Zero();
//                for (int j = 0; j < 3; j++) {
//                    mem.m_kappa(j,j) = kappa;
//                    mem.m_kappa_inv(j,j) = 1.0/kappa;
//                }
//            }
//
//        
//        }
//        
        
        
        
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

//
void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices)
{
    //    {
    //        std::ofstream out("gmeshref.txt");
    //        gmesh->Print(out);
    //    }
    coarseindices.Resize(gmesh->NElements());
    int count = 0;
    for (int64_t el=0; el<gmesh->NElements(); el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()) continue;
        if(gel->Father()) continue;
        coarseindices[count] = el;
        count++;
    }
    coarseindices.Resize(count);
    }
//void  TMRSApproxSpaceGenerator::InsertMaterialObjects(TPZMHMeshControl &control)
//{
//    TPZCompMesh &cmesh = control.CMesh();
//    TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
//    int order = 1;
////    const int typeFlux = 1, typePressure = 0;
//    int dimension =mGeometry->Dimension();
////    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
//    TPZMixedDarcyFlow * volume = nullptr;
//
//    MixedFluxPressureCmesh->SetDefaultOrder(order);
//    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
//    for (int d = 0; d <= dimension; d++) {
//        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
//            std::string material_name = chunk.first;
//            std::cout << "physical name = " << material_name << std::endl;
//            int materia_id = chunk.second;
////            volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,d);
//            volume = new TPZMixedDarcyFlow(materia_id,d);
//
////            volume->SetDataTransfer(mSimData);
//            MixedFluxPressureCmesh->InsertMaterialObject(volume);
//        }
//    }
//    
//    if (!volume) {
//        DebugStop();
//    }
//    
//    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
//    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
//    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
//        int bc_id   = get<0>(chunk);
//        int bc_type = get<1>(chunk);
//        val2(0,0)   = get<2>(chunk);
//        TPZMaterial * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
//        MixedFluxPressureCmesh->InsertMaterialObject(face);
//    }
//    
//}

void TMRSApproxSpaceGenerator::InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
  

        TPZCompMesh &cmesh = control.CMesh();
    
        TPZGeoMesh &gmesh = control.GMesh();
        const int typeFlux = 1, typePressure = 0;
        TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,0.);
    
    
        int dim = gmesh.Dimension();
        cmesh.SetDimModel(dim);
    
        TPZCompMesh *MixedFluxPressureCmesh = &cmesh;
    
        // Material medio poroso
        TMRSDarcyFlowWithMem<TMRSMemory> * mat = new TMRSDarcyFlowWithMem<TMRSMemory>(1,dim);
        mat->SetDataTransfer(mSimData);
//        mat->SetPermeability(1.);
        //    mat->SetForcingFunction(One);
        MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    
    
        // Bc N
        TPZBndCond * bcN = mat->CreateBC(mat, -1, typeFlux, val1, val2Flux);
        //    bcN->SetForcingFunction(0, force);
    
        MixedFluxPressureCmesh->InsertMaterialObject(bcN);
        bcN = mat->CreateBC(mat, -3, typeFlux, val1, val2Flux);
        //    bcN->SetForcingFunction(0, force);
    
        MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
        val2Pressure(0,0)=10;
        TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);
    
        MixedFluxPressureCmesh->InsertMaterialObject(bcS);
        val2Pressure(0,0)=1000;
        bcS = mat->CreateBC(mat, -4, typePressure, val1, val2Pressure);
        MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
//
//    //aqui
//    int dimension = mGeometry->Dimension();
//
//    TPZMixedDarcyFlow * volume = nullptr;
//    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
//    for (int d = 0; d <= dimension; d++) {
//        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
//            std::string material_name = chunk.first;
//            std::cout << "physical name = " << material_name << std::endl;
//            int materia_id = chunk.second;
//            volume = new TPZMixedDarcyFlow(materia_id,d);
//
////            volume->SetDataTransfer(mSimData);
//            MixedFluxPressureCmesh->InsertMaterialObject(volume);
//        }
//    }
//
//    if (!volume) {
//        DebugStop();
//    }
//
//    TPZFMatrix<STATE> val1(1,1,0.0),val2(1,1,0.0);
//    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
//    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
//        int bc_id   = get<0>(chunk);
//        int bc_type = get<1>(chunk);
//        val2(0,0)   = get<2>(chunk);
//        double valimposto = val2(0,0) ;
//        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
//        MixedFluxPressureCmesh->InsertMaterialObject(face);
//    }
//
    
    
}
