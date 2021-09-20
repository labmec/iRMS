//
//  TMRSApproxSpaceGenerator.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#include "TMRSApproxSpaceGenerator.h"
#include "TPZMHMixedMeshWithTransportControl.h"
#include "TPZCompMeshTools.h"
#include "TPZMixedDarcyWithFourSpaces.h"
#include "TPZMHMixedMesh4SpacesControl.h"
#include "TPZFastCondensedElement.h"
#include "TPZReservoirTools.h"
#include "TPZPostProcessResProp.h"
#include "TPZDarcyMemory.h"
#include "TPZDarcyFlowWithMem.h"
#include "TMRSDarcyFractureFlowWithMem.h"
#include "TPZLagrangeMultiplierCS.h"
#include "TPZCompElHDivCollapsed.h"
#include "TMRSDarcyMemory.h"
#include "TMRSTransportMemory.h"
#include "TPZNullMaterialCS.h"
#include "pzcompelwithmem.h"
#include "pzsmanal.h"
#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

using namespace std;

void ComputeCoarseIndices(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);

TMRSApproxSpaceGenerator::TMRSApproxSpaceGenerator()
{
    mGeometry = nullptr;
    mMixedOperator = nullptr;
    mTransportOperator = nullptr;
    mHybridizer = nullptr;
}

TMRSApproxSpaceGenerator &TMRSApproxSpaceGenerator::operator=(const TMRSApproxSpaceGenerator &other){
    DebugStop();
}

TMRSApproxSpaceGenerator::~TMRSApproxSpaceGenerator(){
    if (mHybridizer){
        delete mHybridizer;
    }
    std::cout << __PRETTY_FUNCTION__ << std::endl;
}

void TMRSApproxSpaceGenerator::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TMRSApproxSpaceGenerator::Read(TPZStream &buf, void *context){
    DebugStop();
}

int TMRSApproxSpaceGenerator::ClassId() const{
    DebugStop();
    return -1;
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
void TMRSApproxSpaceGenerator::LoadGeometry(std::string geometry_file,
        TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag){
    
    TPZGmshReader Geometry;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
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
    bool IsQuad= true;
    if (IsQuad) {
        gen.SetElementType(MMeshType::EQuadrilateral);
    }
    else{
         gen.SetElementType(MMeshType::ETriangular);
    }
    
   
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
            int nsides = gel->NSides();
            if(coordinates(2,0)==0){
                if(IsQuad){
                gel->CreateBCGeoEl(20, -1);
                }
                else{
                gel->CreateBCGeoEl(15, -1);
                }
            }
            
            if(coordinates(2,4)==w){
                if(IsQuad){
                    gel->CreateBCGeoEl(25, -1);
                }
                else{
                    gel->CreateBCGeoEl(19, -1);
                }
                
            }
            REAL sum_cords =coordinates(0,0)+coordinates(0,2)+coordinates(0,3)+coordinates(0,5);
            if(coordinates(0,0)==0.0 ){
                if (IsQuad) {
                    gel->CreateBCGeoEl(24, -4);
                }
                else if(sum_cords==0){
                   
                    gel->CreateBCGeoEl(18, -4);
                }
                
            }
            sum_cords =coordinates(1,0)+coordinates(1,1)+coordinates(1,3)+coordinates(1,4);
            if(coordinates(1,0)==0.0 ){
                
                
                if (IsQuad) {
                    gel->CreateBCGeoEl(21, -3);
                }
                else if (sum_cords==0){
                   
                    gel->CreateBCGeoEl(16, -3);
                }
                
            }
            sum_cords =coordinates(0,1)+coordinates(0,2)+coordinates(0,4)+coordinates(0,5);
            if(coordinates(0,1)== L ){
                if (IsQuad) {
                     gel->CreateBCGeoEl(22, -2);
                }
                else if (sum_cords==4*L)
                {
                    
                    gel->CreateBCGeoEl(17, -2);
                }
               
            }
            if(coordinates(1,1)==h){
                if (IsQuad) {
                    gel->CreateBCGeoEl(23, -3);
                }
                else{
                    std::cout<<"Cords: "<<std::endl;
                    coordinates.Print(std::cout);
                    gel->CreateBCGeoEl(17, -3);
                }
                
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

void TMRSApproxSpaceGenerator::ApplyUniformRefinement(int nelref){
    
    
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
void TMRSApproxSpaceGenerator::ApplyUniformRefinement(){
    std::cout << "Applying uniform refinement numref = " << mSimData.mTGeometry.mnref << "\n";

    ApplyUniformRefinement(mSimData.mTGeometry.mnref);
}
void TMRSApproxSpaceGenerator::PrintGeometry(std::string name, bool vtkFile, bool textfile)
{
    if (!mGeometry) {
        DebugStop();
    }
    std::stringstream text_name;
    std::stringstream vtk_name;
    if (vtkFile) {
        vtk_name   << name << "_geometry"  << ".vtk";
        std::ofstream vtkfile(vtk_name.str().c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, vtkfile, true);
    }
    if (textfile) {
        text_name  << name << "_geometry" << ".txt";
        std::ofstream textfile(text_name.str().c_str());
        mGeometry->Print(textfile);
        
    }
    
}


TPZCompMesh * TMRSApproxSpaceGenerator::HdivFluxCmesh(int order){
    
    if (!mGeometry) {
        DebugStop();
    }
    
    // PHIL : nesta malha nao ha necessidade de inserir materiais "de verdade"
    // poderia incluir "NullMaterial"...
    // este observacao vale para todos as malhas atomicas
    
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZNullMaterial<> *volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(order);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    REAL kappa = 1.0;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            int materia_id = chunk.second;
            volume = new TPZNullMaterial<>(materia_id);
            cmesh->InsertMaterialObject(volume);
        }
    }

    if (!volume) {
        DebugStop();
    }
    
    TMRSDarcyFractureFlowWithMem<TMRSMemory> * fracmat = nullptr;
    for(auto chunk : mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag[dimension-1])
    {
        std::string material_name = chunk.first;
        std::cout << "physical name = " << material_name <<
        " material id " << chunk.second << " dimension " << dimension-1 << std::endl;
        int materia_id = chunk.second;
//        matsWithMem.insert(materia_id);
        fracmat = new TMRSDarcyFractureFlowWithMem<TMRSMemory>(materia_id,dimension-1);
        TMRSMemory defaultmem;
        fracmat->SetDataTransfer(mSimData);
        cmesh->InsertMaterialObject(fracmat);

    }
    
    // PHIL : as malhas atomicas nao precisam ser BndCond
    // poderiam ser null material...
    
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]   = get<2>(chunk);
        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
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

#include "pzshapequad.h"

/// create an HDiv mesh for mortar approximation
TPZCompMesh *TMRSApproxSpaceGenerator::HDivMortarFluxCmesh(char fluxmortarlagrange)
{
    std::set<int> matids, bcmatids;
    int dimension = mGeometry->Dimension();
    int nstate = 1;
    GetMaterialIds(dimension, matids, bcmatids);
//    bcmatids.insert(-11);
    
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    cmesh->SetName("FluxMortarMesh");
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    for (auto it:matids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dimension,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    for (auto it:bcmatids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dimension-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId,dimension-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    
    
    std::set<int> buildmatids(matids);
    // create all flux elements as discontinuous elements
    cmesh->AutoBuild(buildmatids);
    // add the Hdiv wrapper elements as boundary of the volumetric elements
    mGeometry->ResetReference();
    {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(mSimData.mTGeometry.m_HdivWrapMatId,dimension-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
#ifdef PZDEBUG
    {
        {
            std::ofstream out("FluxMortarMesh.txt");
//            cmesh->Print(out);
        }
    }
#endif
    int64_t nelflux = cmesh->NElements();
    for (int64_t el = 0; el<nelflux; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != dimension) continue;
        int firstside = 0;
        for(int dim = 0; dim < dimension-1; dim++) firstside+=gel->NSides(dim);
        for (int side=firstside; side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if(neighbour.Element()->MaterialId() != mSimData.mTGeometry.m_HdivWrapMatId) DebugStop();
//            int hdiv_orient = gel->NormalOrientation(side);
            gel->SetReference(cel);
            int64_t index;
            int neighside = neighbour.Side();
            TPZCompEl *celwrap = cmesh->CreateCompEl(neighbour.Element(), index);
//            TPZInterpolationSpace *space = dynamic_cast<TPZInterpolationSpace *>(celwrap);
//            space->SetSideOrient(neighside,hdiv_orient);
            gel->ResetReference();
            neighbour.Element()->ResetReference();
        }
    }

    buildmatids.clear();
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        buildmatids.insert(19);
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(19,2,1);
        cmesh->InsertMaterialObject(nullmat);
    }
    buildmatids.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    buildmatids.insert(bcmatids.begin(),bcmatids.end());
    int borderOrder =mSimData.mTNumerics.m_BorderElementFluxOrder;
    cmesh->SetDefaultOrder(borderOrder);
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    // create all flux elements as discontinuous elements
    cmesh->AutoBuild(buildmatids);

    // vector with computational element pointers for creating fracture elements
    TPZVec<TPZCompEl *> fracsupport(mGeometry->NElements(),0);
    
  
    // set the lagrange level
    {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() < dimension && gel->MaterialId() == mSimData.mTGeometry.m_zeroOrderHdivFluxMatId)
            {
                fracsupport[gel->Index()] = cel;
                cel->Connect(0).SetLagrangeMultiplier(fluxmortarlagrange);
            }
        }
    }
    
   
    
    cmesh->SetDimModel(dimension-1);
    matids.clear();
    bcmatids.clear();
    GetMaterialIds(dimension-1, matids, bcmatids);
    for(auto matid : matids)
    {
        cmesh->InsertMaterialObject(new TPZNullMaterial(matid,dimension-1,1));
    }
    for(auto matid : bcmatids)
    {
        cmesh->InsertMaterialObject(new TPZNullMaterial(matid,dimension-2,1));
    }
    bool complicated = false;
    if(complicated)
    {
        // data structure that stores for each geometric fracture element the two element/side volumetric elements
        // it is facing
        typedef std::pair<TPZGeoElSide,TPZGeoElSide> sidepair;
        // first geometric fracture element
        // second : left/right geoelsides of volumetric elements
        typedef std::pair<TPZGeoEl*,sidepair> gelsideandpair;
        std::list<gelsideandpair> fracElements;
        {
            int64_t nel = mGeometry->NElements();
            for(int64_t el = 0; el<nel; el++)
            {
                TPZGeoEl *gel = mGeometry->Element(el);
                if(!gel || gel->Dimension() != dimension-1) continue;
                int matid = gel->MaterialId();
                // we are looking for fracture elements
                if(matids.find(matid) == matids.end()) continue;
                TPZGeoElSide gelside(gel);
                TPZGeoElSide neighbour = gelside.Neighbour();
                sidepair leftright;
                int nfound = 0;
                // look for volumetric neighbours of the fracture element
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension() == dimension)
                    {
                        if(neighbour.Element()->NormalOrientation(neighbour.Side()) == 1)
                        {
                            leftright.first = neighbour;
                            nfound++;
                        }
                        else
                        {
                            leftright.second = neighbour;
                            nfound++;
                        }
                    }
                    neighbour=neighbour.Neighbour();
                }
                if(nfound != 2) DebugStop();
                gelsideandpair gs;
                gs.first = gel;
                gs.second = leftright;
                fracElements.push_back(gs);
            }
            // now, we have the elements and the H(div) elements connected
            // we will now create the HDiv collapsed elements, link the connects and correct the interface matid
            for(auto it : fracElements)
            {
                TPZGeoEl *fracgel = it.first;
                TPZGeoElSide fracgelside(fracgel);
                TPZGeoElSide leftgelside = it.second.first;
                TPZGeoElSide rightgelside = it.second.second;
                // find the element with material id m_zeroOrderHdivFluxMatId
                TPZGeoElSide zeroflux = fracgelside.HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                TPZCompElSide zerofluxcomp = zeroflux.Reference();
                zerofluxcomp.SetElement(fracsupport[zeroflux.Element()->Index()]);
                if(!zerofluxcomp) DebugStop();
                int zerofluxorder = zerofluxcomp.Element()->Connect(0).Order();
                if(!zeroflux) DebugStop();
                if(!zerofluxcomp) DebugStop();
                TPZGeoElSide intfaceleft;
                TPZGeoElSide intfaceright;
                // based on the order of creation of the geometric wrap elements we reach the interface elements
                {
                    TPZGeoElSide hdivwrap = leftgelside.Neighbour();
                    TPZGeoElSide intface1 = hdivwrap.Neighbour();
                    TPZGeoElSide pressmortar = intface1.Neighbour();
                    intfaceleft = pressmortar.Neighbour();
                }
                {
                    TPZGeoElSide hdivwrap = rightgelside.Neighbour();
                    TPZGeoElSide intface1 = hdivwrap.Neighbour();
                    TPZGeoElSide pressmortar = intface1.Neighbour();
                    intfaceright = pressmortar.Neighbour();
                }
                if(intfaceleft.Element()->MaterialId() != mSimData.mTGeometry.m_negLagrangeMatId) DebugStop();
                if(intfaceright.Element()->MaterialId() != mSimData.mTGeometry.m_posLagrangeMatId) DebugStop();
                // create a second zero flux element and second connect (hybridizing the mesh
                TPZGeoElBC gbc(zeroflux,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                TPZGeoElSide zeroflux2 = gbc.CreatedElement();
                zeroflux.Element()->ResetReference();
                int64_t index;
                cmesh->SetDefaultOrder(zerofluxorder);
                TPZCompEl *celflux2 = cmesh->CreateCompEl(zeroflux2.Element(), index);
                TPZCompElSide zerofluxcomp2 = zeroflux2.Reference();
                
                zeroflux2.Element()->ResetReference();
                cmesh->SetDefaultOrder(1);
                TPZInterpolationSpace *HDivCollapsed = 0;
                if(zeroflux.Element()->Type() == ETriangle)
                {
                    int64_t index;
                    HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeTriang>(*cmesh,fracgel,index);
                }
                else if(zeroflux.Element()->Type() == EQuadrilateral)
                {
                    int64_t index;
                    HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,fracgel,index);
                }
                int nconnects = HDivCollapsed->NConnects();
                // bottom connect index
                int64_t cindex1 = HDivCollapsed->ConnectIndex(nconnects-2);
                // top connect index
                int64_t cindex2 = HDivCollapsed->ConnectIndex(nconnects-1);
                int nsides = fracgel->NSides();
                // changing orientation of top
                HDivCollapsed->SetSideOrient(nsides, -1);
                // linking bottom to first zero flux element - interface left - zeroflux - frac pressure
                HDivCollapsed->SetConnectIndex(nconnects-2, zerofluxcomp.Element()->ConnectIndex(0));
                // linking top to second zero flux element - interface right - zeroflux - frac pressure
                HDivCollapsed->SetConnectIndex(nconnects-1, zerofluxcomp2.Element()->ConnectIndex(0));
                cmesh->ConnectVec()[cindex1].DecrementElConnected();
                cmesh->ConnectVec()[cindex2].DecrementElConnected();
                cmesh->ConnectVec()[cindex1].SetSequenceNumber(-1);
                cmesh->ConnectVec()[cindex2].SetSequenceNumber(-1);
                cmesh->ConnectVec().SetFree(cindex1);
                cmesh->ConnectVec().SetFree(cindex2);
            }
        }
    }
    else // not so complicated
    {
        // create the fracture elements and connect them to the volumetric elements
        int64_t nel = mGeometry->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZGeoEl *gel = mGeometry->Element(el);
            if(!gel || gel->Dimension() != dimension-1) continue;
            int matid = gel->MaterialId();
            // we are looking for fracture elements
            if(matids.find(matid) == matids.end()) continue;
            // now we have a fracture element
            TPZGeoElSide gelside(gel);
            // find the element with material id m_zeroOrderHdivFluxMatId
            TPZGeoElSide zeroflux = gelside.HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
            if(!zeroflux) DebugStop();
            TPZCompElSide zerofluxcomp = zeroflux.Reference();
            zerofluxcomp.SetElement(fracsupport[zeroflux.Element()->Index()]);
            if(!zerofluxcomp) DebugStop();
            int zerofluxorder = zerofluxcomp.Element()->Connect(0).Order();
            
            int domain=0;
            TPZGeoElSide zeroflux2;
            if(mSimData.mTNumerics.m_mhm_mixed_Q){
                if(mSubdomainIndexGel.size() != mGeometry->NElements()) DebugStop();
                domain = mSubdomainIndexGel[zeroflux.Element()->Index()];
                // create a second zero flux element and second connect (hybridizing the mesh
                TPZGeoElBC gbc(zeroflux,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                zeroflux2 = gbc.CreatedElement();
                {
                    auto nel = mGeometry->NElements();
                    mSubdomainIndexGel.resize(nel);
                    mSubdomainIndexGel[gbc.CreatedElement()->Index()] = domain;
                }
            }
            else{
                // create a second zero flux element and second connect (hybridizing the mesh
                TPZGeoElBC gbc(zeroflux,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                zeroflux2 = gbc.CreatedElement();
            }
            
            
           
          
             zeroflux.Element()->ResetReference();
             int64_t index;
             cmesh->SetDefaultOrder(zerofluxorder);
             TPZCompEl *celflux2 = cmesh->CreateCompEl(zeroflux2.Element(), index);
             TPZCompElSide zerofluxcomp2 = zeroflux2.Reference();
             zeroflux2.Element()->ResetReference();
             cmesh->SetDefaultOrder(1);
             TPZInterpolationSpace *HDivCollapsed = 0;
             if(zeroflux.Element()->Type() == ETriangle)
             {
                 int64_t index;
                 HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeTriang>(*cmesh,gel,index);
             }
             else if(zeroflux.Element()->Type() == EQuadrilateral)
             {
                 int64_t index;
                 HDivCollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,gel,index);
             }
             int nconnects = HDivCollapsed->NConnects();
             // bottom connect index (which will be substituted)
             int64_t cindex1 = HDivCollapsed->ConnectIndex(nconnects-2);
             // top connect index (which will be substituted)
             int64_t cindex2 = HDivCollapsed->ConnectIndex(nconnects-1);
             int nsides = gel->NSides();
             // changing orientation of top
             HDivCollapsed->SetSideOrient(nsides, -1);
             // linking bottom to first zero flux element - interface left - zeroflux - frac pressure
             HDivCollapsed->SetConnectIndex(nconnects-2, zerofluxcomp.Element()->ConnectIndex(0));
             // linking top to second zero flux element - interface right - zeroflux - frac pressure
             HDivCollapsed->SetConnectIndex(nconnects-1, zerofluxcomp2.Element()->ConnectIndex(0));
             cmesh->ConnectVec()[cindex1].DecrementElConnected();
             cmesh->ConnectVec()[cindex2].DecrementElConnected();
             cmesh->ConnectVec()[cindex1].SetSequenceNumber(-1);
             cmesh->ConnectVec()[cindex2].SetSequenceNumber(-1);
             cmesh->ConnectVec().SetFree(cindex1);
             cmesh->ConnectVec().SetFree(cindex2);
        }
    }
    
//    TPZNullMaterial* volume = new TPZNullMaterial(-11,1,1);
//    cmesh->InsertMaterialObject(volume);
//    bcmatids.insert(-11);
    
    cmesh->CleanUpUnconnectedNodes();
    // insert the fracture hdiv elements
    cmesh->SetDefaultOrder(1);
    cmesh->SetDimModel(dimension-1);
    cmesh->SetAllCreateFunctionsHDiv();
    // add the boundary HDiv elements
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild(bcmatids);
    cmesh->SetDimModel(dimension);
    cmesh->ExpandSolution();
#ifdef PZDEBUG
    {
        {
//            std::ofstream out("FluxMortarMesh.txt");
//            cmesh->Print(out);
        }
        {
//            std::ofstream out("FluxMortar.vtk");
//            TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
        }
        {
//            std::ofstream out("GMeshAfterFluxMortar.vtk");
//            TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out);
        }

        cmesh->LoadReferences();
        int64_t nelgeo = mGeometry->NElements();
        std::set<int> matidsloc(matids);
        matidsloc.insert(bcmatids.begin(),bcmatids.end());
        for(int64_t el = 0; el<nelgeo; el++)
        {
            TPZGeoEl *gel = mGeometry->Element(el);
            int matid = gel->MaterialId();
            if(matidsloc.find(matid) != matidsloc.end())
            {
                TPZCompEl *cel = gel->Reference();
                if(!cel) DebugStop();
            }
        }
    }
#endif
    return cmesh;
}

/// create a pressure with mortar elements
TPZCompMesh *TMRSApproxSpaceGenerator::PressureMortarCmesh(char firstlagrangepressure,char lagrangepressure, char lagrangemortar)
{
    std::set<int> matids, bcmatids;
    int dimension = mGeometry->Dimension();
    int nstate = 1;
    GetMaterialIds(dimension, matids, bcmatids);
    // GetMaterialIds will accumulate in the matids set!
    GetMaterialIds(dimension-1,matids,bcmatids);
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    cmesh->SetName("PressureWithMortar.txt");
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    for (auto it:matids) {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(it,dimension,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    cmesh->AutoBuild(matids);
    int64_t nconnects = cmesh->NConnects();
    for(int64_t ic = 0; ic<nconnects; ic++)
    {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(lagrangepressure);
    }
    // one connect for each element should be set at lower level
    {
        int64_t nel = cmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            cel->Connect(0).SetLagrangeMultiplier(firstlagrangepressure);
        }
    }
    {
        TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial(mSimData.mTGeometry.m_MortarMatId,dimension,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    // create discontinous elements of dimension-1
    cmesh->SetDimModel(dimension-1);
    std::set<int> mortarids = {mSimData.mTGeometry.m_MortarMatId};
    int borderOrder = mSimData.mTNumerics.m_BorderElementPresOrder;
    cmesh->SetDefaultOrder(borderOrder);
    if(borderOrder == 0)
    {
        cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    cmesh->AutoBuild(mortarids);
    int64_t ncon_new = cmesh->NConnects();
    for (int64_t ic = nconnects; ic<ncon_new; ic++) {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(lagrangemortar);
    }
    cmesh->SetDimModel(dimension);
    
    // insert the fracture hdiv elements
    cmesh->SetDimModel(dimension-1);
    GetMaterialIds(dimension-1, matids, bcmatids);
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild(matids);
    cmesh->SetDimModel(dimension);

#ifdef PZDEBUG
    {
        std::ofstream out("PressureMortarCMesh.txt");
        cmesh->Print(out);
    }
    {
        std::ofstream out("PressureMortarCMesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    }

#endif
    return cmesh;
}


TPZCompMesh * TMRSApproxSpaceGenerator::DiscontinuousCmesh(int order, char lagrange){
 
    if (!mGeometry) {
        DebugStop();
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZNullMaterial<STATE> * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(order);
    for(int dim = 0; dim <= dimension; dim++)
    {
        std::set<int> matids, bcmatids;
        GetMaterialIds(dim, matids, bcmatids);
        if(matids.size() == 0) continue;
        int nstate = 1;
        for (auto material_id : matids) {
            volume = new TPZNullMaterial(material_id,dim,nstate);
            cmesh->InsertMaterialObject(volume);
        }
    
        // PHIL : as malhas de contorno precisam objetos de condicao de contorno?
        cmesh->SetDimModel(dim);
        if(order > 0)
        {
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh->SetAllCreateFunctionsDiscontinuous();
        }
        cmesh->AutoBuild(matids);
        cmesh->InitializeBlock();
    }
    if(lagrange > 0)
    {
        int64_t ncon = cmesh->NConnects();
        //Set Lagrange multiplier
        for(int64_t i=0; i<ncon; i++){
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(lagrange);
        }
    }
    cmesh->SetDimModel(mGeometry->Dimension());
    return cmesh;
}

TPZCompMesh * TMRSApproxSpaceGenerator::TransportCmesh(){
 
    if (!mGeometry) {
        DebugStop();
    }
    
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZNullMaterial<STATE> * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(0);
    
//    for(int dim = 0; dim <= dimension; dim++)
//    {
//        std::set<int> matids, bcmatids;
//        GetMaterialIds(dim, matids, bcmatids);
//
//        int nstate = 1;
//        for(auto material_id : matids)
//        {
//            volume = new TPZNullMaterial(material_id,dim,nstate);
//            cmesh->InsertMaterialObject(volume);
//        }
//        for(auto material_id : bcmatids)
//        {
//            volume = new TPZNullMaterial(material_id,dim-1,nstate);
//            cmesh->InsertMaterialObject(volume);
//        }
//        std::set<int> allmat(matids);
//        allmat.insert(bcmatids.begin(),bcmatids.end());
////        if(allmat.size() == 0) continue;
//        cmesh->SetDimModel(dim);
//        cmesh->SetAllCreateFunctionsDiscontinuous();
//        cmesh->AutoBuild(allmat);
//    }
    //
//    TPZTracerFlow * volume = nullptr;
//    cmesh->SetDefaultOrder(0);
    std::set<int> volIds;
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int material_id = chunk.second;
            //            volume = new TMRSMultiphaseFlow<TMRSMemory>(material_id,d);
            //            volume->SetDataTransfer(mSimData);
            
            volume = new TPZNullMaterial(material_id,d,1);
            volIds.insert(material_id);
            //            volume->SetDataTransfer(mSimData);
            
            cmesh->InsertMaterialObject(volume);
        }
    }
    //            volume->SetDataTransfer(mSimData);
    volume = new TPZNullMaterial(10,2,1);
    cmesh->InsertMaterialObject(volume);
    volIds.insert(10);
    
   
   
    
    
    if (!volume) {
        DebugStop();
    }
    
    std::set<int> boundaryId;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]   = get<2>(chunk);
        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        boundaryId.insert(bc_id);
        cmesh->InsertMaterialObject(face);
    }

    
    //
    cmesh->InitializeBlock();
    cmesh->SetAllCreateFunctionsDiscontinuous();
    int nels = mGeometry->NElements();
    for (auto gel:mGeometry->ElementVec()) {
       int gelId = gel->MaterialId();
        for (auto Mat_id: volIds) {
            if (gelId == Mat_id) {
                CreateTransportElement(0, cmesh, gel, false);
            }
        }
        for (auto Mat_id: boundaryId) {
            if (gelId == Mat_id) {
                CreateTransportElement(0, cmesh, gel, true);
            }
        }
    }
    
    cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int s_order = 0;
    cmesh->SetDefaultOrder(s_order);
    cmesh->ExpandSolution();
//    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "s_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    std::ofstream file("transportAtomic.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file);
    return cmesh;
}

TPZCompMesh * TMRSApproxSpaceGenerator::DiscontinuousCmesh(TPZAlgebraicDataTransfer &Atransfer){
 
    if (!mGeometry) {
        DebugStop();
    }
    int order = 0;
    TPZCompMesh *cmesh = new TPZCompMesh(mGeometry);
    TPZL2Projection<STATE> * volume = nullptr;
    int dimension = mGeometry->Dimension();
    cmesh->SetDefaultOrder(order);
    
    std::vector<std::map<std::string,int>> &DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
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
    
    // PHIL : as malhas de contorno precisam objetos de condicao de contorno?
    
    cmesh->SetAllCreateFunctionsDiscontinuous();

    // The index of the computational volume elements in the transport mesh identified by material id
//    std::map<int,TPZVec<int64_t>> fVolumeElements;
    int64_t vol_index = 0;
    TPZVec<int64_t> &elindices = Atransfer.fVolumeElements;
    for (int64_t i=0; i<elindices.size(); i++) {
        int64_t el = elindices[i];
        TPZCompEl *cel = mTransportOperator->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        int dim = gel->Dimension();
        int matid = gel->MaterialId();
        if(cmesh->FindMaterial(matid) == 0) continue;
        int64_t celindex;
        cmesh->CreateCompEl(gel, celindex);
        vol_index++;
    }
    cmesh->InitializeBlock();
    

    
#ifdef PZDEBUG
//    std::stringstream file_name;
//    if (order == 0) {
//        file_name << "s_cmesh" << ".txt";
//    }
//    else{
//        file_name << "p_cmesh" << ".txt";
//    }
//    std::ofstream sout(file_name.str().c_str());
//    cmesh->Print(sout);
#endif
    
    return cmesh;
}


void TMRSApproxSpaceGenerator::BuildMixedMultiPhysicsCompMesh(int order){
    
    bool cond1 = mSimData.mTNumerics.m_four_approx_spaces_Q;
    bool cond2 = mSimData.mTNumerics.m_mhm_mixed_Q;
    
    switch(mSimData.mTNumerics.m_SpaceType)
    {
        case TMRSDataTransfer::TNumerics::E2Space:
            if(cond1) DebugStop();
            if(cond2) DebugStop();
            BuildMixed2SpacesMultiPhysicsCompMesh(order);
            break;
        case TMRSDataTransfer::TNumerics::E4Space:
//            if(!cond1) DebugStop();
//            if(!cond2) DebugStop();
            BuildMixed4SpacesMultiPhysicsCompMesh(order);
            break;
        case TMRSDataTransfer::TNumerics::E2SpaceMHM:
            if(!cond2) DebugStop();
            if(cond1) DebugStop();
            BuildMHMMixed2SpacesMultiPhysicsCompMesh();
            break;
        case TMRSDataTransfer::TNumerics::E4SpaceMHM:
            if(!cond1) DebugStop();
            if(!cond2) DebugStop();
            BuildMHMMixed4SpacesMultiPhysicsCompMesh();
            break;
        case TMRSDataTransfer::TNumerics::E4SpaceMortar:
            if(!cond1) DebugStop();
//            if(cond2) DebugStop();
            BuildMixed4SpacesMortarMesh();
            break;
        default:
            DebugStop();
    }
    
//    std::string name_ref = "mhm_geo";
//    PrintGeometry(name_ref);
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
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]  = get<2>(chunk);
        TPZBndCondT<REAL> * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        if (HasForcingFunctionBC()) {
            face->SetForcingFunctionBC(mForcingFunctionBC);
        }
        mMixedOperator->InsertMaterialObject(face);
    }

    
    TMRSDarcyFractureFlowWithMem<TMRSMemory> * fracmat = nullptr;
    for(auto chunk : mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag[dimension-1])
    {
        std::string material_name = chunk.first;
        std::cout << "physical name = " << material_name <<
        " material id " << chunk.second << " dimension " << dimension-1 << std::endl;
        int materia_id = chunk.second;
//        matsWithMem.insert(materia_id);
        fracmat = new TMRSDarcyFractureFlowWithMem<TMRSMemory>(materia_id,dimension-1);
        TMRSMemory defaultmem;
        fracmat->SetDataTransfer(mSimData);
        mMixedOperator->InsertMaterialObject(fracmat);

    }
    
    TPZManVector<TPZCompMesh *, 3> mesh_vec(3);
    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order,1);
    mesh_vec[2] = TransportCmesh();
    TPZManVector<int,5> active_approx_spaces(3);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 0;
    mMixedOperator->SetDimModel(dimension);
    
//    DeleteBCsThatAreOnIntersect(mesh_vec[0]);
    
    if (isThereFracIntersection()) {
        mHybridizer = new TPZHybridizeHDiv(mesh_vec);
        HybridizeIntersections(mesh_vec);
        mHybridizer->InsertPeriferalMaterialObjects(mMixedOperator);
    }
    
//    mMixedOperator->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec);

    // Creates interface elements in case there is hybridization for fracture intersection
    if (mHybridizer){
        CreateIntersectionInterfaceElements(mesh_vec);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "mixed_cmesh_two_spaces" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    mMixedOperator->Print(sout);
#endif
    
}

void TMRSApproxSpaceGenerator::DeleteBCsThatAreOnIntersect(TPZCompMesh* hdivcmesh) {
    hdivcmesh->LoadReferences();
    
    for (auto cel : hdivcmesh->ElementVec()) {
        if (!cel)
            continue;
        TPZGeoEl *gel = cel->Reference();
        if (gel->MaterialId() != 5)
            continue;
        
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.HasNeighbour(mMatIDIntersection);
        if (neigh) {
            delete cel;
        }
    
    }
    hdivcmesh->InitializeBlock();
    
}

void TMRSApproxSpaceGenerator::BuildMixed4SpacesMortarMesh(){
    int dimension = mGeometry->Dimension();
    std::cout << __PRETTY_FUNCTION__ << " on input nel geom " << mGeometry->NElements() << std::endl;
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    std::set<int> matsWithMem, matsWithOutMem;
    GetMaterialIds(dimension, matsWithMem, matsWithOutMem);
    
    InsertGeoWrappersForMortar();
    TPZManVector<TPZCompMesh *> meshvec(4,0);
    // hdiv mesh
    char fluxmortar = 5;
    meshvec[0] = HDivMortarFluxCmesh(fluxmortar);
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> fatherAndSons;
        TPZReservoirTools::TakeFatherSonsCorrespondence(meshvec[0], this->mSubdomainIndexGel, fatherAndSons);
        TPZReservoirTools::PutFluxElementsinSubdomain(meshvec[0], this->mSubdomainIndexGel, fatherAndSons);
        TPZReservoirTools::AddDependency(fatherAndSons);
    }
    // pressure mesh
    char firstpressurelagrange = 1;
    char pressurelagrange = 3;
    char pressuremortar = 4;
    meshvec[1] = PressureMortarCmesh(firstpressurelagrange,pressurelagrange,pressuremortar);
    
    std::ofstream file2("PressureCmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(meshvec[1], file2);
    
    // distributed flux mesh
    int porder = 0;
    char distfluxlagrange = 2;
    meshvec[2] = DiscontinuousCmesh(porder,distfluxlagrange);
    // constant pressure mesh
    char avpressurelagrange = 6;
    meshvec[3] = DiscontinuousCmesh(porder,avpressurelagrange);
    // transport mesh
    // @TODO I believe we don't need a transport mesh anymore
//    meshvec[4] = TransportCmesh();
    
    // create the multiphysics mesh
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
    //    TPZMixedDarcyFlow *volume = nullptr;
    mMixedOperator->SetDefaultOrder(1);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    std::cout << "Creating material objects\n";
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name <<
            " material id " << chunk.second << " dimension " << d << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSDarcyFlowWithMem<TMRSMemory>(materia_id,d);
            TMRSMemory defaultmem;
            // neste ponto podemos inserir as propriedades de permeabilidade absoluta
            volume->SetDefaultMem(defaultmem);
//            volume = new TPZMixedDarcyFlow(materia_id, d);
//             volume->SetPermeability(1.0);
            volume->SetDataTransfer(mSimData);
            mMixedOperator->InsertMaterialObject(volume);
        }
    }
    TMRSDarcyFractureFlowWithMem<TMRSMemory> * fracmat = nullptr;
    for(auto chunk : mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag[dimension-1])
    {
        std::string material_name = chunk.first;
        std::cout << "physical name = " << material_name <<
        " material id " << chunk.second << " dimension " << dimension-1 << std::endl;
        int materia_id = chunk.second;
        matsWithMem.insert(materia_id);
        fracmat = new TMRSDarcyFractureFlowWithMem<TMRSMemory>(materia_id,dimension-1);
        TMRSMemory defaultmem;
        // neste ponto podemos inserir as propriedades de permeabilidade absoluta
//        defaultmem.m_kappa = 0.00000001;
//        defaultmem.m_kappa_inv = 1.0/0.00000001;
//        defaultmem.SetPermeability(0.00000001);
//        defaultmem.m_p=100000000.0;
//        fracmat->SetDefaultMem(defaultmem);
//            volume = new TPZMixedDarcyFlow(materia_id, d);
//             volume->SetPermeability(1.0);
        fracmat->SetDataTransfer(mSimData);
        mMixedOperator->InsertMaterialObject(fracmat);

    }
    if (!volume) {
        DebugStop();
    }
    
    {
        TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
        TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
        for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
            int bc_id   = get<0>(chunk);
            int bc_type = get<1>(chunk);
            val2[0]  = get<2>(chunk);
            TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
            mMixedOperator->InsertMaterialObject(face);
        }
    }
    
    {
        TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
        TPZManVector<std::tuple<int, int, REAL>> &BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue;
        for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
            if(!fracmat) DebugStop();
            int bc_id   = get<0>(chunk);
            int bc_type = get<1>(chunk);
            matsWithOutMem.insert(bc_id);
            val2[0]  = get<2>(chunk);
            TPZBndCond * face = fracmat->CreateBC(volume,bc_id,bc_type,val1,val2);
            mMixedOperator->InsertMaterialObject(face);
        }
    }

    {
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<STATE> *nullmat = new TPZNullMaterialCS(mSimData.mTGeometry.m_HdivWrapMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
    {
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<STATE> *nullmat = new TPZNullMaterialCS(mSimData.mTGeometry.m_MortarMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
    {
        int dim = 1;
        int nstate = 1;
        TPZNullMaterialCS<STATE> *nullmat = new TPZNullMaterialCS(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId,dim,nstate);
        mMixedOperator->InsertMaterialObject(nullmat);
    }
//    mGeometry->ResetReference();
//    mMixedOperator->LoadReferences();
    
//    TMRSDarcyFlowWithMem<TMRSMemory>*volume2 = new TMRSDarcyFlowWithMem<TMRSMemory>(19,2);
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        TPZNullMaterialCS<STATE> *volume2 = new TPZNullMaterialCS(mSimData.mTGeometry.m_skeletonMatId,2,1);
        mMixedOperator->InsertMaterialObject(volume2);
    }
   
    
    matsWithOutMem.insert(mSimData.mTGeometry.m_HdivWrapMatId);
    matsWithOutMem.insert(mSimData.mTGeometry.m_MortarMatId);
    matsWithOutMem.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    matsWithOutMem.insert(mSimData.mTGeometry.m_skeletonMatId);

    TPZManVector<int> active_approx_spaces(4,1);
//    active_approx_spaces[4] = 0;
    mMixedOperator->SetDimModel(3);
    gSinglePointMemory = true;
    
    // NS to Jose: Should this method with these arguments be commited in PZ???
    //    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,meshvec,matsWithMem,
//    DebugStop(); // look up
    // NS to Jose: Using this for now. Erase later
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,meshvec);

    //Insert fractures properties
    InitializeFracProperties(mMixedOperator);
    // insert the interface elements
    InsertInterfaceElements();
#ifdef PZDEBUG
    {
        std::stringstream file_name;
        file_name  << "mixed_cmesh_four_space_mortar" << ".txt";
        std::ofstream sout(file_name.str().c_str());
        mMixedOperator->Print(sout);
    }
#endif
    mMixedOperator->ComputeNodElCon();


    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        if(mSimData.mTNumerics.m_UseSubstructures_Q)
        {
            std::cout<<"Num Eq Mixed: "<<mMixedOperator->NEquations()<<std::endl;
            HideTheElements(mMixedOperator);
            int nels = mMixedOperator->NElements();
           
            
            for(int iel =0; iel<nels; iel++){
                TPZCompEl *cel = mMixedOperator->Element(iel);
                if(!cel){continue;}
                TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
                if(subcmesh){
                
                    std::set<int64_t> seed, groups;
                    subcmesh->ComputeNodElCon();
                  
                    
#ifdef LOG4CXX
                    if(logger->isDebugEnabled())
                    {
                        std::stringstream sout;
                        subcmesh->Print(sout);
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    
                    TPZReservoirTools::TakeSeedElements(subcmesh, seed);
                    TPZReservoirTools::GroupNeighbourElements(subcmesh,seed,groups );
                    subcmesh->ComputeNodElCon();
                    std::set<int> volmatId;
                    
                    int nel = subcmesh->NElements();
                    TPZReservoirTools::CondenseElements(subcmesh, pressuremortar, false,volmatId);
                    std::set<int64_t> groups2;
                    // this will act only on volumetric elements
                    TPZReservoirTools::GroupNeighbourElements(subcmesh, groups, groups2);
                    subcmesh->ComputeNodElCon();
                    // this shouldn't affect the fracture elements as they won't have condensable connects
//                    int nel = subcmesh->NElements();
                    for (int64_t el=0; el<nel; el++) {
                        TPZCompEl *cel = subcmesh->Element(el);
                        if (!cel) {
                            continue;
                        }
                        int nconnects = cel->NConnects();
                        for (int icon=0; icon<nconnects; icon++) {
                            TPZConnect &connect = cel->Connect(icon);
                            
                            int lagrangemult = connect.LagrangeMultiplier();
                            //Increment the number of connected elements for the avg pressure in order to not condense them
                            if (lagrangemult==6) {
                                connect.IncrementElConnected();
                            }
                        }
                    }
                    TPZReservoirTools::CondenseElements(subcmesh, fluxmortar, false);
                    subcmesh->ComputeNodElCon();
                    int numThreads =0;
                    int preconditioned =0;
                    TPZAutoPointer<TPZGuiInterface> guiInterface;
                    subcmesh->SetAnalysisSkyline(numThreads, preconditioned, guiInterface);
                    TPZStepSolver<STATE> step;
                    step.SetDirect(ELDLt);
                    subcmesh->Analysis()->SetSolver(step);
                    
                    
                }
            }
            mMixedOperator->ComputeNodElCon();
            std::cout<<"Num Eq Mixed MHM: "<<mMixedOperator->NEquations()<<std::endl;

        }
        else{
            // group and condense the H(div) space (only dimension of the mesh)
            std::set<int64_t> seed, groups;
            int64_t nel = mMixedOperator->NElements();
            int dim = mMixedOperator->Dimension();
            for (int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = mMixedOperator->Element(el);
                if(!cel){
                    continue;
                }
                TPZGeoEl *gel = cel->Reference();
                if(gel->Dimension() == dim) seed.insert(el);
            }
            // this will only group volumetric elements
            TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, seed, groups);
//            mMixedOperator->ComputeNodElCon();
            
            std::set<int> volmatId;
            volmatId.insert(10);
            TPZReservoirTools::CondenseElements(mMixedOperator, pressuremortar, false,volmatId);
            
            
            std::set<int64_t> groups2;
            
            // this will act only on volumetric elements
            TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, groups, groups2);
//            mMixedOperator->ComputeNodElCon();
            // this shouldn't affect the fracture elements as they won't have condensable connects
            TPZReservoirTools::CondenseElements(mMixedOperator, fluxmortar, false);
            mMixedOperator->ComputeNodElCon();
        }

    }
    else
    {
    
        // group and condense the H(div) space (only dimension of the mesh)
        std::set<int64_t> seed, groups;
        int64_t nel = mMixedOperator->NElements();
        int dim = mMixedOperator->Dimension();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = mMixedOperator->Element(el);
            if(!cel){
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() == dim) seed.insert(el);
        }
        // this will only group volumetric elements
        TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, seed, groups);
        mMixedOperator->ComputeNodElCon();
        
        std::set<int> volmatId;
        volmatId.insert(10);
        TPZReservoirTools::CondenseElements(mMixedOperator, pressuremortar, false,volmatId);
        

        std::set<int64_t> groups2;
        
        // this will act only on volumetric elements
        
  
        TPZCompMeshTools::GroupNeighbourElements(mMixedOperator, groups, groups2);
        mMixedOperator->ComputeNodElCon();
        // this shouldn't affect the fracture elements as they won't have condensable connects
        TPZReservoirTools::CondenseElements(mMixedOperator, fluxmortar, false);
    }

    
  
#ifdef PZDEBUG
    {
        std::stringstream file_name;
        file_name  << "mixed_cmesh_four_space_mortar_two_condense" << ".txt";
        std::ofstream sout(file_name.str().c_str());
        mMixedOperator->Print(sout);
    }
#endif
    
}

/// insert the necessary interface elements
void TMRSApproxSpaceGenerator::InsertInterfaceElements()
{
    int dim = mGeometry->Dimension();
    std::set<int> matids, bcmatids, fracmatids, fracbcmatids;
    GetMaterialIds(dim, matids, bcmatids);
    GetMaterialIds(dim-1, fracmatids, fracbcmatids);
//    bcmatids.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    TPZLagrangeMultiplierCS<STATE> *mat1 = new TPZLagrangeMultiplierCS<STATE>(mSimData.mTGeometry.m_posLagrangeMatId,dim,1);
    mMixedOperator->InsertMaterialObject(mat1);
    TPZLagrangeMultiplierCS<STATE> *mat2 = new TPZLagrangeMultiplierCS<STATE>(mSimData.mTGeometry.m_negLagrangeMatId,dim,1);
    mat2->SetMultiplier(-1.);
    mMixedOperator->InsertMaterialObject(mat2);
    mGeometry->ResetReference();
    mMixedOperator->LoadReferences();
    int64_t nel = mMixedOperator->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = mMixedOperator->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        // the interface will be generated between HDiv wrapper and the mortar pressure
        if(matid != mSimData.mTGeometry.m_HdivWrapMatId && matid != mSimData.mTGeometry.m_MortarMatId)
            continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide leftgel(gelside), rightgel;
        TPZGeoElSide IntfaceSide = leftgel.Neighbour();
        if(matid == mSimData.mTGeometry.m_HdivWrapMatId)
        {
            rightgel = IntfaceSide.Neighbour();
        }
        else if(matid == mSimData.mTGeometry.m_MortarMatId)
        {
            // from a mortar element we need to add an interface to either a boundary flux
            // or a zero order flux
            // if the mortar element has a fracture neighbour, then there are two neighbouring
            // zero order flux elements. Depending on the interface matid we need to connect to either
            //  zero order flux element
            int domain =0;
            if(mSimData.mTNumerics.m_mhm_mixed_Q){
                domain = mSubdomainIndexGel[gel->Index()];
                if(domain == -1) DebugStop();
            }
         
            TPZGeoElSide BCGelside = gelside.HasNeighbour(bcmatids);
            TPZGeoElSide FracGelside = gelside.HasNeighbour(fracmatids);
            TPZGeoElSide Zerofluxside = gelside.HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
            
            int intfacematid = IntfaceSide.Element()->MaterialId();
            if(BCGelside)
            {
                rightgel = BCGelside;
            }
            else if(Zerofluxside)
            {
                if(mSimData.mTNumerics.m_mhm_mixed_Q){
                    int64_t zerodomain = mSubdomainIndexGel[Zerofluxside.Element()->Index()];
                    if(zerodomain != domain)
                    {
                        Zerofluxside = Zerofluxside.Neighbour().HasNeighbour(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
                        if(!Zerofluxside) DebugStop();
                        zerodomain = mSubdomainIndexGel[Zerofluxside.Element()->Index()];
                        if(zerodomain != domain) DebugStop();
                    }
                    // if the lagrange multiplier is positive, the this is the top lagrange multiplier
                    // the top flux element is the second zero flux element
                    if(FracGelside && intfacematid == mSimData.mTGeometry.m_posLagrangeMatId)
                    {
                        Zerofluxside = Zerofluxside.Neighbour();
                        if(Zerofluxside.Element()->MaterialId() != mSimData.mTGeometry.m_zeroOrderHdivFluxMatId)
                        {
                            DebugStop();
                        }
                    }
                    rightgel = Zerofluxside;
                }
                else{
                    
                    if(FracGelside && intfacematid == mSimData.mTGeometry.m_posLagrangeMatId)
                    {
                        Zerofluxside = Zerofluxside.Neighbour();
                        if(Zerofluxside.Element()->MaterialId() != mSimData.mTGeometry.m_zeroOrderHdivFluxMatId)
                        {
                            DebugStop();
                        }
                    }
                    rightgel = Zerofluxside;
                }
               
            }
            else
            {
                DebugStop();
            }
        }
        TPZCompElSide leftcel = leftgel.Reference();
        TPZCompElSide rightcel = rightgel.Reference();
        int64_t index;
        TPZMultiphysicsInterfaceElement *intface =
        new TPZMultiphysicsInterfaceElement(*mMixedOperator,IntfaceSide.Element(),index,leftcel,rightcel);
    }
}

/// insert wrapper elements necessary for creating the (hybridized) mortar spaces
void TMRSApproxSpaceGenerator::InsertGeoWrappersForMortar()
{
    int dimension = mGeometry->Dimension();
    std::set<int> matids, bcmatids;
    GetMaterialIds(dimension, matids, bcmatids);
    bcmatids.insert(mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
    // around each volumetric element create an
    // - hdiv wrap element
    // - interface geo element
    // - pressure element
    // - second interface element
    // - a flux hdiv boundary element (conditionally)
    int64_t nel = mGeometry->NElements();
    int nElSubDomain = mSubdomainIndexGel.size();
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        if(nel != nElSubDomain){
            DebugStop();
        }
    }
    
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = mGeometry->Element(el);
        if(!gel || gel->Dimension()!= dimension || gel->HasSubElement()) continue;
        int firstside = 0;
        for(int dim = 0; dim < dimension-1; dim++)
        {
            firstside+=gel->NSides(dim);
        }
        for (int side=firstside; side < gel->NSides()-1; side++) {
            TPZGeoElSide gelside(gel,side);
            GeoWrappersForMortarGelSide(gelside,bcmatids);
        }
    }
    
    nel = mGeometry->NElements();
    nElSubDomain = mSubdomainIndexGel.size();
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        if(nel != nElSubDomain){
            DebugStop();
        }
        std::ofstream out("gmesh_withwrap.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out, mSubdomainIndexGel);
    }
    else{
        std::ofstream out("gmesh_withwrap.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, out);
    }
    

}
void TMRSApproxSpaceGenerator::GeoWrappersForMortarGelSide(TPZGeoElSide &gelside, std::set<int> bcmatids){
   
    TPZGeoEl *gel = gelside.Element();
    int side = gelside.Side();
    int gelindex = gel->Index();
    int hdiv_orient = gel->NormalOrientation(side);
    int subDomainIndex =-1;
    int subDomainIndexNeig =-1;
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
         subDomainIndex = mSubdomainIndexGel[gelindex];
         subDomainIndexNeig = FindNeighSubDomain(gelside);
        if(subDomainIndex<0){
            DebugStop();
        }
    }
    int first_lagrange = mSimData.mTGeometry.m_posLagrangeMatId;
    int second_lagrange = mSimData.mTGeometry.m_negLagrangeMatId;
    bool cond1 = hdiv_orient < 0 ;
    bool cond2 = false;
    
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        cond2 =subDomainIndexNeig!=-1 && (subDomainIndexNeig > subDomainIndex);
    }
    
    if(cond1 || cond2)
    {
        first_lagrange = mSimData.mTGeometry.m_negLagrangeMatId;
        second_lagrange = mSimData.mTGeometry.m_posLagrangeMatId;
    }

    
    TPZGeoElBC gbc1(gelside,mSimData.mTGeometry.m_HdivWrapMatId);
    TPZGeoElSide gelwrapside(gbc1.CreatedElement());
    int nBCCreated = 4;
    int index1 =gbc1.CreatedElement()->Index();
    
    TPZGeoElBC gbc2(gelwrapside,first_lagrange);
    TPZGeoElSide gelintface1(gbc2.CreatedElement());
    int index2 =gbc2.CreatedElement()->Index();
    
    
    TPZGeoElBC gbc3(gelintface1,mSimData.mTGeometry.m_MortarMatId);
    TPZGeoElSide gelmortar(gbc3.CreatedElement());
    int index3 =gbc3.CreatedElement()->Index();
    
    TPZGeoElBC gbc4(gelmortar,second_lagrange);
    TPZGeoElSide gelinterface2(gbc4.CreatedElement());
    int index4 =gbc4.CreatedElement()->Index();
   
    int size =mSubdomainIndexGel.size();
    if(mSimData.mTNumerics.m_mhm_mixed_Q){
        mSubdomainIndexGel.resize(size + nBCCreated);
        if(index1 != size || index2!=(size+1) || index3!=(size+2) || index4!=(size+3)){
            DebugStop();
        }
        mSubdomainIndexGel[index1]=subDomainIndex;
        mSubdomainIndexGel[index2]=subDomainIndex;
        mSubdomainIndexGel[index3]=subDomainIndex;
        mSubdomainIndexGel[index4]=subDomainIndex;
    }
    
    
    // this method has to be adjusted if we create MHM meshes
    if(!gelinterface2.HasNeighbour(bcmatids)){
      
        TPZGeoElBC gbc5(gelinterface2,mSimData.mTGeometry.m_zeroOrderHdivFluxMatId);
        if(mSimData.mTNumerics.m_mhm_mixed_Q){
            mSubdomainIndexGel.resize(size + nBCCreated +1);
            mSubdomainIndexGel[size + nBCCreated]=subDomainIndex;
        }
        
    }
 
}

int TMRSApproxSpaceGenerator::FindNeighSubDomain(TPZGeoElSide &gelside){
    TPZGeoEl *gel = gelside.Element();
    int indexGel = gel->Index();
    int indexSubDomain = mSubdomainIndexGel[indexGel];
    int subdomainIndex =-1;
    TPZGeoElSide neigh = gelside.Neighbour();
    while(neigh!= gelside){
        int neighIndex=neigh.Element()->Index();
        if(neigh.Element()->Dimension() != gel->Dimension()){
            neigh = neigh.Neighbour();
            continue;
        }
        subdomainIndex = mSubdomainIndexGel[neighIndex];
        break;
    }
    
    return subdomainIndex;
}

/// return the material ids and boundary condition material ids
void TMRSApproxSpaceGenerator::GetMaterialIds(int dim, std::set<int> &matids, std::set<int> &bcmatids)
{
#ifdef PZDEBUG
    std::cout << "GetMaterialIds - Identifying material objects for dimension " << dim << std::endl;
#endif
    if(dim == mGeometry->Dimension())
    {
        std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
        for (auto chunk : DomainDimNameAndPhysicalTag[dim]) {
#ifdef PZDEBUG
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << " matid " << chunk.second<< std::endl;
#endif
            matids.insert(chunk.second);
        }
        TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
        for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
            int bc_id   = get<0>(chunk);
#ifdef PZDEBUG
            std::cout << "boundary condition matid " << bc_id << std::endl;
#endif
            bcmatids.insert(bc_id);
        }
    } else
    {
        std::vector<std::map<std::string,int>> &DomainDimNameAndPhysicalTag =
        mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag;
        for (auto chunk : DomainDimNameAndPhysicalTag[dim]) {
#ifdef PZDEBUG
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << " matid " << chunk.second<< std::endl;
#endif
            matids.insert(chunk.second);
        }
        if(dim == mGeometry->Dimension()-1)
        {
            TPZManVector<std::tuple<int, int, REAL>> &BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue;
            for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
                int bc_id   = get<0>(chunk);
#ifdef PZDEBUG
                std::cout << "boundary condition matid " << bc_id << std::endl;
#endif
                bcmatids.insert(bc_id);
            }
        }

    }
}


void TMRSApproxSpaceGenerator::BuildMixed4SpacesMultiPhysicsCompMesh(int order){
    
    int dimension = mGeometry->Dimension();
    mMixedOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
//    TPZDarcyFlowWithMem *volume = nullptr;
    mMixedOperator->SetDefaultOrder(order);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int material_id = chunk.second;
            volume = new TMRSDarcyFlowWithMem<TMRSMemory>(material_id,d);
            volume->SetDataTransfer(mSimData);
            
//                volume = new TPZDarcyFlowWithMem(material_id, d);
//                volume->SetPermeability(1.0);
            mMixedOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        REAL val = get<2>(chunk);
        val2[0] =val;
        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mMixedOperator->InsertMaterialObject(face);
    }
    
    // PHIL : Vamos poder criar apenas 4 espacos...
  
    
    TPZManVector<TPZCompMesh *, 5> mesh_vec(5);
    mesh_vec[0] = HdivFluxCmesh(order);
    mesh_vec[1] = DiscontinuousCmesh(order,0);
    mesh_vec[2] = DiscontinuousCmesh(0,0);
    mesh_vec[3] = DiscontinuousCmesh(0,0);
    mesh_vec[4] = TransportCmesh();
    
    int ncon = mesh_vec[1]->NConnects();
    //Set Lagrange multiplier
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = mesh_vec[1]->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
     ncon = mesh_vec[2]->NConnects();
    //Set Lagrange multiplier
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = mesh_vec[2]->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(2);
    }
    ncon = mesh_vec[3]->NConnects();
    //Set Lagrange multiplier
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = mesh_vec[3]->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(3);
        newnod.IncrementElConnected();
    }
    
    
    TPZManVector<int,5> active_approx_spaces(5);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    active_approx_spaces[3] = 1;
    active_approx_spaces[4] = 0;
    mMixedOperator->SetDimModel(dimension);
    mMixedOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,mesh_vec);
//    mMixedOperator->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "mixed_cmesh_four_spaces" << ".txt";
    std::ofstream sout(file_name.str().c_str());
//    mMixedOperator->Print(sout);
#endif
    TPZReservoirTools::CreatedCondensedElements(mMixedOperator, false, true);
    
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
            
            TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
            TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
            for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
                int bc_id   = get<0>(chunk);
                int bc_type = get<1>(chunk);
                val2[0]  = get<2>(chunk);
                matids.insert(bc_id);
                TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
                mhm->CMesh()->InsertMaterialObject(face);
            }
            mhm->fMaterialBCIds = matids;
        }

        mhm->SetInternalPOrder(1);
        mhm->SetSkeletonPOrder(1);
        
        mhm->DivideSkeletonElements(mSimData.mTGeometry.mSkeletonDiv);
        mhm->DivideBoundarySkeletonElements();
    
        bool substructure = true;
    
        mhm->SetApproxSpaceGenerator(this);
        mhm->BuildComputationalMesh(substructure);
        
        std::cout << "MHM Hdiv Computational meshes created\n";
        
        std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;
    
    TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
    
    TPZMultiphysicsCompMesh *cmeshtest = dynamic_cast<TPZMultiphysicsCompMesh*>(MixedMesh);
    
    mMixedOperator = cmeshtest;
    
    std::cout << "WHAT WILL HAPPEN WITH THE mhm VARIABLE????\n";

}

void TMRSApproxSpaceGenerator::BuildMHMMixed4SpacesMultiPhysicsCompMesh(){
   
    int dimension = mGeometry->Dimension();
    TPZMHMixedMesh4SpacesControl *mhm = new TPZMHMixedMesh4SpacesControl(mGeometry);
    TPZVec<int64_t> coarseindices;
    ComputeCoarseIndices(mGeometry, coarseindices);
    
    mhm->DefinePartitionbyCoarseIndices(coarseindices);
    // Create geometric elements
    
    {
        std::set<int> matids;
        for (auto omId:mSimData.mTGeometry.mDomainDimNameAndPhysicalTag[dimension] ) {      //Materials
            int valor = omId.second;
            matids.insert(omId.second);
        } ;
        mhm->fMaterialIds = matids;
        matids.clear();
        
        for (auto omId:mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue ) {      //BC
         
            int MatId = get<0>(omId);
            matids.insert(MatId);
        } ;
        mhm->fMaterialBCIds = matids;
    }
    
    InsertMaterialObjects(*mhm);
    
    mhm->SetInternalPOrder(1);
    mhm->SetSkeletonPOrder(1);
    
    mhm->DivideSkeletonElements(mSimData.mTGeometry.mSkeletonDiv);
    mhm->DivideBoundarySkeletonElements();
    if (0) {
        std::ofstream file_geo("geometry.txt");
        mhm->CMesh()->Reference()->Print(file_geo);
    }
    
    bool substructure = true;
    
    mhm->BuildComputationalMesh(substructure);
    std::cout << "Number of equations MHMixed " << mhm->CMesh()->NEquations() << std::endl;
    
 
    TPZCompMesh *MixedMesh = mhm->CMesh().operator->();
    mMixedOperator = dynamic_cast<TPZMultiphysicsCompMesh *>(MixedMesh);
    std::cout << "Changing the analysis to sparse\n";
    {
        int64_t nel = MixedMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = MixedMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
//                                TPZSymetricSpStructMatrixEigen matrix(sub);
////                                TPZSymetricSpStructMatrix matrix(sub);
//                                int numinternal = sub->NumInternalEquations();
//                                matrix.EquationFilter().SetMinMaxEq(0, numinternal);
//                                TPZAutoPointer<TPZMatrix<STATE> > mat = matrix.Create();
//                                matrix.EquationFilter().Reset();
//                                matrix.SetNumThreads(0);
//                                sub->Analysis()->SetStructuralMatrix(matrix);
//                                TPZStepSolver<STATE> step;
//                                step.SetDirect(ELDLt);
//                                sub->Analysis()->SetSolver(step);
            }
        }
    }
    
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
    transport_meshvec[2] = TransportCmesh();
    
    
    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
//    TPZTracerFlow * volume = nullptr;
    
    mTransportOperator->SetDefaultOrder(0);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int materia_id = chunk.second;
            volume = new TMRSMultiphaseFlow<TMRSMemory>(materia_id,d);
            volume->SetDataTransfer(mSimData);
//               volume = new  TPZTracerFlow(materia_id,d);
            mTransportOperator->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]  = get<2>(chunk);
        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mTransportOperator->InsertMaterialObject(face);
        
    }
    
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    {
        TMRSMultiphaseFlow<TMRSMemory> * interface = new TMRSMultiphaseFlow<TMRSMemory>(transport_matid,dimension-1);
        interface->SetDataTransfer(mSimData);
        mTransportOperator->InsertMaterialObject(interface);
        
//        TPZTracerFlow * interface = new TPZTracerFlow(transport_matid,dimension-1);
////        interface->SetDataTransfer(mSimData);
//        mTransportOperator->InsertMaterialObject(interface);
    }
    
    mTransportOperator->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 1;
    mTransportOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,transport_meshvec);
//    mTransportOperator->BuildMultiphysicsSpace(active_approx_spaces,transport_meshvec);
    

    
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
    transport_meshvec[4] = TransportCmesh();
    
    
    int dimension = mGeometry->Dimension();
    mTransportOperator = new TPZMultiphysicsCompMesh(mGeometry);
    
//    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
    TPZTracerFlow * volume = nullptr;
    mTransportOperator->SetDefaultOrder(0);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int material_id = chunk.second;
//            volume = new TMRSMultiphaseFlow<TMRSMemory>(material_id,d);
//            volume->SetDataTransfer(mSimData);
            
            volume = new TPZTracerFlow(material_id,d);
//            volume->SetDataTransfer(mSimData);
            
            mTransportOperator->InsertMaterialObject(volume);
        }
    }
    volume = new TPZTracerFlow(10,2);
    mTransportOperator->InsertMaterialObject(volume);
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]  = get<2>(chunk);
        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        mTransportOperator->InsertMaterialObject(face);
        
    }
    //crear controle no simdata
    TPZBndCond * face3 = volume->CreateBC(volume,-11,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face3);
    
    int fracvol1ID = mSimData.mTGeometry.mInterface_material_idFracInf;
    int fracvol2ID = mSimData.mTGeometry.mInterface_material_idFracSup;
    int fracFracID = mSimData.mTGeometry.mInterface_material_idFracFrac;
    int fracbounId = mSimData.mTGeometry.mIterface_material_idFracBound;
    
    TPZBndCond * face = volume->CreateBC(volume,fracvol1ID,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face);

    TPZBndCond * face4 = volume->CreateBC(volume,fracvol2ID,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face4);
    
    TPZBndCond * face2 = volume->CreateBC(volume,fracFracID,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face2);
//
    TPZBndCond * face5 = volume->CreateBC(volume,fracbounId,1,val1,val2);
    mTransportOperator->InsertMaterialObject(face5);
    
    
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    {
//        TMRSMultiphaseFlow<TMRSMemory> * interface = new TMRSMultiphaseFlow<TMRSMemory>(transport_matid,dimension-1);
         TPZTracerFlow * interface = new TPZTracerFlow (transport_matid,dimension-1);
//        interface->SetDataTransfer(mSimData);
        mTransportOperator->InsertMaterialObject(interface);
  
    }
    
    mTransportOperator->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(5); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 0;
    active_approx_spaces[3] = 0;
    active_approx_spaces[4] = 1;
//    mTransportOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,transport_meshvec);
    mTransportOperator->BuildMultiphysicsSpace(active_approx_spaces,transport_meshvec);

    
    {
        mTransportOperator->Reference()->ResetReference();
        mTransportOperator->LoadReferences();
        CreateInterfaces(mTransportOperator);
    }
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    mTransportOperator->Print(transport);
#endif
    
}
TPZMultiphysicsCompMesh *TMRSApproxSpaceGenerator::BuildAuxPosProcessCmesh(TPZAlgebraicDataTransfer &Atransfer){
    
    if (!mGeometry) {
        DebugStop();
    }
    
//    TPZManVector<TPZCompMesh *,5> mixed_meshvec = mMixedOperator->MeshVector();
    TPZManVector<TPZCompMesh *,5> transport_meshvec(5);
    
    transport_meshvec[0] = DiscontinuousCmesh(Atransfer);
    transport_meshvec[1] = DiscontinuousCmesh(Atransfer);
    transport_meshvec[2] = DiscontinuousCmesh(Atransfer);
    transport_meshvec[3] = DiscontinuousCmesh(Atransfer);
    transport_meshvec[4] = DiscontinuousCmesh(Atransfer);
    
    
    int dimension = mGeometry->Dimension();
    TPZMultiphysicsCompMesh *auxmesh = new TPZMultiphysicsCompMesh(mGeometry);
    
    //    TMRSMultiphaseFlow<TMRSMemory> * volume = nullptr;
    TPZPostProcessResProp * volume = nullptr;
    auxmesh->SetDefaultOrder(0);
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int material_id = chunk.second;
            //            volume = new TMRSMultiphaseFlow<TMRSMemory>(material_id,d);
            //            volume->SetDataTransfer(mSimData);
            
            volume = new TPZPostProcessResProp(material_id,d);
            //            volume->SetDataTransfer(mSimData);
            
            auxmesh->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]  = get<2>(chunk);
        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        auxmesh->InsertMaterialObject(face);
    }
    
    auxmesh->SetDimModel(dimension);
    TPZManVector<int,5> active_approx_spaces(5); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 0;
    active_approx_spaces[3] = 0;
    active_approx_spaces[4] = 1;
    //    mTransportOperator->BuildMultiphysicsSpaceWithMemory(active_approx_spaces,transport_meshvec);
    auxmesh->BuildMultiphysicsSpace(active_approx_spaces,transport_meshvec);
    return auxmesh;
}

void TMRSApproxSpaceGenerator::SetDataTransfer(TMRSDataTransfer & SimData){
    mSimData = SimData;
    
    if (mSimData.mTGeometry.mGmeshFileName=="") {
        std::string geoname = "PreProcess/meshes/"+mSimData.mSimulationName + "_nLayers_"+ std::to_string(mSimData.mTGeometry.mnLayers)  +"_nRef_"+std::to_string(mSimData.mTGeometry.mnref)+".txt" ;
        mSimData.mTGeometry.mGmeshFileName = geoname;
    }
   
    std::ifstream file(mSimData.mTGeometry.mGmeshFileName);
    
#ifdef USING_BOOST
    boost::posix_time::ptime LoadMesh1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    if(!mGeometry){
        if (file) {
            std::cout<<"The geometric mesh will be loaded from the " + mSimData.mTGeometry.mGmeshFileName + " file."<<std::endl;
            std::string filename = mSimData.mTGeometry.mGmeshFileName;
            TPZPersistenceManager::OpenRead(filename);
            TPZSavable *restore = TPZPersistenceManager::ReadFromFile();
            mGeometry = dynamic_cast<TPZGeoMesh *>(restore);
            PrintGeometry(mSimData.mSimulationName,1,0);
            std::cout<<"The geometric mesh has been loaded successfully."<<std::endl;
#ifdef USING_BOOST
            boost::posix_time::ptime LoadMesh2 = boost::posix_time::microsec_clock::local_time();
            auto delta_LoadMesh = LoadMesh2-LoadMesh1;
            std::cout << "Load GMesh time " << delta_LoadMesh << std::endl;;
#endif
        }
        else{
           std::cout<<" Geometric mesh information has not been entered. Please enter the mesh in a text file or set it in the approximation space object."<<std::endl;
            DebugStop();
        }
        
    }
 
    
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
    TPZMatWithMem<TMRSMemory> * mat_with_memory_d = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material_d);
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
            mem.m_phi = 1.0;
            
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
            mem.m_phi = 1.0;
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

void TMRSApproxSpaceGenerator::FillMaterialMemoryDarcy(int material_id, TPZMultiphysicsCompMesh * MixedOperator, TPZAlgebraicTransport *algebraicTranspor){
    
    if (!MixedOperator || !algebraicTranspor) {
        DebugStop();
    }
    
    TPZCompMesh * cmesh = MixedOperator;
    TPZMaterial * material = cmesh->FindMaterial(material_id);
    if (!material) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZDarcyMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TPZDarcyMemory> * >(material);
    if (!mat_with_memory) {
        DebugStop();
    }
    int nels = cmesh->NElements();
    for (int iel =0; iel<=nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if (!cel) {
            continue;
        }
        TPZManVector<int64_t,20> indices;
        cel->GetMemoryIndices(indices);
        
        
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
    
    TPZMatWithMem<TMRSMemory> * mat_with_memory_trans = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
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
//    const int typeFlux = 1, typePressure = 0;
//    int dimension =mGeometry->Dimension();
//    TMRSDarcyFlowWithMem<TMRSMemory> * volume = nullptr;
//    TPZMixedDarcyFlow * volume = nullptr;
//
//    MixedFluxPressureCmesh->SetDefaultOrder(order);
//    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag =   mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
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
//    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
//    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
//    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
//        int bc_id   = get<0>(chunk);
//        int bc_type = get<1>(chunk);
//        val2[0]  = get<2>(chunk);
//        TPZBndCond * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
//        MixedFluxPressureCmesh->InsertMaterialObject(face);
//    }
//    
//}

void TMRSApproxSpaceGenerator::InsertMaterialObjects(TPZMHMixedMeshControl &control)
{
  
    TPZCompMesh &cmesh = control.CMesh();
    TPZGeoMesh &gmesh = control.GMesh();
    TPZCompMesh *MixedFluxPressureCmesh =  &cmesh;
    
    int dim = gmesh.Dimension();
    MixedFluxPressureCmesh->SetDimModel(dim);
    
    
    int dimension = mGeometry->Dimension();
    
//    TPZMixedDarcyWithFourSpaces *volume = nullptr;
    TPZDarcyFlowWithMem *volume = nullptr;
    
    std::vector<std::map<std::string,int>> DomainDimNameAndPhysicalTag = mSimData.mTGeometry.mDomainDimNameAndPhysicalTag;
    for (int d = 0; d <= dimension; d++) {
        for (auto chunk : DomainDimNameAndPhysicalTag[d]) {
            std::string material_name = chunk.first;
            std::cout << "physical name = " << material_name << std::endl;
            int material_id = chunk.second;
//            volume = new TPZMixedDarcyWithFourSpaces(material_id, d);
            volume = new TPZDarcyFlowWithMem(material_id, d);
//            volume->SetPermeability(1.0);
            volume->SetGravity(mSimData.mTNumerics.m_gravity);
            volume->mSimData = mSimData;
            MixedFluxPressureCmesh->InsertMaterialObject(volume);
        }
    }
    
    if (!volume) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> val1(1,1,0.0); TPZVec<STATE> val2(1,0.0);
    TPZManVector<std::tuple<int, int, REAL>> BCPhysicalTagTypeValue =  mSimData.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue;
    for (std::tuple<int, int, REAL> chunk : BCPhysicalTagTypeValue) {
        int bc_id   = get<0>(chunk);
        int bc_type = get<1>(chunk);
        val2[0]  = get<2>(chunk);
//        std::cout<<"val: "<<val2(0,0)<<std::endl;
        TPZBndCondT<STATE> * face = volume->CreateBC(volume,bc_id,bc_type,val1,val2);
        MixedFluxPressureCmesh->InsertMaterialObject(face);
    }
    
}
void TMRSApproxSpaceGenerator::InitializeFracProperties(TPZMultiphysicsCompMesh * MixedOperator)
{
    //
    TPZCompMesh * cmesh = MixedOperator;
    if (!MixedOperator) {
        DebugStop();
    }
    //
    int n_mats = mSimData.mTReservoirProperties.m_permeabilitiesbyId.size();
    for(int imat=0; imat<n_mats; imat++){
        int matId = mSimData.mTReservoirProperties.m_permeabilitiesbyId[imat].first;
        REAL kappa = mSimData.mTReservoirProperties.m_permeabilitiesbyId[imat].second;
        
        
        TPZMaterial * material1 = cmesh->FindMaterial(matId); ;
        TPZMatWithMem<TMRSMemory> * mat_with_memory1 = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material1);

        
        // Set initial porosity, permeability, saturations, etc ...
        {
            std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector1 = mat_with_memory1->GetMemory();
            int ndata1 = memory_vector1->NElements();
            
            for (int i = 0; i < ndata1; i++) {
                TMRSMemory &mem = memory_vector1.get()->operator [](i);
                mem.m_kappa.Resize(3, 3);
                mem.m_kappa.Zero();
                mem.m_kappa_inv.Resize(3, 3);
                mem.m_kappa_inv.Zero();
                for (int j = 0; j < 3; j++) {
                    mem.m_kappa(j,j) = kappa;
                    mem.m_kappa_inv(j,j) = 1.0/kappa;
                }
            }
        }
    }
    
  
    TPZMaterial * material = cmesh->FindMaterial(10); //matIdFractures;
    if (!material) {
        return;
    }
    TPZMatWithMem<TMRSMemory> * mat_with_memory = dynamic_cast<TPZMatWithMem<TMRSMemory> * >(material);
    if (!mat_with_memory) {
        DebugStop();
    }
    
    // Set initial porosity, permeability, saturations, etc ...
    {
        std::shared_ptr<TPZAdmChunkVector<TMRSMemory>> & memory_vector = mat_with_memory->GetMemory();
        int ndata = memory_vector->NElements();
        
        for (int i = 0; i < ndata; i++) {
            TMRSMemory &mem = memory_vector.get()->operator [](i);
            mem.m_sw = 0.0;
            mem.m_phi = 0.1;
            REAL kappa = mSimData.mTFracProperties.m_Permeability;
            mem.m_kappa.Resize(3, 3);
            mem.m_kappa.Zero();
            mem.m_kappa_inv.Resize(3, 3);
            mem.m_kappa_inv.Zero();
            for (int j = 0; j < 3; j++) {
                mem.m_kappa(j,j) = kappa;
                mem.m_kappa_inv(j,j) = 1.0/kappa;
            }
        }
    }
}


void TMRSApproxSpaceGenerator::findNeighElementbyMatId(TPZGeoElSide &gelside, std::vector<TPZGeoElSide> &neihside, std::set<int> VolMatIds){
    
    int verif =0;
    TPZGeoElSide NeihSideAux = gelside.Neighbour();
    while (gelside != NeihSideAux) {
        for (auto Mat_id: VolMatIds) {
            if (NeihSideAux.Element()->MaterialId() == Mat_id) {
                neihside.push_back(NeihSideAux);
                verif=1;
                break;
            }
        }
        NeihSideAux = NeihSideAux.Neighbour();
    }

}
void TMRSApproxSpaceGenerator::CreateElementInterfaces(TPZGeoEl *gel){
    
    int dimension = gel->Dimension();
    int nsides = gel->NSides();
    int ncoord = gel->NCornerNodes();
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    auto frac = mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag[2];
    int fracmat = frac["Fractures"];
    std::set<int> FracMatID;
    FracMatID.insert(fracmat);
    int val =0;
    TPZManVector<int64_t,3> left_mesh_indexes(2,0);
    left_mesh_indexes[0] = 0;
    left_mesh_indexes[1] = 4;
    TPZManVector<int64_t,3> right_mesh_indexes(1,0);
    right_mesh_indexes[0] = 4;
    
    if (dimension==3) {
        val = gel->NSides(1);
    }
    for (int iside = ncoord+val; iside < nsides-1; iside++) {
        TPZGeoElSide gelside(gel,iside);
        std::vector<TPZGeoElSide> gelneigVec;
        //Volumetric Elements
        TPZCompElSide celside_l=gelside.Reference();
        std::set<int> VolMatIds;
        VolMatIds.insert(1);
        VolMatIds.insert(2);
        findNeighElementbyMatId(gelside,gelneigVec,VolMatIds);
        if(gelneigVec.size()==1){
            TPZGeoElSide gelneig =gelneigVec[0];
            std::vector <TPZGeoElSide >gelfracVec;
            findNeighElementbyMatId(gelside, gelfracVec, FracMatID);
            if (gelfracVec.size()==0) {
                if (gel->Id() < gelneig.Element()->Id()) {
                    TPZCompElSide celside_r = gelneig.Reference();
                    TPZGeoElBC gbc(gelside,transport_matid);
                    int64_t index;
                    TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                }
            }
            else{
//                std::cout<<"Not created Vol-Vol (Frac-Vol -> ok)"<<std::endl;
            }
           
        }
        //BoundaryElements
        std::set<int> Boundaries;
        Boundaries.insert(-1);
        Boundaries.insert(-2);
        Boundaries.insert(-4);
        if(gelneigVec.size()==0){
            findNeighElementbyMatId(gelside,gelneigVec,Boundaries);
            TPZGeoElSide gelneig =gelneigVec[0];
                TPZCompElSide celside_r = gelneig.Reference();
            int matId = gelneig.Element()->MaterialId();
                TPZGeoElBC gbc(gelside,matId);
                int64_t index;
//            if (gel->Id() < gelneig.Element()->Id()) {
                TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
//            }
//            else
//            {
//                TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index,celside_r, celside_l);
//                mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
//            }
            
            
            
        }
    }
}

void TMRSApproxSpaceGenerator::CreateInterfaces(TPZMultiphysicsCompMesh *cmesh){
    
    int nels = cmesh->NElements();
    int dim = cmesh->Dimension();
    int fracId = mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"];
    TPZManVector<std::vector<int64_t>,4> cel_indexes(4);
    for (int64_t el = 0; el < nels; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) DebugStop();
        TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!celmp) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        
        int gel_dim = gel->Dimension();
        if (gel_dim == 2) {
            if (gel->MaterialId() != fracId) {
                continue;
            }
        }
        cel_indexes[gel_dim].push_back(el);
    }
    //Creating Interfaces vol-vol and vol-bound without vol-fractures.
    for(auto elindex: cel_indexes[dim]){
        TPZCompEl *cel = cmesh->Element(elindex);
        TPZGeoEl *gel = cel->Reference();
        CreateElementInterfaces(gel);
    }
    //Creating Interfaces frac-vol frac-frac and frac-boundary
    for(auto elindex: cel_indexes[dim-1]){
        TPZCompEl *cel = cmesh->Element(elindex);
        TPZGeoEl *gel = cel->Reference();
        CreateFracInterfaces(gel);
    }
    
//    std::ofstream file("NewInterfaces.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, file);
    
}
void TMRSApproxSpaceGenerator::CreateFracInterfaces(TPZGeoEl *gel){
   
    int dimension = gel->Dimension();
    int nsides = gel->NSides();
    int ncoord = gel->NCornerNodes();
    int transport_matid = mSimData.mTGeometry.mInterface_material_id;
    auto frac = mSimData.mTGeometry.mDomainFracDimNameAndPhysicalTag[2];
    int fracmat = frac["Fractures"];
    std::set<int> FracMatID;
    FracMatID.insert(fracmat);
    int val =0;
    TPZManVector<int64_t,3> left_mesh_indexes(2,0);
    left_mesh_indexes[0] = 0;
    left_mesh_indexes[1] = 4;
    TPZManVector<int64_t,3> right_mesh_indexes(1,0);
    right_mesh_indexes[0] = 4;
    
    if (dimension==3) {
        val = gel->NSides(1);
    }
    //
    for (int iside = ncoord; iside < nsides; iside++) {
        TPZGeoElSide gelside(gel,iside);
        std::vector<TPZGeoElSide> gelneigVec;
        TPZCompElSide celside_l=gelside.Reference();
        //FracVol Elements
        if(iside == nsides-1 ){
            std::set<int> VolMatIds;
            VolMatIds.insert(1);
            VolMatIds.insert(2);
            findNeighElementbyMatId(gelside,gelneigVec,VolMatIds);
            if(gelneigVec.size()==2){
                for (int ivol = 0; ivol<2; ivol++) {
                    TPZGeoElSide gelneig =gelneigVec[ivol];
                    TPZCompElSide celside_r = gelneig.Reference();
                    int matid=mSimData.mTGeometry.mInterface_material_idFracSup;
                    if (ivol==0) {
                        matid=mSimData.mTGeometry.mInterface_material_idFracInf;
                    }
                    TPZGeoElBC gbc(gelside,matid);
                    int64_t index;
                    if(gel->Id() < gelneig.Element()->Id()){
                        TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                        mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                        
                    }
                    else{
                        TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index,celside_r, celside_l);
                        mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    }
                   
                }
            }
        }
        else{
            std::set<int> FracNeihVec;
            FracNeihVec.insert(10);
            findNeighElementbyMatId(gelside,gelneigVec,FracNeihVec);
            int nneih=gelneigVec.size();
            if (nneih>1){
                DebugStop();
            }
            //Frac-Frac Interfaces
            if (nneih==1) {
                TPZGeoElSide gelneig =gelneigVec[0];
                if (gel->Id() < gelneig.Element()->Id()) {
                TPZCompElSide celside_r = gelneig.Reference();
                int matid=103;
                TPZGeoElBC gbc(gelside,matid);
                int64_t index;
                    if(gel->Id() < gelneig.Element()->Id()){
                        TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                        mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    }
                    else{
                        TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index,celside_r, celside_l);
                        mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);

                    }
            }
            }
            //Frac - Boound Interfaces
            if(nneih==0){
                std::set<int> FracBoundary;
                FracBoundary.insert(-11);
                std::vector<TPZGeoElSide > boundaries;
                findNeighElementbyMatId(gelside,boundaries,FracBoundary);
                int nneihbound=boundaries.size();
                if(nneihbound==1){
                    TPZGeoElSide gelneig =boundaries[0];
                        TPZCompElSide celside_r = gelneig.Reference();
                        int matid=104;
                        TPZGeoElBC gbc(gelside,matid);
                        int64_t index;
                        if(gel->Id() < gelneig.Element()->Id()){
                            TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index, celside_l,celside_r);
                            mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                            
                        }
                        else{
                            TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*mTransportOperator, gbc.CreatedElement(), index,celside_r, celside_l);
                            mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                        
                        }
                }
                else{
                    std::cout<<"Error"<<std::endl;
                }
             
            }
        }
    }
    
}
void TMRSApproxSpaceGenerator::CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC){
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
void TMRSApproxSpaceGenerator::HideTheElements(TPZCompMesh *cmesh){
    int KeepOneLagrangian = 3;
//    if (fHybridize) {
//        KeepOneLagrangian = false;
//    }
//    if (true) {
//        KeepOneLagrangian = false;
//    }
//
    
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
//    TPZGeoMesh *gmesh = fCMesh->Reference();
//    gmesh->ResetReference();
    int64_t nel = mSubdomainIndexGel.size();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel){
            continue;
        }
        TPZCompEl *cel =gel->Reference();
        if(!cel){continue;}
        int indexel = cel->Index();
     
//        TPZCompEl *cel2 = cmesh->ElementVec()[indexel];
//        if(!cel2){
//            DebugStop();
//        }
        int64_t domain = mSubdomainIndexGel[el];
        if (domain == -1) {
            std::cout<<"matId: "<<gel->MaterialId()<<std::endl;
            continue;
        }
        ElementGroups[domain].insert(indexel);
    }
    if (ElementGroups.size() <= 5)
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<int64_t,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << cmesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    
    
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(cmesh, ElementGroups, submeshindices, KeepOneLagrangian);
    
    
    std::cout << "After putting in substructures\n";
//    fMHMtoSubCMesh = submeshindices;
    cmesh->ComputeNodElCon();
    {
        int64_t nc = cmesh->NConnects();
        for (int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = cmesh->ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency())
            {
                c.RemoveDepend();
            }
        }
    }
//    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    
//    GroupandCondenseElements();
    //    GroupandCondenseElementsEigen();
    
    std::cout << "Finished substructuring\n";
}

const bool TMRSApproxSpaceGenerator::isThereFracIntersection() const {
        
    for (auto gel : mGeometry->ElementVec()){
        const int matid = gel->MaterialId();
        if (matid == mMatIDIntersection || matid == mMatIDIntersection+1) { // intersection matid
            return true;
        }
    }
    
    return false;
}

void TMRSApproxSpaceGenerator::HybridizeIntersections(TPZManVector<TPZCompMesh *, 3>& meshvec_Hybrid) {
    
    if (!mHybridizer){
        DebugStop();
    }

    const int matidfrac = 1;
    
    TPZCompMesh* fluxmesh = meshvec_Hybrid[0];
    TPZGeoMesh* gmesh = fluxmesh->Reference();
    fluxmesh->LoadReferences();
    mHybridizer->InsertPeriferalMaterialObjects(meshvec_Hybrid);
        
    int dim = fluxmesh->Dimension();
    for (auto gel : gmesh->ElementVec()) {
        const int gelmatid = gel->MaterialId();
        if (gelmatid != mMatIDIntersection && gelmatid != mMatIDIntersection+1) {
            continue;
        }
        if (gel->Dimension() != 1) {
            DebugStop();
        }
        const bool isIntersectEnd = gelmatid == mMatIDIntersection+1 ? true : false;
        // Search for first neighbor that that is domain
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        while(neigh != gelside){
            TPZGeoEl* gelneigh = neigh.Element();
            int neighmatid = gelneigh->MaterialId();
            int neighdim = gelneigh->Dimension();
            
            if (neighmatid == matidfrac && neighdim == dim) {
                cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " is an intersection element" << endl;
                cout << "===> Trying to split the connects of the flux mesh and create pressure element..." << endl;
                TPZCompEl* celneigh = gelneigh->Reference();
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (celneigh);
                if (!intel)
                    DebugStop();
                
                const int side = neigh.Side();
                TPZCompElSide celsideleft(intel, side);
                bool isNewInterface = mHybridizer->HybridizeInterface(celsideleft,intel,side,meshvec_Hybrid,isIntersectEnd);
                if (isNewInterface) {
                    cout << "=====> Connects splitted succesfuly!" << endl;
                    break;
                }
                else{
                    DebugStop();
                }

            }
            neigh = neigh.Neighbour();
        } // while
        
    }
}

void TMRSApproxSpaceGenerator::CreateIntersectionInterfaceElements(TPZManVector<TPZCompMesh *, 3>& meshvec_Hybrid) {
    TPZCompMesh* cmeshpressure = mMixedOperator->MeshVector()[1];
    mMixedOperator->Reference()->ResetReference();
    cmeshpressure->LoadReferences();
    const int lagrangematid = mHybridizer->lagrangeInterfaceMatId();
    const int lagrangematidend = mHybridizer->lagrangeInterfaceEndMatId();
    for (auto cel : cmeshpressure->ElementVec()) {
        const int celmatid = cel->Material()->Id();
        if (!cel && celmatid != lagrangematid && celmatid != lagrangematidend) {
            continue;
        }
        TPZGeoEl* gel = cel->Reference();
        
        mHybridizer->CreateInterfaceElementsForGeoEl(mMixedOperator, meshvec_Hybrid, gel);
        
    }
}
