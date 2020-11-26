#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "SBFemHdiv.h"
#include <stack> 

#ifdef _AUTODIFF
TLaplaceExample1 LaplaceExact;
#endif

bool gHybrid, gSbfemhdiv, gHdivcollapsed;

TPZGeoMesh * SetupSquareMesh(int nelx)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh, Emat1);

    // Setting boundary conditions
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);

    gmesh->BuildConnectivity();

    if (gHdivcollapsed)
    {
        TPZGeoMesh *gmeshcollapsed = new TPZGeoMesh();
        gmeshcollapsed->SetDimension(1);
        auto nnodes = gmesh->NNodes();
        for (auto node : gmesh->NodeVec())
        {
            auto newnode0 = gmeshcollapsed->NodeVec().AllocateNewElement();
            TPZManVector<REAL,3> xco(3,0.);
            node.GetCoordinates(xco);
            gmeshcollapsed->NodeVec()[newnode0].Initialize(xco, *gmeshcollapsed);
        }
        
        for (auto gel : gmesh->ElementVec())
        {
            if (!gel) continue;
            if (gel->MaterialId() != Emat1) continue;

            for (auto isides = 0; isides < gel->NSides(); isides++)
            {
                TPZGeoElSide gelside = gel->Neighbour(isides);
                if (!gelside || gelside.Dimension() != 1) continue;
                if (gelside.Element()->Index() < gel->Index()) continue; // if the neighbour is a previous element, continue (the interface element was already created)
                
                TPZManVector<int64_t> nodes(2,0);
                nodes[0] = gelside.SideNodeIndex(0);
                nodes[1] = gelside.SideNodeIndex(1);
                int64_t index;
                gmeshcollapsed->CreateGeoElement(EOned, nodes, ESkeleton, index);
            }
        }
        
        gmeshcollapsed->BuildConnectivity();
        
        for (auto gel : gmeshcollapsed->ElementVec())
        {
            if(!gel) continue;
            if (gel->MaterialId() != ESkeleton) continue;
            
            TPZManVector<REAL,3> xco1(3), xco2(3);
            gel->Node(0).GetCoordinates(xco1);
            gel->Node(1).GetCoordinates(xco2);
            if (IsZero(xco1[1] + 1) && IsZero(xco2[1] + 1))
            {
                TPZGeoElBC(gel,2,Ebc1);
            }
            else if(IsZero(xco1[0] - 1) && IsZero(xco2[0] - 1))
            {
                TPZGeoElBC(gel,2,Ebc2);
            }
            else if(IsZero(xco1[1] - 1) && IsZero(xco2[1] - 1))
            {
                TPZGeoElBC(gel,2,Ebc3);
            }
            else if(IsZero(xco1[0] + 1) && IsZero(xco2[0] + 1))
            {
                TPZGeoElBC(gel,2,Ebc4);
            }
        }
        std::ofstream out("gmesh.txt");
        gmeshcollapsed->Print(out);
        return gmeshcollapsed;
    }
    
    return gmesh;
}

TPZGeoMesh * SetupCollapsedMesh(int nelx)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh, Emat1);

    // Setting boundary conditions
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);

    gmesh->BuildConnectivity();

    TPZGeoMesh *gmeshcollapsed = new TPZGeoMesh();
    gmeshcollapsed->SetDimension(2);
    auto nnodes = gmesh->NNodes();
    for (auto node : gmesh->NodeVec())
    {
        auto newnode0 = gmeshcollapsed->NodeVec().AllocateNewElement();
        TPZManVector<REAL,3> xco(3,0.);
        node.GetCoordinates(xco);
        gmeshcollapsed->NodeVec()[newnode0].Initialize(xco, *gmeshcollapsed);
    }
    
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != Emat1) continue;

        // creating node in the element's center
        TPZManVector<REAL,3> xco(3,0.);
        auto nnodes = gel->NNodes();
        TPZManVector<int64_t,4> nodeindices(nnodes,0);
        gel->GetNodeIndices(nodeindices);

        for (auto inode : nodeindices)
        {
            auto node = gel->Node(inode);
            TPZManVector<REAL,3> xnode(3,0.);
            node.GetCoordinates(xnode);
            xco[0] += xnode[0]; xco[1] += xnode[1]; xco[2] += xnode[2];
        }
        xco[0] = xco[0]/nnodes; xco[1] = xco[1]/nnodes; xco[2] = xco[2]/nnodes;

        auto newnode0 = gmeshcollapsed->NodeVec().AllocateNewElement();
        gmeshcollapsed->NodeVec()[newnode0].Initialize(xco, *gmeshcollapsed);

        for (auto isides = 0; isides < gel->NSides(); isides++)
        {
            TPZGeoElSide gelside = gel->Neighbour(isides);
            if (!gelside || gelside.Dimension() != 1) continue;
            if (gelside.Element()->Index() < gel->Index()) continue; // if the neighbour is a previous element, continue (the interface element was already created)
            
            TPZManVector<int64_t> nodes(2,0);
            nodes[0] = gelside.SideNodeIndex(0);
            nodes[1] = gelside.SideNodeIndex(1);
            int64_t index;
            gmeshcollapsed->CreateGeoElement(EOned, nodes, ESkeleton, index);
            nodes.Resize(4);
            nodes[2] = newnode0; nodes[3] = newnode0;
            gmeshcollapsed->CreateGeoElement(EQuadrilateral, nodes, Emat0, index);
        }
    }
    
    gmeshcollapsed->BuildConnectivity();
    
    // Boundary conditions
    for (auto gel : gmeshcollapsed->ElementVec())
    {
        if(!gel) continue;
        if (gel->MaterialId() != ESkeleton) continue;
        TPZGeoElBC(gel,2,Ebc1);
        
        TPZManVector<REAL,3> xco1(3), xco2(3);
        gel->Node(0).GetCoordinates(xco1);
        gel->Node(1).GetCoordinates(xco2);
        if (IsZero(xco1[1] + 1) && IsZero(xco2[1] + 1))
        {
            TPZGeoElBC(gel,2,Ebc1);
        }
        else if(IsZero(xco1[0] - 1) && IsZero(xco2[0] - 1))
        {
            TPZGeoElBC(gel,2,Ebc2);
        }
        else if(IsZero(xco1[1] - 1) && IsZero(xco2[1] - 1))
        {
            TPZGeoElBC(gel,2,Ebc3);
        }
        else if(IsZero(xco1[0] + 1) && IsZero(xco2[0] + 1))
        {
            TPZGeoElBC(gel,2,Ebc4);
        }
    }

    gmeshcollapsed->BuildConnectivity();

    if (0)
    {
        for (auto gel : gmeshcollapsed->ElementVec())
        {
            if (!gel || gel->MaterialId() != ESkeleton)
            {
                continue;
            }
            TPZGeoElSide gelside(gel, 2);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (gelside != neighbour)
            {
                switch (neighbour.Element()->MaterialId())
                {
                case Emat0:
                    TPZGeoElBC(gelside, Eintleft);
                    break;
                case Ebc1:
                    TPZGeoElBC(gelside, Eintright);
                    break;
                case Ebc2:
                    TPZGeoElBC(gelside, Eintright);
                    break;
                case Ebc3:
                    TPZGeoElBC(gelside, Eintright);
                    break;
                case Ebc4:
                    TPZGeoElBC(gelside, Eintright);
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, Eint);
            }
        }
    }
    
    std::ofstream out("gmesh.txt");
    gmeshcollapsed->Print(out);
    return gmeshcollapsed;
}

TPZCompMesh * flux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    
    TPZVecL2 *mat = new TPZVecL2(Emat1);
    mat->SetDimension(dim);
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    
    cmesh->InsertMaterialObject(mat); //Insere material na malha
    cmesh->SetDimModel(dim);

    int nstate = 1;
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

TPZCompMesh * fluxhdivcollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = 2;
    TPZCompMesh * cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->SetDefaultOrder(POrder);
    cmeshcollapsed->CleanUp();
    cmeshcollapsed->SetDimModel(dim);

    auto mat = new TPZNullMaterial(ESkeleton);
    mat->SetNStateVariables(1);
    mat->SetDimension(dim);
    cmeshcollapsed->InsertMaterialObject(mat); //Insere material na malha
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->Dimension() != 1) continue;
        if (gel->MaterialId() != ESkeleton) continue;
        int64_t index;
        auto celhdivc = new TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>(*cmeshcollapsed, gel, index);
        if (1)
        {
            celhdivc->Print(std::cout);
        }
    }

    cmeshcollapsed->SetDimModel(2);
    cmeshcollapsed->SetAllCreateFunctionsHDiv();
    
    int nstate = 1;
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    std::set<int> matids = {Ebc1, Ebc2, Ebc3, Ebc4};
    cmeshcollapsed->Reference()->ResetReference();
    cmeshcollapsed->AutoBuild(matids);
    cmeshcollapsed->AdjustBoundaryElements();
    cmeshcollapsed->CleanUpUnconnectedNodes();
    
    cmeshcollapsed->ExpandSolution();
    return cmeshcollapsed;
}

using namespace std;
TPZCompMesh * fluxhdivsbfem(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = 1;
    auto cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->SetDefaultOrder(POrder);
    cmeshcollapsed->SetDimModel(dim);
    cmeshcollapsed->CleanUp();

    auto mat = new TPZNullMaterial(ESkeleton);
    mat->SetNStateVariables(1);
    mat->SetDimension(dim);
    // mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    // mat->SetForcingFunctionExact(LaplaceExact.Exact());
    
    cmeshcollapsed->InsertMaterialObject(mat); //Insere material na malha
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != Emat0) continue;
        int64_t index;
        TPZGeoElSide gelside(gel,4);
        TPZGeoEl * gel1d = gelside.Neighbour().Element();
        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(*cmeshcollapsed, gel1d, gelside, index);
    }
    cmeshcollapsed->SetDimModel(dim);
    cmeshcollapsed->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshcollapsed->AutoBuild();

    // Creating HdivBound elements for the external fluxes
    // TPZManVector<int64_t,4> indexleft(cmeshcollapsed->NElements(),-1);
    // TPZManVector<int64_t,4> indexright(cmeshcollapsed->NElements(),-1);
    // int count = 0;
    for (auto cel : cmeshcollapsed->ElementVec())
    {
        if (!cel) continue;

        TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdiv = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
        if (!celhdiv) continue;

        auto gel = cel->Reference();
        if (!gel)
        {
            DebugStop();
        }
        
        if (1)
        {
            int64_t index;
            auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(*cmeshcollapsed,gel,index);
            hdivboundleft->SetConnectIndex(0,celhdiv->ConnectIndex(3));
            // indexleft[count] = index;

            auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(*cmeshcollapsed,gel,index);
            hdivboundright->SetConnectIndex(0,celhdiv->ConnectIndex(4));
            // indexright[count] = index;
        }
        // count++;
    }
    cmeshcollapsed->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    // cmeshcollapsed->AutoBuild();

    // for (auto cel : cmeshcollapsed->ElementVec())
    // {
    //     if (!cel) continue;

    //     TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdiv = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
    //     if (!celhdiv) continue;

    //     auto gel = celhdiv->Reference();
    //     if (!gel)
    //     {
    //         DebugStop();
    //     }
    //     {
    //         auto index = indexleft[cont];
    //         auto celbound = cmeshcollapsed->ElementVec()[index];
    //         auto celboundhdiv = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeLinear> * >(celbound);
    //         celhdiv->SetConnectIndex(3,celboundhdiv->ConnectIndex(0));
    //     }
    //     {
    //         auto index = indexright[cont];
    //         auto celbound = cmeshcollapsed->ElementVec()[index];
    //         auto celboundhdiv = dynamic_cast<TPZCompElHDivBound2<pzshape::TPZShapeLinear> * >(celbound);
    //         celhdiv->SetConnectIndex(4,celboundhdiv->ConnectIndex(0));
    //     }
    // }
    cmeshcollapsed->LoadReferences();
    // cmeshcollapsed->AdjustBoundaryElements();
    cmeshcollapsed->CleanUpUnconnectedNodes();
    // cmeshcollapsed->ExpandSolution();

    int nstate = 1;
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }

    std::set<int> matids = {Ebc1, Ebc2, Ebc3, Ebc4};
    cmeshcollapsed->Reference()->ResetReference();
    cmeshcollapsed->AutoBuild(matids);
    cmeshcollapsed->AdjustBoundaryElements();
    cmeshcollapsed->CleanUpUnconnectedNodes();
    
    cmeshcollapsed->ExpandSolution();
    
    return cmeshcollapsed;
}


TPZCompMesh * pressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{        
    auto fmat = Emat1;
    auto dim = 2;
    if (gHdivcollapsed || gSbfemhdiv)
    {
        fmat = ESkeleton;
        dim = 1;
    }

    auto cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    TPZMaterial *mat;
    
    if (gHdivcollapsed || gSbfemhdiv)
    {
        mat = new TPZHybridPoissonCollapsed(ESkeleton,dim);
    }
    else
    {
        mat = new TPZMixedPoisson(fmat, dim);
    }
    cmesh->InsertMaterialObject(mat); //Insere material na malha

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();

    for(auto newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZMultiphysicsCompMesh *  multiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder)
{
    gmesh->ResetReference();
    auto dim = gmesh->Dimension();
    auto fmat = Emat1;

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    auto mat = new TPZMixedPoisson(fmat, dim);
    
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    std::cout << "Creating multiphysics mesh\n";
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshf;
    meshvector[1] = cmeshp;

    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    
    return cmesh;
}

TPZMultiphysicsCompMesh *  multiphysicscollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder)
{
    // gmesh->ResetReference();
    auto fmat = ESkeleton;
    int dim = 2;

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    auto mat = new TPZHybridPoissonCollapsed(ESkeleton,dim);
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }

    std::cout << "Creating multiphysics mesh\n";
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshf;
    meshvector[1] = cmeshp;

    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    
    return cmesh;
}