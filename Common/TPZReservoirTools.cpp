#include "TPZReservoirTools.h"


void TPZReservoirTools::CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix)
{
    //    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        if (KeepOneLagrangian) {
            int count = 0;
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.LagrangeMultiplier() > 0) {
                    c.IncrementElConnected();
                    count++;
                    if(count == 1 && c.NState() == 1)
                    {
                        break;
                    } else if(count == 2 && c.NState() == 2)
                    {
                        break;
                    } else if(count == 3 && c.NState() == 3)
                    {
                        break;
                    }
                    
                }
            }
        }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
        {
            TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
            cond->SetLambda(1.0);
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}
/// create a condensed element and do not condense the connect with a given lagrange level
// the method does the same procedure as CreatedCondensedElements, but has different policy for
// keeping a connect out the condensation loop
void TPZReservoirTools::CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix)
{
    //    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        //        std::cout << "Element " << el << std::endl;
        TPZCompEl *cel = cmesh->Element(el);
     
        
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int nc = cel->NConnects();
        bool found = false;
        if(LagrangeLevelNotCondensed >=0)
        {
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                //                std::cout << "ic ";
                //                c.Print(*cmesh,std::cout);
        
                if((c.LagrangeMultiplier() >= LagrangeLevelNotCondensed && c.NElConnected() == 1) )
                {
                    c.IncrementElConnected();
                    found = true;
                }
            }
        }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
        {
            if(LagrangeLevelNotCondensed >= 0 && !found) DebugStop();
            TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
            cond->SetLambda(1.0);
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}

//Condense by MATID

void TPZReservoirTools::CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix, std::set<int> matids)
{
    
    //    cmesh->ComputeNodElCon();
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        //        std::cout << "Element " << el << std::endl;
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
     
        int nc = cel->NConnects();
        bool found = false;
        if(LagrangeLevelNotCondensed >=0)
        {
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                char lag =c.LagrangeMultiplier();
                int ncon =c.NElConnected();
                if((c.LagrangeMultiplier() >= LagrangeLevelNotCondensed && c.NElConnected() == 1))
                {
                    c.IncrementElConnected();
                    found = true;
                }
            }
        }
        int ic;
        for (ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() || c.NElConnected() > 1) {
                continue;
            }
            break;
        }
        bool cancondense = (ic != nc);
        if(cancondense)
        {
            int gel_matId = 0;
            if (gel) {
                    gel_matId = gel->MaterialId();
            }
            
            if(LagrangeLevelNotCondensed >= 0 && !found) DebugStop();
            int verif = 0;
            
            for (auto matId:matids) {
                if (gel_matId==matId) {
                    verif=1;
                    break;
                }
            }
            if (verif==1) {
//                TPZFastCondensedElement *cond = new TPZFastCondensedElement(cel, keepmatrix);
//                cond->SetLambda(1.0);
            }
            else{
//                TPZCondensedCompEl *cond = new TPZCondensedCompEl(cel, keepmatrix);
            }
        }
        
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}

void TPZReservoirTools::FindCondensed(TPZCompEl *cel, TPZStack<TPZFastCondensedElement *> &condensedelements)
{
    TPZFastCondensedElement *cond = dynamic_cast<TPZFastCondensedElement *>(cel);
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
    if(cond)
    {
        condensedelements.Push(cond);
        FindCondensed(cond->ReferenceCompEl(), condensedelements);
        return;
    }
    if(elgr)
    {
        const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
        int nel = elvec.size();
        for (int el = 0; el<nel; el++) {
            TPZCompEl *loccel = elvec[el];
            FindCondensed(loccel, condensedelements);
        }
    }
}

void TPZReservoirTools::AddDependency( std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons){
    
    for(auto pairfatherson:fatherAndSons){
        TPZCompEl *celfat = pairfatherson.first->Reference();
        TPZInterpolatedElement *intelFather = dynamic_cast<TPZInterpolatedElement*>(celfat);
        if(celfat->NConnects() !=1){
            DebugStop();
        }
        for(auto son:pairfatherson.second){
            TPZCompEl *celSon= son->Reference();
            if(celSon->NConnects() !=1){
                DebugStop();
            }
            TPZInterpolatedElement *intelSon = dynamic_cast<TPZInterpolatedElement*>(celSon);
            intelSon->RestrainSide(son->NSides()-1, intelFather, pairfatherson.first->NSides()-1);
        }
    }
}
void TPZReservoirTools::TakeFatherSonsCorrespondence(TPZCompMesh *fluxCmesh,  std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons){
    
    int nels = fluxCmesh->NElements();
    TPZGeoMesh *gmesh = fluxCmesh->Reference();
//    gmesh->ResetReference();
    fluxCmesh->LoadReferences();
    std::map<int, std::vector<TPZGeoEl* >> interfaces;
    std::vector<int> matIds;
    int skeletonCoarseId = 19;
    matIds.push_back(skeletonCoarseId);
    TakeElementsbyID(gmesh, interfaces, matIds);
    int count=-1;
    for(auto gel:interfaces[skeletonCoarseId] ){
        TPZStack<TPZGeoEl *> Sons;
        gel->YoungestChildren(Sons);
        TPZCompEl *celFat = gel->Reference();
        std::vector<TPZGeoEl*> sons;
        TPZGeoElSide fatside(gel, gel->NSides()-1);
        TPZGeoElSide neigfatside = fatside.Neighbour();
        while(fatside!=neigfatside){
            TPZCompEl *cel = neigfatside.Element()->Reference();
            int matId =neigfatside.Element()->MaterialId();
            int nconnects =-1;
            if(cel){
                nconnects=cel->NConnects();
            }
            if(matId==19 && nconnects>0){
                std::pair fatson = std::make_pair(neigfatside.Element(),sons);
                fatherAndSons.push_back(fatson);
                count++;
            }
            neigfatside = neigfatside.Neighbour();
        }
        for(auto songel: Sons){
            TPZGeoElSide gelside(songel, songel->NSides()-1);
            TPZGeoElSide neig = gelside.Neighbour();
            std::vector<TPZCompEl *> fluxEls;
            while(neig!=gelside){
                TPZGeoEl *neigGel = neig.Element();
                int matId = neigGel->MaterialId();
                TPZCompEl *celneih = neigGel->Reference();
                int nconnects = -1;
                if(celneih){
                    nconnects=celneih->NConnects();
                }
                
//                TPZCompEl *neigCel = neigGel->Reference();
                if( matId==40){
                    fatherAndSons[count].second.push_back(neigGel);
//                    TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *>(neigCel);
                }
                neig=neig.Neighbour();
            }
        }
    }
}
void TPZReservoirTools::TakeElementsbyID(TPZGeoMesh *mGeometry, std::map<int, std::vector<TPZGeoEl* >> &interfaces, std::vector<int> &matIds){
    int nels = mGeometry->NElements();
    for(int iel =0; iel <nels; iel++){
        TPZGeoEl *gel = mGeometry->Element(iel);
        int matId = gel->MaterialId();
        for(auto mat: matIds){
            if(matId==mat){
                interfaces[matId].push_back(gel);
                break;
            }
        }
    }
}
