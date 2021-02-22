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
        int nc = cel->NConnects();
        bool found = false;
        if(LagrangeLevelNotCondensed >=0)
        {
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                //                std::cout << "ic ";
                //                c.Print(*cmesh,std::cout);
                if(c.LagrangeMultiplier() >= LagrangeLevelNotCondensed && c.NElConnected() == 1)
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


