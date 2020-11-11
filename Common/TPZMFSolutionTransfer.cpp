//
//  TPZMFSolutionTransfer.cpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#include "TPZMFSolutionTransfer.h"


void TPZMFSolutionTransfer::Match::TransferFromMultiphysics(TPZCompMesh * cmesh){
    
    TPZBlock<STATE> &blockMF = cmesh->Block();
    TPZBlock<STATE>* blockAto = fblockTarget.first;
    int seqMF = fblockTarget.second;
    int seqAto = fblocknumber;
    int blocksizeAto = blockAto->Size(seqAto);
    int blocksizeTarg = blockMF.Size(seqMF);
    if(blocksizeAto!=blocksizeTarg){
        DebugStop();
    }
    for (int idf=0; idf<blocksizeAto; idf++) {
        blockAto->Put(seqAto, idf, 0, blockMF.Get(seqMF, idf, 0));
    }
}
void TPZMFSolutionTransfer::Match::TransferToMultiphysics(TPZCompMesh * cmesh){
//    cmesh->InitializeBlock();
    TPZBlock<STATE> &blockToTransfer = cmesh->Block();
    TPZBlock<STATE>* blockAto = fblockTarget.first;
    int seqtoTrans = fblockTarget.second;
    int seqAto = fblocknumber;
    
    int blocksizeAto = blockAto->Size(seqAto);
    int blocksizeTarg = blockToTransfer.Size(seqtoTrans);
    if(blocksizeAto!=blocksizeTarg){
        DebugStop();
    }
    for (int idf=0; idf<blocksizeAto; idf++) {
        blockToTransfer.Put(seqtoTrans, idf, 0, blockAto->Get(seqAto, idf, 0));
    }
}

void TPZMFSolutionTransfer::MeshTransferData::TransferToMultiphysics(){
   
    int nmatch = fconnecttransfer.size();
    fcmesh_ori->InitializeBlock();
    for (int imatch=0; imatch<nmatch; imatch++) {
        fconnecttransfer[imatch].TransferToMultiphysics(fcmesh_ori);
    }
}
void TPZMFSolutionTransfer::MeshTransferData::TransferFromMultiphysics(){
    int nmatch = fconnecttransfer.size();
    for (int imatch=0; imatch<nmatch; imatch++) {
        fconnecttransfer[imatch].TransferFromMultiphysics(fcmesh_ori);
    }
}

void TPZMFSolutionTransfer::MeshTransferData::BuildTransferData(TPZCompMesh* cmesh){
    TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    
    TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh);
    TPZCompMesh * targetmesh = 0;
    
    if(subcmesh){
        targetmesh =subcmesh;
    }
    
    if(multcmesh){
        targetmesh =multcmesh;
    }

    int nels = targetmesh->NElements();
    for (int iel =0; iel <nels; iel++){
        TPZCompEl *celtarget = cmesh->Element(iel);
        TPZMultiphysicsElement *mul = dynamic_cast<TPZMultiphysicsElement *>(celtarget);
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(celtarget);
        TPZElementGroup *elgroup;
        
        if(cond){
            mul =dynamic_cast<TPZMultiphysicsElement *>(cond->ReferenceCompEl());
            
            elgroup = dynamic_cast<TPZElementGroup *>(cond->ReferenceCompEl());
            
        }
        
        bool condition1 = !mul && !cond;
        bool condition2 = cond && (!mul && !elgroup);
        
        if((condition1||condition2) )  {
            continue;
        }
        int nsubel =0;
        TPZVec<TPZCompEl*> subels;
        
        if(((elgroup) && (!mul)) || mul){
            
            if((elgroup) && (!mul)){
                subels = elgroup->GetElGroup();
            }
            else{
                subels.resize(1);
                subels[0] = mul;
            }
            nsubel = subels.size();
            for(int isub =0; isub<nsubel; isub++){
                TPZCompEl *subcel = subels[isub];
                mul =dynamic_cast<TPZMultiphysicsElement *>(subcel);
                TPZVec<int> act_spacesEl = mul->GetActiveApproxSpaces();
                int initialtest =0;
                for (int iespacetes=0; iespacetes<act_spacesEl.size(); iespacetes++) {
                    if (act_spacesEl[iespacetes]==0) {
//                        initialtest += mul->Element(iespacetes)->NConnects();
                        continue;
                    }
                    TPZCompEl *celfrom = mul->Element(iespacetes);
                    if (!celfrom) {
                        continue;
                    }
                    
                    int nconnects = celfrom->NConnects();
                    for (int icon = 0 ; icon < nconnects; icon++){
                        TPZConnect &conectFrom = celfrom->Connect(icon);
                        int seqnumberAto = conectFrom.SequenceNumber();
                        TPZConnect &conectTarget = mul->Connect(initialtest + icon);
                        int seqnumberMF = conectTarget.SequenceNumber();
                        Match currentmatch;
                        currentmatch.fblocknumber = seqnumberAto;
                        TPZBlock<STATE> *blockAto = &celfrom->Mesh()->Block();
                        std::pair<TPZBlock<STATE> *, int> target = std::make_pair(blockAto,seqnumberMF);
                        currentmatch.fblockTarget = target;
                        fconnecttransfer.push_back(currentmatch);
                    }
                    initialtest += nconnects;
                }
            }
        }
    }
}


void TPZMFSolutionTransfer::TransferFromMultiphysics(){
    
    int nsolutionstransfers = fmeshTransfers.size();
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        fmeshTransfers[isoltrans].TransferFromMultiphysics();
    }
    
}
void TPZMFSolutionTransfer::TransferToMultiphysics(){

    int nsolutionstransfers = fmeshTransfers.size();
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        fmeshTransfers[isoltrans].TransferToMultiphysics();
    }
}

void TPZMFSolutionTransfer::BuildTransferData(TPZCompMesh* cmesh){
    
    MeshTransferData transdata;
    transdata.fcmesh_ori =cmesh;
    transdata.BuildTransferData(cmesh);
    fmeshTransfers.push_back(transdata);
    int nels = cmesh->NElements();
    for (int iel=0; iel<nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        TPZSubCompMesh * subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            continue;
        }
        MeshTransferData transdatasub;
        transdatasub.fcmesh_ori =subcmesh;
        transdatasub.BuildTransferData(subcmesh);
        fmeshTransfers.push_back(transdatasub);
    }
}
