//
//  TPZMFSolutionTransfer.cpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#include "TPZMFSolutionTransfer.h"

void TPZMFSolutionTransfer::TransferFromMultiphysics(){
    
}
void TPZMFSolutionTransfer::TransferToMultiphysics(){

    int nsolutionstransfers = fmeshTransfers.size();
    
    for (int isoltrans=0; isoltrans<nsolutionstransfers ; isoltrans++) {
        
        TPZCompMesh *cmesh = fmeshTransfers[isoltrans].fmultiphysicsmesh;
        cmesh->InitializeBlock();
        TPZBlock<STATE> &blockMF = cmesh->Block();
        int nblocks = blockMF.NBlocks();
        for (int iblock =0; iblock<nblocks; iblock++) {
            int nsizeblock =blockMF.Size(iblock);
            std::cout<<"nsizeblock: "<<nsizeblock<<std::endl;
        }
        
        int nmatch = fmeshTransfers[isoltrans].fconnecttransfer.size();
        for (int imatch=0; imatch<nmatch; imatch++) {
            
            std::pair<TPZCompMesh*, int64_t> match =fmeshTransfers[isoltrans].fconnecttransfer[imatch].fblockTarget;
            
            int seqMF = std::get<1>(match);
            int seqAto =fmeshTransfers[isoltrans].fconnecttransfer[imatch].fblocknumber;
            TPZCompMesh *cmeshAto = std::get<0>(match);
            TPZBlock<STATE> &blockAto = cmeshAto->Block();
            int blocksize = blockAto.Size(seqAto);
            
            int64_t posMF = cmesh->Block().Position(seqMF);
            int64_t posAto = cmeshAto->Block().Position(seqAto);
            for (int idf=0; idf<blocksize; idf++) {
//                blockMF.Put(seqMF, idf, 0, blockAto.Get(seqAto, idf, 0));
                

                cmesh->Solution()(posMF+idf,0) = cmeshAto->Solution()(posAto+idf,0);
            }
        }
    }
}
void TPZMFSolutionTransfer::BuildTransferData(TPZCompMesh* cmesh){
    
    TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    MeshTransferData transdata;
    transdata.fmultiphysicsmesh =multcmesh;
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
        transdatasub.fmultiphysicsmesh = multcmesh;
        transdatasub.BuildTransferData(subcmesh);
        fmeshTransfers.push_back(transdatasub);
        
    }
}

