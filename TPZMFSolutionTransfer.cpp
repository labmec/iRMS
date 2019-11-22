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
        int aka=0;
    }
}
