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
    TPZVec<int> activspaces = multcmesh->GetActiveApproximationSpaces();
    TPZVec<TPZCompMesh*> cmesh_vec = multcmesh->MeshVector();
    for (int ispace =0; ispace< activspaces.size(); ispace++) {
        if (activspaces[ispace]==0) {
            continue;
        }
        int nel = cmesh_vec[ispace]->NElements();
        for (int iel =0; iel<nel; iel++) {
            TPZCompEl *cel = cmesh_vec[ispace]->Element(iel);
            
            
        }
        
    }
    
    
    int nels = cmesh->NElements();
    for (int iel = 0 ; iel<nels; iel++) {
        TPZCompEl* cel = cmesh->Element(iel);
        TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(cel);
        mult->GetActiveApproxSpaces();
        MeshTransferData meshdata;
        //
        //  calcular meshdata//
        //
        fmeshTransfers.push_back(meshdata);

    }
}
