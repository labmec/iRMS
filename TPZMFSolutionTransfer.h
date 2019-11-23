//  TPZMFSolutionTransfer.hpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#ifndef TPZMFSolutionTransfer_hpp
#define TPZMFSolutionTransfer_hpp

#include <stdio.h>
#include "TPZMultiphysicsCompMesh.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


class TPZMFSolutionTransfer
{
    struct Match{
        
        int64_t fblocknumber; // sequence number
        std::pair<TPZBlock<STATE> *, int64_t> fblockTarget;// block target
        
        Match(){
            
        }
        Match(const Match & copy){
            fblocknumber=copy.fblocknumber;
            fblockTarget = copy.fblockTarget;
        }
        Match &operator=(const Match &other){
            fblocknumber=other.fblocknumber;
            fblockTarget = other.fblockTarget;
            return *this;
        }
        ~Match(){
            
        }
        
        void TransferFromMultiphysics(TPZCompMesh * mfmesh);
        void TransferToMultiphysics(TPZCompMesh * cmesh);
        
    };
    struct MeshTransferData{
        TPZCompMesh * fcmesh_ori;
        TPZStack<Match> fconnecttransfer;
        
        MeshTransferData(){
            
        }
        
        MeshTransferData(const MeshTransferData & copy){
            fcmesh_ori=copy.fcmesh_ori;
            fconnecttransfer=copy.fconnecttransfer;
           
        }
        MeshTransferData &operator=(const MeshTransferData &other){
            fcmesh_ori=other.fcmesh_ori;
            fconnecttransfer=other.fconnecttransfer;
            return *this;
        }
        ~MeshTransferData(){
            
        }
        
        void TransferFromMultiphysics();
        void TransferToMultiphysics();
        void BuildTransferData(TPZCompMesh* cmesh);
        
    };
    public:
    
        TPZStack<MeshTransferData> fmeshTransfers;
    
        void TransferFromMultiphysics();
        void TransferToMultiphysics();
        void BuildTransferData(TPZCompMesh*);
        
};

#endif /* TPZMFSolutionTransfer_hpp */
