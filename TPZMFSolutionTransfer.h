//
//  TPZMFSolutionTransfer.hpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#ifndef TPZMFSolutionTransfer_hpp
#define TPZMFSolutionTransfer_hpp

#include <stdio.h>
#include "TPZMultiphysicsCompMesh.h"

class TPZMFSolutionTransfer
{
    struct Match{
        
        int64_t fconnectindex; //index del connect de la malla multifisica
        std::pair<TPZCompMesh*, int64_t> fconnectTarget; //index del connect de la atomica
        
        Match(){
            
        }
        
        Match(const Match & copy){
            fconnectindex=copy.fconnectindex;
            
        }
        Match &operator=(const Match &other){
            fconnectindex=other.fconnectindex;
            return *this;
        }
        ~Match(){
            
        }
        
        void TransferFromMultiphysics(TPZCompMesh * mfmesh);
        void TransferToMultiphysics(TPZCompMesh * cmesh);
        
    };
    struct MeshTransferData{
        TPZCompMesh * fmultiphysicsmesh;
        TPZStack<Match> fconnecttransfer;
        
        
        MeshTransferData(){
            
        }
        
        MeshTransferData(const MeshTransferData & copy){
            fmultiphysicsmesh=copy.fmultiphysicsmesh;
            fconnecttransfer=copy.fconnecttransfer;
           
        }
        MeshTransferData &operator=(const MeshTransferData &other){
            fmultiphysicsmesh=other.fmultiphysicsmesh;
            fconnecttransfer=other.fconnecttransfer;
            return *this;
        }
        ~MeshTransferData(){
            
        }
        
        void TransferFromMultiphysics();
        void TransferToMultiphysics();
        static void BuildTransferData(TPZCompMesh*, TPZStack<MeshTransferData>&){
            
        }
        void BuildConnectTransfers();
        
    };
    public:
    
        TPZStack<MeshTransferData> fmeshTransfers;
        void TransferFromMultiphysics();
        void TransferToMultiphysics();
        void BuildTransferData(TPZCompMesh*);
        
};

#endif /* TPZMFSolutionTransfer_hpp */
