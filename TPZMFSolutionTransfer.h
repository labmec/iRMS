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
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


class TPZMFSolutionTransfer
{
    struct Match{
        
        int64_t fblocknumber; // sequence number
        std::pair<TPZCompMesh*, int64_t> fblockTarget;// block target
        
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
        TPZMultiphysicsCompMesh * fmultiphysicsmesh;
        
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
        
//        void BuildTransferData(TPZCompMesh* cmesh, TPZStack<MeshTransferData>& ftansdata){
        void BuildTransferData(TPZCompMesh* cmesh){
            
            TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
            
//            if(multcmesh){
//                fmultiphysicsmesh = cmesh;
//            }
            TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh);
            
            TPZCompMesh * targetmesh = 0;
            
            if(subcmesh){
                targetmesh =subcmesh;
            }
            if(multcmesh){
                targetmesh =multcmesh;
            }
            TPZVec<TPZCompMesh*> mesh_vec = fmultiphysicsmesh->MeshVector();
            TPZVec<int> actspaces= fmultiphysicsmesh->GetActiveApproximationSpaces();
            
            int nspaces =actspaces.size();
                for (int ispace =0; ispace <nspaces; ispace++){
                if(actspaces[ispace]==0){       //Computes only in actives spaces
                    continue;
                }
                
                    
                int nels = targetmesh->NElements();
                    
                fmultiphysicsmesh = multcmesh;
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
                    if(condition1||condition2)  {
                        continue;
                    }
                    int nsubel =0;
                    TPZVec<TPZCompEl*> subels;
                    if(elgroup){
                        
                         subels = elgroup->GetElGroup();
                         nsubel = subels.size();
                         for(int isub =0; isub<nsubel; isub++){
                            TPZCompEl *subcel = subels[isub];
                            mul =dynamic_cast<TPZMultiphysicsElement *>(subcel);
                            
                             //
                             Vector<int> initials(actspaces.size()+1,0);
                             for (int ini=0; ini<ispace; ini++) {
                                 if(mul->Element(ini)){
                                     initials[ini+1] = initials[ini]+mul->Element(ini)->NConnects();
                                 }
                             }
                             
                             int initial = initials[ispace];
                             TPZCompEl *celfrom = mul->Element(ispace);
                             if (!celfrom) {
                                 continue;
                             }
                             int nconnects = celfrom->NConnects();
                             for (int icon = 0 ; icon < nconnects; icon++){
                                 TPZConnect &conectFrom = celfrom->Connect(icon);
                                 int seqnumberAto = conectFrom.SequenceNumber();
                                 int nconectsmult = mul->NConnects();
                                 TPZConnect &conectTarget = mul->Connect(initial + icon);
                                 int seqnumberMF = conectTarget.SequenceNumber();
                                 Match currentmatch;
                                 currentmatch.fblocknumber = seqnumberAto;
                                 std::pair<TPZCompMesh*, int> target = std::make_pair(mesh_vec[ispace],seqnumberMF);
                                 currentmatch.fblockTarget =target;
                                 fconnecttransfer.push_back(currentmatch);
                             }
                             //
                         }
                    }
    
                    
                    Vector<int> initials(actspaces.size()+1,0);
                    for (int ini=0; ini<ispace; ini++) {
                        if(mul->Element(ini)){
                            initials[ini+1] = initials[ini]+mul->Element(ini)->NConnects();
                        }
                    }

                    int initial = initials[ispace];
                    TPZCompEl *celfrom = mul->Element(ispace);
                    if (!celfrom) {
                        continue;
                    }
                    int nconnects = celfrom->NConnects();
                    for (int icon = 0 ; icon < nconnects; icon++){
                        TPZConnect &conectFrom = celfrom->Connect(icon);
                        int seqnumberAto = conectFrom.SequenceNumber();
                        int nconectsmult = mul->NConnects();
                        TPZConnect &conectTarget = mul->Connect(initial + icon);
                        int seqnumberMF = conectTarget.SequenceNumber();
                        Match currentmatch;
                        currentmatch.fblocknumber = seqnumberAto;
                        std::pair<TPZCompMesh*, int> target = std::make_pair(mesh_vec[ispace],seqnumberMF);
                        currentmatch.fblockTarget =target;
                        fconnecttransfer.push_back(currentmatch);
                    }
                
                }
            }
            
            
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
