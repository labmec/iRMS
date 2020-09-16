//
//  TPZSymetricSpStructMatrixEigenEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#include "TPZSymetricSpStructMatrixEigen.h"


#include "pzstrmatrix.h"

#include "pzgeoelbc.h"
#include "pzgmesh.h"
#include "pzcmesh.h"

#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"

#include "pzdxmesh.h"
#include <fstream>

#include "pzelmat.h"

#include "pzsysmp.h"
#include "pzmetis.h"
#include "pzbndcond.h"
#include "TPZTimer.h"
#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

using namespace std;

TPZStructMatrix * TPZSymetricSpStructMatrixEigen::Clone(){
    return new TPZSymetricSpStructMatrixEigen(*this);
}
TPZMatrix<STATE> * TPZSymetricSpStructMatrixEigen::CreateAssemble(TPZFMatrix<STATE> &rhs,
                                                             TPZAutoPointer<TPZGuiInterface> guiInterface){
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        LOGPZ_DEBUG(logger,"TPZSymetricSpStructMatrixEigen::CreateAssemble starting");
    }
#endif
    
    int64_t neq = fMesh->NEquations();
    
    if(fMesh->FatherMesh()) {
        cout << "TPZSymetricSpStructMatrixEigen should not be called with CreateAssemble for a substructure mesh\n";
        return new TPZSYsmpMatrixEigen<STATE>(0,0);
    }
    //    std::cout << "Creating\n";
    TPZMatrix<STATE> *stiff = Create();//new TPZFYsmpMatrix(neq,neq);
    
    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (stiff);
    rhs.Redim(neq,1);
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
    TPZTimer before("Assembly of a sparse matrix");
    //    std::cout << "Assembling\n";
    before.start();
#ifdef LOG4CXX
    if(logger->isDebugEnabled()) LOGPZ_DEBUG(logger,"TPZSymetricSpStructMatrixEigen::CreateAssemble calling Assemble()");
#endif
    Assemble(*stiff,rhs,guiInterface);
    mat->ComputeDiagonal();
    
    //    std::cout << "Rhs norm " << Norm(rhs) << std::endl;
    
    before.stop();
    //std::cout << __PRETTY_FUNCTION__ << " " << before << std::endl;
    //    mat->ComputeDiagonal();
    //stiff->Print("Stiffness TPZFYsmpMatrix :: CreateAssemble()");
#ifdef LOG4CXX
    if(logger->isDebugEnabled()) LOGPZ_DEBUG(logger,"TPZSymetricSpStructMatrixEigen::CreateAssemble exiting");
#endif
    return stiff;
}
TPZMatrix<STATE> * TPZSymetricSpStructMatrixEigen::Create(){
    
    /**
     *Longhin implementation
     */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    fMesh->ComputeElGraph(elgraph,elgraphindex);
    TPZMatrix<STATE> * mat = SetupMatrixData(elgraph, elgraphindex);
    return mat;
}

TPZMatrix<STATE> * TPZSymetricSpStructMatrixEigen::SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex){
    
    int64_t neq = fEquationFilter.NActiveEquations();
    TPZSYsmpMatrixEigen<STATE> * mat = new TPZSYsmpMatrixEigen<STATE>(neq,neq);
    
    /**Creates a element graph*/
    TPZMetis metis;
    metis.SetElementsNodes(elgraphindex.NElements() -1 ,fMesh->NIndependentConnects());
    metis.SetElementGraph(elgraph,elgraphindex);
    
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    metis.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    int64_t i;
    int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(i=0;i<nblock;i++){
        int64_t iblsize = fMesh->Block().Size(i);
        int64_t iblpos = fMesh->Block().Position(i);
        int64_t numactive = fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        totaleq += iblsize;
        int64_t icfirst = nodegraphindex[i];
        int64_t iclast = nodegraphindex[i+1];
        int64_t j;
        //longhin
        totalvar+=(iblsize*(iblsize+1))/2;
        for(j=icfirst;j<iclast;j++) {
            int64_t col = nodegraph[j];
            if (col < i) {
                continue;
            }
            
            if (col == i) {
                DebugStop();
            }
            
            int64_t colsize = fMesh->Block().Size(col);
            int64_t colpos = fMesh->Block().Position(col);
            int64_t numactive = fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
            totalvar += iblsize*colsize;
        }
    }
    
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
    
    nblock=fMesh->NIndependentConnects();
    
    // number of non-zeros
    int64_t nnzeros = 0;
    TPZVec<int64_t> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<STATE> EqValue(totalvar,0.);
    for(i=0;i<nblock;i++){
        int64_t iblsize = fMesh->Block().Size(i);
        int64_t iblpos = fMesh->Block().Position(i);
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t i=0; i<iblsize; i++) {
            rowdestindices[i] = iblpos+i;
        }
        fEquationFilter.Filter(rowdestindices);
        
        int64_t ibleq;
        // working equation by equation
        for(ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
            Eq[ieq] = pos;
            int64_t colsize,colpos,jbleq;
            int64_t diagonalinsert = 0;
            int64_t icfirst = nodegraphindex[i];
            int64_t iclast = nodegraphindex[i+1];
            int64_t j;
            for(j=icfirst;j<iclast;j++)
            {
                int64_t col = nodegraph[j];
                if (col < i) {
                    continue;
                }
                // force the diagonal block to be inserted
                // the nodegraph does not contain the pointer to itself
                if(!diagonalinsert && col > i)
                {
                    diagonalinsert = 1;
                    int64_t colsize = fMesh->Block().Size(i);
                    int64_t colpos = fMesh->Block().Position(i);
                    TPZManVector<int64_t> destindices(colsize);
                    for (int64_t i=0; i<colsize; i++) {
                        destindices[i] = colpos+i;
                    }
                    fEquationFilter.Filter(destindices);
                    int64_t jbleq;
                    for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                        //             if(colpos+jbleq == ieq) continue;
                        int64_t jeq = destindices[jbleq];
                        if (jeq < ieq) {
                            nnzeros++;
                            continue;
                        }
                        EqCol[pos] = destindices[jbleq];
                        EqValue[pos] = 0.;
                        //            colpos++;
                        pos++;
                        nnzeros++;
                    }
                }
                colsize = fMesh->Block().Size(col);
                colpos = fMesh->Block().Position(col);
                if (fEquationFilter.NumActive(colpos, colpos+colsize) == 0) {
                    continue;
                }
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
                for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                    int64_t jeq = destindices[jbleq];
                    if (jeq < ieq) {
                        nnzeros++;
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    colpos++;
                    pos++;
                    nnzeros++;
                }
            }
            // all elements are below (last block certainly)
            if(!diagonalinsert)
            {
                diagonalinsert = 1;
                int64_t colsize = fMesh->Block().Size(i);
                int64_t colpos = fMesh->Block().Position(i);
                TPZManVector<int64_t> destindices(colsize);
                for (int64_t i=0; i<colsize; i++) {
                    destindices[i] = colpos+i;
                }
                fEquationFilter.Filter(destindices);
                int64_t jbleq;
                for(jbleq=0; jbleq<destindices.size(); jbleq++) {
                    //             if(colpos+jbleq == ieq) continue;
                    int64_t jeq = destindices[jbleq];
                    if (jeq < ieq) {
                        nnzeros++;
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    //            colpos++;
                    nnzeros++;
                    pos++;
                }
            }
            ieq++;
        }
    }
    
    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
//    m_triplets.resize(nnzeros);
    m_triplets.clear();
    
    return mat;
}

TPZSymetricSpStructMatrixEigen::TPZSymetricSpStructMatrixEigen() : TPZStructMatrix(){
    
}

TPZSymetricSpStructMatrixEigen::TPZSymetricSpStructMatrixEigen(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{}

void TPZSymetricSpStructMatrixEigen::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    
    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matRed = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (&stiffness);
    if (mat) {
        Serial_AssembleGlob(stiffness,rhs, guiInterface);
    }
    if (matRed) {
        Serial_AssembleSub(stiffness,rhs, guiInterface);
    }
    
}
void TPZSymetricSpStructMatrixEigen::Serial_AssembleSub(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matRed = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (&stiffness);
    TPZMatrix<STATE> * matpzmat=matRed->K00().operator->();
    TPZSYsmpMatrixEigen<STATE> *matk00eigen =dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmat);
    
    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
    
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
    int64_t count = 0;
    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = fMaterialIds.size();
        if(matidsize){
            if(!el->NeedsComputing(fMaterialIds)) continue;
        }
        
        count++;
        if (!(count % 1000)) {
            std::cout << '*';
            std::cout.flush();
        }
        if (!(count % 20000)) {
            std::cout << "\n";
        }
        calcstiff.start();
        ek.Reset();
        ef.Reset();
        el->CalcStiff(ek, ef);
        if (guiInterface) if (guiInterface->AmIKilled()) {
            return;
        }
        
        
        calcstiff.stop();
        assemble.start();
        
        if (!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            
        }
        
        assemble.stop();
    }//fim for iel
    matk00eigen->SetFromTriplets();
    
}
void TPZSymetricSpStructMatrixEigen::Serial_AssembleGlob(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
    
        TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
    
        int64_t nelem = fMesh->NElements();
    
    #ifdef LOG4CXX
        bool globalresult = true;
        bool writereadresult = true;
    #endif
        TPZTimer calcstiff("Computing the stiffness matrices");
        TPZTimer assemble("Assembling the stiffness matrices");
        TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
        // scan all the mesh to compute nnzeros.
        int64_t nnzeros = 0; // counter for nonzeros;
        std::vector<int64_t> iel_to_l_index(nelem);
    
        TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
    
        for (int64_t iel = 0; iel < nelem; iel++) {
            TPZCompEl *el = elementvec[iel];
            if (!el) continue;
            int matid = 0;
            TPZGeoEl *gel = el->Reference();
            if (gel) {
                matid = gel->MaterialId();
            }
            int matidsize = fMaterialIds.size();
            if(matidsize){
                if(!el->NeedsComputing(fMaterialIds)) continue;
            }
    
            ek.Reset();
            ef.Reset();
            el->CalcStiff(ek, ef);
            if (guiInterface) if (guiInterface->AmIKilled()) {
                return;
            }
    
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            iel_to_l_index[iel] = (nnzeros);
            nnzeros += ek.fSourceIndex.NElements() * ek.fSourceIndex.NElements();
            // needs to separate rhs
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
        }
        m_triplets.resize(nnzeros);
    
        #ifdef USING_TBB
            tbb::parallel_for(size_t(0), size_t(nelem), size_t(1),
                              [this,&elementvec,&iel_to_l_index] (size_t & iel){
    
                TPZCompEl *el = elementvec[iel];
                if (el) {
                    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
    
                    ek.Reset();
                    ef.Reset();
                    el->CalcStiff(ek, ef);
    
                    ek.ComputeDestinationIndices();
                    fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
                    int64_t i,j = 0;
                    REAL value=0.;
                    int64_t ipos,jpos;
    
                    int64_t l=0;
                    for(i=0;i<ek.fSourceIndex.NElements();i++){
                        for(j=0;j<ek.fSourceIndex.NElements();j++){
                            ipos=ek.fDestinationIndex[i];
                            jpos=ek.fDestinationIndex[j];
                            value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
                            Triplet3<REAL> trip(ipos, jpos, value);
                            m_triplets[iel_to_l_index[iel] + l] = trip;
                            l++;
                        }
                    }
                }
            }
        );
        #else
            for (iel = 0; iel < nelem; iel++) {
    
                TPZCompEl *el = elementvec[iel];
                if (!el) continue;
                int matid = 0;
                TPZGeoEl *gel = el->Reference();
                if (gel) {
                    matid = gel->MaterialId();
                }
                int matidsize = fMaterialIds.size();
                if(matidsize){
                    if(!el->NeedsComputing(fMaterialIds)) continue;
                }
    
                ek.Reset();
                ef.Reset();
                el->CalcStiff(ek, ef);
                if (guiInterface) if (guiInterface->AmIKilled()) {
                    return;
                }
    
                ek.ComputeDestinationIndices();
                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
                int64_t i,j = 0;
                REAL value=0.;
                int64_t ipos,jpos;
    
                int64_t l=0;
                for(i=0;i<ek.fSourceIndex.NElements();i++){
                    for(j=0;j<ek.fSourceIndex.NElements();j++){
                        ipos=ek.fDestinationIndex[i];
                        jpos=ek.fDestinationIndex[j];
                        value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
                        Triplet3<REAL> trip(ipos, jpos, value);
                        m_triplets[iel_to_l_index[iel] + l] = trip;
                        l++;
                    }
                }
            }
        #endif
    
        mat->fsparse_eigen.setFromTriplets(m_triplets.begin(), m_triplets.end());
        m_triplets.clear();
}
//void TPZSymetricSpStructMatrixEigen::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
//
//
//    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
//
//    int64_t nelem = fMesh->NElements();
//
//#ifdef LOG4CXX
//    bool globalresult = true;
//    bool writereadresult = true;
//#endif
//    TPZTimer calcstiff("Computing the stiffness matrices");
//    TPZTimer assemble("Assembling the stiffness matrices");
//    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
//
//    // scan all the mesh to compute nnzeros.
//    int64_t nnzeros = 0; // counter for nonzeros;
//    std::vector<int64_t> iel_to_l_index(nelem);
//
//    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
//    if (mat) {
//
//    for (int64_t iel = 0; iel < nelem; iel++) {
//        TPZCompEl *el = elementvec[iel];
//        if (!el) continue;
//        int matid = 0;
//        TPZGeoEl *gel = el->Reference();
//        if (gel) {
//            matid = gel->MaterialId();
//        }
//        int matidsize = fMaterialIds.size();
//        if(matidsize){
//            if(!el->NeedsComputing(fMaterialIds)) continue;
//        }
//
//        ek.Reset();
//        ef.Reset();
//        el->CalcStiff(ek, ef);
//        if (guiInterface) if (guiInterface->AmIKilled()) {
//            return;
//        }
//
//        ek.ComputeDestinationIndices();
//        fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//        iel_to_l_index[iel] = (nnzeros);
//        nnzeros += ek.fSourceIndex.NElements() * ek.fSourceIndex.NElements();
//        // needs to separate rhs
//        rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//    }
//    m_triplets.resize(nnzeros);
//
//    #ifdef USING_TBB
//        tbb::parallel_for(size_t(0), size_t(nelem), size_t(1),
//                          [this,&elementvec,&iel_to_l_index] (size_t & iel){
//
//            TPZCompEl *el = elementvec[iel];
//            if (el) {
//                TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
//
//                ek.Reset();
//                ef.Reset();
//                el->CalcStiff(ek, ef);
//
//                ek.ComputeDestinationIndices();
//                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                int64_t i,j = 0;
//                REAL value=0.;
//                int64_t ipos,jpos;
//
//                int64_t l=0;
//                for(i=0;i<ek.fSourceIndex.NElements();i++){
//                    for(j=0;j<ek.fSourceIndex.NElements();j++){
//                        ipos=ek.fDestinationIndex[i];
//                        jpos=ek.fDestinationIndex[j];
//                        value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
//                        Triplet3<REAL> trip(ipos, jpos, value);
//                        m_triplets[iel_to_l_index[iel] + l] = trip;
//                        l++;
//                    }
//                }
//            }
//        }
//    );
//    #else
//        for (iel = 0; iel < nelem; iel++) {
//
//            TPZCompEl *el = elementvec[iel];
//            if (!el) continue;
//            int matid = 0;
//            TPZGeoEl *gel = el->Reference();
//            if (gel) {
//                matid = gel->MaterialId();
//            }
//            int matidsize = fMaterialIds.size();
//            if(matidsize){
//                if(!el->NeedsComputing(fMaterialIds)) continue;
//            }
//
//            ek.Reset();
//            ef.Reset();
//            el->CalcStiff(ek, ef);
//            if (guiInterface) if (guiInterface->AmIKilled()) {
//                return;
//            }
//
//            ek.ComputeDestinationIndices();
//            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//            int64_t i,j = 0;
//            REAL value=0.;
//            int64_t ipos,jpos;
//
//            int64_t l=0;
//            for(i=0;i<ek.fSourceIndex.NElements();i++){
//                for(j=0;j<ek.fSourceIndex.NElements();j++){
//                    ipos=ek.fDestinationIndex[i];
//                    jpos=ek.fDestinationIndex[j];
//                    value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
//                    Triplet3<REAL> trip(ipos, jpos, value);
//                    m_triplets[iel_to_l_index[iel] + l] = trip;
//                    l++;
//                }
//            }
//        }
//    #endif
//
//    mat->fsparse_eigen.setFromTriplets(m_triplets.begin(), m_triplets.end());
//    m_triplets.clear();
//    }
//    else{
//            int64_t iel;
//            int64_t nelem = fMesh->NElements();
//            TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
//
//            TPZTimer calcstiff("Computing the stiffness matrices");
//            TPZTimer assemble("Assembling the stiffness matrices");
//            TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
//
//            int64_t count = 0;
//            for (iel = 0; iel < nelem; iel++) {
//                TPZCompEl *el = elementvec[iel];
//                if (!el) continue;
//                int matid = 0;
//                TPZGeoEl *gel = el->Reference();
//                if (gel) {
//                    matid = gel->MaterialId();
//                }
//                int matidsize = fMaterialIds.size();
//                if(matidsize){
//                    if(!el->NeedsComputing(fMaterialIds)) continue;
//                }
//
//                count++;
//                if (!(count % 1000)) {
//                    std::cout << '*';
//                    std::cout.flush();
//                }
//                if (!(count % 20000)) {
//                    std::cout << "\n";
//                }
//                calcstiff.start();
//                ek.Reset();
//                ef.Reset();
//                el->CalcStiff(ek, ef);
//                if (guiInterface) if (guiInterface->AmIKilled()) {
//                    return;
//                }
//
//
//                calcstiff.stop();
//                assemble.start();
//
//                if (!ek.HasDependency()) {
//                    ek.ComputeDestinationIndices();
//                    fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                    stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//
//                    rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//                } else {
//                    // the element has dependent nodes
//                    ek.ApplyConstraints();
//                    ef.ApplyConstraints();
//                    ek.ComputeDestinationIndices();
//                    fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                    stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//                    rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//
//                }
//
//                assemble.stop();
//            }//fim for iel
//            if (count > 1000) std::cout << std::endl;
//                TPZMatRed<STATE, TPZFMatrix<STATE> > *mat = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (&stiffness);
//                if (mat) {
//                    int val=1;
//                    TPZMatrix<STATE> * matpzmat = mat->K00().operator->();
//                    TPZSYsmpMatrixEigen<STATE> *matk00eigen = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmat);
//
//                    TPZMatrix<STATE> * matpzmatk11 = &mat->K11();
//                    TPZSYsmpMatrixEigen<STATE> *matk11eigen = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmatk11);
//
//                    TPZMatrix<STATE> * matpzmatk01 = &mat->K01();
//                    TPZSYsmpMatrixEigen<STATE> *matk01eigen = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmatk01);
//
//                    TPZMatrix<STATE> * matpzmatk10 = &mat->K11();
//                    TPZSYsmpMatrixEigen<STATE> *matk10eigen = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmatk10);
//                    if(matk00eigen){
//                       matk00eigen->SetFromTriplets();
//                    }
//                    if (matk11eigen) {
//                        matk11eigen->SetFromTriplets();
//                    }
//                    if (matk01eigen) {
//                        matk01eigen->SetFromTriplets();
//                    }
//                    if (matk10eigen) {
//                        matk10eigen->SetFromTriplets();
//                    }
//                }
//    }
//    //
//}
//

//
//void  TPZSymetricSpStructMatrixEigen::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
//
//    int64_t iel;
//    int64_t nelem = fMesh->NElements();
//    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
//
//    TPZTimer calcstiff("Computing the stiffness matrices");
//    TPZTimer assemble("Assembling the stiffness matrices");
//    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
//
//    int64_t count = 0;
//    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
//    TPZMatRed<STATE, TPZFMatrix<STATE> > *matRed = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (&stiffness);
//    TPZMatrix<STATE> * matpzmat=0;
//    TPZSYsmpMatrixEigen<STATE> *matk00eigen =0;
//    if (matRed) {
//        matpzmat = matRed->K00().operator->();
//        matk00eigen = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *>(matpzmat);
//    }
//   
//    
//    if (mat) {
//        for (iel = 0; iel < nelem; iel++) {
//            TPZCompEl *el = elementvec[iel];
//            if (!el) continue;
//            int matid = 0;
//            TPZGeoEl *gel = el->Reference();
//            if (gel) {
//                matid = gel->MaterialId();
//            }
//            int matidsize = fMaterialIds.size();
//            if(matidsize){
//                if(!el->NeedsComputing(fMaterialIds)) continue;
//            }
//            
//            count++;
//            if (!(count % 1000)) {
//                std::cout << '*';
//                std::cout.flush();
//            }
//            if (!(count % 20000)) {
//                std::cout << "\n";
//            }
//            calcstiff.start();
//            ek.Reset();
//            ef.Reset();
//            el->CalcStiff(ek, ef);
//            if (guiInterface) if (guiInterface->AmIKilled()) {
//                return;
//            }
//            
//            
//            calcstiff.stop();
//            assemble.start();
//            
//            if (!ek.HasDependency()) {
//                ek.ComputeDestinationIndices();
//                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//            
//                int64_t i,j = 0;
//                REAL value=0.;
//                int64_t ipos,jpos;
//                for(i=0;i<ek.fSourceIndex.NElements();i++){
//                    for(j=0;j<ek.fSourceIndex.NElements();j++){
//                        ipos=ek.fDestinationIndex[i];
//                        jpos=ek.fDestinationIndex[j];
//                        value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
//                        Triplet3<REAL> trip(ipos, jpos, value);
//                        mat->m_triplets.push_back(trip);
//                    }
//                }
//            } else {
//                // the element has dependent nodes
//                ek.ApplyConstraints();
//                ef.ApplyConstraints();
//                ek.ComputeDestinationIndices();
//                fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                int64_t i,j = 0;
//                REAL value=0.;
//                int64_t ipos,jpos;
//                for(i=0;i<ek.fSourceIndex.NElements();i++){
//                    for(j=0;j<ek.fSourceIndex.NElements();j++){
//                        ipos=ek.fDestinationIndex[i];
//                        jpos=ek.fDestinationIndex[j];
//                        value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
//                        Triplet3<REAL> trip(ipos, jpos, value);
//                        mat->m_triplets.push_back(trip);
//                    }
//                }
//                
//            }
//            
//            assemble.stop();
//        }//fim for iel
//        if (count > 1000) std::cout << std::endl;
//        mat->hastriplets =1;
//        mat->SetFromTriplets();
//    }
//    else{
//                    int64_t iel;
//                    int64_t nelem = fMesh->NElements();
//                    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
//        
//                    TPZTimer calcstiff("Computing the stiffness matrices");
//                    TPZTimer assemble("Assembling the stiffness matrices");
//                    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
//        
//                    int64_t count = 0;
//                    for (iel = 0; iel < nelem; iel++) {
//                        TPZCompEl *el = elementvec[iel];
//                        if (!el) continue;
//                        int matid = 0;
//                        TPZGeoEl *gel = el->Reference();
//                        if (gel) {
//                            matid = gel->MaterialId();
//                        }
//                        int matidsize = fMaterialIds.size();
//                        if(matidsize){
//                            if(!el->NeedsComputing(fMaterialIds)) continue;
//                        }
//        
//                        count++;
//                        if (!(count % 1000)) {
//                            std::cout << '*';
//                            std::cout.flush();
//                        }
//                        if (!(count % 20000)) {
//                            std::cout << "\n";
//                        }
//                        calcstiff.start();
//                        ek.Reset();
//                        ef.Reset();
//                        el->CalcStiff(ek, ef);
//                        if (guiInterface) if (guiInterface->AmIKilled()) {
//                            return;
//                        }
//        
//        
//                        calcstiff.stop();
//                        assemble.start();
//        
//                        if (!ek.HasDependency()) {
//                            ek.ComputeDestinationIndices();
//                            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//        
//                            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
//                        } else {
//                            // the element has dependent nodes
//                            ek.ApplyConstraints();
//                            ef.ApplyConstraints();
//                            ek.ComputeDestinationIndices();
//                            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
//                            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//                            rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//        
//                        }
//        
//                        assemble.stop();
//                    }//fim for iel
//                    matk00eigen->SetFromTriplets();
//    }
//    
//}


#ifndef STATE_COMPLEX
#include "pzmat2dlin.h"

int TPZSymetricSpStructMatrixEigen::main() {
    
    int refine=5;
    int order=5;
    
    TPZGeoMesh gmesh;
    TPZCompMesh cmesh(&gmesh);
    double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},
        {0.,1.,0.}};
    
    int i,j;
    TPZVec<REAL> coord(3,0.);
    for(i=0; i<4; i++) {
        // initializar as coordenadas do no em um vetor
        for (j=0; j<3; j++) coord[j] = coordstore[i][j];
        
        // identificar um espa� no vetor onde podemos armazenar
        // este vetor
        gmesh.NodeVec ().AllocateNewElement ();
        
        // initializar os dados do n�       gmesh.NodeVec ()[nodeindex].Initialize (i,coord,gmesh);
    }
    int el;
    TPZGeoEl *gel;
    for(el=0; el<1; el++) {
        
        // initializar os indices dos nos
        TPZVec<int64_t> indices(4);
        for(i=0; i<4; i++) indices[i] = i;
        // O proprio construtor vai inserir o elemento na malha
        //       gel = new TPZGeoElQ2d(el,indices,1,gmesh);
        int64_t index;
        gel = gmesh.CreateGeoElement(EQuadrilateral,indices,1,index);
    }
    gmesh.BuildConnectivity ();
    
    TPZVec<TPZGeoEl *> subel;
    
    cout << "Refinement ";
    cin >> refine;
    cout << endl;
    DebugStop();
    //    UniformRefine(refine,gmesh);
    
    TPZGeoElBC gelbc(gel,4,-4);
    TPZMat2dLin *meumat = new TPZMat2dLin(1);
    TPZFMatrix<STATE> xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
    meumat->SetMaterial (xk,xc,xf);
    TPZMaterial * meumatptr(meumat);
    cmesh.InsertMaterialObject(meumatptr);
    
    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
    TPZMaterial * bnd = meumat->CreateBC (meumatptr,-4,0,val1,val2);
    cmesh.InsertMaterialObject(bnd);
    
    cout << "Interpolation order ";
    cin >> order;
    cout << endl;
    
    //TPZCompEl::gOrder = order;
    cmesh.SetDefaultOrder(order);
    
    cmesh.AutoBuild();
    //    cmesh.AdjustBoundaryElements();
    cmesh.InitializeBlock();
    
    ofstream output("outputPar.dat");
    TPZAnalysis an(&cmesh,true,output);
    
    //TPZVec<int> numelconnected(cmesh.NEquations(),0);
    TPZSymetricSpStructMatrixEigen mat(&cmesh);
    
    an.SetStructuralMatrix(mat);
    
    TPZStepSolver<STATE> sol;
    sol.SetJacobi(100,1.e-5,0);
    
    
    an.SetSolver(sol);
    an.Run(output);
    output.flush();
    cout.flush();
    return 0;
    
}
#endif
