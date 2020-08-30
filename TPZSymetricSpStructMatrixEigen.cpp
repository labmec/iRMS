//
//  TPZSymetricSpStructMatrixEigenEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/27/20.
//

#include "TPZSymetricSpStructMatrixEigen.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpMatrixEigen.h"

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
                            continue;
                        }
                        EqCol[pos] = destindices[jbleq];
                        EqValue[pos] = 0.;
                        //            colpos++;
                        pos++;
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
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    colpos++;
                    pos++;
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
                        continue;
                    }
                    EqCol[pos] = jeq;
                    EqValue[pos] = 0.;
                    //            colpos++;
                    pos++;
                }
            }
            ieq++;
        }
    }
    
    Eq[ieq] = pos;
    mat->SetData(Eq,EqCol,EqValue);
    return mat;
}

TPZSymetricSpStructMatrixEigen::TPZSymetricSpStructMatrixEigen() : TPZStructMatrix(){
    
}

TPZSymetricSpStructMatrixEigen::TPZSymetricSpStructMatrixEigen(TPZCompMesh *mesh) : TPZStructMatrix(mesh)
{}


void TPZSymetricSpStructMatrixEigen::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) {
#ifdef PZDEBUG
    TExceptionManager activateExceptions;
#endif
    if (!fMesh) {
        LOGPZ_ERROR(logger, "Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef LOG4CXX
    if (loggerelmat->isDebugEnabled()) {
        if (dynamic_cast<TPZSubCompMesh *> (fMesh)) {
            std::stringstream sout;
            sout << "AllEig = {};";
            LOGPZ_DEBUG(loggerelmat, sout.str())
        }
    }
#endif
    
#ifdef PZDEBUG
    if (rhs.Rows() != fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif
    
    std::vector<Triplet3<REAL> > triplets;
    TPZSYsmpMatrixEigen<STATE> *mat = dynamic_cast<TPZSYsmpMatrixEigen<STATE> *> (&stiffness);
    if (!mat) {
        DebugStop();
    }
    else{
        mat->Zero();
    }
//    int nnonzeros = mat->fsparse_eigen.nonZeros();
//    triplets.reserve(nnonzeros);
    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK), ef(fMesh, TPZElementMatrix::EF);
#ifdef LOG4CXX
    bool globalresult = true;
    bool writereadresult = true;
#endif
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
        
#ifdef LOG4CXX
        if (loggerelmat->isDebugEnabled()) {
            if (dynamic_cast<TPZSubCompMesh *> (fMesh)) {
                std::stringstream objname;
                objname << "Element" << iel;
                std::string name = objname.str();
                objname << " = ";
                std::stringstream sout;
                ek.fMat.Print(objname.str().c_str(), sout, EMathematicaInput);
                sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";
                
                LOGPZ_DEBUG(loggerelmat, sout.str())
           
            }
        }
#endif
        
#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if (!result) {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck, sout.str())
        }
        
#endif
        
        calcstiff.stop();
        assemble.start();
        
        if (!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            int64_t i,j = 0;
            REAL value=0.;
            int64_t ipos,jpos;
   
            for(i=0;i<ek.fSourceIndex.NElements();i++){
                for(j=0;j<ek.fSourceIndex.NElements();j++){
                    ipos=ek.fDestinationIndex[i];
                    jpos=ek.fDestinationIndex[j];
                    value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
                    Triplet3<REAL> trip(ipos, jpos, value);
                    triplets.push_back(trip);
                }
            }
            
//            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            
#ifdef PZDEBUG
            REAL rhsnorm = Norm(ef.fMat);
            REAL eknorm = Norm(ek.fMat);
            if (rhsnorm > 10e15) {
                DebugStop();
            }
            if (std::isnan(rhsnorm) || std::isnan(eknorm)) {
                std::cout << "element " << iel << " has norm " << rhsnorm << std::endl;
                el->Print();
                ek.fMat.Print("ek",std::cout);
                ef.fMat.Print("ef",std::cout);
            }
#endif
            rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            
#ifdef LOG4CXX
            if (loggerel->isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element without associated geometric element index " << el->Index() << "\n";
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            int64_t nelem = ek.fSourceIndex.NElements();
            int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
            REAL prevval;
//            if(IsSimetric()) {
//                for(icoef=0; icoef<nelem; icoef++) {
//                    ieq = ek.fDestinationIndex[icoef];
//                    ieqs = ek.fSourceIndex[icoef];
//                    for(jcoef=icoef; jcoef<nelem; jcoef++) {
//                        jeq = ek.fDestinationIndex[jcoef];
//                        jeqs = ek.fSourceIndex[jcoef];
//                        prevval = ek.fMat.GetVal(ieq,jeq);
//                        prevval += ek.fMat(ieqs,jeqs);
//                        Triplet3<REAL> trip(ieq,jeq,prevval);
//                        triplets.push_back(trip);
////                        PutVal(ieq,jeq,prevval);
//                    }
//                }
//            }
            
            
//            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
//            int64_t i,j = 0;
//            REAL value=0.;
//            int64_t ipos,jpos;
//            ek.fMat.Print("ek= ",std::cout, EMathematicaInput);
//            std::cout<<std::endl;
//            ek.fSourceIndex.Print(std::cout);
//            std::cout<<std::endl;
//            ek.fDestinationIndex.Print(std::cout);
//            for(i=0;i<ek.fSourceIndex.NElements();i++){
//                for(j=0;j<ek.fSourceIndex.NElements();j++){
//                    ipos=ek.fDestinationIndex[i];
//                    jpos=ek.fDestinationIndex[j];
//                    value=ek.fMat.GetVal(ek.fSourceIndex[i],ek.fSourceIndex[j]);
//                    Triplet3<REAL> trip(ipos, jpos, value);
//                    triplets.push_back(trip);
//                }
//            }
            rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            
#ifdef LOG4CXX
            if (loggerel->isDebugEnabled()) {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                
                if (gel) {
                    TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                    gel->CenterPoint(gel->NSides() - 1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                } else {
                    sout << "Stiffness for computational element index " << iel << std::endl;
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel, sout.str())
            }
#endif
        }
       
        
       
        assemble.stop();
        
    }//fim for iel
//    std::cout<<mat->fsparse_eigen;
     mat->fsparse_eigen.setFromTriplets(triplets.begin(), triplets.end());
    if (count > 1000) std::cout << std::endl;
    
#ifdef LOG4CXX
    if (loggerCheck->isDebugEnabled()) {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
    if (loggerGlobStiff->isDebugEnabled())
    {
        std::stringstream sout;
        stiffness.Print("GK = ",sout,EMathematicaInput);
        rhs.Print("GR = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
    
#endif
    
}

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
