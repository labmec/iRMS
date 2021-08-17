//
//  TPZSparceMatrixEigen.cpp
//  ALL_BUILD
//
//  Created by Jose on 7/17/20.
//
#include "TPZAnalysisAuxEigen.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif

void TPZAnalysisAuxEigen::AssembleMass(){
    int n_cells = fAlgebraicTransport->fCellsData.fVolume.size();
    m_mass =  Eigen::SparseMatrix<REAL>( n_cells, n_cells );
    
    m_mass_triplets.resize(n_cells);
    //Volumetric Elements
    
#ifdef USING_TBB2
    tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                      [this] (size_t & ivol){
                          int eqindex = fAlgebraicTransport->fCellsData.fEqNumber[ivol];
                          TPZFMatrix<double> elmat, ef;
                          elmat.Resize(1, 1);
                          ef.Resize(1,1);
                          fAlgebraicTransport->Contribute(ivol, elmat, ef);
                          m_mass_triplets[ivol] = Triplet2<REAL>(eqindex,eqindex, elmat(0,0));
                      }
                      );
#else
    for (int ivol = 0; ivol < n_cells; ivol++) {
        int eqindex = fAlgebraicTransport->fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1,1);
        fAlgebraicTransport->Contribute(ivol, elmat, ef);
//        std::cout<<"mass: "<<elmat(0,0)<<std::endl;
        m_mass_triplets[ivol] = Triplet2<REAL>(eqindex,eqindex, elmat(0,0));
    }
#endif
    
    m_mass.setFromTriplets( m_mass_triplets.begin(), m_mass_triplets.end() );
    m_mass_triplets.clear();
}

void TPZAnalysisAuxEigen::Assemble(){
    int n_cells = fAlgebraicTransport->fCellsData.fVolume.size();
    m_transmissibility =  Eigen::SparseMatrix<REAL>( n_cells, n_cells );
    m_rhs =  Eigen::SparseMatrix<REAL>( n_cells, 1 );
//    if(!this->fSolver){
//        DebugStop();
//    }
    
    m_rhs.setZero();
    m_transmissibility.setZero();
    
    int internal_faces_id = 100;
    int internal_faces_id1 = 101;
    int internal_faces_id2 = 102;
    int internal_faces_id3 = 103;
    int outletfrac_id = 104;
    int inlet_faces_id = -2;
    int outlet_faces_id = -4;
    
    int n_internal_faces = fAlgebraicTransport->fInterfaceData[internal_faces_id].fFluxSign.size();
    int n_internal_faces1 = fAlgebraicTransport->fInterfaceData[internal_faces_id1].fFluxSign.size();
    int n_internal_faces2 = fAlgebraicTransport->fInterfaceData[internal_faces_id2].fFluxSign.size();
    int n_internal_faces3 = fAlgebraicTransport->fInterfaceData[internal_faces_id3].fFluxSign.size();
    int n_outletfrac = fAlgebraicTransport->fInterfaceData[outletfrac_id].fFluxSign.size();
    
    int n_inlet_faces = fAlgebraicTransport->fInterfaceData[inlet_faces_id].fFluxSign.size();
    int n_outlet_faces = fAlgebraicTransport->fInterfaceData[outlet_faces_id].fFluxSign.size();
    
    size_t n_nzeros_trans =( (n_internal_faces + n_internal_faces1 +n_internal_faces2 + n_internal_faces3) * 4) + n_outlet_faces + n_outletfrac;
    size_t n_nzeros_res = n_cells + 2*(n_internal_faces + n_internal_faces1  + n_internal_faces2 + n_internal_faces3) +n_inlet_faces + n_outlet_faces + n_outletfrac ;
    m_trans_triplets.resize(n_nzeros_trans);
    m_rhs_triplets.resize(n_nzeros_res);


#ifdef USING_TBB
    
    tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                      [this] (size_t & ivol){
                          int eqindex = fAlgebraicTransport->fCellsData.fEqNumber[ivol];
                          TPZFMatrix<double> ef;
                          ef.Resize(1,1);
                          fAlgebraicTransport->ContributeResidual(ivol, ef);
                          m_rhs_triplets[ivol] = Triplet2<REAL>(eqindex,0, ef(0,0));
                      }
                      );
    
    
    tbb::parallel_for(size_t(0), size_t(n_internal_faces), size_t(1),
                      [this,&internal_faces_id,&n_cells] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int right = lrindex.second;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
                          
                          TPZVec<int64_t> indexes(2);
                          indexes[0]=lefteq;
                          indexes[1]=righteq;
                          TPZFMatrix<double> elmat, ef;
                          elmat.Resize(2, 2);
                          ef.Resize(2, 1);
                          fAlgebraicTransport->ContributeInterface(iface,elmat, ef);
#ifdef PZDEBUG
                          {
                              STATE normek = Norm(elmat);
                              STATE normef = Norm(ef);
                              if(std::isnan(normek) || std::isnan(normef))
                              {
                                  std::cout << "matriz ou residuo nan";
                                  elmat.Print("elmat");
                                  ef.Print("ef");
                                  DebugStop();
                              }
                          }
#endif
                          size_t i_begin = 2*2*(iface);
                          m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
                          m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
                          m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
                          m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));
                          
                          size_t i_rhs_begin = 2*(iface) + n_cells;
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                          m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
                      }
                      );
    // Frac-Vol1 Interfaces
    
    tbb::parallel_for(size_t(0), size_t(n_internal_faces1), size_t(1),
                      [this,&internal_faces_id1,&n_internal_faces ,&n_cells] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id1].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int right = lrindex.second;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
                          
                          TPZVec<int64_t> indexes(2);
                          indexes[0]=lefteq;
                          indexes[1]=righteq;
                          TPZFMatrix<double> elmat, ef;
                          elmat.Resize(2, 2);
                          ef.Resize(2, 1);
                          fAlgebraicTransport->ContributeInterface(iface,elmat, ef);
#ifdef PZDEBUG
                          {
                              STATE normek = Norm(elmat);
                              STATE normef = Norm(ef);
                              if(std::isnan(normek) || std::isnan(normef))
                              {
                                  std::cout << "matriz ou residuo nan";
                                  elmat.Print("elmat");
                                  ef.Print("ef");
                                  DebugStop();
                              }
                          }
#endif
                          size_t i_begin = 2*2*(iface) + 4*n_internal_faces;
                          m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
                          m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
                          m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
                          m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));
                          
                          size_t i_rhs_begin = 2*(iface) + n_cells;
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                          m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
                      }
                      );
    
     // Frac-Vol2 Interfaces
    tbb::parallel_for(size_t(0), size_t(n_internal_faces2), size_t(1),
                      [this,&internal_faces_id2,&n_cells, &n_internal_faces , &n_internal_faces1] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id2].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int right = lrindex.second;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
                          
                          TPZVec<int64_t> indexes(2);
                          indexes[0]=lefteq;
                          indexes[1]=righteq;
                          TPZFMatrix<double> elmat, ef;
                          elmat.Resize(2, 2);
                          ef.Resize(2, 1);
                          fAlgebraicTransport->ContributeInterface(iface,elmat, ef);
#ifdef PZDEBUG
                          {
                              STATE normek = Norm(elmat);
                              STATE normef = Norm(ef);
                              if(std::isnan(normek) || std::isnan(normef))
                              {
                                  std::cout << "matriz ou residuo nan";
                                  elmat.Print("elmat");
                                  ef.Print("ef");
                                  DebugStop();
                              }
                          }
#endif
                          size_t i_begin = 2*2*(iface) + 4*n_internal_faces +4*n_internal_faces1 ;
                          m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
                          m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
                          m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
                          m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));
                          
                          size_t i_rhs_begin = 2*(iface) + n_cells + 2*(n_internal_faces);
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                          m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
                      }
                      );
    //Frac-Frac
    
    tbb::parallel_for(size_t(0), size_t(n_internal_faces3), size_t(1),
                      [this,&internal_faces_id3,&n_cells,&n_internal_faces ,&n_internal_faces1, &n_internal_faces2] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id3].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int right = lrindex.second;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
                          
                          TPZVec<int64_t> indexes(2);
                          indexes[0]=lefteq;
                          indexes[1]=righteq;
                          TPZFMatrix<double> elmat, ef;
                          elmat.Resize(2, 2);
                          ef.Resize(2, 1);
                          fAlgebraicTransport->ContributeInterface(iface,elmat, ef);
#ifdef PZDEBUG
                          {
                              STATE normek = Norm(elmat);
                              STATE normef = Norm(ef);
                              if(std::isnan(normek) || std::isnan(normef))
                              {
                                  std::cout << "matriz ou residuo nan";
                                  elmat.Print("elmat");
                                  ef.Print("ef");
                                  DebugStop();
                              }
                          }
#endif
                          size_t i_begin = 2*2*(iface) + 4*(n_internal_faces +n_internal_faces1 + n_internal_faces2) ;
                          m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
                          m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
                          m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
                          m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));
                          
                          size_t i_rhs_begin = 2*(iface) + n_cells + 2*(n_internal_faces + n_internal_faces1 + n_internal_faces2);
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                          m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
                      }
                      );
    
    tbb::parallel_for(size_t(0), size_t(n_inlet_faces), size_t(1),
                      [this,&inlet_faces_id,&n_internal_faces,&n_internal_faces1 , &n_internal_faces2, &n_internal_faces3,&n_cells] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          
                          TPZVec<int64_t> indexes(1);
                          indexes[0]=lefteq;
                          TPZFMatrix<double> ef;
                          ef.Resize(1, 1);
                          fAlgebraicTransport->ContributeBCInletInterface(iface,ef);
                          size_t i_rhs_begin = (iface) + n_cells + 2*(n_internal_faces + n_internal_faces1 + n_internal_faces2 + n_internal_faces3);
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                      }
                      );
    
    tbb::parallel_for(size_t(0), size_t(n_outlet_faces), size_t(1),
                      [this,&outlet_faces_id,&n_internal_faces,&n_inlet_faces,&n_cells, &n_internal_faces1, &n_internal_faces2, &n_internal_faces3] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport->fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          
                          TPZVec<int64_t> indexes(1);
                          indexes[0]=lefteq;
                          TPZFMatrix<double> elmat, ef;
                          elmat.Resize(1, 1);
                          ef.Resize(1, 1);
                          fAlgebraicTransport->ContributeBCOutletInterface(iface,elmat, ef);
                          size_t i_begin = iface +   + 4*(n_internal_faces + n_internal_faces1 + n_internal_faces2 + n_internal_faces3) ;
                          m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
                          
                          size_t i_rhs_begin = (iface) + n_cells + 2*(n_internal_faces + n_internal_faces1 + n_internal_faces2 + n_internal_faces3) + n_inlet_faces;
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                      }
                      );
    
#else
    
    for (int ivol = 0; ivol < n_cells; ivol++) {
        int eqindex = fAlgebraicTransport->fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> ef;
        ef.Resize(1,1);
        fAlgebraicTransport->ContributeResidual(ivol, ef);
        m_rhs_triplets[ivol] = Triplet2<REAL>(eqindex,0, ef(0,0));
   
    }
    
    for (int iface = 0; iface < n_internal_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterface(iface,elmat, ef);
        size_t i_begin = 2*2*(iface);
        m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
        m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
        m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
        m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));
        
        size_t i_rhs_begin = 2*(iface) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));

    }
    //FracInterfaces
    
    for (int iface = 0; iface < n_internal_faces1; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id1].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        
//        indexes[0]=righteq;
//        indexes[1]=lefteq;
        
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterface(iface,elmat, ef, internal_faces_id1);
        size_t i_begin = 2*2*(iface) + 4*n_internal_faces;
        m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
        m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
        m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
        m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));
        
        
        
        size_t i_rhs_begin = 2*(iface) + n_cells + 2*(n_internal_faces);
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
        
//        std::cout<<"leftEq= "<<lefteq<<std::endl;
//        std::cout<<"rightEq= "<<righteq<<std::endl;
//        elmat.Print(std::cout);
//        ef.Print(std::cout);
    }
    
    for (int iface = 0; iface < n_internal_faces2; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id2].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        
//        indexes[0]=righteq;
//        indexes[1]=lefteq;
        
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterface(iface,elmat, ef, internal_faces_id2);
        size_t i_begin = 2*2*(iface) + 4*n_internal_faces +4*n_internal_faces1 ;
        m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
        m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
        m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
        m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));

        
        size_t i_rhs_begin = 2*(iface) + n_cells + 2*n_internal_faces + 2*n_internal_faces1;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
        
//        std::cout<<"leftEq= "<<lefteq<<std::endl;
//        std::cout<<"rightEq= "<<righteq<<std::endl;
//        elmat.Print(std::cout);
//        ef.Print(std::cout);
        
    }
    for (int iface = 0; iface < n_internal_faces3; iface++) {

        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id3].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];

        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        
//        indexes[0]=righteq;
//        indexes[1]=lefteq;

        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterface(iface,elmat, ef, internal_faces_id3);
        size_t i_begin = 2*2*(iface) + 4*(n_internal_faces +n_internal_faces1 + n_internal_faces2) ;
        m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
        m_trans_triplets[i_begin+1] = (Triplet2<REAL>(indexes[0],indexes[1], elmat(0,1)));
        m_trans_triplets[i_begin+2] = (Triplet2<REAL>(indexes[1],indexes[0], elmat(1,0)));
        m_trans_triplets[i_begin+3] = (Triplet2<REAL>(indexes[1],indexes[1], elmat(1,1)));

        size_t i_rhs_begin = 2*(iface) + n_cells + 2*(n_internal_faces + n_internal_faces1 + n_internal_faces2);
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
        
//        std::cout<<"leftEq= "<<lefteq<<std::endl;
//        std::cout<<"rightEq= "<<righteq<<std::endl;
//        elmat.Print(std::cout);
//        ef.Print(std::cout);

    }

    for (int iface = 0; iface < n_inlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
      
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport->ContributeBCInletInterface(iface,ef);
        size_t i_rhs_begin = (iface) + n_cells + 2*(n_internal_faces + n_internal_faces1 + n_internal_faces2 + n_internal_faces3);
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
    
        
    }
    
    for (int iface = 0; iface < n_outlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport->fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        
        
    
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1, 1);
        fAlgebraicTransport->ContributeBCOutletInterface(iface,elmat, ef);
        size_t i_begin = iface +   + 4*(n_internal_faces + n_internal_faces1 + n_internal_faces2 + n_internal_faces3) ;
        m_trans_triplets[i_begin] = (Triplet2<REAL>(indexes[0],indexes[0], elmat(0,0)));
        
        size_t i_rhs_begin = (iface) + n_cells + 2*(n_internal_faces+n_internal_faces1+n_internal_faces2+n_internal_faces3) + n_inlet_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));

    }
    
    
    
#endif
    
    m_rhs.setFromTriplets( m_rhs_triplets.begin(), m_rhs_triplets.end() );
    m_rhs_triplets.clear();
    
    m_transmissibility.setFromTriplets( m_trans_triplets.begin(), m_trans_triplets.end() );
    m_trans_triplets.clear();
    
    
}
void TPZAnalysisAuxEigen::AssembleResidual(){
    int n_cells = fAlgebraicTransport->fCellsData.fVolume.size();
    m_rhs =  Eigen::SparseMatrix<REAL>( n_cells, 1 );
    
    //
    int internal_faces_id = 100;
    int internal_faces_id1 = 101;
    int internal_faces_id2 = 102;
    int internal_faces_id3 = 103;
    int inlet_faces_id = -2;
    int outlet_faces_id = -4;
    
    int n_internal_faces = fAlgebraicTransport->fInterfaceData[internal_faces_id].fFluxSign.size();
    int n_internal_faces1 = fAlgebraicTransport->fInterfaceData[internal_faces_id1].fFluxSign.size();
    int n_internal_faces2 = fAlgebraicTransport->fInterfaceData[internal_faces_id2].fFluxSign.size();
    int n_internal_faces3 = fAlgebraicTransport->fInterfaceData[internal_faces_id3].fFluxSign.size();
    int n_inlet_faces = fAlgebraicTransport->fInterfaceData[inlet_faces_id].fFluxSign.size();
    int n_outlet_faces = fAlgebraicTransport->fInterfaceData[outlet_faces_id].fFluxSign.size();
    
   
    size_t n_nzeros_res = n_cells + (n_internal_faces * 2) + (n_internal_faces1 * 2) + (n_internal_faces2 * 2) + (n_internal_faces3 * 2) +n_inlet_faces + n_outlet_faces;
    m_rhs_triplets.resize(n_nzeros_res);
    m_rhs.setZero();
    //
    
#ifdef USING_TBB2
    
    tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1),
                      [this] (size_t & ivol){
                          int eqindex = fAlgebraicTransport->fCellsData.fEqNumber[ivol];
                          TPZFMatrix<double> ef;
                          ef.Resize(1,1);
                          fAlgebraicTransport->ContributeResidual(ivol, ef);
                          m_rhs_triplets[ivol] = Triplet2<REAL>(eqindex,0, ef(0,0));
                      }
                      );
    
    
    tbb::parallel_for(size_t(0), size_t(n_internal_faces), size_t(1),
                      [this,&internal_faces_id,&n_cells] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int right = lrindex.second;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
                          
                          TPZVec<int64_t> indexes(2);
                          indexes[0]=lefteq;
                          indexes[1]=righteq;
                          TPZFMatrix<double> ef;
                          ef.Resize(2, 1);
                          fAlgebraicTransport->ContributeInterfaceResidual(iface, ef);
                          
                          size_t i_rhs_begin = 2*(iface) + n_cells;
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                          m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
                      }
                      );
    
    tbb::parallel_for(size_t(0), size_t(n_inlet_faces), size_t(1),
                      [this,&inlet_faces_id,&n_internal_faces,&n_cells] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          
                          TPZVec<int64_t> indexes(1);
                          indexes[0]=lefteq;
                          TPZFMatrix<double> ef;
                          ef.Resize(1, 1);
                          fAlgebraicTransport->ContributeBCInletInterface(iface,ef);
                          size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces;
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                      }
                      );
    
    tbb::parallel_for(size_t(0), size_t(n_outlet_faces), size_t(1),
                      [this,&outlet_faces_id,&n_internal_faces,&n_inlet_faces,&n_cells] (size_t & iface){
                          std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport->fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
                          int left = lrindex.first;
                          int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
                          
                          TPZVec<int64_t> indexes(1);
                          indexes[0]=lefteq;
                          TPZFMatrix<double>  ef;
                          ef.Resize(1, 1);
                          fAlgebraicTransport->ContributeBCOutletInterfaceResidual(iface, ef);
                          
                          size_t i_rhs_begin = (iface) + n_cells + 2*n_internal_faces + n_inlet_faces;
                          m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
                      }
                      );
    
#else
    
    for (int ivol = 0; ivol < n_cells; ivol++) {
        int eqindex = fAlgebraicTransport->fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> ef;
        ef.Resize(1,1);
        fAlgebraicTransport->ContributeResidual(ivol, ef);
        m_rhs_triplets[ivol] = Triplet2<REAL>(eqindex,0, ef(0,0));
//        std::cout<<"pos: "<<eqindex<<" val: "<<ef(0,0)<<std::endl;
    }
    
    for (int iface = 0; iface < n_internal_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> ef;
        
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterfaceResidual(iface, ef, internal_faces_id);
        
        size_t i_rhs_begin = 2*(iface) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
    }
    
    //
    for (int iface = 0; iface < n_internal_faces1; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id1].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> ef;
        
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterfaceResidual(iface, ef,internal_faces_id1);
        
        size_t i_rhs_begin = 2*(iface) + 2*(n_internal_faces) + n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
    
        
    }
    
    for (int iface = 0; iface < n_internal_faces2; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id2].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> ef;
        
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterfaceResidual(iface, ef,internal_faces_id2);
        
        size_t i_rhs_begin = 2*(iface) + 2*(n_internal_faces+ n_internal_faces1) +  n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
        
    }
    //
    
    for (int iface = 0; iface < n_internal_faces3; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[internal_faces_id3].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport->fCellsData.fEqNumber[right];
        
        TPZVec<int64_t> indexes(2);
        indexes[0]=lefteq;
        indexes[1]=righteq;
        TPZFMatrix<double> ef;
        
        ef.Resize(2, 1);
        fAlgebraicTransport->ContributeInterfaceResidual(iface, ef,internal_faces_id3);
        
        size_t i_rhs_begin = 2*(iface) + 2*(n_internal_faces+ n_internal_faces1 + n_internal_faces2) +  n_cells;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        m_rhs_triplets[i_rhs_begin+1] = Triplet2<REAL>(indexes[1],0, ef(1,0));
        
        //        std::cout<<"pos: "<<indexes[0]<<" val: "<<ef(0,0)<<std::endl;
        //        std::cout<<"pos: "<<indexes[1]<<" val: "<<ef(1,0)<<std::endl;
    }
    //
    
    for (int iface = 0; iface < n_inlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport->fInterfaceData[inlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport->ContributeBCInletInterface(iface,ef);
        size_t i_rhs_begin = (iface) + n_cells + 2*(n_internal_faces + n_internal_faces1+n_internal_faces2+n_internal_faces3);
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));
        
//        std::cout<<"pos: "<<indexes[0]<<" val: "<<ef(0,0)<<std::endl;
    }
    
    for (int iface = 0; iface < n_outlet_faces; iface++) {
        
        std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport->fInterfaceData[outlet_faces_id].fLeftRightVolIndex[iface];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport->fCellsData.fEqNumber[left];
        
        TPZVec<int64_t> indexes(1);
        indexes[0]=lefteq;
        TPZFMatrix<double> ef;
        ef.Resize(1, 1);
        fAlgebraicTransport->ContributeBCOutletInterfaceResidual(iface, ef);
        
        size_t i_rhs_begin = (iface) + n_cells + 2*(n_internal_faces+n_internal_faces1+n_internal_faces2+n_internal_faces3) + n_inlet_faces;
        m_rhs_triplets[i_rhs_begin] = Triplet2<REAL>(indexes[0],0, ef(0,0));

    }
    
    
    
#endif
    
    m_rhs.setFromTriplets( m_rhs_triplets.begin(), m_rhs_triplets.end() );
    m_rhs_triplets.clear();
}

void TPZAnalysisAuxEigen::AnalyzePattern(){
    AssembleMass();
    Assemble();
    m_transmissibility += m_mass;
//    std::cout<<m_transmissibility.toDense()<<std::endl;
    m_analysis.analyzePattern(m_transmissibility);
    
}

void TPZAnalysisAuxEigen::Solve(){
    TPZFMatrix<REAL> gmatrixmass(m_transmissibility.rows(),m_transmissibility.rows(),0.0);
    TPZFMatrix<REAL> gmatrix(m_transmissibility.rows(),m_transmissibility.rows(),0.0);
    TPZFMatrix<REAL> rhs(m_transmissibility.rows(),1,0.0);
    
//        for(int i=0; i< m_mass.rows(); i++){
//            for(int j=0; j< m_mass.rows(); j++){
//                gmatrixmass(i,j) =m_transmissibility.coeffRef(i, j);
//            }
//        }
////        gmatrixmass.Print("Ek= ",std::cout, EMathematicaInput);
    
    m_transmissibility += m_mass;
    
    
    
     m_rhs *= -1.0;

//        for(int j=0; j< m_mass.rows(); j++){
//            rhs(j,0) =m_rhs.coeffRef(j, 0);
//        }
//
//    rhs.Print("rhs= ",std::cout, EMathematicaInput);
    
//    std::cout<<"Matrix: "<<std::endl;
//    std::cout<<m_transmissibility.toDense()<<std::endl;
   
 
//    gmatrixmass.Print("Ek= ",std::cout, EMathematicaInput);
//    rhs.Print("RHS= ",std::cout, EMathematicaInput);
    
//    std::cout<<"RHS: "<<std::endl;
//    std::cout<<m_rhs.toDense()<<std::endl;
    

    m_analysis.factorize(m_transmissibility);
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> ds = m_analysis.solve(m_rhs);
    m_ds=ds;
    
   
//    TPZFMatrix<REAL> vectorm(m_transmissibility.rows(),1,0.0);
//    for(int i=0; i< m_transmissibility.rows(); i++){
//        vectorm(i,0) = ds.coeffRef(i, 0);
//    }
//    vectorm.Print("Sol= ",std::cout, EMathematicaInput);
//    std::cout<<"MatrixTransportMortar"<<std::endl;
//    std::ofstream filemortar("matrixmortar.txt");
//    gmatrixmass.Print("MassMotar=", filemortar, EMathematicaInput);
//    gmatrix.Print("KgMotar=", filemortar, EMathematicaInput);
//    vectorm.Print("vector=", filemortar, EMathematicaInput);
    
    
//    std::cout<<"Solution"<<std::endl;
//    std::cout<<ds<<std::endl;
}
Eigen::SparseMatrix<REAL> TPZAnalysisAuxEigen::Rhs(){
    return m_rhs;
}
int TPZAnalysisAuxEigen::NRhsRows(){
    return m_rhs.rows();
}
