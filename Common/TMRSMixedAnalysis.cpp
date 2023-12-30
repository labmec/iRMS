//
//  TMRSMixedAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.

#include "TMRSMixedAnalysis.h"
#include "TPZSpStructMatrix_Eigen.h"
#include "TPZSpMatrixEigen.h"
#include <pzshapequad.h>
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "TPZCompElHDivCollapsed.h"

// Uses the new vtk function developed by Fran
#define USENEWVTK

#ifdef USENEWVTK
#include "TPZVTKGenerator.h"
#endif

using namespace std;

TMRSMixedAnalysis::TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::~TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::TMRSMixedAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,
                                     const RenumType& renumtype) : TPZLinearAnalysis(cmesh_mult, renumtype){
    fsoltransfer.BuildTransferData(cmesh_mult);
}

void TMRSMixedAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    
}

TMRSDataTransfer * TMRSMixedAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSMixedAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSMixedAnalysis::Configure(int n_threads, bool UsePardiso_Q,bool UsePZ){
    
    if (UsePardiso_Q) {
        if(UsePZ){
            TPZSSpStructMatrix<STATE> matrix(Mesh());
            matrix.SetNumThreads(n_threads);
            SetStructuralMatrix(matrix);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            SetSolver(step);
        }
        else{
            TPZSymetricSpStructMatrixEigen matrix(Mesh());
            matrix.SetNumThreads(n_threads);
            SetStructuralMatrix(matrix);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            SetSolver(step);
        }


    }else{
        TPZSkylineStructMatrix<STATE> matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
    
   
    //    Assemble();
}

void TMRSMixedAnalysis::RunTimeStep(){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh)
        DebugStop();
    
    int n = m_sim_data->mTNumerics.m_max_iter_mixed;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_mixed;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_mixed;
    
    TPZFMatrix<STATE> dx,x(Solution());
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
       
        NewtonIteration();
        dx = Solution();
        corr_norm = Norm(dx);
        res_norm = Norm(Rhs());
        x +=dx;
//        cmesh->UpdatePreviousState(-1.);
        fsoltransfer.TransferFromMultiphysics();

        Assemble();
      
        res_norm = Norm(Rhs());
        REAL normsol = Norm(Solution());
        

#ifdef PZDEBUG
        {
            if(std::isnan(corr_norm) || std::isnan(res_norm))
                DebugStop();
        }
#endif

        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_Q) {
            std::cout << "\n\n\t================================================" << std::endl;
            std::cout << "Mixed operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            break;
        }
    }
}


void TMRSMixedAnalysis::NewtonIteration(){
    
    if(mIsFirstAssembleQ == true)
    {
        fStructMatrix->SetNumThreads(m_sim_data->mTNumerics.m_nThreadsMixedProblem);
        int64_t nel = fCompMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCompMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
//            int numthreads = 0;
//                sub->SetAnalysisSparse(0); sub->Analysis()->StructMatrix()->SetNumThreads(m_sim_data->mTNumerics.m_nThreadsMixedProblem);
//                TPZSymetricSpStructMatrixEigen matrix(sub);
//                //matrix.SetNumThreads(n_threads);
//                sub->Analysis()->SetStructuralMatrix(matrix);
            }
        }
        mIsFirstAssembleQ=false;
    }

    Assemble();
    Solve();

}

void TMRSMixedAnalysis::PostProcessTimeStep(int dimToPost){
    
    const int dim = Mesh()->Dimension();
    auto start_time_pp = std::chrono::steady_clock::now();
    cout << "\n--------------------- Post process dim = " << dimToPost << " ---------------------\n" << endl;
        
    TPZStack<std::string,10> scalnames, vecnames;
    
    scalnames = m_sim_data->mTPostProcess.m_scalnamesDarcy;
    vecnames = m_sim_data->mTPostProcess.m_vecnamesDarcy;
    
    int div = 0;
    if (dimToPost < 0){
        dimToPost = Mesh()->Reference()->Dimension();
    }
//    dim = 2;
//    std::set<int> mat_id_2D;
//    mat_id_2D.insert(10);
//    std::string file_frac("fractureFlux_s.vtk");
//    DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
//    PostProcess(div,2);
    std::string file = m_sim_data->mTPostProcess.m_file_name_mixed;

    constexpr int vtkRes{0}; //resolucao do vtk
    if (dimToPost == dim-1){
        std::set<int> matids;
        map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
        for (it = m_sim_data->mTFracProperties.m_fracprops.begin(); it != m_sim_data->mTFracProperties.m_fracprops.end(); it++)
        {
            int matfracid = it->first;
            matids.insert(matfracid);
        }

#ifdef USENEWVTK
        const std::string plotfile = file.substr(0, file.find(".")) + "_frac";
        for (auto nm : vecnames) {
            scalnames.Push(nm);
        }
        auto vtk = TPZVTKGenerator(fCompMesh, matids, scalnames, plotfile, vtkRes);
        vtk.SetNThreads(8);
        vtk.Do();
#else
        file = file.substr(0, file.find(".")) + "_frac.vtk";
        DefineGraphMesh(dimToPost, matids, scalnames, vecnames, file);
        PostProcess(div,dimToPost);
#endif
    }
    else{
#ifdef USENEWVTK
        const std::string plotfile = file.substr(0, file.find(".")); //sem o .vtk no final
        for (auto nm : vecnames) {
            scalnames.Push(nm);
        }
        
        auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dimToPost);
        vtk.SetNThreads(8);
        vtk.Do();
#else
        DefineGraphMesh(dimToPost,scalnames,vecnames,file);
        PostProcess(div,dimToPost);
#endif
    }
    
    
    
    auto total_time_pp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_pp).count()/1000.;
    cout << "Total time post process = " << total_time_pp << " seconds" << endl;

}

void TMRSMixedAnalysis::Assemble(){

    auto start_time_ass = std::chrono::steady_clock::now();

    cout << "\n---------------------- Assemble Flux Problem ----------------------" << endl;
    cout << "Number of equations: " << fCompMesh->NEquations() << endl;
    cout << "Number of elements: " << fCompMesh->NElements() << endl;
    TPZLinearAnalysis::Assemble();


    auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count()/1000.;
    cout << "\nTotal time assemble = " << total_time_ass << " seconds" << endl;
}

void TMRSMixedAnalysis::Solve(){

    auto start_time_solve = std::chrono::steady_clock::now();
    
    cout << "\n---------------------- Solve Flux Problem ----------------------" << endl;
    TPZLinearAnalysis::Solve();
    
    auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
    cout << "Total time solve = " << total_time_solve << " seconds" << endl;
}

void TMRSMixedAnalysis::VerifyElementFluxes(){
    const REAL tol = 1.e-10;
	TPZMultiphysicsCompMesh *mixedmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh()) ;
	TPZCompMesh *cmesh = mixedmesh->MeshVector()[0];
//	std::ofstream file("fuxmesh.txt");
	const TPZFMatrix<STATE> &meshSol = cmesh->Solution();
//	mixedmesh->MeshVector()[0]->Print(file);
	int nels =mixedmesh->MeshVector()[0]->NElements();
	for (int iel =0; iel<nels-1; iel++) {
		TPZCompEl *cel = mixedmesh->MeshVector()[0]->Element(iel);
		if (!cel) continue;
//            if(cel->Dimension() != cmesh->Dimension()) continue;
		TPZCompElHDiv<pzshape::TPZShapeCube> *hdivel = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(cel);
		TPZCompElHDiv<pzshape::TPZShapeTetra> *hdiveltet = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(cel);
		TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivelq = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(cel);
		TPZCompElHDiv<pzshape::TPZShapeTriang> *hdivelqTri = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTriang> *>(cel);
		TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollaps = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivelq);
		TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *hdivCollapsTri = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *>(hdivelqTri);
		if (!hdivel && !hdiveltet && !hdivCollaps && !hdivCollapsTri) {
			continue;
		}
		TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
		
		int ncon = cel->NConnects();
		int nCorners = cel->Reference()->NCornerNodes();
		int nsides1 = cel->Reference()->NSides(1);
		int nsides = cel->Reference()->NSides();
		REAL sumel=0.0;
		for (int icon=0; icon<ncon-1; icon++) {
			TPZConnect &con = cel->Connect(icon);
			int sideOrient =0;
			const bool isColaps = hdivCollaps || hdivCollapsTri;
			const int iconNum = hdivCollaps ? 5 : 4;
			const int sideForEl = hdivCollaps ? 8 : 6;
			if(!isColaps){
				 sideOrient = intel->GetSideOrient(nCorners+nsides1+icon);
			}
			else{
				if (isColaps && icon < iconNum) {
					sideOrient = intel->GetSideOrient(nCorners+icon);
				}
				 
			}
			if (isColaps && icon == iconNum-1) {
				ncon++;
				continue;
			}
			if (isColaps && icon == iconNum) {
				sideOrient = intel->GetSideOrient(sideForEl);
			}
			if (isColaps && icon == iconNum+1) {
				sideOrient = intel->GetSideOrient(sideForEl+1);
			}
			
		   
			int sequence =con.fSequenceNumber;
			int64_t pos = cmesh->Block().Position(sequence);
			if(sequence==-1) continue;
			int blocksize = cmesh->Block().Size(sequence);
			for(int ieq=0; ieq< blocksize; ieq++)
			{
                sumel += sideOrient*meshSol.GetVal(cmesh->Block().Index(sequence,ieq),0);
			}
		}
		if(std::abs(sumel)> tol ){
            std::cout << "\n\nERROR! Conservation of element index " << cel->Reference()->Index() << " is " << sumel << std::endl;
			DebugStop();
		}
	}
    std::cout << "\n\n===> Nice! All flux elements satisfy conservation up to tolerance " << tol << std::endl;
}
