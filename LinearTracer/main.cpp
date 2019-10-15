
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "pzgengrid.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pznoderep.h"
#include "pzshapetetra.h"
#include "RSimulatorConfiguration.h"
#include "pzshapepiram.h"

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"

TPZAnalysis * CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,  bool must_opt_band_width_Q, int n_threads, bool UsePardiso_Q);

TPZAnalysis * CreateTransportAnalysis(TPZMultiphysicsCompMesh * cmesh, bool must_opt_band_width_Q, int n_threads, bool UsePardiso_Q);

TPZFMatrix<STATE> TimeForward(TMRSDataTransfer & sim_data, TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag);

TMRSDataTransfer Setting2D();

TMRSDataTransfer Setting3D();

int main(){
 
    bool is_3D_Q = false;
    bool is_2D_Coarse_Q = false;
    TMRSDataTransfer sim_data;
    std::string geometry_file, name;
    if (is_3D_Q) {
        geometry_file = "gmsh/reservoir_3d.msh";
        name = "reservoir_3d";
        sim_data = Setting3D();
    }else{
        if (is_2D_Coarse_Q) {
            geometry_file = "gmsh/reservoir_2d_coarse.msh";
        }
        else{
            geometry_file = "gmsh/reservoir_2d_fine.msh";
        }
        name = "reservoir_2d";
        sim_data = Setting2D();
    }
    
    
    TMRSApproxSpaceGenerator aspace;
    aspace.LoadGeometry(geometry_file);
    aspace.PrintGeometry(name);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 24;
    bool UsePardiso_Q = true;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();

    aspace.BuildTransportMultiPhysicsCompMesh();
    TPZMultiphysicsCompMesh * transport_operator = aspace.GetTransportOperator();
    
    
    // Linking the memory between the operators
    {
        TMRSApproxSpaceGenerator::AdjustMemory(mixed_operator, transport_operator);
        for (auto item : sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]) {
            int material_id = item.second;
            TMRSApproxSpaceGenerator::UnifyMaterialMemory(material_id, mixed_operator, transport_operator);
            TMRSApproxSpaceGenerator::FillMaterialMemory(material_id, mixed_operator, transport_operator);
        }
    }
    
    // Solving global darcy problem
    {
        TPZAnalysis * analysis =  CreateAnalysis(mixed_operator,  must_opt_band_width_Q, n_threads, UsePardiso_Q);
        
        analysis->Assemble();
        analysis->Solve();
        mixed_operator->LoadSolutionFromMultiPhysics();
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        int dim = aspace.GetGeometry()->Dimension();
        std::string fileresult("flux_and_p.vtk");
        analysis->DefineGraphMesh(dim,scalnames,vecnames,fileresult);
        analysis->PostProcess(div,dim);
    }
    
    TPZAnalysis *tracer_an = CreateTransportAnalysis(transport_operator, must_opt_band_width_Q, n_threads, UsePardiso_Q);
    
    // Solving transport problem
    int n_steps = 50;
    REAL dt     = 1000.0;
    TPZFMatrix<STATE> M_diag;
    TPZFMatrix<STATE> saturations = TimeForward(sim_data, tracer_an, n_steps, dt, M_diag);
    
    return 0;
    
}

TMRSDataTransfer Setting2D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_rock"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_wbregion_p"] = 2;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["d_wbregion_i"] = 3;
    
    //zero flux boundaries
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 30.0;
    REAL pressure_out = 20.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(6);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(4,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(5,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(6,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(7,N_Type,zero_flux);
    
    //well pressure
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[4] = std::make_tuple(10,D_Type,pressure_out);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[5] = std::make_tuple(11,D_Type,pressure_in);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(6);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(4,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(5,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(6,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(7,bc_inlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[4] = std::make_tuple(10,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[5] = std::make_tuple(11,bc_inlet,sat_in);
    
    //Relative permermeabilities
    
    sim_data.mTPetroPhysics.mLayer_Krw_RelPerModel.Resize(1);
    TRSLinearInterpolator krw, krow ;
    std::string name_krw("PetroPhysics/krw.txt");
    std::string name_krow("PetroPhysics/krow.txt");
    krw.ReadData(name_krw);
    krow.ReadData(name_krow);
    sim_data.mTPetroPhysics.mLayer_Krow_RelPerModel[0] ;
    
    
    return sim_data;
}

TMRSDataTransfer Setting3D(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["d_rock"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["d_wbregion_p"] = 2;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["d_wbregion_i"] = 3;
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(10);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(4,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(5,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(6,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(7,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[4] = std::make_tuple(8,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[5] = std::make_tuple(9,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[6] = std::make_tuple(10,1,0.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[7] = std::make_tuple(11,1,0.0);
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[8] = std::make_tuple(12,0,20.0);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[9] = std::make_tuple(13,0,30.0);
    
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(10);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(4,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(5,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(6,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(7,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[4] = std::make_tuple(8,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[5] = std::make_tuple(9,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[6] = std::make_tuple(10,0,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[7] = std::make_tuple(11,0,0.0);
    
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[8] = std::make_tuple(12,1,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[9] = std::make_tuple(13,0,1.0);
    
    return sim_data;
}


TPZAnalysis * CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_mult,  bool must_opt_band_width_Q, int n_threads, bool UsePardiso_Q){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh_mult, must_opt_band_width_Q);
    
    if (UsePardiso_Q) {
        
        TPZSymetricSpStructMatrix matrix(cmesh_mult);
        matrix.SetNumThreads(n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);

    }else{
        
        TPZSkylineStructMatrix matrix(cmesh_mult);
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);

    }
    
    return analysis;
    
}

TPZAnalysis * CreateTransportAnalysis(TPZMultiphysicsCompMesh * cmesh, bool must_opt_band_width_Q, int n_threads, bool UsePardiso_Q){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, must_opt_band_width_Q);
    
    if (UsePardiso_Q) {
        
        TPZSpStructMatrix matrix(cmesh);
        matrix.SetNumThreads(n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        
        return analysis;
    }else{
        
        TPZSkylineNSymStructMatrix matrix(cmesh);
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
        
    }
    
}

#define NewTimeForward_Q

TPZFMatrix<STATE> TimeForward(TMRSDataTransfer & sim_data, TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag){
    
    TPZMultiphysicsCompMesh * cmesh_transport = dynamic_cast<TPZMultiphysicsCompMesh *>(tracer_analysis->Mesh());
    
    if (!cmesh_transport) {
        DebugStop();
    }
    
    std::set<int> volumetric_mat_ids = {1,2,3,100};
    TPZManVector<TPZCompMesh *,3> meshtrvec = cmesh_transport->MeshVector();
    
    int64_t n_eq = tracer_analysis->Mesh()->NEquations();
    TPZFMatrix<STATE> saturations(n_eq,n_steps);
    
#ifdef NewTimeForward_Q
    
    /// Time evolution
    std::string file_reservoir("transport.vtk");
    
    {
        int div = 0;
        int n_iterations = 1;
        bool stop_criterion_Q = false;
        REAL res_norm = 1.0;
        REAL res_tol = 1.0e-8;
        TPZFMatrix<REAL> ds,s;
        s = tracer_analysis->Solution();
        for (int it = 0; it < n_steps; it++) {
            
            // Newton like process
            for (int k = 0; k < n_iterations; k++) {
                
                /// Newton correction
                tracer_analysis->Assemble();
//                tracer_analysis->Rhs().Print("r = ",std::cout,EMathematicaInput);
                tracer_analysis->Rhs() *= -1.0;
                tracer_analysis->Solve(); /// (LU decomposition)
                ds = tracer_analysis->Solution();
//                ds.Print("ds = ",std::cout,EMathematicaInput);
                s += ds;
                tracer_analysis->LoadSolution(s);
                cmesh_transport->LoadSolutionFromMultiPhysics();
                tracer_analysis->AssembleResidual();
                
                res_norm = Norm(tracer_analysis->Rhs());
                stop_criterion_Q = res_norm <= res_tol;
                if (stop_criterion_Q) {
                    std::cout << "Residual norm = " << res_norm << std::endl;
                    // Updating memory
                    TMRSApproxSpaceGenerator::SetUpdateMemory(2, sim_data, cmesh_transport, true);
                    tracer_analysis->AssembleResidual();
                    TMRSApproxSpaceGenerator::SetUpdateMemory(2, sim_data, cmesh_transport, false);
                }else{
                    DebugStop();
                }
                

                
                /// postprocess ...
                TPZStack<std::string,10> scalnames, vecnames;
                scalnames.Push("Sw");
                scalnames.Push("So");
                
                std::map<int,int> volumetric_ids;
                volumetric_ids.insert(std::make_pair(1, 2));
                
                std::map<int,int> fracture_ids;
                fracture_ids.insert(std::make_pair(6, 2));
                
                std::map<int,int> fracture_intersections_ids;
                fracture_intersections_ids.insert(std::make_pair(7, 1));
                
                for (auto data: volumetric_ids) {
                    TPZMaterial * mat = cmesh_transport->FindMaterial(data.first);
                    TMRSMultiphaseFlow<TMRSMemory> * volume = dynamic_cast<TMRSMultiphaseFlow<TMRSMemory> * >(mat);
                    if (!volume) {
                        DebugStop();
                    }
                    volume->SetDimension(data.second);
                }
                
                int dimension = tracer_analysis->Mesh()->Reference()->Dimension();
                std::set<int> mat_id_3D = volumetric_mat_ids;
                
                std::string file_reservoir("saturation.vtk");
                tracer_analysis->DefineGraphMesh(dimension,mat_id_3D,scalnames,vecnames,file_reservoir);
                tracer_analysis->PostProcess(div,dimension);
                
            }
            
        }
        
        
    }
    
    return saturations;
    
    
#else
    
    /// Compute mass matrix M.
    TPZAutoPointer<TPZMatrix<STATE> > M;
    TPZFMatrix<REAL> F_inlet;
    {
        
        bool mass_matrix_Q = true;
        
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing Mass Matrix." << std::endl;
        tracer_analysis->Assemble();
        std::cout << "Mass Matrix is computed." << std::endl;
        M = tracer_analysis->Solver().Matrix()->Clone();
    }
    
    //    M->Print("M = ",std::cout,EMathematicaInput);
    
    int n_rows = M->Rows();
    M_diag.Resize(n_rows,M->Cols());
    
    for (int64_t i = 0; i < n_rows; i++) {
        M_diag(i,0) = M->Get(i, i);
    }
    
    
    
    
    int64_t n_eq = tracer_analysis->Mesh()->NEquations();
    TPZFMatrix<STATE> saturations(n_eq,n_steps);
    
    {
        bool mass_matrix_Q = false;
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                continue;
            }
            volume->SetTimeStep(dt);
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing transport operator K = M + T, and F_inlet " << std::endl;
        tracer_analysis->Assemble();
        //        tracer_analysis->Rhs().Print("f = ",std::cout,EMathematicaInput);
        
        F_inlet = tracer_analysis->Rhs();
        
    }
    
    /// Time evolution
    std::string file_reservoir("transport.vtk");
    
    {
        int div = 0;
        TPZFMatrix<REAL> s_n(n_eq,1,0.0);
        TPZFMatrix<REAL> last_state_mass(n_eq,1,0.0);
        TPZFMatrix<REAL> s_np1;
        
        for (int it = 0; it < n_steps; it++) {
            
            for (int64_t i = 0; i < n_eq; i++) {
                last_state_mass(i,0) = M_diag(i,0)*s_n(i,0);
            }
            
            tracer_analysis->Rhs() = F_inlet - last_state_mass;
            tracer_analysis->Rhs() *= -1.0;
            
            tracer_analysis->Solve(); /// (LU decomposition)
            
            s_np1 = tracer_analysis->Solution();
            tracer_analysis->LoadSolution(s_np1);
            
            /// postprocess ...
            TPZStack<std::string,10> scalnames, vecnames;
            scalnames.Push("Sw");
            scalnames.Push("So");
            
            std::map<int,int> volumetric_ids;
            volumetric_ids.insert(std::make_pair(1, 2));
            
            std::map<int,int> fracture_ids;
            fracture_ids.insert(std::make_pair(6, 2));
            
            std::map<int,int> fracture_intersections_ids;
            fracture_intersections_ids.insert(std::make_pair(7, 1));
            
            for (auto data: volumetric_ids) {
                TPZMaterial * mat = cmesh_transport->FindMaterial(data.first);
                TMRSMultiphaseFlow<TMRSMemory> * volume = dynamic_cast<TMRSMultiphaseFlow<TMRSMemory> * >(mat);
                if (!volume) {
                    DebugStop();
                }
                volume->SetDimension(data.second);
            }
            
            int div = 0;
            int dimension = tracer_analysis->Mesh()->Reference()->Dimension();
            std::set<int> mat_id_3D = volumetric_mat_ids;
            
            std::string file_reservoir("saturation.vtk");
            tracer_analysis->DefineGraphMesh(dimension,mat_id_3D,scalnames,vecnames,file_reservoir);
            tracer_analysis->PostProcess(div,dimension);
            
            //                std::set<int> mat_id_2D;
            //                mat_id_2D.insert(6);
            //                std::string file_frac("fracture_s.vtk");
            //                tracer_analysis->DefineGraphMesh(1,mat_id_2D,scalnames,vecnames,file_frac);
            //                tracer_analysis->PostProcess(div,1);
            
            
            
            // configuring next time step
            s_n = s_np1;
            for (int64_t i = 0; i < n_eq; i++) {
                saturations(i,it) = cmesh_transport->Solution()(i,0);
            }
            
        }
        
        
    }
    
    return saturations;
    
#endif


}
