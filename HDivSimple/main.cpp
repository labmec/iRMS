
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void SimpleTest2DHDiv();

void LearningReadFracMesh();

TMRSDataTransfer SettingPaper3D();


TPZGeoMesh *ReadFractureMesh();

TMRSDataTransfer SettingSimple2DHdiv();
void  ForcingFunction (const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

int main(){
    InitializePZLOG();
    LearningReadFracMesh();
//    SimpleTest2DHDiv();
}

void SimpleTest2DHDiv(){
    
    //sim_data contains the simulation configuration, including contour conditions and muneric controls
    TMRSDataTransfer sim_data  = SettingSimple2DHdiv();
    TMRSApproxSpaceGenerator aspace;
    
    //The problem geometry is a rectangle with LxY dimentions and elements number in x y and nxXny
    int nx=2;
    int ny=2;
    REAL L =10.0;
    REAL Y =10.0;
    aspace.CreateUniformMesh(nx, L, ny, Y);
    

    std::string name = "2D_geo";
//    aspace.PrintGeometry(name);
    aspace.ApplyUniformRefinement(0);
    std::string name_ref = "2D_ref_geo";
    aspace.PrintGeometry(name_ref);
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    
    //This parameter should be always "true"
    bool UsePardiso_Q = true;
    
    //Multiphysic mesh creation of the mixed problem
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    {
        std::ofstream mout("mphysics.txt");
        mixed_operator->Print(mout);
    }
    TMRSMixedAnalysis *mixedAnal = new TMRSMixedAnalysis(mixed_operator,must_opt_band_width_Q);
    
    //If the parameter "UsingPzSparse" is true, it uses the pz sparse matrix, otherwise it uses eigen sparse matrix
    bool UsingPzSparse = false;

    //The parallelism is just implemented for the "UsingPzSparse=True" case, with eigen for now is running in serial (the next task to do)
    mixedAnal->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    mixedAnal->SetDataTransfer(&sim_data);
    mixedAnal->Assemble();
    size_t n_dof = mixedAnal->Solver().Matrix()->Rows();

    #ifdef USING_BOOST
        boost::posix_time::ptime tsim1 = boost::posix_time::microsec_clock::local_time();
    #endif
    mixedAnal->Solve();
    #ifdef USING_BOOST
        boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
        auto deltat = tsim2-tsim1;
        std::cout << "Overal solve calling time " << deltat << std::endl;
    #endif
    std::cout << "Number of dof = " << n_dof << std::endl;
    mixed_operator->UpdatePreviousState(-1);
    mixedAnal->fsoltransfer.TransferFromMultiphysics();
    mixedAnal->PostProcessTimeStep();
  {
      std::ofstream mout("mphysics.txt");
      mixed_operator->Print(mout);
  }

}

TMRSDataTransfer SettingSimple2DHdiv(){
    
    TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[2]["RockMatrix"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressurein = 25;
    REAL pressureout = 10;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressurein);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,D_Type,pressureout);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-3,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_outlet,1.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mReferencePressure =1.013e5;
    sim_data.mTFluidProperties.mWaterViscosity = 0.8;
    sim_data.mTFluidProperties.mOilViscosity = 0.7;
    sim_data.mTPetroPhysics.mWaterViscosity = sim_data.mTFluidProperties.mWaterViscosity;
    sim_data.mTPetroPhysics.mOilViscosity = sim_data.mTFluidProperties.mOilViscosity;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilCompressibility = 1.0e-9;
    sim_data.mTFluidProperties.mWaterCompressibility = 1.0e-9;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 2;
    sim_data.mTNumerics.m_max_iter_transport = 10;
    sim_data.mTNumerics.m_max_iter_sfi = 20;
    sim_data.mTNumerics.m_sfi_tol = 0.001;
    sim_data.mTNumerics.m_res_tol_transport = 0.00001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.00001;
    sim_data.mTNumerics.m_n_steps = 1;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.03 ;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
    std::vector<REAL> grav(3,0.0);
    grav[1] = -0.0;//-9.81; //Needs to be multiplied by a constant
    sim_data.mTNumerics.m_gravity = grav;
    
    //Update to m_ISLinearKrModelQ=true if is required a linear relative permeabilities model
    sim_data.mTNumerics.m_ISLinearKrModelQ = false;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    return sim_data;
}
void  ForcingFunction (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    double fx= 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    disp[0]=fx;
}

void MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh);

TPZGeoMesh *ReadFractureMesh()
{
    std::string fileCoarse("../../FracMeshes/flem_case1_Coarse_BC.msh");
    std::string fileFine("../../FracMeshes/flem_case1_Submesh_Fractures.msh");
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
/*
    2 4 "inlet"
    2 5 "outlet"
    2 6 "noflux"
    3 3 "k33"
    3 10 "k31"
*/
    dim_name_and_physical_tagCoarse[3]["k33"] = 1;
    dim_name_and_physical_tagCoarse[3]["k31"] = 2;
    dim_name_and_physical_tagCoarse[2]["inlet"] = -2;
    dim_name_and_physical_tagCoarse[2]["outlet"] = -4;
    dim_name_and_physical_tagCoarse[2]["noflux"] = -1;

/*
 2 2 "Fractures"
 3 1 "c1"
 */
    dim_name_and_physical_tagFine[2]["Fractures"] = 10;
    for(int i=1; i<=100; i++)
    {
        std::stringstream sout;
        sout << "c" << i;
        dim_name_and_physical_tagFine[3][sout.str()] = i+10;
    }
    TPZGmshReader GeometryCoarse, GeometryFine;
    TPZGeoMesh *gmeshCoarse, *gmeshFine;
    {
        REAL l = 1.0;
        GeometryCoarse.SetCharacteristiclength(l);
        GeometryCoarse.SetFormatVersion("4.1");
        GeometryCoarse.SetDimNamePhysical(dim_name_and_physical_tagCoarse);
        gmeshCoarse = GeometryCoarse.GeometricGmshMesh(fileCoarse);
        GeometryCoarse.PrintPartitionSummary(std::cout);
    }
    {
        REAL l = 1.0;
        GeometryFine.SetCharacteristiclength(l);
        GeometryFine.SetFormatVersion("4.1");
        GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
        gmeshFine = GeometryFine.GeometricGmshMesh(fileFine);
        GeometryCoarse.PrintPartitionSummary(std::cout);
    }
    {
        std::ofstream fileCoarse("mesh3dCoarse.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshCoarse, fileCoarse);
    }
    {
        std::ofstream fileFine("mesh3dFine.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFine);
    }
    MergeMeshes(gmeshFine, gmeshCoarse);
    {
        std::ofstream fileFine("mesh3dFineMerge.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFine);
    }
    delete gmeshCoarse;
    return gmeshFine;
}

void MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh)
{
    std::map<int,int64_t> MatFinetoCoarseElIndex;
    std::map<int64_t,int64_t> NodeCoarseToNodeFine;
    std::map<int64_t,int64_t> ElCoarseToElFine;
    int temp_bc_mat = -10;
    // create boundary elements for elements without neighbour
    {
        std::map<int, int> num_created;
        int64_t nel_fine = finemesh->NElements();
        for (int64_t el = 0; el<nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            int dim = gel->Dimension();
            int nsides = gel->NSides();
            int firstside = nsides-gel->NSides(dim-1)-1;
            for (int side = firstside; side<nsides-1; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour.Element()->Dimension() != dim)
                {
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == gelside)
                {
                    TPZGeoElBC gelbc(gelside,temp_bc_mat);
                    num_created[gel->MaterialId()]++;
                }
            }
        }
        for (auto it : num_created) {
            std::cout << "For matid " << it.first << " number of elements created " << it.second << std::endl;
        }
    }
    // find the correspondence between coarse nodes and fine nodes
    {
        int64_t nnode_coarse = coarsemesh->NNodes();
        for (int64_t n = 0; n<nnode_coarse; n++) {
            TPZGeoNode &no = coarsemesh->NodeVec()[n];
            if(no.Id() == -1) continue;
            TPZManVector<REAL,3> co(3);
            no.GetCoordinates(co);
            int64_t fineindex;
            TPZGeoNode *finenode = finemesh->FindNode(co,fineindex);
            NodeCoarseToNodeFine[n] = fineindex;
        }
    }
    // identify the correspondence between the material id of the fine mesh and the element index
    // of the coarse mesh
    {
        int64_t nel_fine = finemesh->NElements();
        int dim = finemesh->Dimension();
        for (int64_t el = 0; el<nel_fine; el++) {
            auto *gel = finemesh->Element(el);
            if(gel->Dimension() != dim) continue;
            int matid = gel->MaterialId();
            if(MatFinetoCoarseElIndex.find(matid) == MatFinetoCoarseElIndex.end())
            {
                TPZManVector<REAL,3> xcenter(3);
                TPZGeoElSide gelside(gel);
                gelside.CenterX(xcenter);
                TPZManVector<REAL,3> qsi(dim,0.);
                int64_t coarse_index = 0;
                TPZGeoEl *gelcoarse = coarsemesh->FindElementCaju(xcenter, qsi, coarse_index, dim);
                MatFinetoCoarseElIndex[matid] = coarse_index;
            }
        }
    }
    // modify the material id of the boundary elements of the fine mesh
    {
        int64_t nel_fine = finemesh->NElements();
        int meshdim = finemesh->Dimension();
        std::map<int,int> created_by_mat;
        for (int64_t el = 0; el < nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            if(gel->MaterialId() == temp_bc_mat)
            {
                int dim = gel->Dimension();
                if(dim != meshdim-1) continue;
                TPZGeoElSide gelside(gel);
                TPZManVector<REAL,3> xcenter(3);
                gelside.CenterX(xcenter);
                int64_t elindex3D = 0;
                TPZManVector<REAL, 3> qsi3D(3,0.);
                coarsemesh->FindElementCaju(xcenter, qsi3D, elindex3D, meshdim);
                TPZGeoEl *coarsegel3D = coarsemesh->Element(elindex3D);
                int coarseside3D = coarsegel3D->WhichSide(qsi3D);
                if(coarsegel3D->SideDimension(coarseside3D) != dim) DebugStop();
                TPZGeoElSide BCSide(coarsegel3D,coarseside3D);
                TPZGeoElSide neighbour = BCSide.Neighbour();
                while(neighbour != BCSide)
                {
                    if(neighbour.Element()->Dimension() == dim) break;
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == BCSide) DebugStop();
                int bc_id = neighbour.Element()->MaterialId();
                created_by_mat[bc_id]++;
                gel->SetMaterialId(bc_id);
            }
        }
        for (auto it : created_by_mat) {
            std::cout << "For matid " << it.first << " number of elements created " << it.second << std::endl;
        }
    }
    // duplicate the elements of the coarse mesh within the fine mesh
    {
        int64_t nelcoarse = coarsemesh->NElements();
        int meshdim = coarsemesh->Dimension();
        for (int64_t el = 0; el<nelcoarse; el++) {
            auto gel = coarsemesh->Element(el);
            if(gel->Dimension() != meshdim) continue;
            int nnode = gel->NNodes();
            TPZManVector<int64_t, 8> nodeindices(nnode);
            for(int n=0; n<nnode; n++)
            {
                int64_t node_index_coarse = gel->NodeIndex(n);
                int64_t node_index_fine = NodeCoarseToNodeFine[node_index_coarse];
                nodeindices[n] = node_index_fine;
            }
            int matid = gel->MaterialId();
            auto eltype = gel->Type();
            int64_t fine_index;
            finemesh->CreateGeoElement(eltype, nodeindices, matid, fine_index);
            ElCoarseToElFine[el] = fine_index;
        }
    }
    finemesh->BuildConnectivity();
    // set the father index of the fine elements to the indices of the coarse elements
    // modify the material id of the boundary elements of the fine mesh
    {
        int64_t nel_fine = finemesh->NElements();
        int dim = finemesh->Dimension();
        for (int64_t el = 0; el < nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            if(gel->Dimension() != dim) continue;
            int matid = gel->MaterialId();
            int64_t coarse_index = MatFinetoCoarseElIndex[matid];
            int64_t fine_index = ElCoarseToElFine[coarse_index];
            gel->SetFatherIndex(fine_index);
        }
    }
}

TMRSDataTransfer SettingPaper3D(){
    
   TMRSDataTransfer sim_data;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["Reservoir"] = 1;
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 15.0;
    REAL pressure_out = 10.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-3,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-4,D_Type,pressure_out);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-3,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[3] = std::make_tuple(-4,bc_outlet,0.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.001;
    sim_data.mTFluidProperties.mOilViscosity = 0.001;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 500.0;
    
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 50;
    sim_data.mTNumerics.m_max_iter_sfi = 50;
    
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0000001;
    sim_data.mTNumerics.m_n_steps = 10;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 0.01;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = 2*sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    return sim_data;
}

void LearningReadFracMesh()
{
    TPZGeoMesh *gmesh = ReadFractureMesh();
    
    TMRSApproxSpaceGenerator aspace;
    TMRSDataTransfer sim_data  = SettingPaper3D();
    sim_data.mTGeometry.mSkeletonDiv =0;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
    aspace.SetGeometry(gmesh);
    
    aspace.SetDataTransfer(sim_data);

    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);


    delete gmesh;
}

