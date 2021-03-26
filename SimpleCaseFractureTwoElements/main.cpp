
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TPZRefPatternTools.h"
#include "TPZReservoirTools.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void LearningReadFracMesh();
void FracSimpleCase();

TMRSDataTransfer SettingSimpleFracCase();
TPZGeoMesh *ReadFractureMesh();
TMRSDataTransfer SettingSimple2DHdiv();

int main(){
    InitializePZLOG();
    LearningReadFracMesh();

}

void LearningReadFracMesh()
{
    /*
      the different lagrange levels for this mesh layout
     char fluxmortar = 5;
     char firstpressurelagrange = 1;
     char pressurelagrange = 3;
     char pressuremortar = 4;
     char distfluxlagrange = 2;
     char avpressurelagrange = 6;

    */
    
    TPZGeoMesh *gmesh = ReadFractureMesh();
//    TPZVec<int64_t> subdomain;
//    TPZGeoMesh *gmesh =ReadFractureMesh(subdomain);
    
    TMRSApproxSpaceGenerator aspace;
    TMRSDataTransfer sim_data  = SettingSimpleFracCase();
    sim_data.mTGeometry.mSkeletonDiv =0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
 
    aspace.SetGeometry(gmesh);
    std::string name("fractureTest.vtk");
    aspace.PrintGeometry(name);
    
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();

    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    
    //This parameter should be always "true"
    bool UsePardiso_Q = true;
    
    //
    {
      
        TMRSMixedAnalysis *mixedAnal = new TMRSMixedAnalysis(mixed_operator,must_opt_band_width_Q);
        
        
        //If the parameter "UsingPzSparse" is true, it uses the pz sparse matrix, otherwise it uses eigen sparse matrix
        bool UsingPzSparse = false;
        
        //The parallelism is just implemented for the "UsingPzSparse=True" case, with eigen for now is running in serial (the next task to do)
        mixedAnal->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
        mixedAnal->SetDataTransfer(&sim_data);
        
        mixedAnal->Assemble();
        mixedAnal->Solve();
        mixed_operator->UpdatePreviousState(-1);
        
        TPZFastCondensedElement::fSkipLoadSolution = false;
        mixed_operator->LoadSolution(mixed_operator->Solution());
        mixedAnal->fsoltransfer.TransferFromMultiphysics();
        mixedAnal->PostProcessTimeStep();
    }
    
    
    
}

TMRSDataTransfer SettingSimpleFracCase(){
    
    TMRSDataTransfer sim_data;
    
    //    dim_name_and_physical_tagCoarse[3]["k33"] = 1;
    //    dim_name_and_physical_tagCoarse[3]["k31"] = 2;
    //    dim_name_and_physical_tagCoarse[2]["inlet"] = -2;
    //    dim_name_and_physical_tagCoarse[2]["outlet"] = -4;
    //    dim_name_and_physical_tagCoarse[2]["noflux"] = -1;
    
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["k33"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["k31"] = 2;
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = 10;
    sim_data.mTGeometry.mInterface_material_id = 100;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 4.0;
    REAL pressure_out = 1.0;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-4,D_Type,pressure_out);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-5,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(1);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] =
    std::make_tuple(-11,N_Type,zero_flux);
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-4,bc_outlet,0.0);
    
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
    
    //FracProperties
    sim_data.mTFracProperties.m_Permeability = 1.0e-3;
    
    
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



TPZGeoMesh *ReadFractureMesh(){
    std::string fileFine("../../FracMeshes/jose_simple2.msh");
    
    //    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
    /*
     2 4 "inlet"
     2 5 "outlet"
     2 6 "noflux"
     3 3 "k33"
     3 10 "k31"
     */
    
    dim_name_and_physical_tagFine[3]["c1"] = 1;
    dim_name_and_physical_tagFine[2]["inlet"] = -2;
    dim_name_and_physical_tagFine[2]["outlet"] = -4;
    dim_name_and_physical_tagFine[2]["noflux"] = -1;
    
    dim_name_and_physical_tagFine[2]["Fracture2"] = 10;
    dim_name_and_physical_tagFine[1]["BCfrac0"] = -11;
    /*
     2 2 "Fractures"
     3 1 "c1"
     */
    
    //    for(int i=1; i<=100; i++)
    //    {
    //        std::stringstream sout;
    //        sout << "c" << i;
    //        dim_name_and_physical_tagFine[3][sout.str()] = i+9;
    //    }
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    
    {
        REAL l = 1.0;
        GeometryFine.SetCharacteristiclength(l);
        GeometryFine.SetFormatVersion("4.1");
        
        GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
        gmeshFine = GeometryFine.GeometricGmshMesh(fileFine);
        GeometryFine.PrintPartitionSummary(std::cout);
    }
    
    {
        std::ofstream fileFine("mesh3dFine.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFine);
    }
    
    
    return gmeshFine;
}
