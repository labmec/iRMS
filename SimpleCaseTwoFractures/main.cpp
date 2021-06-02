
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
#include "pzlog.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "pzlog.h"

void CaseSimple2Frac();
TMRSDataTransfer Setting2Fractures();
TPZGeoMesh *ReadFractureMesh();



int main(){
    TPZLogger::InitializePZLOG();
    CaseSimple2Frac();

}

void CaseSimple2Frac()
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
    // vector with subdomain index of the geometric elements
    TPZVec<int64_t> subdomain;
    //    TPZGeoMesh *gmesh = ReadFractureMesh(subdomain);
    TPZGeoMesh *gmesh = ReadFractureMesh();
    
    TMRSApproxSpaceGenerator aspace;
    TMRSDataTransfer sim_data  = Setting2Fractures();
    sim_data.mTGeometry.mSkeletonDiv =0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
    //mSimData.mTGeometry.mDomainDimNameAndPhysicalTag
    aspace.SetGeometry(gmesh);
    aspace.SetSubdomainIndexes(subdomain);
    //    std::ofstream name("fractureTest.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name, subdomain);
    
    
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    
    //This parameter should be always "true"
    bool UsingPzSparse = false;
    bool UsePardiso_Q = true;
    
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    mixAnalisys->PostProcessTimeStep();
    delete gmesh;
}

TMRSDataTransfer Setting2Fractures(){
    
    TMRSDataTransfer sim_data;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["k33"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["k31"] = 2;
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = 10;
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mIterface_material_idFracBound = 104;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 4.0 ;
    REAL pressure_out = 1.0 ;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-4,D_Type,pressure_out);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-5,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] =
    std::make_tuple(-11,N_Type,zero_flux);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[1] =
    std::make_tuple(-12,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[2] =
    std::make_tuple(-13,D_Type,pressure_out);
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-4,bc_outlet,0.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.1;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 1;
    sim_data.mTNumerics.m_max_iter_sfi = 1;
    //BorderElementOrder
    sim_data.mTNumerics.m_BorderElementPresOrder=1;
    sim_data.mTNumerics.m_BorderElementFluxOrder=1;
    
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    //FracAndReservoirProperties
        sim_data.mTFracProperties.m_Permeability = 0.00001;
        REAL kappa=1.0;
        int  id1=1;
        int  id2=2;
        std::vector<std::pair<int, REAL>> idPerm(2);
        idPerm[0]= std::make_pair(id1,kappa);
        idPerm[1]= std::make_pair(id2,kappa);
        sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    
    
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
    return sim_data;
}




TPZGeoMesh *ReadFractureMesh(){
    std::string fileFine("../../FracMeshes/Case2Frac.msh");
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
    /*
     2 4 "inlet"
     2 5 "outlet"
     2 6 "noflux"
     3 3 "k33"
     3 10 "k31"
     */
    
    //    dim_name_and_physical_tagFine[3]["c1"] = 1;
    dim_name_and_physical_tagFine[3]["k33"] = 2;
    dim_name_and_physical_tagFine[3]["k31"] = 1;
    dim_name_and_physical_tagFine[2]["inlet"] = -2;
    dim_name_and_physical_tagFine[2]["outlet"] = -4;
    dim_name_and_physical_tagFine[2]["noflux"] = -1;
    
    dim_name_and_physical_tagFine[2]["Fractures"] = 10;
    dim_name_and_physical_tagFine[2]["Fracture2"] = 10;
    dim_name_and_physical_tagFine[1]["BCfrac0"] = -11;
    
    dim_name_and_physical_tagFine[1]["BCFracInlet"] = -12;
    dim_name_and_physical_tagFine[1]["BCFracOutlet"] = -13;
    dim_name_and_physical_tagFine[1]["FracIntersection"] = -14;
    
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    GeometryFine.SetFormatVersion("4.1");
    
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(fileFine);
    GeometryFine.PrintPartitionSummary(std::cout);
    std::ofstream fileFinevtk("mesh3dFine.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFinevtk);
    return gmeshFine;
}
