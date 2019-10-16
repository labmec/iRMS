
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
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"

TMRSDataTransfer Setting2D();

TMRSDataTransfer Setting3D();

int main(){
 
    bool is_3D_Q = false;
    bool is_2D_Coarse_Q = true;
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
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->Configure(n_threads, UsePardiso_Q);
    sfi_analysis->SetDataTransfer(&sim_data);
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    for (int it = 0; it < n_steps; it++) {
        sfi_analysis->RunTimeStep();
        sfi_analysis->PostProcessTimeStep();
    }
    
    return 0;
    
    // Solving global darcy problem
    TMRSMixedAnalysis * mixed_analysis = new TMRSMixedAnalysis(mixed_operator,  must_opt_band_width_Q);
    mixed_analysis->Configure(n_threads, UsePardiso_Q);
    mixed_analysis->SetDataTransfer(&sim_data);
    mixed_analysis->RunTimeStep();
    mixed_analysis->PostProcessTimeStep();
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, sim_data, mixed_operator, true);
    mixed_analysis->AssembleResidual();
    TMRSApproxSpaceGenerator::SetUpdateMemory(2, sim_data, mixed_operator, false);
    
    // Solving transport problem
    TMRSTransportAnalysis * transport_analysis = new  TMRSTransportAnalysis(transport_operator,  must_opt_band_width_Q);
    transport_analysis->Configure(n_threads, UsePardiso_Q);
    transport_analysis->SetDataTransfer(&sim_data);
    
    int n_steps1 = sim_data.mTNumerics.m_n_steps;
    for (int it = 0; it < n_steps1; it++) {
        transport_analysis->RunTimeStep();
        transport_analysis->PostProcessTimeStep();
        // Updating memory
        TMRSApproxSpaceGenerator::SetUpdateMemory(2, sim_data, transport_operator, true);
        transport_analysis->AssembleResidual();
        TMRSApproxSpaceGenerator::SetUpdateMemory(2, sim_data, transport_operator, false);
    }
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
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 2;
    sim_data.mTNumerics.m_max_iter_transport = 2;
    sim_data.mTNumerics.m_max_iter_sfi = 2;
    sim_data.mTNumerics.m_res_tol_mixed = 1.0e-8;
    sim_data.mTNumerics.m_res_tol_transport = 1.0e-8;
    sim_data.mTNumerics.m_n_steps = 50;
    sim_data.mTNumerics.m_dt      = 1000.0;
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    
    
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
