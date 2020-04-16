//
//  RSimulatorConfiguration.hpp
//  MonophasicTest
//  Created by Jose on 27/8/19.
//

#ifndef RSimulatorConfiguration_hpp
#define RSimulatorConfiguration_hpp

#include <stdio.h>
#include "includes.h"

class RSimulatorConfiguration{
    
private:
    TPZGeoMesh *fGmesh;
    TPZMultiphysicsCompMesh *multmesh;
    SimulationCase fsim_case;
    
public:
    RSimulatorConfiguration();
    RSimulatorConfiguration(SimulationCase sim_case);
    
    RSimulatorConfiguration(TPZGeoMesh *gmesh);
    RSimulatorConfiguration(RSimulatorConfiguration *conf);
    
    void CreateGeomesh(int nx, int ny, double l, double h, MMeshType type);
    TPZGeoMesh * GetGeomesh();
    
    TPZCompMesh * CreateFluxCmesh(TPZGeoMesh * gmesh, int order);
    TPZCompMesh * CreatePressureCmesh(TPZGeoMesh * gmesh, int order);
    TPZCompMesh * CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh, int ref);
    void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);
    TPZMultiphysicsCompMesh *CreateMultiPhysicsCompMesh(TPZGeoMesh * gmesh);
       TPZMultiphysicsCompMesh * MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZManVector<TPZCompMesh *> meshvec);
    
    TPZAnalysis * CreateAnalysis(TPZMultiphysicsCompMesh * cmesh_c,  bool must_opt_band_width_Q, int number_threads, bool UsePardiso_Q);
    
    void PrintCmesh(int mesh_index, std::ofstream &file_name);
    void PrintGmesh( std::ofstream &file_name);
    void PosProcess();
    void Run();
    void InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh);
    void InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes);
    TPZAnalysis * CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
    TPZFMatrix<STATE> TimeForward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt, TPZFMatrix<STATE> & M_diag);
};

#endif /* RSimulatorConfiguration_hpp */
