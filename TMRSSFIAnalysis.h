//
//  TMRSSFIAnalysis.h
//
//  Created by Omar Durán on 10/15/19.
//

#ifndef TMRSSFIAnalysis_h
#define TMRSSFIAnalysis_h

#include <stdio.h>
#include "TPZMultiphysicsCompMesh.h"
#include "TMRSDataTransfer.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TMRSApproxSpaceGenerator.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZAlgebraicDataTransfer.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

class TMRSSFIAnalysis {
    
private:
    
    /// Data transfer object
    TMRSDataTransfer * m_sim_data;

    TPZFMatrix<STATE> m_x_mixed;
    TPZFMatrix<STATE> m_x_transport;
    
public:
    
    /// Number of iterations
    int m_k_iteration = 0;
    
    /// Mixed module
    TMRSMixedAnalysis * m_mixed_module;
    
    /// Transport module
    TMRSTransportAnalysis * m_transport_module;
    TPZAlgebraicDataTransfer fAlgebraicDataTransfer;
    
    /// Default constructor
    TMRSSFIAnalysis();
    
    /// Default destructor
    ~TMRSSFIAnalysis();
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q);
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<REAL(const TPZVec<REAL> & )> & kx, std::function<REAL(const TPZVec<REAL> & )> & ky, std::function<REAL(const TPZVec<REAL> & )> & kz, std::function<REAL(const TPZVec<REAL> & )> & phi, std::function<REAL(const TPZVec<REAL> & )> & s0);
    
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q, std::function<std::vector<REAL>(const TPZVec<REAL> & )> & kappa_phi, std::function<REAL(const TPZVec<REAL> & )> & s0);
   
    void FillMaterialMemoryDarcy(int material_id);
    
    void FillProperties(std::string fileprops, TPZAlgebraicTransport *algebraicTransport);
    
    static  void ReadProperties(std::string name, bool print_table_Q, std::vector<REAL> &Kx, std::vector<REAL> &Ky, std::vector<REAL> &Kz, std::vector<REAL> &Phi);
    /// Configurates iternal members
    void Configure(int n_threads, bool UsePardiso_Q);
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer * sim_data);
    
    /// Get data transfer object
    TMRSDataTransfer * GetDataTransfer();
    
    /// Get the number of iterations
    int GetNumberOfIterations();
    
    /// Run a time step
    void RunTimeStep();
    
    /// Render a vtk file with requested variables for a time step
    void PostProcessTimeStep(int val=0);
    
    /// Perform a SFI iteration
    void SFIIteration();
    
    void TransferToTransportModule();
    
    void TransferToMixedModule();
    
    void UpdateMemoryMixedModule();
    
    void UpdateMemoryTransportModule();
    
    void UpdateMemoryInModules();
    void PostProcessResevoirProp();
    
    // transfer the permeability and lambda to the element solution for post processing
    void SetMixedMeshElementSolution(TPZCompMesh *cmesh);
    
};

#endif /* TMRSSFIAnalysis_h */
