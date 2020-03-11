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

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

class TMRSSFIAnalysis {
    
private:
    
    /// Data transfer object
    TMRSDataTransfer * m_sim_data;
    
    /// Number of iterations
    int m_k_iteration;
    
    /// Mixed module
    TMRSMixedAnalysis * m_mixed_module;
    
    TPZFMatrix<STATE> m_x_mixed;
    
    TPZFMatrix<STATE> m_x_transport;
    
   
    
public:
    /// Transport module
    TMRSTransportAnalysis * m_transport_module;
    
    /// Default constructor
    TMRSSFIAnalysis();
    
    /// Default destructor
    ~TMRSSFIAnalysis();
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSSFIAnalysis(TPZMultiphysicsCompMesh * cmesh_mixed, TPZMultiphysicsCompMesh * cmesh_transport, bool must_opt_band_width_Q);
    
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
    
};

#endif /* TMRSSFIAnalysis_h */
