//
//  TMRSMixedAnalysis.hpp
//
//  Created by Omar Durán on 10/15/19.
//

#ifndef TMRSMixedAnalysis_h
#define TMRSMixedAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TMRSDataTransfer.h"
#include "TPZMFSolutionTransfer.h"


class TMRSMixedAnalysis : public TPZAnalysis {
    
private:
    
    /// Data transfer object
    TMRSDataTransfer * m_sim_data;
    

    /// Number of iterations
    int m_k_iteration;
    
public:
    
    TPZMFSolutionTransfer m_soltransportTransfer;

    
    /// Default constructor
    TMRSMixedAnalysis();
    
    /// Default destructor
    ~TMRSMixedAnalysis();
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSMixedAnalysis(TPZMultiphysicsCompMesh * cmesh_mult, bool must_opt_band_width_Q);
    
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
    void PostProcessTimeStep();
    
    /// Perform a Newton iteration
    void NewtonIteration();
    
};

#endif /* TMRSMixedAnalysis_h */
