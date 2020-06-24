//
//  TRMSpatialPropertiesMap.h
//  PZ
//
//  Created by omar duran on 21/06/2020.
//
//

#ifndef __PZ__TRMSpatialPropertiesMap__
#define __PZ__TRMSpatialPropertiesMap__

#include <stdio.h>
#include "pzmanvector.h"
#include "pzfmatrix.h"

#include "pzgmesh.h"
#include "tpzhierarquicalgrid.h"
#include "pzcmesh.h"
#include "pzl2projection.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "pzanalysis.h"
#include "pzinterpolationspace.h"

class TRMSpatialPropertiesMap{
    
    
public:
    
    /// Enumerate defining the map type
    enum EMapType { ENone = -1, ECartesianGrid = 0, ECornerPointGrid = 1};
    
    EMapType m_map_type;
    
    std::vector<std::vector<REAL>> m_grid_coordinates;
    
    std::vector<std::vector<REAL>> m_properties;
    
    REAL m_kappa_default;
    
    REAL m_phi_default;
    
    // Structured Auxiliary Mesh with Element Storage (SAMe)
    std::vector<REAL> m_n_SAMe_blocks;
    
    std::vector<REAL> m_size_SAMe_blocks;
    
    // SPE case 10
    std::vector<size_t> m_n_spe_blocks;
       
    std::vector<REAL> m_size_spe_blocks;
    
    std::vector<REAL> m_spe_translation;
    
    std::map<std::vector<unsigned int>,std::vector<std::vector<unsigned int>>> m_SAMe_to_SPE;
    
        
public:
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap();
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap(TRMSpatialPropertiesMap &other)
    {
        m_grid_coordinates  = other.m_grid_coordinates;
        m_properties        = other.m_properties;
    }
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap &operator=(TRMSpatialPropertiesMap &other)
    {
        // check for self-assignment
        if(&other == this){
            return *this;
        }
        m_grid_coordinates  = other.m_grid_coordinates;
        m_properties        =   other.m_properties;
        
        return *this;
    }
    
    /** @brief default destructor */
    ~TRMSpatialPropertiesMap();

    
    /** @brief Evaluate the permeabilty as diagonal tensor and porosity */
    void SampleKappaAndPhi(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi);

    void SampleKappaAndPhiCartesianGrid(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi);
    
    void SampleKappaAndPhiCornerGrid(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi);
    
    /** @brief Set Carterian mesh data  */
    void SetCartesianMeshData(std::vector<size_t> n_blocks, std::vector<REAL> size_blocks, std::string perm_data_name, std::string phi_data_name, std::vector<REAL> translation = {0,0,0});
    
    void BuildSAMe();

    /** @brief Set Corner grid mesh data  */
    void SetCornerGridMeshData(size_t n_cells, std::string corner_data_name, std::string props_data_name);
    
};

#endif /* defined(__PZ__TRMSpatialPropertiesMap__) */
