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
    
    std::vector<std::vector<REAL>> m_grid_coordinates;
    
    std::vector<std::vector<REAL>> m_properties;
    
    REAL m_kappa_default;
    
    REAL m_phi_default;
        
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

    
    /** @brief Set Carterian mesh data  */
    void SetCartesianMeshData(){
        
    }

    /** @brief Set Corner grid mesh data  */
    void SetCornerGridMeshData(size_t n_cells, std::string corner_data_name, std::string props_data_name);
    
};

#endif /* defined(__PZ__TRMSpatialPropertiesMap__) */
