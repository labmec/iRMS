//
//  TRMSpatialPropertiesMap.cpp
//  PZ
//
//  Created by omar duran on 21/06/2020.
//
//

#include "TRMSpatialPropertiesMap.h"


TRMSpatialPropertiesMap::TRMSpatialPropertiesMap(){
    m_grid_coordinates.clear();
    m_properties.clear();
    m_kappa_default = 1.0e-7;
    m_phi_default = 0.1;
}

TRMSpatialPropertiesMap::~TRMSpatialPropertiesMap(){
    
}

void TRMSpatialPropertiesMap::SampleKappaAndPhi(TPZManVector<REAL,3> &x, std::vector<REAL> &kappa_and_phi){

    std::vector<REAL> pt = {x[0],x[1],x[2]};
    
    auto substract = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {a[0] - b[0],a[1] - b[1],a[2] - b[2]};
    };
    
    auto cross = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> std::vector<REAL> {
            return {-(a[2]*b[1]) + a[1]*b[2],a[2]*b[0] - a[0]*b[2],-(a[1]*b[0]) + a[0]*b[1]};
    };
    
    auto dot = [] (const std::vector<REAL> & a, const std::vector<REAL> & b) -> REAL {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };
    
    auto sign = [] (const REAL a) -> int {
            return (REAL(0.0) < a) - (a < REAL(0.0));
    };
    
    auto SameSide = [&substract,&cross,&dot,&sign] (const std::vector<REAL> & v1, const std::vector<REAL> & v2, const std::vector<REAL> & v3, const std::vector<REAL> & v4, const std::vector<REAL> & pt) -> bool {
        std::vector<REAL> normal = cross(substract(v2,v1),substract(v3,v1));
        REAL dot_v4 = dot(normal,substract(v4,v1));
        REAL dot_pt = dot(normal,substract(pt,v1));
        bool checkQ = sign(dot_v4) == sign(dot_pt);
        return checkQ;
    };
    
    auto IsTetrahedronMemberQ = [&SameSide] (const std::vector<REAL> & v1, const std::vector<REAL> & v2, const std::vector<REAL> & v3, const std::vector<REAL> & v4, const std::vector<REAL> & pt) -> bool {
        bool check1Q = SameSide(v1, v2, v3, v4, pt);
        bool check2Q = SameSide(v2, v3, v4, v1, pt);
        bool check3Q = SameSide(v3, v4, v1, v2, pt);
        bool check4Q = SameSide(v4, v1, v2, v3, pt);
        return check1Q && check2Q && check3Q && check4Q;
    };
    
    unsigned int pos;
    bool IsCellMemberQ = false;
    for (unsigned int i = 0; i < m_grid_coordinates.size(); i++) {
        std::vector<REAL> & cell_data = m_grid_coordinates[i];
        
        std::vector<std::vector<REAL>> cell(8);
        unsigned int c = 0;
        for (unsigned int j = 0; j < 24; j += 3) {
            cell[c] = {cell_data[j],cell_data[j+1],cell_data[j+2]};
            c++;
        }
        
        std::vector<REAL> normal = cross(substract(cell[1],cell[0]),substract(cell[2],cell[0]));
        bool isCollapsedQ = dot(normal,normal) < 1.0e-12;
        if (isCollapsedQ) {
            continue;
        }
        
        bool checktet1Q = IsTetrahedronMemberQ(cell[0],cell[1],cell[3],cell[4],pt);
        bool checktet2Q = IsTetrahedronMemberQ(cell[3],cell[7],cell[4],cell[1],pt);
        bool checktet3Q = IsTetrahedronMemberQ(cell[4],cell[5],cell[7],cell[1],pt);
        bool checktet4Q = IsTetrahedronMemberQ(cell[1],cell[2],cell[3],cell[6],pt);
        bool checktet5Q = IsTetrahedronMemberQ(cell[3],cell[6],cell[7],cell[1],pt);
        bool checktet6Q = IsTetrahedronMemberQ(cell[5],cell[6],cell[7],cell[1],pt);
        
        IsCellMemberQ = checktet1Q || checktet2Q || checktet3Q || checktet4Q || checktet5Q || checktet6Q;
        
        if (IsCellMemberQ) {
            pos = i;
            break;
        }
        
    }
    
    if (!IsCellMemberQ) {
        kappa_and_phi = {m_kappa_default,m_kappa_default,m_kappa_default,m_phi_default};
        return ;
    }else{
        std::vector<REAL> chunk = m_properties[pos];
        kappa_and_phi = chunk;
    }
}

void TRMSpatialPropertiesMap::SetCornerGridMeshData(size_t n_cells, std::string corner_data_name, std::string props_data_name){
    
    std::ifstream stream_corners (corner_data_name.c_str());
    std::ifstream stream_props (props_data_name.c_str());
    REAL x,y,z,kx,ky,kz,phi;
    m_grid_coordinates.resize(n_cells);
    m_properties.resize(n_cells);
    for (unsigned int i = 0; i < n_cells; i++) {
        
        for (unsigned int j = 0; j < 24; j += 3) {
            stream_corners >> x;
            stream_corners >> y;
            stream_corners >> z;
            
            m_grid_coordinates[i].push_back(x);
            m_grid_coordinates[i].push_back(y);
            m_grid_coordinates[i].push_back(z);
        }
        
        stream_props >> kx;
        stream_props >> ky;
        stream_props >> kz;
        stream_props >> phi;
        
        m_properties[i].push_back(kx);
        m_properties[i].push_back(ky);
        m_properties[i].push_back(kz);
        m_properties[i].push_back(phi);
    }
}
