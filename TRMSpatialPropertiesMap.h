//
//  TRMSpatialPropertiesMap.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
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
#include "TRMSpatialMap.h"
#include "pzanalysis.h"
#include "pzinterpolationspace.h"

class TRMSpatialPropertiesMap{
    
    
private:
    
    /** @brief L2 computational mesh for compute spatial properties of SPE10 */
    TPZCompMesh * fSPE10Cmesh;

    /** @brief Gmsh grid file */
    std::string fGridName;
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > fPermPorFields;
    
    /** @brief number of blocks i, j and k  */
    TPZStack<int> fNBlocks;
    
    /** @brief size of blocks dx, dy and dz  */
    TPZStack<REAL> fBlocks_sizes;
    
    /** @brief size of blocks dx, dy and dz  */
    std::map<int64_t,int64_t> fm_data;
        
public:
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap();
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap(TRMSpatialPropertiesMap &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMSpatialPropertiesMap &operator=(TRMSpatialPropertiesMap &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief default destructor */
    ~TRMSpatialPropertiesMap();


    /** @brief Compute permeabilty diagonal tensor */
    void Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa);
    
    /** @brief compute porosity */
    void phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi);
    
    /** @brief Set SPE10 fields file */
    void SetSpatialFields(TPZStack<int> NBlocks, TPZStack<REAL> Blocks_sizes, std::pair< std::string , std::string > PermPorFields){
        fNBlocks = NBlocks;
        fBlocks_sizes = Blocks_sizes;
        fPermPorFields = PermPorFields;
    }
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > & SpatialFields(){
        return fPermPorFields;
    }
    
    /** @brief Get SPE10 fields file */
    TPZStack<int> & NBlocks(){
        return fNBlocks;
    }
    
    /** @brief Get SPE10 fields file */
    TPZStack<REAL> & Blocks_sizes(){
        return fBlocks_sizes;
    }
    
    /** @brief Load spatial properties from SPE10 cartesian intact fields Kx, ky, kz and phi */
    void LoadSPE10Map(bool PrintMapQ);
    
    /** @brief Load spatial properties from SPE10 cartesian intact fields Kx, ky, kz and phi */
    bool ComputePropertieSPE10Map(long & index, TPZVec<STATE> &x, TPZFMatrix<STATE> &kappa, TPZFMatrix<STATE> &inv_kappa, REAL & phi);

    /** @brief Insert spatial properties from SPE10 on pz mesh of order zero */
    bool Insert_Inside_Map(int n_data);
    
    /** @brief Get dof for spatial properties from SPE10 on pz mesh with connect solution (kx,ky,kz,phi) */
    void ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes);
    
    /** @brief Create a reservoir-box geometry */
    TPZGeoMesh * CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz);
    
    static void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    static void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    static void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    void ExpandGeomesh(TPZGeoMesh *gmesh, REAL sx, REAL  sy, REAL  sz);
    
    void TraslateGeomesh(TPZGeoMesh *gmesh, TPZVec<REAL> t_vec);
    
    void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
    
};

#endif /* defined(__PZ__TRMSpatialPropertiesMap__) */
