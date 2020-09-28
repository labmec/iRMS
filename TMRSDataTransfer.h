//
//  TMRSDataTransfer.hpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#ifndef TMRSDataTransfer_hpp
#define TMRSDataTransfer_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include "TMRSSavable.h"
#include "pzmanvector.h"
#include <map> // for map
#include <tuple> // for tuple
#include "TRSLinearInterpolator.h"


/// Object that represents GUI state and store all the required input/output data
class TMRSDataTransfer : public TMRSSavable {
    
public:
    
    /// Default constructor
    TMRSDataTransfer();
    
    std::string mSimulationName = "";
    
    /// Copy constructor
    TMRSDataTransfer(const TMRSDataTransfer &other);
    
    // Copy assignment operator
    TMRSDataTransfer &operator=(const TMRSDataTransfer &other);
    
    /// Destructor
    ~TMRSDataTransfer();
    
    /// Write object state
    void Write(TPZStream &buf, int withclassid) const;
    
    /// Read object state
    void Read(TPZStream &buf, void *context);
    
    /// Read object state
    virtual int ClassId() const;
    
    /**
     * @brief Class that stores geometric information
     */
    class TGeometry : public TMRSSavable {
        
    public:
        
        /** @brief
         Contains the dimension, name and physical tag of the domain. */
        std::vector<std::map<std::string,int>> mDomainDimNameAndPhysicalTag;
        
        /** @brief
         Contains the dimension and physical tag of the fractures. */
        std::vector<std::map<std::string,int>> mDomainFracDimNameAndPhysicalTag;
        
        /** @brief
          MaterialID of the interface element that will be inserted in the transport mesh
         */
        int mInterface_material_id = 100;
        int mSkeletonDiv =0;
        int mnLayers = 1;
        int mnref =0;
        std::string mGmeshFileName="";
        
        
         /** @brief Default constructor */
        TGeometry(){
            
            mDomainDimNameAndPhysicalTag.resize(4);
            mDomainFracDimNameAndPhysicalTag.resize(3);
        }
        
        /** @brief Destructor */
        ~TGeometry(){
            
        }
        
        /** @brief Copy constructor */
        TGeometry(const TGeometry &other){
            mDomainDimNameAndPhysicalTag = other.mDomainDimNameAndPhysicalTag;
            mDomainFracDimNameAndPhysicalTag = other.mDomainFracDimNameAndPhysicalTag;
            mInterface_material_id=other.mInterface_material_id;
            mSkeletonDiv = other.mSkeletonDiv;
            mnLayers = other.mnLayers;
            mnref = other.mnref;
            mGmeshFileName = other.mGmeshFileName;
            
        }
        /** @brief Copy assignment operator*/
        TGeometry &operator=(const TGeometry &other){
            if (this != & other) // prevent self-assignment
            {
                mDomainDimNameAndPhysicalTag = other.mDomainDimNameAndPhysicalTag;
                mDomainFracDimNameAndPhysicalTag = other.mDomainFracDimNameAndPhysicalTag;
                mInterface_material_id=other.mInterface_material_id;
                mSkeletonDiv = other.mSkeletonDiv;
                mnLayers = other.mnLayers;
                mnref = other.mnref;
                mGmeshFileName = other.mGmeshFileName;
               
            }
            return *this;
        }
        
    };
    
    /**
     * @brief Class that stores PetroPhysics information
     */
   
    class TPetroPhysics : public TMRSSavable {
        
    public:
        REAL mOilViscosity;
        REAL mWaterViscosity;
        
        /** @brief Contains the water relative permeability model for each layer */
        std::vector<TRSLinearInterpolator > mLayer_Krw_RelPerModel;
        /** @brief Contains the oil relative permeability model for each layer */
        std::vector<TRSLinearInterpolator > mLayer_Kro_RelPerModel;
        
        std::function<std::tuple<REAL, REAL>(REAL &)> mKro;
        std::function<std::tuple<REAL, REAL>(REAL &)> mKrw;
        std::function<std::tuple<REAL, REAL>(REAL &)> mFo;
        std::function<std::tuple<REAL, REAL>(REAL &)> mFw;
        std::function<std::tuple<REAL, REAL>(REAL &)> mLambdaW;
        std::function<std::tuple<REAL, REAL>(REAL &)> mLambdaO;
        std::function<std::tuple<REAL, REAL>(REAL &)> mLambdaTotal;
        
        /** @brief Default constructor */
        TPetroPhysics(){
            mOilViscosity=1.0;
            mWaterViscosity=1.0;
            CreateLinearKrModel();
            mLayer_Krw_RelPerModel.clear();
            mLayer_Kro_RelPerModel.clear();
        }
        
        /** @brief Destructor */
        ~TPetroPhysics(){
            
        }
        
        /** @brief Copy constructor */
        TPetroPhysics(const TPetroPhysics &other){
            mOilViscosity=other.mOilViscosity;
            mWaterViscosity=other.mWaterViscosity;
            mKro = other.mKro;
            mKrw = other.mKrw;
            mFo  = other.mFo;
            mFw  = other.mFw;
            mLambdaW = other.mLambdaW;
            mLambdaO = other.mLambdaO;
            mLambdaTotal = other.mLambdaTotal;
            mLayer_Krw_RelPerModel=other.mLayer_Krw_RelPerModel;
            mLayer_Kro_RelPerModel=other.mLayer_Kro_RelPerModel;
        }
        
        /** @brief Copy assignment operator*/
        TPetroPhysics &operator=(const TPetroPhysics &other){
            if (this != & other) // prevent self-assignment
            {
                mOilViscosity=other.mOilViscosity;
                mWaterViscosity=other.mWaterViscosity;
                mKro = other.mKro;
                mKrw = other.mKrw;
                mFo  = other.mFo;
                mFw  = other.mFw;
                mLambdaW = other.mLambdaW;
                mLambdaO = other.mLambdaO;
                 mLambdaTotal = other.mLambdaTotal;
                mLayer_Krw_RelPerModel=other.mLayer_Krw_RelPerModel;
                mLayer_Kro_RelPerModel=other.mLayer_Kro_RelPerModel;
               
            }
            return *this;
        }
        void CreateLinearKrModel();
        void CreateQuadraticKrModel();
        void UpdateLambdasAndFracFlows();
    };
    
    /**
     * @brief Class that stores Fluid properties
     */
    class TFluidProperties : public TMRSSavable {
        
    public:
        REAL mOilViscosity;
        REAL mWaterViscosity;
        REAL mWaterCompressibility;
        REAL mOilCompressibility;
        REAL mOilDensityRef;
        REAL mWaterDensityRef;
        REAL mReferencePressure;
        std::function<std::tuple<REAL, REAL>(REAL &)> mOilDensityF;
        std::function<std::tuple<REAL, REAL>(REAL &)> mWaterDensityF;
        
        /** @brief Default constructor */
        TFluidProperties(){
            mOilViscosity= 1.0;
            mWaterViscosity=1.0;
            mOilDensityRef = 800;
            mWaterDensityRef = 1000;
            mWaterCompressibility = 1.0e-8;
            mOilCompressibility = 1.0e-7;
            mReferencePressure = 1.013e5;
            mOilDensityF = 0;
            mWaterDensityF =0;
            CreateLinearDensityFunction();
        }
        
        /** @brief Destructor */
        ~TFluidProperties(){
            
        }
        
        /** @brief Copy constructor */
        TFluidProperties(const TFluidProperties &other){
            mOilViscosity = other.mOilViscosity;
            mWaterViscosity = other.mWaterViscosity;
            mWaterCompressibility = other.mWaterCompressibility;
            mOilCompressibility = other.mOilCompressibility;

            mReferencePressure = other.mReferencePressure;
            mOilDensityRef = other.mOilDensityRef;
            mWaterDensityRef = other.mWaterDensityRef;
            mOilDensityF = other.mOilDensityF;
            mWaterDensityF =other.mWaterDensityF;

        }
        
        /** @brief Copy assignment operator*/
        TFluidProperties &operator=(const TFluidProperties &other){
            mOilViscosity = other.mOilViscosity;
            mWaterViscosity = other.mWaterViscosity;
            mWaterCompressibility = other.mWaterCompressibility;
            mOilCompressibility = other.mOilCompressibility;
            mOilDensityRef = other.mOilDensityRef;
            mWaterDensityRef = other.mWaterDensityRef;
             mReferencePressure = other.mReferencePressure;
            mOilDensityF = other.mOilDensityF;
            mWaterDensityF =other.mWaterDensityF;
            return *this;
        }
        void CreateLinearDensityFunction();
        void CreateExponentialDensityFunction();
        
    };
    class TReservoirProperties : public TMRSSavable {
        
    public:
        REAL mPorosity;
        bool fPropsFromPreProcess = false;
        std::function<std::vector<REAL>(const TPZVec<REAL> & )> kappa_phi;
        std::function<REAL(const TPZVec<REAL> & )> s0;
        std::string mPropsFileName="";
        /** @brief Default constructor */
        TReservoirProperties(){
            mPorosity=1.0;
            fPropsFromPreProcess = false;
        }
        
        /** @brief Destructor */
        ~TReservoirProperties(){
            
        }
        
        /** @brief Copy constructor */
        TReservoirProperties(const TReservoirProperties &other){
            mPorosity = other.mPorosity;
            fPropsFromPreProcess = other.fPropsFromPreProcess;
            kappa_phi = other.kappa_phi;
            s0=other.s0;
            mPropsFileName=other.mPropsFileName;
        }
        
        /** @brief Copy assignment operator*/
        TReservoirProperties &operator=(const TReservoirProperties &other){
            mPorosity = other.mPorosity;
            fPropsFromPreProcess = other.fPropsFromPreProcess;
            kappa_phi = other.kappa_phi;
            s0=other.s0;
            mPropsFileName=other.mPropsFileName;
            return *this;
        }
        
        
    };
    
    /**
     * @brief Class that stores multiphase functions
     */
    class TMultiphaseFunctions : public TMRSSavable {
        
    public:
        
        /**
         * contains the fractional flow of water for each layer
         */
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_fw;
        /**
         * contains the fractional flow of oil for each layer
         */
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_fo;
        
        /**
         * contains the mobility (lambda) for each layer
         */
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_lambda;
        
        /**
         * contains the Gravitational term (Glambda) for each layer
         */
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_Glambda;
        
        /** @brief Default constructor */
        TMultiphaseFunctions(){
            mLayer_fw.clear();
            mLayer_fo.clear();
            mLayer_lambda.clear();
            mLayer_Glambda.clear();
        }
        
        /** @brief Destructor */
        ~TMultiphaseFunctions(){
            
        }
        
        /** @brief Copy constructor */
        TMultiphaseFunctions(const TMultiphaseFunctions &other){
            mLayer_fw = other.mLayer_fw;
            mLayer_fo = other.mLayer_fo;
            mLayer_lambda = other.mLayer_lambda;
            mLayer_Glambda = other.mLayer_Glambda;
            
        }
        
        /** @brief Copy assignment operator*/
        TMultiphaseFunctions &operator=(const TMultiphaseFunctions &other){
            if (this != & other) // prevent self-assignment
            {
                mLayer_fw = other.mLayer_fw;
                mLayer_fo = other.mLayer_fo;
                mLayer_lambda = other.mLayer_lambda;
                mLayer_Glambda = other.mLayer_Glambda;
            }
            return *this;
        }
        
        
    };
    
    /**
     * @brief Class that stores the boundary conditions of the problem
     */
    class TBoundaryConditions : public TMRSSavable {
        
    public:
        
        /**
         * @brief Contains the boundary conditions (material_id), condition type and value of the mixed problem
         */
        TPZManVector<std::tuple<int, int, REAL>> mBCMixedPhysicalTagTypeValue;
        /**
         * @brief Contains the boundary conditions (material_id), condition type and value of the fractures
         */
        TPZManVector<std::tuple<int, int, REAL>> mBCMixedFracPhysicalTagTypeValue;
        /**
         * @brief Contains the boundary conditions (material_id), condition type and value of the transport problem
         */
        TPZManVector<std::tuple<int, int, REAL>> mBCTransportPhysicalTagTypeValue;
        
        /**
         * @brief Contains the boundary conditions (material_id), condition type and value of fractures in the transport problem
         */
        TPZManVector<std::tuple<int, int, REAL>> mBCTransportFracPhysicalTagTypeValue;
       
        /** @brief Default constructor */
        TBoundaryConditions(){
            
            mBCMixedPhysicalTagTypeValue.Resize(0);
            
            mBCMixedFracPhysicalTagTypeValue.Resize(0);
            
            mBCTransportPhysicalTagTypeValue.Resize(0);
            
            mBCTransportFracPhysicalTagTypeValue.Resize(0);
        }
        
        /** @brief Destructor */
        ~TBoundaryConditions(){
            
        }
        
        /** @brief Copy constructor */
        TBoundaryConditions(const TBoundaryConditions &other){
            mBCMixedPhysicalTagTypeValue = other.mBCMixedPhysicalTagTypeValue;
            mBCMixedFracPhysicalTagTypeValue = other.mBCMixedFracPhysicalTagTypeValue;
            mBCTransportPhysicalTagTypeValue = other.mBCTransportPhysicalTagTypeValue;
            mBCTransportFracPhysicalTagTypeValue = other.mBCTransportFracPhysicalTagTypeValue;
        }
        
        /** @brief Copy assignment operator*/
        TBoundaryConditions &operator=(const TBoundaryConditions &other){
            if (this != & other) // prevent self-assignment
            {
                mBCMixedPhysicalTagTypeValue = other.mBCMixedPhysicalTagTypeValue;
                mBCMixedFracPhysicalTagTypeValue = other.mBCMixedFracPhysicalTagTypeValue;
                mBCTransportPhysicalTagTypeValue = other.mBCTransportPhysicalTagTypeValue;
                mBCTransportFracPhysicalTagTypeValue = other.mBCTransportFracPhysicalTagTypeValue;
            }
            return *this;
        }
        
    };
    
    /**
     * @brief Class that stores the numerical parameters of the simulation
     */
    class TNumerics : public TMRSSavable {
        
    public:
        
        /**
         * @brief time step size
         */
        REAL m_dt;
     
        
        /**
         * @brief Residual tolerance for mixed operator
         */
        REAL m_res_tol_mixed;
        
        /**
         * @brief Residual tolerance for transport operator
         */
        REAL m_res_tol_transport;
        
        /**
         * @brief Correction tolerance for mixed operator
         */
        REAL m_corr_tol_mixed;
        
        /**
         * @brief Correction tolerance for transport operator
         */
        REAL m_corr_tol_transport;
        
        REAL m_sfi_tol;
        /**
         * @brief Maximum number of iterations per time step for mixed operator
         */
        int m_max_iter_mixed;
        
        /**
         * @brief Maximum number of iterations per time step for transport operator
         */
        int m_max_iter_transport;
        
        /**
         * @brief Maximum number of Sequential Fully Implicit (SFI) iterations per time step
         */
        int m_max_iter_sfi;
        
        /**
         * @brief Number of time steps
         */
        int m_n_steps;
        
        /**
         * @brief Directive for the use of four spaces
         */
        bool m_four_approx_spaces_Q;
        
        /**
         * @brief Directive MHM mixed approximation
         */
        bool m_mhm_mixed_Q;
        
        std::vector<REAL> m_gravity;
        
        bool m_ISLinearKrModelQ;
        
        int m_nThreadsMixedProblem = 0;
        /** @brief Default constructor */
        TNumerics(){
            
            m_dt                    = 0.0;
            m_res_tol_mixed         = 1.0e-4;
            m_res_tol_transport     = 1.0e-7;
            m_corr_tol_mixed        = 1.0e-4;
            m_corr_tol_transport    = 1.0e-7;
            m_max_iter_mixed        = 0;
            m_max_iter_transport    = 0;
            m_max_iter_sfi          = 0;
            m_sfi_tol               = 1.0e-3;
            m_n_steps               = 0;
            m_four_approx_spaces_Q  = false;
            m_mhm_mixed_Q           = false;
            m_gravity.resize(3,0.0);
            m_gravity[2] = -10.0;
            m_ISLinearKrModelQ = true;
            m_nThreadsMixedProblem  =0;
            
        }
         /** @brief Destructor */
        ~TNumerics(){
            
        }
        
        /** @brief Copy constructor */
        TNumerics(const TNumerics & other){
            
            m_dt                    = other.m_dt;
            m_res_tol_mixed         = other.m_res_tol_mixed;
            m_res_tol_transport     = other.m_res_tol_transport;
            m_corr_tol_mixed        = other.m_corr_tol_mixed;
            m_corr_tol_transport    = other.m_corr_tol_transport;
            m_max_iter_mixed        = other.m_max_iter_mixed;
            m_max_iter_transport    = other.m_max_iter_transport;
            m_max_iter_sfi          = other.m_max_iter_sfi;
            m_sfi_tol               = other.m_sfi_tol;
            m_n_steps               = other.m_n_steps;
            m_four_approx_spaces_Q  = other.m_four_approx_spaces_Q;
            m_mhm_mixed_Q           = other.m_mhm_mixed_Q;
            m_gravity               = other.m_gravity;
            m_ISLinearKrModelQ      = other.m_ISLinearKrModelQ;
            m_nThreadsMixedProblem  = other.m_nThreadsMixedProblem;
        }
        
        /** @brief Copy assignment operator*/
        TNumerics & operator=(const TNumerics &other){
            
            // check for self-assignment
            if(&other == this){
                return *this;
            }
            
            m_dt                    = other.m_dt;
            m_res_tol_mixed         = other.m_res_tol_mixed;
            m_res_tol_transport     = other.m_res_tol_transport;
            m_corr_tol_mixed        = other.m_corr_tol_mixed;
            m_corr_tol_transport    = other.m_corr_tol_transport;
            m_max_iter_mixed        = other.m_max_iter_mixed;
            m_max_iter_transport    = other.m_max_iter_transport;
            m_max_iter_sfi          = other.m_max_iter_sfi;
            m_sfi_tol               = other.m_sfi_tol;
            m_n_steps               = other.m_n_steps;
            m_four_approx_spaces_Q  = other.m_four_approx_spaces_Q;
            m_mhm_mixed_Q           = other.m_mhm_mixed_Q;
            m_gravity               = other.m_gravity;
            m_ISLinearKrModelQ      = other.m_ISLinearKrModelQ;
            m_nThreadsMixedProblem  = other.m_nThreadsMixedProblem;
            return *this;
        }
        
        bool operator==(const TNumerics &other){
            
            // check for self-assignment
            if(&other == this){
                return true;
            }
            
            return
            m_dt                    == other.m_dt &&
            m_res_tol_mixed         == other.m_res_tol_mixed &&
            m_res_tol_transport     == other.m_res_tol_transport &&
            m_corr_tol_mixed        == other.m_corr_tol_mixed &&
            m_corr_tol_transport    == other.m_corr_tol_transport &&
            m_max_iter_mixed        == other.m_max_iter_mixed &&
            m_max_iter_transport    == other.m_max_iter_transport &&
            m_max_iter_sfi          == other.m_max_iter_sfi &&
            m_sfi_tol               ==  other.m_sfi_tol &&
            m_n_steps               == other.m_n_steps &&
            m_four_approx_spaces_Q  == other.m_four_approx_spaces_Q &&
            m_mhm_mixed_Q           == other.m_mhm_mixed_Q&&
            m_gravity               == other.m_gravity&&
            m_ISLinearKrModelQ      == other.m_ISLinearKrModelQ&&
            m_nThreadsMixedProblem  == other.m_nThreadsMixedProblem;
        }
        
        void Write(TPZStream &buf, int withclassid) const{ //ok
            buf.Write(&m_dt);
            buf.Write(&m_res_tol_mixed);
            buf.Write(&m_res_tol_transport);
            buf.Write(&m_corr_tol_mixed);
            buf.Write(&m_max_iter_mixed);
            buf.Write(&m_max_iter_transport);
            buf.Write(&m_max_iter_sfi);
            buf.Write(&m_sfi_tol);
            buf.Write(&m_n_steps);
            int temp = m_four_approx_spaces_Q;
            buf.Write(&temp);
            temp = m_mhm_mixed_Q;
            buf.Write(&temp);
            buf.Write(m_gravity);
            buf.Write(m_ISLinearKrModelQ);
            buf.Write(m_nThreadsMixedProblem);
        }
        
        void Read(TPZStream &buf, void *context){ //ok
            buf.Read(&m_dt);
            buf.Read(&m_res_tol_mixed);
            buf.Read(&m_res_tol_transport);
            buf.Read(&m_corr_tol_mixed);
            buf.Read(&m_corr_tol_transport);
            buf.Read(&m_max_iter_mixed);
            buf.Read(&m_max_iter_transport);
            buf.Read(&m_max_iter_sfi);
            buf.Read(&m_sfi_tol);
            buf.Read(&m_n_steps);
            int temp;
            buf.Read(&temp);
            m_four_approx_spaces_Q = temp;
            buf.Read(&temp);
            m_mhm_mixed_Q = temp;
            buf.Read(m_ISLinearKrModelQ);
            buf.Read(&m_nThreadsMixedProblem);
           
        }
        
        virtual int ClassId() const {
            return Hash("TMRSDataTransfer::TNumerics");
        }
        
        void Print() const {
            std::cout << m_dt << std::endl;
            std::cout << m_res_tol_mixed << std::endl;
            std::cout << m_res_tol_transport << std::endl;
            std::cout << m_corr_tol_mixed << std::endl;
            std::cout << m_corr_tol_transport << std::endl;
            std::cout << m_max_iter_mixed << std::endl;
            std::cout << m_max_iter_transport << std::endl;
            std::cout << m_max_iter_sfi << std::endl;
            std::cout << m_n_steps << std::endl;
            std::cout << m_four_approx_spaces_Q << std::endl;
            std::cout << m_mhm_mixed_Q << std::endl;
            std::cout << m_ISLinearKrModelQ << std::endl;
            std::cout << m_nThreadsMixedProblem << std::endl;
        }
        
    };
    
    /**
     * @brief Class that stores the PostProcess information
     */
    class TPostProcess : public TMRSSavable {
        
    public:
        
        /**
         * @brief Mixed operator vtk file name
         */
        std::string m_file_name_mixed;
        
        /**
         * @brief Transpor operator vtk file name
          */
        std::string m_file_name_transport;
        
        /**
         * @brief Contains scalar variables that will be postprocessed
         */
        TPZStack<std::string,10> m_scalnames;
        /**
         * @brief Contains vector variables that will be postprocessed
         */
        TPZStack<std::string,10> m_vecnames;
        
        /**
         * @brief Period of time post-processed data is printed
         */
        REAL m_file_time_step;
        
        /**
         * @brief Contains the times at which post-processed data is printed
         */
        TPZStack<REAL,100> m_vec_reporting_times;
        
        /**
         * @brief Default constructor
         */
        TPostProcess(){
            
            m_file_name_mixed       = "";
            m_file_name_transport   = "";
            m_scalnames.Resize(0);
            m_vecnames.Resize(0);
            m_file_time_step = 0.0;
            m_vec_reporting_times.Resize(0);
            
        }
        /**
         * @brief Destructor
         */
        ~TPostProcess(){
            
        }
        
        /**
         * @brief Copy constructor
         */
        TPostProcess(const TPostProcess & other){
            m_file_name_mixed       = other.m_file_name_mixed;
            m_file_name_transport   = other.m_file_name_transport;
            m_vecnames              = other.m_vecnames;
            m_scalnames             = other.m_scalnames;
            m_file_time_step        = other.m_file_time_step;
            m_vec_reporting_times   = other.m_vec_reporting_times;
            
        }
        /**
         * @brief Copy assignment operator
         */
        TPostProcess & operator=(const TPostProcess &other){
            
            // check for self-assignment
            if(&other == this){
                return *this;
            }
            
            m_file_name_mixed       = other.m_file_name_mixed;
            m_file_name_transport   = other.m_file_name_transport;
            m_vecnames              = other.m_vecnames;
            m_scalnames             = other.m_scalnames;
            m_file_time_step        = other.m_file_time_step;
            m_vec_reporting_times   = other.m_vec_reporting_times;
            
            return *this;
        }
        
        bool operator==(const TPostProcess &other){
            
            // check for self-assignment
            if(&other == this){
                return true;
            }
            
            return
            m_file_name_mixed       == other.m_file_name_mixed &&
            m_file_name_transport   == other.m_file_name_transport&&
            m_vecnames              == other.m_vecnames&&
            m_scalnames             == other.m_scalnames&&
            m_file_time_step        == other.m_file_time_step&&
            m_vec_reporting_times   == other.m_vec_reporting_times;
            
        }
        
    
        void Write(TPZStream &buf, int withclassid) const{ //ok
            buf.Write(&m_file_name_mixed);
            buf.Write(&m_file_name_transport);
            buf.Write(m_vecnames);
            buf.Write(m_scalnames);
            buf.Write(&m_file_time_step);
            buf.Write(m_vec_reporting_times);
        }
        
        void Read(TPZStream &buf, void *context){ //ok
            buf.Read(&m_file_name_mixed);
            buf.Read(&m_file_name_transport);
            buf.Read(m_vecnames);
            buf.Read(m_scalnames);
            buf.Read(&m_file_time_step);
            buf.Read(m_vec_reporting_times);
        }
        
        virtual int ClassId() const {
            return Hash("TMRSDataTransfer::TPostProcess");
        }
        
        
        
        void Print() const {
            std::cout << m_file_name_mixed << std::endl;
            std::cout << m_file_name_transport << std::endl;
            std::cout << m_file_time_step << std::endl;
            std::cout << m_vec_reporting_times << std::endl;
            //scalnames and vecnames
        }
        
    };
    
    TGeometry mTGeometry;
    TPetroPhysics mTPetroPhysics;
    TFluidProperties mTFluidProperties;
    TReservoirProperties mTReservoirProperties;
    TMultiphaseFunctions mTMultiphaseFunctions;
    TBoundaryConditions mTBoundaryConditions;
    TNumerics mTNumerics;
    TPostProcess mTPostProcess;
    
};

#endif /* TMRSDataTransfer_h */
