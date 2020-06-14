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
    
    
    class TGeometry : public TMRSSavable {
        
    public:
        
        std::vector<std::map<std::string,int>> mDomainDimNameAndPhysicalTag;
        
        std::vector<std::map<std::string,int>> mDomainFracDimNameAndPhysicalTag;

        
        TGeometry(){
            
            mDomainDimNameAndPhysicalTag.resize(4);
            mDomainFracDimNameAndPhysicalTag.resize(3);
        }
        
        ~TGeometry(){
            
        }
        
    
        TGeometry(const TGeometry &other){
            mDomainDimNameAndPhysicalTag = other.mDomainDimNameAndPhysicalTag;
            mDomainFracDimNameAndPhysicalTag = other.mDomainFracDimNameAndPhysicalTag;
        }
        
        TGeometry &operator=(const TGeometry &other){
            if (this != & other) // prevent self-assignment
            {
                mDomainDimNameAndPhysicalTag = other.mDomainDimNameAndPhysicalTag;
                mDomainFracDimNameAndPhysicalTag = other.mDomainFracDimNameAndPhysicalTag;
            }
            return *this;
        }
        
    };
    
    class TPetroPhysics : public TMRSSavable {
        
    public:
        
        std::vector<TRSLinearInterpolator > mLayer_Krw_RelPerModel;
        
        std::vector<TRSLinearInterpolator > mLayer_Kro_RelPerModel;
        
        
        TPetroPhysics(){
            mLayer_Krw_RelPerModel.clear();
            mLayer_Kro_RelPerModel.clear();
        }
        
        ~TPetroPhysics(){
            
        }
    
        TPetroPhysics(const TPetroPhysics &other){
            mLayer_Krw_RelPerModel = other.mLayer_Krw_RelPerModel;
            mLayer_Kro_RelPerModel = other.mLayer_Kro_RelPerModel;
           
        }
        
        TPetroPhysics &operator=(const TPetroPhysics &other){
            if (this != & other) // prevent self-assignment
            {
                mLayer_Krw_RelPerModel = other.mLayer_Krw_RelPerModel;
                mLayer_Kro_RelPerModel = other.mLayer_Kro_RelPerModel;
               
            }
            return *this;
        }
    
        
    };
    
    class TFluidProperties : public TMRSSavable {
        
    public:
        
    };
    
    class TMultiphaseFunctions : public TMRSSavable {
        
    public:
        
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_fw;
        
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_fo;
        
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_lambda;
        
        std::map<int, std::function<std::tuple<double, double, double> (TRSLinearInterpolator &, TRSLinearInterpolator &, double, double)> > mLayer_Glambda;
        
        
        TMultiphaseFunctions(){
            mLayer_fw.clear();
            mLayer_fo.clear();
            mLayer_lambda.clear();
            mLayer_Glambda.clear();
        }
        
        ~TMultiphaseFunctions(){
            
        }
        
        TMultiphaseFunctions(const TMultiphaseFunctions &other){
            mLayer_fw = other.mLayer_fw;
            mLayer_fo = other.mLayer_fo;
            mLayer_lambda = other.mLayer_lambda;
            mLayer_Glambda = other.mLayer_Glambda;
            
        }
        
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
    
    class TBoundaryConditions : public TMRSSavable {
        
    public:
        
        TPZManVector<std::tuple<int, int, REAL>> mBCMixedPhysicalTagTypeValue;
        
        TPZManVector<std::tuple<int, int, REAL>> mBCMixedFracPhysicalTagTypeValue;
        
        TPZManVector<std::tuple<int, int, REAL>> mBCTransportPhysicalTagTypeValue;
        
        TPZManVector<std::tuple<int, int, REAL>> mBCTransportFracPhysicalTagTypeValue;
        
        TBoundaryConditions(){
            
            mBCMixedPhysicalTagTypeValue.Resize(0);
            
            mBCMixedFracPhysicalTagTypeValue.Resize(0);
            
            mBCTransportPhysicalTagTypeValue.Resize(0);
            
            mBCTransportFracPhysicalTagTypeValue.Resize(0);
        }
        
        ~TBoundaryConditions(){
            
        }
        
        TBoundaryConditions(const TBoundaryConditions &other){
            mBCMixedPhysicalTagTypeValue = other.mBCMixedPhysicalTagTypeValue;
            mBCMixedFracPhysicalTagTypeValue = other.mBCMixedFracPhysicalTagTypeValue;
            mBCTransportPhysicalTagTypeValue = other.mBCTransportPhysicalTagTypeValue;
            mBCTransportFracPhysicalTagTypeValue = other.mBCTransportFracPhysicalTagTypeValue;
        }
        
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
    
    class TNumerics : public TMRSSavable {
        
    public:
        
        /// time step size
        REAL m_dt;
     
        
        /// Residual tolerance for mixed operator
        REAL m_res_tol_mixed;
        
        /// Residual tolerance for transport operator
        REAL m_res_tol_transport;
        
        /// Correction tolerance for mixed operator
        REAL m_corr_tol_mixed;
        
        /// Correction tolerance for transport operator
        REAL m_corr_tol_transport;
        
        /// Maximum number of iterations per time step for mixed operator
        int m_max_iter_mixed;
        
        /// Maximum number of iterations per time step for transport operator
        int m_max_iter_transport;
        
        /// Maximum number of Sequential Fully Implicit (SFI) iterations per time step
        int m_max_iter_sfi;
        
        /// Number of time steps
        int m_n_steps;
        
        /// Directive for the use of four spaces
        bool m_four_approx_spaces_Q;
        
        /// Directive MHM mixed approximation
        bool m_mhm_mixed_Q;
        
        
        TNumerics(){
            
            m_dt                    = 0.0;
            m_res_tol_mixed         = 1.0e-7;
            m_res_tol_transport     = 1.0e-7;
            m_corr_tol_mixed        = 1.0e-7;
            m_corr_tol_transport    = 1.0e-7;
            m_max_iter_mixed        = 0;
            m_max_iter_transport    = 0;
            m_max_iter_sfi          = 0;
            m_n_steps               = 0;
            m_four_approx_spaces_Q  = false;
            m_mhm_mixed_Q           = false;
            
        }
        
        ~TNumerics(){
            
        }
        
        TNumerics(const TNumerics & other){
            
            m_dt                    = other.m_dt;
            m_res_tol_mixed         = other.m_res_tol_mixed;
            m_res_tol_transport     = other.m_res_tol_transport;
            m_corr_tol_mixed        = other.m_corr_tol_mixed;
            m_corr_tol_transport    = other.m_corr_tol_transport;
            m_max_iter_mixed        = other.m_max_iter_mixed;
            m_max_iter_transport    = other.m_max_iter_transport;
            m_max_iter_sfi          = other.m_max_iter_sfi;
            m_n_steps               = other.m_n_steps;
            m_four_approx_spaces_Q  = other.m_four_approx_spaces_Q;
            m_mhm_mixed_Q           = other.m_mhm_mixed_Q;
            
            
        }
        
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
            m_n_steps               = other.m_n_steps;
            m_four_approx_spaces_Q  = other.m_four_approx_spaces_Q;
            m_mhm_mixed_Q           = other.m_mhm_mixed_Q;
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
            m_n_steps               == other.m_n_steps &&
            m_four_approx_spaces_Q  == other.m_four_approx_spaces_Q &&
            m_mhm_mixed_Q           == other.m_mhm_mixed_Q;
            
        }
        
        void Write(TPZStream &buf, int withclassid) const{ //ok
            buf.Write(&m_dt);
            buf.Write(&m_res_tol_mixed);
            buf.Write(&m_res_tol_transport);
            buf.Write(&m_corr_tol_mixed);
            buf.Write(&m_max_iter_mixed);
            buf.Write(&m_max_iter_transport);
            buf.Write(&m_max_iter_sfi);
            buf.Write(&m_n_steps);
            int temp = m_four_approx_spaces_Q;
            buf.Write(&temp);
            temp = m_mhm_mixed_Q;
            buf.Write(&temp);
          
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
            buf.Read(&m_n_steps);
            int temp;
            buf.Read(&temp);
            m_four_approx_spaces_Q = temp;
            buf.Read(&temp);
            m_mhm_mixed_Q = temp;
           
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
        }
        
    };
    
    class TPostProcess : public TMRSSavable {
        
    public:
        
        /// Mixed operator vtk file name
        std::string m_file_name_mixed;
        
        /// Transpor operator vtk file name
        std::string m_file_name_transport;
        
        TPZStack<std::string,10> m_scalnames;
        TPZStack<std::string,10> m_vecnames;
        REAL m_file_time_step;
        TPZStack<REAL,100> m_vec_reporting_times;
        
        TPostProcess(){
            
            m_file_name_mixed       = "";
            m_file_name_transport   = "";
            m_scalnames.Resize(0);
            m_vecnames.Resize(0);
            m_file_time_step = 0.0;
            m_vec_reporting_times.Resize(0);
            
        }
        
        ~TPostProcess(){
            
        }
        
        TPostProcess(const TPostProcess & other){
            m_file_name_mixed       = other.m_file_name_mixed;
            m_file_name_transport   = other.m_file_name_transport;
            m_vecnames              = other.m_vecnames;
            m_scalnames             = other.m_scalnames;
            m_file_time_step        = other.m_file_time_step;
            m_vec_reporting_times   = other.m_vec_reporting_times;
            
        }
        
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
    
    TMultiphaseFunctions mTMultiphaseFunctions;
    
    TBoundaryConditions mTBoundaryConditions;
    
    TNumerics mTNumerics;
    
    TPostProcess mTPostProcess;
    
};

#endif /* TMRSDataTransfer_h */
