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
#include<tuple> // for tuple
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
           TPZManVector<std::tuple<int, TRSLinearInterpolator>> mLayer_Krw_RelPerModel;
        TPZManVector<std::tuple<int, TRSLinearInterpolator>> mLayer_Krow_RelPerModel;
        
        TPetroPhysics(){
            mLayer_Krw_RelPerModel.Resize(1);
            mLayer_Krow_RelPerModel.Resize(1);
        }
        
        ~TPetroPhysics(){
            
        }
    
        TPetroPhysics(const TPetroPhysics &other){
            mLayer_Krw_RelPerModel = other.mLayer_Krw_RelPerModel;
            mLayer_Krow_RelPerModel = other.mLayer_Krow_RelPerModel;
           
        }
        
        TPetroPhysics &operator=(const TPetroPhysics &other){
            if (this != & other) // prevent self-assignment
            {
                mLayer_Krw_RelPerModel = other.mLayer_Krw_RelPerModel;
                mLayer_Krow_RelPerModel = other.mLayer_Krow_RelPerModel;
               
            }
            return *this;
        }
    
        
    };
    
    class TFluidProperties : public TMRSSavable {
        
    public:
        
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
        
        /// time step size
        int m_report_time;
        
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
            m_report_time           = 3;
            
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
            m_report_time           = other.m_report_time;
            
            
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
            m_report_time           = other.m_report_time;
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
            m_report_time           == other.m_report_time;
            
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
            buf.Write(&m_report_time);
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
            buf.Read(&m_report_time);
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
            std::cout << m_report_time << std::endl;
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
        
        TPostProcess(){
            
            m_file_name_mixed       = "";
            m_file_name_transport   = "";
            m_scalnames.Resize(0);
            m_vecnames.Resize(0);
            
        }
        
        ~TPostProcess(){
            
        }
        
        TPostProcess(const TPostProcess & other){
            m_file_name_mixed       = other.m_file_name_mixed;
            m_file_name_transport   = other.m_file_name_transport;
            m_vecnames              = other.m_vecnames;
            m_scalnames             = other.m_scalnames;
        }
        
        TPostProcess & operator=(const TPostProcess &other){
            
            // check for self-assignment
            if(&other == this){
                return *this;
            }
            
            m_file_name_mixed       = other.m_file_name_mixed;
            m_file_name_transport   = other.m_file_name_transport;
            m_vecnames              = other.m_vecnames;
            m_scalnames              = other.m_scalnames;
            
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
            m_scalnames              == other.m_scalnames;
        }
        
        void Write(TPZStream &buf, int withclassid) const{ //ok
            buf.Write(&m_file_name_mixed);
            buf.Write(&m_file_name_transport);
            buf.Write(&m_file_name_transport);
            buf.Write(&m_vecnames);
            buf.Write(&m_scalnames);
        }
        
        void Read(TPZStream &buf, void *context){ //ok
            buf.Read(&m_file_name_mixed);
            buf.Read(&m_file_name_transport);
//            buf.Read(&m_vecnames);
//            buf.Read(&m_scalnames);
        }
        
        virtual int ClassId() const {
            return Hash("TMRSDataTransfer::TPostProcess");
        }
        
        
        
        void Print() const {
            std::cout << m_file_name_mixed << std::endl;
            std::cout << m_file_name_transport << std::endl;
            //scalnames and vecnames
        }
        
    };
    
    TGeometry mTGeometry;
    
    TPetroPhysics mTPetroPhysics;
    
    TFluidProperties mTFluidProperties;
    
    TBoundaryConditions mTBoundaryConditions;
    
    TNumerics mTNumerics;
    
    TPostProcess mTPostProcess;
    
};

#endif /* TMRSDataTransfer_h */
