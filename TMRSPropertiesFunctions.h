//
//  TMRSPropertiesFunctions.h
//  Buckley_Levertt
//
//  Created by Omar Durán on 6/16/20.
//

#ifndef TMRSPropertiesFunctions_h
#define TMRSPropertiesFunctions_h

#include "pzreal.h"

class TMRSPropertiesFunctions
{
    public:
    
    /// Enumerate defining the function type
    enum EFunctionType { EConstantFunction = 0, EPiecewiseFunction = 1, ECircleLevelSetFunction = 2, EUNISIMFunction = 3, ESPECase10Function = 4};
    
    
    TMRSPropertiesFunctions(){
        m_function_type_kappa = EConstantFunction;
        m_function_type_phi = EConstantFunction;
        m_function_type_s0 = EConstantFunction;
    }
    
    ~TMRSPropertiesFunctions(){
        
    }
    
    void set_function_type_kappa(EFunctionType function_type_kappa){
        m_function_type_kappa   = function_type_kappa;
    }
    
    void set_function_type_phi(EFunctionType function_type_phi){
        m_function_type_phi     = function_type_phi;
    }
    
    void set_function_type_s0(EFunctionType function_type_s0){
        m_function_type_s0      = function_type_s0;
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Kx(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kx;
                            x = pt[0];
                            y = pt[1];
                            kx = 1.0e-7;;
                            return kx;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kx,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                kx = 1.0e-10;
                            }else{
                                kx = 1.0e-7;
                            }
                            return kx;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Ky(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,ky;
                            x = pt[0];
                            y = 1.0e-5;
                            ky = y;
                            return ky +1;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,ky,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                 ky = 1.0e-10;
                             }else{
                                 ky = 1.0e-7;
                             }
                             return ky;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_Kz(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kz;
                            x = pt[0];
                            y = pt[1];
                            kz = 1.0e-7;
                            return kz ;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kz,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                 kz = 1.0e-10;
                             }else{
                                 kz = 1.0e-7;
                             }
                             return kz;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_phi(){
        
        switch (m_function_type_phi) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,phi;
                            x = pt[0];
                            y = pt[1];
                            phi = 0.1;
                            return phi;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,phi,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                phi = 0.001;
                            }else{
                                phi = 0.1;
                            }
                            return phi;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_s0(){
        
        switch (m_function_type_s0) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s;
                            x = pt[0];
                            y = pt[1];
                            s = 0.25;
                            return s;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f>0){
                                s = 0.0;
                            }else{
                                s = 1.0;
                            }
                            return s;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s;
                            if(y<5){
                                 s = 0.0;
                             }else{
                                 s = 1.0;
                             }
                             return s;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    private:
    
    EFunctionType m_function_type_kappa;
    EFunctionType m_function_type_phi;
    EFunctionType m_function_type_s0;
  
};

#endif /* TMRSPropertiesFunctions_h */
