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
    enum EFunctionType { EConstantFunction = 0, EPiecewiseFunction = 1, EUNISIMFunction = 2, ESPECase10Function = 3};
    
    
    TMRSPropertiesFunctions(){
        m_function_type = EConstantFunction;
    }
    
    ~TMRSPropertiesFunctions(){
        
    }
    
    void set_function_type(EFunctionType function_type){
        m_function_type = function_type;
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Kx(){
        
        switch (m_function_type) {
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
        
        switch (m_function_type) {
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
        
        switch (m_function_type) {
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
        
        switch (m_function_type) {
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
            case EPiecewiseFunction:
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
    
    private:
    
    EFunctionType m_function_type;
  
};

#endif /* TMRSPropertiesFunctions_h */
