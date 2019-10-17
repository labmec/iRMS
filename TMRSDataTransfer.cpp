//
//  TMRSDataTransfer.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#include "TMRSDataTransfer.h"


TMRSDataTransfer::TMRSDataTransfer() : mTGeometry(), mTPetroPhysics(), mTFluidProperties(), mTBoundaryConditions(), mTNumerics(), mTPostProcess(){
    
}


TMRSDataTransfer::TMRSDataTransfer(const TMRSDataTransfer &other){
    mTGeometry = other.mTGeometry;
    mTPetroPhysics = other.mTPetroPhysics;
    mTFluidProperties = other.mTFluidProperties;
    mTMultiphaseFunctions = other.mTMultiphaseFunctions;
    mTBoundaryConditions = other.mTBoundaryConditions;
    mTNumerics = other.mTNumerics;
    mTPostProcess = other.mTPostProcess;
}


TMRSDataTransfer & TMRSDataTransfer::operator=(const TMRSDataTransfer &other){
    
    if (this != & other) // prevent self-assignment
    {
        mTGeometry = other.mTGeometry;
        mTPetroPhysics = other.mTPetroPhysics;
        mTFluidProperties = other.mTFluidProperties;
        mTMultiphaseFunctions = other.mTMultiphaseFunctions;
        mTBoundaryConditions = other.mTBoundaryConditions;
        mTNumerics = other.mTNumerics;
        mTPostProcess = other.mTPostProcess;
    }
    return *this;
}

TMRSDataTransfer::~TMRSDataTransfer(){

}

void TMRSDataTransfer::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TMRSDataTransfer::Read(TPZStream &buf, void *context){
    DebugStop();
}

int TMRSDataTransfer::ClassId() const{
    DebugStop();
}
