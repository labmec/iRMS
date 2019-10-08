//
//  TMRSDataTransfer.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/7/19.
//

#include "TMRSDataTransfer.h"


TMRSDataTransfer::TMRSDataTransfer(){
    DebugStop();
}


TMRSDataTransfer::TMRSDataTransfer(const TMRSDataTransfer &other){
    DebugStop();
}


TMRSDataTransfer & TMRSDataTransfer::operator=(const TMRSDataTransfer &other){
    DebugStop();
}

TMRSDataTransfer::~TMRSDataTransfer(){
    DebugStop();
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
