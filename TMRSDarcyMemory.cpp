//
//  TMRSDarcyMemory.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#include "TMRSDarcyMemory.h"

TMRSDarcyMemory::TMRSDarcyMemory(){
    DebugStop();
}

TMRSDarcyMemory::TMRSDarcyMemory(const TMRSDarcyMemory & other){
    DebugStop();
}

const TMRSDarcyMemory & TMRSDarcyMemory::operator=(const TMRSDarcyMemory & other){
    DebugStop();
}

TMRSDarcyMemory::~TMRSDarcyMemory(){
    DebugStop();
}

const std::string TMRSDarcyMemory::Name() const{
    DebugStop();
}

void TMRSDarcyMemory::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TMRSDarcyMemory::Read(TPZStream &buf, void *context){
    DebugStop();
}

void TMRSDarcyMemory::Print(std::ostream &out) const{
    DebugStop();
}

int TMRSDarcyMemory::ClassId() const{
    DebugStop();
}
