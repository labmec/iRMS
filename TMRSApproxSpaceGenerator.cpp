//
//  TMRSApproxSpaceGenerator.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/9/19.
//

#include "TMRSApproxSpaceGenerator.h"


TMRSApproxSpaceGenerator::TMRSApproxSpaceGenerator(){
    mGeometry = nullptr;
}

TMRSApproxSpaceGenerator & TMRSApproxSpaceGenerator::operator=(const TMRSApproxSpaceGenerator &other){
    DebugStop();
}

TMRSApproxSpaceGenerator::~TMRSApproxSpaceGenerator(){
    DebugStop();
}

void TMRSApproxSpaceGenerator::Write(TPZStream &buf, int withclassid) const{
    DebugStop();
}

void TMRSApproxSpaceGenerator::Read(TPZStream &buf, void *context){
    DebugStop();
}

int TMRSApproxSpaceGenerator::ClassId() const{
    DebugStop();
}

void TMRSApproxSpaceGenerator::SetGeometry(TPZGeoMesh * geometry){
    mGeometry = geometry;
}

TPZGeoMesh * TMRSApproxSpaceGenerator::GetGeometry(){
    return mGeometry;
}


void TMRSApproxSpaceGenerator::LoadGeometry(std::string geometry_file){
    
    TPZGmshReader Geometry;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    mGeometry = Geometry.GeometricGmshMesh(geometry_file);
    Geometry.PrintPartitionSummary(std::cout);
#ifdef PZDEBUG
    if (!mGeometry)
    {
        std::cout << "The geometrical mesh was not generated." << std::endl;
        DebugStop();
    }
#endif
}

void TMRSApproxSpaceGenerator::PrintGeometry(std::string name)
{
    if (!mGeometry) {
        DebugStop();
    }
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name  << name << "_geometry" << ".txt";
    vtk_name   << name << "_geometry"  << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    mGeometry->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(mGeometry, vtkfile, true);
    
}
