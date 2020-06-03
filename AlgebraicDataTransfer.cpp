//
//  AlgebraicDataTransfer.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "AlgebraicDataTransfer.h"


/// Default constructor
AlgebraicDataTransfer::AlgebraicDataTransfer() : fFluxMesh(0), fTransportMesh(0)
{
    
}

/// Copy constructor
AlgebraicDataTransfer::AlgebraicDataTransfer(const AlgebraicDataTransfer & other){
}

/// Assignement constructor
const AlgebraicDataTransfer & AlgebraicDataTransfer::operator=(const AlgebraicDataTransfer & other){
    return *this;
    
}

/// Default desconstructor
 AlgebraicDataTransfer::~AlgebraicDataTransfer(){
    
}

// compute the data transfer data structures between the fluxmesh and transport class
void AlgebraicDataTransfer::BuildTransportDataStructure(AlgebraicTransport &transport)
{
    IdentifyInterfaceGeometricElements();
    Print();
}

// Identify the geometric elements corresponding to interface elements. Order them as
// a function of the number of corner nodes
void AlgebraicDataTransfer::IdentifyInterfaceGeometricElements()
{
    // look for the geometric elements corresponding to interface elements
    // order them as a function of the number of corner nodes
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    std::pair<int, int64_t> defpair(100,-1);
    int nel_gmesh= gmesh->NElements();
    TPZVec<std::pair<int,int64_t>> interfaces(nel_gmesh,defpair);
    
    int64_t neltr = fTransportMesh->NElements();
    for(int64_t el=0; el<neltr; el++)
    {
        TPZCompEl *cel = fTransportMesh->Element(el);
        TPZMultiphysicsInterfaceElement *interf = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if(!interf)
        {
            continue;
        }
        TPZGeoEl *gel = interf->Reference();
        int ncorner = gel->NCornerNodes();
        interfaces[gel->Index()] = std::pair<int,int64_t>(ncorner,gel->Index());
    }
    std::sort(&interfaces[0],(&interfaces[0]+nel_gmesh));
    std::map<int,int64_t> numinterfaces;
    std::map<int,int64_t> total_internal_interfaces;
    for(int64_t el = 0; el<nel_gmesh; el++)
    {
        if(interfaces[el].first > 4) break;
        int64_t gelindex = interfaces[el].second;
        int matid = gmesh->Element(gelindex)->MaterialId();
        total_internal_interfaces[matid]++;
        numinterfaces[interfaces[el].first]++;
    }
    
    // then the number of interfaces elements will be the size of this data structure
    for(auto it = total_internal_interfaces.begin(); it != total_internal_interfaces.end(); it++)
    {
        fInterfaceGelIndexes[it->first].Resize(it->second);
    }
    std::map<int,int64_t> count;
    int64_t prev_count = 0;
    for(auto it = numinterfaces.begin(); it != numinterfaces.end(); it++)
    {
        count[it->first] = prev_count;
        prev_count += it->second;
    }
    for(int64_t el = 0; el<nel_gmesh; el++)
    {
        if(interfaces[el].first > 4) break;
        int64_t gelindex = interfaces[el].second;
        int matid = gmesh->Element(gelindex)->MaterialId();
        TInterfaceWithVolume &intface = fInterfaceGelIndexes[matid][count[matid]];
        intface.fInterface = gelindex;
        count[matid]++;
    }

}


// print the datastructure
void AlgebraicDataTransfer::Print(std::ostream &out)
{
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    out << "Number of interface materials " << fInterfaceGelIndexes.size();
    for(auto it = fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++)
    {
        out << "Element indexes for material id " << it->first << std::endl;
        TPZVec<TInterfaceWithVolume> &gelindex = it->second;
        int64_t nel = gelindex.NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZGeoEl *gel = gmesh->Element(gelindex[el].fInterface);
            int matid = gel->MaterialId();
            int ncorner = gel->NCornerNodes();
            out << "gel index " << gelindex[el].fInterface << " ncorner " << ncorner << " matid " << matid << std::endl;
        }
    }
}
