//
//  AlgebraicDataTransfer.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "TPZAlgebraicDataTransfer.h"


/// Default constructor
TPZAlgebraicDataTransfer::TPZAlgebraicDataTransfer() : fFluxMesh(0), fTransportMesh(0)
{
    
}

/// Copy constructor
TPZAlgebraicDataTransfer::TPZAlgebraicDataTransfer(const TPZAlgebraicDataTransfer & other){
}

/// Assignement constructor
const TPZAlgebraicDataTransfer & TPZAlgebraicDataTransfer::operator=(const TPZAlgebraicDataTransfer & other){
    return *this;
    
}

/// Default desconstructor
 TPZAlgebraicDataTransfer::~TPZAlgebraicDataTransfer(){
    
}

// compute the data transfer data structures between the fluxmesh and transport class
void TPZAlgebraicDataTransfer::BuildTransportDataStructure(TPZAlgebraicTransport &transport)
{
    IdentifyInterfaceGeometricElements();
    Print();
}

// Identify the geometric elements corresponding to interface elements. Order them as
// a function of the number of corner nodes
void TPZAlgebraicDataTransfer::IdentifyInterfaceGeometricElements()
{
    // look for the geometric elements corresponding to interface elements
    // order them as a function of the number of corner nodes
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    std::pair<int, int64_t> defpair(100,-1);
    int nel_gmesh= gmesh->NElements();
    
    int64_t neltr = fTransportMesh->NElements();
    TPZVec<std::pair<int,int64_t>> interfaces(neltr,defpair);
    int64_t num_interfaces = 0;
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
        interfaces[num_interfaces++] = std::pair<int,int64_t>(ncorner,el);
    }
    // this vector has all interface elements
    interfaces.Resize(num_interfaces);
    std::sort(&interfaces[0],(&interfaces[0]+num_interfaces));
    
    // number of interface elements by matid and number of corner nodes
    std::map<int,std::map<int,int64_t>> numinterfaces_matid;
    // number of interface elements by matid
    std::map<int,int64_t> total_internal_interfaces_matid;
    for(int64_t el = 0; el<num_interfaces; el++)
    {
        if(interfaces[el].first > 4) DebugStop();
        int64_t celindex = interfaces[el].second;
        TPZCompEl *cel = fTransportMesh->Element(celindex);
        int64_t gelindex = cel->Reference()->Index();
        int matid = gmesh->Element(gelindex)->MaterialId();
        total_internal_interfaces_matid[matid]++;
        numinterfaces_matid[matid][interfaces[el].first]++;
    }
    
    // then the number of interfaces elements will be the size of this data structure
    for(auto it = total_internal_interfaces_matid.begin(); it != total_internal_interfaces_matid.end(); it++)
    {
        fInterfaceGelIndexes[it->first].Resize(it->second);
    }
    // compute an index within the vector of interfaces
    std::map<int,std::map<int,int64_t>> count_interfaces;
    for(auto matit=count_interfaces.begin(); matit != count_interfaces.end(); matit++)
    {
        int64_t prev_count = 0;
        for(auto it = matit->second.begin(); it != matit->second.end(); it++)
        {
            count_interfaces[matit->first][it->first] = prev_count;
            prev_count += it->second;
        }
    }
    
    
    for(int64_t el = 0; el<num_interfaces; el++)
    {
        if(interfaces[el].first > 4) break;
        int64_t celindex = interfaces[el].second;
        TPZCompEl *cel = fTransportMesh->Element(celindex);
        TPZGeoEl *gel = cel->Reference();
        int64_t gelindex = gel->Index();
        int matid = gmesh->Element(gelindex)->MaterialId();
        int ncornernodes = gel->NCornerNodes();
        TInterfaceWithVolume &intface = fInterfaceGelIndexes[matid][count_interfaces[matid][ncornernodes]];
        intface.fInterface_gelindex = gelindex;
        intface.fInterface_celindex = celindex;
        count_interfaces[matid][ncornernodes]++;
    }

}

// Identify volume information to the interface data structure (TInterfaceWithVolume)
void TPZAlgebraicDataTransfer::IdentifyVolumeGeometricElements()
{
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    int64_t nel = gmesh->NElements();
    TPZVec<int64_t> geometricvolume(nel,0);
    for(auto it = fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++)
    {
        TPZVec<TInterfaceWithVolume> &facevec = it->second;
        
    }
}


// print the datastructure
void TPZAlgebraicDataTransfer::Print(std::ostream &out)
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
            TPZGeoEl *gel = gmesh->Element(gelindex[el].fInterface_gelindex);
            int matid = gel->MaterialId();
            int ncorner = gel->NCornerNodes();
            out << "gel index " << gelindex[el].fInterface_gelindex << " ncorner " << ncorner << " matid " << matid << std::endl;
        }
    }
}
