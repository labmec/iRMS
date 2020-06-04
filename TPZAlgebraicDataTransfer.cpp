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
TPZAlgebraicDataTransfer::TPZAlgebraicDataTransfer(const TPZAlgebraicDataTransfer & other) : fFluxMesh(other.fFluxMesh), fTransportMesh(other.fTransportMesh), fInterfaceGelIndexes(other.fInterfaceGelIndexes), fVolumeElements(other.fVolumeElements)
{
}

/// Assignement constructor
const TPZAlgebraicDataTransfer & TPZAlgebraicDataTransfer::operator=(const TPZAlgebraicDataTransfer & other){
    fFluxMesh = other.fFluxMesh;
    fTransportMesh = other.fTransportMesh;
    fInterfaceGelIndexes = other.fInterfaceGelIndexes;
    fVolumeElements = other.fVolumeElements;
    return *this;
    
}

/// Default desconstructor
 TPZAlgebraicDataTransfer::~TPZAlgebraicDataTransfer(){
    
}

// compute the data transfer data structures between the fluxmesh and transport class
void TPZAlgebraicDataTransfer::BuildTransportDataStructure(TPZAlgebraicTransport &transport)
{
    IdentifyInterfaceGeometricElements();
    this->IdentifyVolumeGeometricElements();
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
    std::map<int,int64_t> volumecount_matid;
    TPZVec<int64_t> VolumeElementIndex(fTransportMesh->NElements(),-1);
    
    for(auto it = fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++)
    {
        TPZVec<TInterfaceWithVolume> &facevec = it->second;
        int64_t nfaces = facevec.size();
        for(int64_t iface = 0; iface<nfaces; iface++)
        {
            int64_t gelindex = facevec[iface].fInterface_gelindex;
            int64_t celindex = facevec[iface].fInterface_celindex;
            TPZCompEl *cel = fTransportMesh->Element(celindex);
            TPZMultiphysicsInterfaceElement *intface = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
            TPZCompEl* left = intface->LeftElement();
            TPZGeoEl *leftgel = left->Reference();
            int leftmatid = left->Reference()->MaterialId();
            int64_t leftindex = left->Index();
            if(VolumeElementIndex[leftindex] == -1)
            {
                VolumeElementIndex[leftindex] = volumecount_matid[leftmatid];
                volumecount_matid[leftmatid]++;
            }
            TPZCompEl *right = intface->RightElement();
            TPZGeoEl *rightgel = right->Reference();
            int rightmatid = right->Reference()->MaterialId();
            int64_t rightindex = right->Index();
            if(VolumeElementIndex[rightindex] == -1)
            {
                VolumeElementIndex[rightindex] = volumecount_matid[rightmatid];
                volumecount_matid[rightmatid]++;
            }
            facevec[iface].fLeftRightGelIndex = {leftgel->Index(),rightgel->Index()};
            facevec[iface].fLeftRightVolIndex = {VolumeElementIndex[leftindex], VolumeElementIndex[rightindex]};
        }
    }
    
    for(auto it = volumecount_matid.begin(); it != volumecount_matid.end(); it++)
    {
        fVolumeElements[it->first].Resize(it->second, -1);
    }
    
    int64_t nelcomp = fTransportMesh->NElements();
    for(int64_t el = 0; el<nelcomp; el++)
    {
        if(VolumeElementIndex[el] == -1) continue;
        TPZCompEl *cel = fTransportMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        int64_t volindex = VolumeElementIndex[el];
        fVolumeElements[matid][volindex] = el;
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
            TInterfaceWithVolume &intface = gelindex[el];
            out << "el = " << el << " gel index " << intface.fInterface_gelindex << " ncorner " << ncorner << " matid " << matid << std::endl;
            out << "         cel index " << intface.fInterface_celindex << " leftright gelindex " <<
            intface.fLeftRightGelIndex.first << " " << intface.fLeftRightGelIndex.second << std::endl;
            out << "         leftright alg vol index " << intface.fLeftRightVolIndex << std::endl;
        }
    }
    out << "Number of volume materials " << fVolumeElements.size() << std::endl;
    
    for(auto it = fVolumeElements.begin(); it != fVolumeElements.end(); it++)
    {
        out << "Volume elements for material id = " << it->first << std::endl;
        TPZVec<int64_t> &volel_indexes = it->second;
        int64_t nel = volel_indexes.size();
        for(int64_t el = 0; el<nel; el++)
        {
            int64_t celindex = volel_indexes[el];
            TPZCompEl *cel = fTransportMesh->Element(celindex);
            TPZGeoEl *gel = cel->Reference();
            int64_t gelindex = gel->Index();
            out << "el = " << el << " gel index " << gelindex << " matid " << gel->MaterialId() << " dim " << gel->Dimension() << std::endl;
        }
    }
}
