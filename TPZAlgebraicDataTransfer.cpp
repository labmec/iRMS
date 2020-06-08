//
//  AlgebraicDataTransfer.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "TPZAlgebraicDataTransfer.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzmultiphysicselement.h"

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
    BuildMixedToTransportDataStructures(fFluxMesh);
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
        if(ncorner > 4)
        {
            // an interface element should not have more than 4 corner nodes
            DebugStop();
        }
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
        if(interfaces[el].first > 4)
        {
            // an interface element should not have more than 4 corner nodes
            DebugStop();
        }
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
        if(interfaces[el].first > 4) DebugStop();
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

int SideLowerIndex(TPZGeoEl *gel, int side)
{
    int dim = gel->Dimension();
    if(dim == 1) return side;
    if(dim == 2) return side-gel->NCornerNodes();
    return side - gel->NCornerNodes() - gel->NSides(1);
}

int SideOriginalIndex(TPZGeoEl *gel, int side)
{
    int dim = gel->Dimension();
    if(dim == 1) return side;
    if(dim == 2) return side+gel->NCornerNodes();
    return side + gel->NCornerNodes() + gel->NSides(1);
}

// Identify volume information to the interface data structure (TInterfaceWithVolume)
void TPZAlgebraicDataTransfer::IdentifyVolumeGeometricElements()
{
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    int64_t nel = gmesh->NElements();
    fInterfaceByGeom.Redim(nel,6);
    for(int64_t el = 0; el<nel; el++)
    {
        for(int f=0; f<6; f++)
        {
            fInterfaceByGeom(el,f) = -1;
        }
    }
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
            TPZCompElSide leftside = intface->Left();
            TPZCompEl *left = leftside.Element();
            TPZGeoEl *leftgel = leftside.Element()->Reference();
            TPZGeoElSideIndex leftgelside(leftgel->Index(),leftside.Side());
            int lowerindex = SideLowerIndex(leftgel,leftside.Side());
            fInterfaceByGeom(leftgel->Index(),lowerindex) = iface;
            int leftmatid = left->Reference()->MaterialId();
            int64_t leftindex = left->Index();
            if(VolumeElementIndex[leftindex] == -1)
            {
                VolumeElementIndex[leftindex] = volumecount_matid[leftmatid];
                volumecount_matid[leftmatid]++;
            }
            TPZCompElSide rightside = intface->Right();
            TPZCompEl *right = rightside.Element();
            TPZGeoEl *rightgel = right->Reference();
            TPZGeoElSideIndex rightgelside(rightgel->Index(),rightside.Side());
            int rightmatid = rightside.Element()->Reference()->MaterialId();
            int64_t rightindex = right->Index();
            lowerindex = SideLowerIndex(rightgel,rightside.Side());
            fInterfaceByGeom(rightgel->Index(),lowerindex) = iface;
            if(VolumeElementIndex[rightindex] == -1)
            {
                VolumeElementIndex[rightindex] = volumecount_matid[rightmatid];
                volumecount_matid[rightmatid]++;
            }
            facevec[iface].fLeftRightGelSideIndex = {leftgelside,rightgelside};
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

// identify material of a face which is connected to a given geometric element
int TPZAlgebraicDataTransfer::IdentifyMaterial(TPZGeoElSideIndex gelindex, int64_t faceindex)
{
    for(auto it=fInterfaceGelIndexes.begin(); it != fInterfaceGelIndexes.end(); it++)
    {
        TPZVec<TInterfaceWithVolume> &faces = it->second;
        if(faceindex >= faces.size()) continue;
        TInterfaceWithVolume &intface = faces[faceindex];
        if(intface.fLeftRightGelSideIndex.first == gelindex || intface.fLeftRightGelSideIndex.second == gelindex)
        {
            return it->first;
        }
    }
    std::cout << __PRETTY_FUNCTION__ << "couldnt find a material corresponding to gelindex " << gelindex <<
            " and faceindex " << faceindex << std::endl;
    DebugStop();
    return -10000;
}

// identify material of a face which is connected to a given geometric element
int TPZAlgebraicDataTransfer::IdentifyMaterial(const TPZGeoElSide &gelside)
{
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    int64_t el = gelside.Element()->Index();
    int i = SideLowerIndex(gelside.Element(), gelside.Side());
    int64_t interface = fInterfaceByGeom(el,i);
    TPZGeoElSide refside(gelside);
    if(interface < 0)
    {
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            i = SideLowerIndex(neighbour.Element(), neighbour.Side());
            interface = fInterfaceByGeom(neighbour.Element()->Index(),i);
            if(interface >= 0)
            {
                refside = neighbour;
                el = refside.Element()->Index();
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        if(interface < 0) DebugStop();
    }
    TPZGeoElSideIndex gelsideindex(el,refside.Side());
    return IdentifyMaterial(gelsideindex, interface);
}

// find the neighbouring interface element
TPZGeoElSide TPZAlgebraicDataTransfer::IdentifyInterfaceElement(const TPZGeoElSide &gelside)
{
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    int64_t el = gelside.Element()->Index();
    int i = SideLowerIndex(gelside.Element(), gelside.Side());
    int64_t interface = fInterfaceByGeom(el,i);
    TPZGeoElSide refside(gelside);
    if(interface < 0)
    {
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            i = SideLowerIndex(neighbour.Element(), neighbour.Side());
            interface = fInterfaceByGeom(neighbour.Element()->Index(),i);
            if(interface >= 0)
            {
                refside = neighbour;
                return refside;
            }
            neighbour = neighbour.Neighbour();
        }
        if(interface < 0) DebugStop();
    }
    return gelside;
}



static void ExtractElement(TPZCompEl *cel, std::list<TPZCompEl *> &ellist)
{
    TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
    TPZGeoEl *gel = cel->Reference();
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
    if(cond)
    {
        TPZCompEl *celref = cond->ReferenceCompEl();
        ExtractElement(celref,ellist);
    }
    else if(gel)
    {
        ellist.push_back(cel);
    }
    else if(elgr)
    {
        const TPZVec<TPZCompEl *> &elvec = elgr->GetElGroup();
        int64_t nel = elvec.size();
        for(int64_t el=0; el<nel; el++)
        {
            ExtractElement(elvec[el],ellist);
        }
    }
    else
    {
        DebugStop();
    }
}
/// extract the list of computational element from a substructured computational mesh
// this method also searches for elements in element groups and condensed elements
// each computational element has an associated geometric element
void TPZAlgebraicDataTransfer::GetElementAndSubmeshPointers(TPZCompMesh &mixedmesh, std::list<TPZCompEl *> &elpointers, std::list<TPZSubCompMesh *> &submeshes)
{
    int64_t nel = mixedmesh.NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZCompEl *cel = mixedmesh.Element(el);
        if(!cel) continue;
        TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if(submesh)
        {
            submeshes.push_back(submesh);
        }
        else
        {
            ExtractElement(cel,elpointers);
        }
    }
}

/// build the data structure from mixed to transport
void TPZAlgebraicDataTransfer::BuildMixedToTransportDataStructures(TPZCompMesh *fluxmesh)
{
    if(!fluxmesh) DebugStop();
    std::list<TPZCompEl *> cellist;
    std::list<TPZSubCompMesh *> submeshes;
    GetElementAndSubmeshPointers(*fluxmesh,cellist,submeshes);
    TPZVec<int64_t> shouldtransfer(fluxmesh->NConnects(),0);
    TPZVec<int64_t> targetindex(fluxmesh->NConnects(),0);
    if(cellist.size())
    {
        // compute the number of connects that will transfer information
        std::map<int,TPZManVector<int64_t,4>> ncontransfer;
        std::map<int,TPZStack<int64_t>> connectindexes;
        // build the connect list to be transferred
        // identify the face index in the algebraic data structure
        for(auto it = cellist.begin(); it != cellist.end(); it++)
        {
            TPZCompEl *cel = *it;
            TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
            TPZGeoEl *gel = cel->Reference();
            int64_t gelindex = gel->Index();
            TPZCompEl *hdiv = mphys->ReferredElement(0);
            int nc = hdiv->NConnects();
            if(nc == 1)
            {
                // a boundary condition element
                // boundary elements are not considered
                continue;
                int64_t cindex = cel->ConnectIndex(0);
                if(shouldtransfer[cindex] != 0) continue;
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                int matid = IdentifyMaterial(gelside);
                connectindexes[matid].Push(cindex);
                if(ncontransfer.find(matid) == ncontransfer.end())
                {
                    ncontransfer[matid].Resize(4, 0);
                }
                int nshape = cel->Connect(0).NShape();
                if(nshape > 4) DebugStop();
                for(int is = 0; is<nshape; is++)
                {
                    ncontransfer[matid][is]++;
                }
                shouldtransfer[cindex] = 1;
                TPZGeoElSide refside = IdentifyInterfaceElement(gelside);
                int side = refside.Side();
                int sidelower = SideLowerIndex(refside.Element(),side);
                targetindex[cindex] = fInterfaceByGeom(refside.Element()->Index(),sidelower);
                if(targetindex[cindex] < 0) DebugStop();
            }
            else
            {
                int nsides = gel->NSides();
                int dim = gel->Dimension();
                int numfaces = gel->NSides(dim-1);
                int firstside = nsides-numfaces-1;
                if(nc != numfaces+1) DebugStop();
                for(int ic=0; ic<nc-1; ic++)
                {
                    int64_t cindex = cel->ConnectIndex(ic);
                    if(shouldtransfer[cindex] != 0) continue;
                    int side = firstside+ic;
                    TPZGeoElSide gelside(gel,side);
                    int matid = IdentifyMaterial(gelside);
                    if(ncontransfer.find(matid) == ncontransfer.end())
                    {
                        ncontransfer[matid].Resize(4, 0);
                    }

                    connectindexes[matid].Push(cindex);
                    int nshape = cel->Connect(ic).NShape();
                    for(int is=0; is<nshape; is++)
                    {
                        ncontransfer[matid][is]++;
                    }
                    shouldtransfer[cindex] = 1;
                    TPZGeoElSide refside = IdentifyInterfaceElement(gelside);
                    int64_t elindex = refside.Element()->Index();
                    int sidelower = SideLowerIndex(refside.Element(),refside.Side());
                    targetindex[cindex] = fInterfaceByGeom(elindex,sidelower);
                    if(targetindex[cindex] < 0) DebugStop();
                }
            }
        }
        for(auto it : ncontransfer)
        {
            int matid = it.first;
            int nc = connectindexes[matid].size();
            for(int is=0; is<4; is++)
            {
                if(ncontransfer[matid][is] == 0) continue;
                TFromMixedToTransport transport;
                transport.fMatid = matid;
                transport.fGather.resize(ncontransfer[matid][is]);
                transport.fScatter.resize(ncontransfer[matid][is]);
                transport.fMixedMesh = fluxmesh;
                transport.flux_sequence = is;
                int64_t count = 0;
                for(int64_t ic = 0; ic<nc; ic++)
                {
                    int64_t cindex = connectindexes[matid][ic];
                    int64_t seqnum = fluxmesh->ConnectVec()[cindex].SequenceNumber();
                    int64_t position = fluxmesh->Block().Position(seqnum);
                    int64_t blsize = fluxmesh->Block().Size(seqnum);
                    if(is>= blsize) continue;
                    if(count >= ncontransfer[matid][is]) DebugStop();
                    transport.fGather[count] = position+is;
                    transport.fScatter[count] = targetindex[cindex];
                    count++;
                }
                if(count != ncontransfer[matid][is]) DebugStop();
                fTransferMixedToTransport[matid].push_back(transport);
            }
        }
    }
    if(submeshes.size())
    {
        for(auto it = submeshes.begin(); it != submeshes.end(); it++)
        {
            BuildMixedToTransportDataStructures(*it);
        }
    }
}


// print the datastructure
void TPZAlgebraicDataTransfer::Print(std::ostream &out)
{
    TPZGeoMesh *gmesh = fTransportMesh->Reference();
    out << "Number of interface materials " << fInterfaceGelIndexes.size() << std::endl;
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
            intface.fLeftRightGelSideIndex.first << " " << intface.fLeftRightGelSideIndex.second << std::endl;
            out << "         leftright alg vol index " << intface.fLeftRightVolIndex << std::endl;
        }
    }
    out << "fVolumeElements Number of volume materials " << fVolumeElements.size() << std::endl;
    
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
    
    out << "fInterfaceByGeom For each geometric element, which are the algebraic faces connected to it\n";
    int64_t nel_geo = gmesh->NElements();
    for(int64_t el = 0; el<nel_geo; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        bool hasface = false;
        for(int i=0; i<6; i++) if(fInterfaceByGeom(el,i) != -1)
        {
            hasface=true;
        }
        if(hasface)
        {
            out << "gel index = " << el << std::endl;
            for(int i=0; i<6; i++) if(fInterfaceByGeom(el,i) != -1)
            {
                int side = SideOriginalIndex(gel,i);
                TPZGeoElSideIndex gelside(el,side);
                out << "i = " << i <<  " side " << side << " face index " << fInterfaceByGeom(el,i) << " matid "
                << IdentifyMaterial(gelside,fInterfaceByGeom(el,i)) << std::endl;
            }
        }
    }
    
    if(fTransferMixedToTransport.size())
    {
        out << "Gather scatter to bring the flux data to the transport mesh\n";
        for(auto it : fTransferMixedToTransport)
        {
            int matid = it.first;
            out << "Flux data transfer for material id : " << matid << std::endl;
            std::list<TFromMixedToTransport> &list = it.second;
            for(auto itlist : list)
            {
                TFromMixedToTransport &transport = itlist;
                transport.Print(out);
            }
        }
    }
}


TPZAlgebraicDataTransfer::TFromMixedToTransport::TFromMixedToTransport() : fMatid(-1),
    flux_sequence(-1), fFrom(0), fTarget(0), fMixedMesh(0)
{
    
}

TPZAlgebraicDataTransfer::TFromMixedToTransport::TFromMixedToTransport(const TFromMixedToTransport &copy) : fMatid(copy.fMatid), flux_sequence(copy.flux_sequence),
    fGather(copy.fGather), fScatter(copy.fScatter), fFrom(copy.fFrom),
    fTarget(copy.fTarget), fMixedMesh(copy.fMixedMesh)
{
    fMixedMesh = copy.fMixedMesh;
}

TPZAlgebraicDataTransfer::TFromMixedToTransport &TPZAlgebraicDataTransfer::TFromMixedToTransport::operator=(const TFromMixedToTransport &copy)
{
    fMatid = copy.fMatid;
    flux_sequence = copy.flux_sequence;
    fFrom = copy.fFrom;
    fTarget = copy.fTarget;
    fMixedMesh = copy.fMixedMesh;
    return *this;
}



void TPZAlgebraicDataTransfer::TFromMixedToTransport::Print(std::ostream &out)
{
    out << "FromMixedToTransport mesh = " << (void *) fMixedMesh << " flux index " << flux_sequence << " matid " << fMatid << std::endl;
    out << "Gather vector ";
    for (auto const& value : fGather) out << value << ' ';
    out << std::endl;
    out << "Scatter vector ";
    for (auto const& value : fScatter) out << value << ' ';
    out << std::endl;
}
