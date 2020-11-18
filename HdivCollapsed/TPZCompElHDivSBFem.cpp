/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivSBFem methods.
 */


#include "TPZCompElHDivSBFem.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"
#include "pzelchdiv.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDivSBFem"));
#endif

// Initialize with the geometry of the SBFemVolume
template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &skeleton, int64_t &index) :
fGeoElVolume(gel), fSkeleton(skeleton)
{
    TPZCompEl * celskeleton = mesh.Element(skeleton);
    if (!celskeleton)
    {
        std::cout << "TPZCompElHDivSBFem<TSHAPE>:: Cannot find computational element for skeleton\n";
        DebugStop();
    }
    TPZGeoEl * gelskeleton = celskeleton->Reference();
    if (!gelskeleton)
    {
        std::cout << "TPZCompElHDivSBFem<TSHAPE>:: Cannot find geometric element for skeleton\n";
        DebugStop();
    }
    
	TPZCompElHDivCollapsed<TSHAPE>(mesh, gelskeleton, index);
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem(TPZCompMesh &mesh, const TPZCompElHDivSBFem<TSHAPE> &copy) :
TPZCompElHDivCollapsed<TSHAPE>(mesh,copy)
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::TPZCompElHDivSBFem() : TPZCompElHDivCollapsed<TSHAPE>()
{
}

template<class TSHAPE>
TPZCompElHDivSBFem<TSHAPE>::~TPZCompElHDivSBFem()
{
}

template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi)
{
    // Compute data for the 1d functions first
    bool needssol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredData(data, qsi);
    data.fNeedsSol = needssol;

    ComputeDeformedDirections(data);

    data.ComputeFunctionDivergence();
    if (data.fNeedsSol)
    {
        // ComputeSolution(qsi, data);
    }
}

template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::ComputeDeformedDirections(TPZMaterialData &data)
{
    HDivCollapsedDirections(data);
    int cont = 0;
    int firstface = TSHAPE::NSides - TSHAPE::NFacets - 1;
    int lastface = TSHAPE::NSides - 1;

    for(int side = firstface; side < lastface; side++)
    {
        int nvec = TSHAPE::NContainedSides(side);
        for (int ivet = 0; ivet<nvec; ivet++)
        {
            for (int il = 0; il<3; il++)
            {
                data.fDeformedDirections(il,ivet+cont) *= fSideOrient[side-firstface];
            }
        }
        cont += nvec;
    }
}

// This function will compute the directions for the HDiv collapsed based on the information of neighbourhood
template<class TSHAPE>
void TPZCompElHDivSBFem<TSHAPE>::HDivCollapsedDirections(TPZMaterialData &data)
{
    // Computing the deformed directions for the 2d functions using the information of the neighbourhood
    // Inspired in TPZSBFemVolume::ComputeKMatrices

    // The Reference element will be the skeleton
    TPZGeoEl *Ref1D = this->Reference();
    int dim1 = TSHAPE::Dimension;

    // The Volume element will be the one passed as an argument by the TPZSBFemVolumeHDiv
    // The SBFemVolumeHdiv will compute the Contribute for the SBFemHdiv element
    // Before that, the material data must be computed (It will call ComputeRequiredData from TPZCompElHDivSBFem)
    int dim2 = fGeoElVolume->Dimension();

    // Find the first face side of the volumetric element
    int nsides = fGeoElVolume->NSides();
    int faceside;
    for (faceside = 0; faceside < nsides; faceside++) {
        if (fGeoElVolume->SideDimension(faceside) == dim1) {
            break;
        }
    }
    TPZGeoElSide thisside(fGeoElVolume, faceside);

    // Find the higher side of the skeleton element
    TPZGeoElSide SkeletonSide(Ref1D, Ref1D->NSides() - 1);

    // Create a transformation between these sides - Why? I didn't understand this part:
    TPZTransform<REAL> tr(dim2, dim1);
    tr = SkeletonSide.NeighbourSideTransform(thisside);
    TPZTransform<REAL> t2 = fGeoElVolume->SideToSideTransform(thisside.Side(), fGeoElVolume->NSides() - 1);
    tr = t2.Multiply(tr);

    // Applying the transformation
    TPZManVector<REAL, 3> xi(dim1), xiquad(dim2), xivol(dim2);
    xi = data.xParametric;
    tr.Apply(xi, xiquad);
    xivol = xiquad;
    xivol[dim2 - 1] = -0.5;

    // Computing the data for the volumetric element
    fGeoElVolume->Jacobian(xiquad, data.jacobian, data.axes, data.detjac, data.jacinv);
    fGeoElVolume->X(xivol, data.x);

    if (dim2 == 3) {
        // AdjustAxes3D(axes, data.axes, data.jacobian, data.jacinv, data.detjac);
    }

    //
}

#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivSBFem<TPZShapeQuad>>;


template class TPZCompElHDivSBFem<TPZShapeTriang>;
template class TPZCompElHDivSBFem<TPZShapeLinear>;
template class TPZCompElHDivSBFem<TPZShapeQuad>;