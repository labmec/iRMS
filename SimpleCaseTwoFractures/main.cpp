
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TPZRefPatternTools.h"
#include "TPZReservoirTools.h"
#include "pzlog.h"
#include "imrs_config.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "TPZHybridizeHDiv.h"

#include "pzlog.h"

using namespace std;

void CaseSimple2Frac();
TMRSDataTransfer Setting2Fractures();
TPZGeoMesh *ReadFractureMesh();

void identifyElementSidesOnIntersection(TPZGeoMesh *gmesh);
void HybridizeIntersection(TPZMultiphysicsCompMesh*& mmesh);
TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side);

std::tuple<int64_t, int> SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid);


int main(){
    TPZLogger::InitializePZLOG();
    CaseSimple2Frac();

}

void CaseSimple2Frac()
{
    /*
     the different lagrange levels for this mesh layout
     char fluxmortar = 5;
     char firstpressurelagrange = 1;
     char pressurelagrange = 3;
     char pressuremortar = 4;
     char distfluxlagrange = 2;
     char avpressurelagrange = 6;
     
     */
    // vector with subdomain index of the geometric elements
    TPZVec<int64_t> subdomain;
    //    TPZGeoMesh *gmesh = ReadFractureMesh(subdomain);
    TPZGeoMesh *gmesh = ReadFractureMesh();
  
    identifyElementSidesOnIntersection(gmesh);
    
    
    TMRSApproxSpaceGenerator aspace;
    TMRSDataTransfer sim_data  = Setting2Fractures();
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
//    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E2Space;
    //mSimData.mTGeometry.mDomainDimNameAndPhysicalTag
    aspace.SetGeometry(gmesh);
    aspace.SetSubdomainIndexes(subdomain);
    //    std::ofstream name("fractureTest.vtk");
    //    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name, subdomain);
    
    
    aspace.SetDataTransfer(sim_data);
    
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
    
    // ------ Hybridizing fracture intersection -------
//    TPZHybridizeHDiv hybridizer(mixed_operator->MeshVector());
//    hybridizer.DimToHybridize() = 2; // fracture dimension
//    hybridizer.MatIDToHybridize() = -14; // hardcoded for now
//    hybridizer.Hybridize(mixed_operator);
    // ------ End of hybridizing fracture intersection -------
    
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    
    //This parameter should be always "true"
    bool UsingPzSparse = false;
    bool UsePardiso_Q = true;
    
    TMRSMixedAnalysis *mixAnalisys = new TMRSMixedAnalysis(mixed_operator, must_opt_band_width_Q);
    mixAnalisys->SetDataTransfer(&sim_data);
    mixAnalisys->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
    mixAnalisys->Assemble();
    mixAnalisys->Solve();
    mixed_operator->UpdatePreviousState(-1.);
    TPZFastCondensedElement::fSkipLoadSolution = false;
    mixed_operator->LoadSolution(mixed_operator->Solution());
    mixAnalisys->fsoltransfer.TransferFromMultiphysics();
    mixAnalisys->PostProcessTimeStep();
    delete gmesh;
}

void HybridizeIntersection(TPZMultiphysicsCompMesh*& mmesh){
    TPZCompMesh* fluxmesh = mmesh->MeshVector()[0];
    TPZGeoMesh* gmesh = fluxmesh->Reference();
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    int dim = gmesh->Dimension();
    
    int64_t nel = fluxmesh->NElements();
    std::list<std::tuple<int64_t, int> > pressures;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        
        // We only want 2d fracture, thus dim-1
        if (!intel || intel->Reference()->Dimension() != dim-1) {
            continue;
        }
        // loop over the side of dimension dim-2 (edges of fracture)
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side < gel->NSides() - 1; side++) {
            if (gel->SideDimension(side) != dim - 2) {
                continue;
            }
            
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neigh = gelside.Neighbour();
            while(neigh != gelside){
                int neighmatid = neigh.Element()->MaterialId();
                int neighdim = neigh.Dimension();
                
                if (neighmatid == -14 && neighdim == 1) {
                    cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " has side number " << side << " with dim = " << neigh.Dimension() << " touching the intersection" << endl;
                    cout << "===> Splitting its connects" << endl;
                    TPZCompElSide celside(intel, side);
                    TPZCompElSide neighcomp = RightElement(intel, side);
                    if (neighcomp) {
                        // SplitConnects returns the geometric element index and interpolation order
                        pressures.push_back(SplitConnects(celside, neighcomp, mmesh->MeshVector()));
                    }
                }
                neigh = neigh.Neighbour();
            } // while
            
        }
    }
    fluxmesh->InitializeBlock();
    fluxmesh->ComputeNodElCon();
    
    TPZCompMesh *pressuremesh = mmesh->MeshVector()[1];
    gmesh->ResetReference();
    pressuremesh->SetDimModel(gmesh->Dimension()-1);
    for (auto pindex : pressures) {
        int64_t elindex;
        int order;
        std::tie(elindex, order) = pindex;
        TPZGeoEl *gel = gmesh->Element(elindex);
        int64_t celindex;
        TPZCompEl *cel = pressuremesh->ApproxSpace().CreateCompEl(gel, *pressuremesh, celindex);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
        if (intel){
            intel->PRefine(order);
            //            intel->SetSideOrder(gel->NSides() - 1, order);
        } else if (intelDisc) {
            intelDisc->SetDegree(order);
            intelDisc->SetTrueUseQsiEta();
        } else {
            DebugStop();
        }
        int n_connects = cel->NConnects();
        for (int i = 0; i < n_connects; ++i) {
            cel->Connect(i).SetLagrangeMultiplier(2);
        }
        gel->ResetReference();
    }
    pressuremesh->InitializeBlock();
    pressuremesh->SetDimModel(gmesh->Dimension());
}

std::tuple<int64_t, int> SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid) {

    const int matIdDivWrap = -100, matIdLagrage = -101;
    
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    //TPZCompMesh *pressuremesh = meshvec[1];
    //TPZGeoMesh *gmesh = fluxmesh->Reference();
    TPZGeoElSide gleft(left.Reference());
    TPZGeoElSide gright(right.Reference());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
    TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *> (right.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    intelright->SetSideOrient(right.Side(), 1);
    TPZStack<TPZCompElSide> equalright;
    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());

    if (cleft.HasDependency()) {
        // check whether the wrap of the large element was already created
        gright.EqualLevelCompElementList(equalright,1,0);
        // only one wrap element should exist
#ifdef PZDEBUG
        if(equalright.size() > 1)
        {
            DebugStop();
        }
        if(equalright.size()==1)
        {
            TPZGeoEl *equalgel = equalright[0].Element()->Reference();
            if(equalgel->Dimension() != fluxmesh->Dimension()-1)
            {
                DebugStop();
            }
        }
#endif
        // reset the reference of the wrap element
        if(equalright.size()) equalright[0].Element()->Reference()->ResetReference();
        cleft.RemoveDepend();
    }
    else
    {
        int64_t index = fluxmesh->AllocateNewConnect(cleft);
        TPZConnect &newcon = fluxmesh->ConnectVec()[index];
        cleft.DecrementElConnected();
        newcon.ResetElConnected();
        newcon.IncrementElConnected();
        newcon.SetSequenceNumber(fluxmesh->NConnects() - 1);

        int rightlocindex = intelright->SideConnectLocId(0, right.Side());
        intelright->SetConnectIndex(rightlocindex, index);
    }
    int sideorder = cleft.Order();
    fluxmesh->SetDefaultOrder(sideorder);
    // create HDivBound on the sides of the elements
    TPZCompEl *wrap1, *wrap2;
    {
        intelright->Reference()->ResetReference();
        intelleft->LoadElementReference();
        intelleft->SetPreferredOrder(sideorder);
        TPZGeoElBC gbc(gleft, matIdDivWrap);
        int64_t index;
        wrap1 = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh, index);
        if(cleft.Order() != sideorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap1);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrap1->Reference()->ResetReference();
    }
    // if the wrap of the large element was not created...
    if(equalright.size() == 0)
    {
        intelleft->Reference()->ResetReference();
        intelright->LoadElementReference();
        TPZConnect &cright = intelright->SideConnect(0,right.Side());
        int rightprevorder = cright.Order();
        intelright->SetPreferredOrder(cright.Order());
        TPZGeoElBC gbc(gright, matIdDivWrap);
        int64_t index;
        wrap2 = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh, index);
        if(cright.Order() != rightprevorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap2);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelright->Reference()->ResetReference();
        wrap2->Reference()->ResetReference();
    }
    else
    {
        wrap2 = equalright[0].Element();
    }
    wrap1->LoadElementReference();
    wrap2->LoadElementReference();
    int64_t pressureindex;
    int pressureorder;
    {
        TPZGeoElBC gbc(gleft, matIdLagrage);
        pressureindex = gbc.CreatedElement()->Index();
        pressureorder = sideorder;
    }
    intelleft->LoadElementReference();
    intelright->LoadElementReference();
    return std::make_tuple(pressureindex, pressureorder);
}

TPZCompElSide RightElement(TPZInterpolatedElement *intel, int side) {
    bool isrestrained = false;
    {
        TPZConnect &c = intel->SideConnect(0, side);
        if (c.HasDependency()) {
            isrestrained = true;
        }
    }
    TPZGeoEl *gel = intel->Reference();
    TPZGeoElSide gelside(gel, side);

    if (isrestrained == true) {
        /// if the side is restrained we will hybridize between the element and the larger element
        TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
        if (!celside) DebugStop();
        TPZGeoEl *neigh = celside.Element()->Reference();
        /// we assume that a larger element should not be a boundary element
        if (neigh->Dimension() != gel->Dimension()) {
            DebugStop();
        }
        return celside;
    } else {
        // if the connect is not restrained
        // - there should be only one neighbour
        //   if there is more than one neighbour the element side has already been hybridized
        // - the neighbour should be of the same dimension
        //   if the neighbour is of lower dimension it is a boundary element
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        if (celstack.size() == 1) {
            TPZGeoEl *neigh = celstack[0].Element()->Reference();
            if (neigh->Dimension() == gel->Dimension()) {
                return celstack[0];
            }
        }
    }
    return TPZCompElSide();
}


void identifyElementSidesOnIntersection(TPZGeoMesh *gmesh) {
  
  TPZAdmChunkVector<TPZGeoEl *>::iterator it = gmesh->ElementVec().begin();
  for (; it != gmesh->ElementVec().end(); it++) {
    TPZGeoEl* gel = (*it);
    if (!gel)
      continue;
    
    int matid = gel->MaterialId();
    if (matid == -14) {
      continue; // this is the intersection
    }
    
    const int nsides = gel->NSides();
    for (int i = 0; i < nsides; i++) {
      TPZGeoElSide gelside(gel,i);
      TPZGeoElSide neigh = gelside.Neighbour();
      while(neigh!= gelside){
        int neighmatid = neigh.Element()->MaterialId();
        if (neighmatid == -14) {
          if (neigh.Dimension() == 1){
            cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " has side number " << i << " with dim = " << neigh.Dimension() << " touching the intersection" << endl;
          }
        }
        neigh = neigh.Neighbour();
      } // while
    } // for i
  } // for it
}

TMRSDataTransfer Setting2Fractures(){
    
    TMRSDataTransfer sim_data;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["k33"] = 1;
    sim_data.mTGeometry.mDomainDimNameAndPhysicalTag[3]["k31"] = 2;
    sim_data.mTGeometry.mDomainFracDimNameAndPhysicalTag[2]["Fractures"] = 10;
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mIterface_material_idFracBound = 104;
    
    int D_Type = 0;
    int N_Type = 1;
    REAL zero_flux = 0.0;
    REAL pressure_in = 1.0 ;
    REAL pressure_out = 1.0 ;
    
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue.Resize(4);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[0] = std::make_tuple(-1,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[1] = std::make_tuple(-2,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[2] = std::make_tuple(-4,D_Type,pressure_out);
    sim_data.mTBoundaryConditions.mBCMixedPhysicalTagTypeValue[3] = std::make_tuple(-5,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[0] =
    std::make_tuple(-11,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[1] =
    std::make_tuple(-12,D_Type,pressure_in);
    sim_data.mTBoundaryConditions.mBCMixedFracPhysicalTagTypeValue[2] =
    std::make_tuple(-13,D_Type,pressure_out);
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 1.0;
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue.Resize(3);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[0] = std::make_tuple(-1,bc_outlet,0.0);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[1] = std::make_tuple(-2,bc_inlet,sat_in);
    sim_data.mTBoundaryConditions.mBCTransportPhysicalTagTypeValue[2] = std::make_tuple(-4,bc_outlet,0.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.1;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 1;
    sim_data.mTNumerics.m_max_iter_transport = 1;
    sim_data.mTNumerics.m_max_iter_sfi = 1;
    //BorderElementOrder
    sim_data.mTNumerics.m_BorderElementPresOrder=1;
    sim_data.mTNumerics.m_BorderElementFluxOrder=1;
    
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 1 ;
    sim_data.mTNumerics.m_dt      = 1.0;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    //FracAndReservoirProperties
    sim_data.mTFracProperties.m_Permeability = 1.e-5;
    REAL kappa=1.0;
    int  id1=1;
    int  id2=2;
    std::vector<std::pair<int, REAL>> idPerm(2);
    idPerm[0]= std::make_pair(id1,kappa);
    idPerm[1]= std::make_pair(id2,kappa);
    sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;
    
    
    
    // PostProcess controls
    sim_data.mTPostProcess.m_file_name_mixed = "mixed_operator.vtk";
    sim_data.mTPostProcess.m_file_name_transport = "transport_operator.vtk";
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    if (sim_data.mTNumerics.m_four_approx_spaces_Q) {
        scalnames.Push("g_average");
        scalnames.Push("p_average");
    }
    sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
    sim_data.mTPostProcess.m_vecnames = vecnames;
    sim_data.mTPostProcess.m_scalnames = scalnames;
    return sim_data;
}




TPZGeoMesh *ReadFractureMesh(){
    string basemeshpath(FRACMESHES);
    std::string fileFine = basemeshpath + "/Case2Frac.msh";
//    std::string fileFine("../../FracMeshes/Case2Frac.msh");
//    std::string fileFine("../../FracMeshes/dfnExport.msh");
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
    /*
     2 4 "inlet"
     2 5 "outlet"
     2 6 "noflux"
     3 3 "k33"
     3 10 "k31"
     */
    
    //    dim_name_and_physical_tagFine[3]["c1"] = 1;
    dim_name_and_physical_tagFine[3]["k33"] = 2;
    dim_name_and_physical_tagFine[3]["k31"] = 1;
    dim_name_and_physical_tagFine[2]["inlet"] = -2;
    dim_name_and_physical_tagFine[2]["outlet"] = -4;
    dim_name_and_physical_tagFine[2]["noflux"] = -1;
    
    dim_name_and_physical_tagFine[2]["Fractures"] = 10;
    dim_name_and_physical_tagFine[2]["Fracture2"] = 10;
    dim_name_and_physical_tagFine[1]["BCfrac0"] = -11;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = -15;
    
    dim_name_and_physical_tagFine[1]["BCFracInlet"] = -12;
    dim_name_and_physical_tagFine[1]["BCFracOutlet"] = -13;
    dim_name_and_physical_tagFine[1]["FracIntersection"] = -14;
    
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    GeometryFine.SetFormatVersion("4.1");
    
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(fileFine);
    GeometryFine.PrintPartitionSummary(std::cout);
    std::ofstream fileFinevtk("mesh3dFine.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFinevtk);
    return gmeshFine;
}
