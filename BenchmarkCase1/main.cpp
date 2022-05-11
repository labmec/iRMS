
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
#include "imrs_config.h"
#include "pzlog.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void LearningReadFracMesh();
void FracSimpleCase();
void TransferLamdasToCondensedCompel(TPZMultiphysicsCompMesh *mixed_operator);



TPZGeoMesh *ReadFractureMesh(TPZVec<int64_t> &subdomain);
TPZGeoMesh *ReadFractureMesh();
void BenchmarkCase1();
TMRSDataTransfer SettingBenchmarkCase1();

void  ForcingFunction (const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

int main(){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    BenchmarkCase1();
}


void MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh, TPZVec<int64_t> &subdomain);

TPZGeoMesh *ReadFractureMesh(TPZVec<int64_t> &subdomain)
{
      std::string basemeshpath(FRACMESHES);
//        std::string fileCoarse=  basemeshpath + "/flem_case1_Coarse_BC.msh";
//        std::string fileFine=  basemeshpath + "/flem_case1_Submesh_Fractures.msh";
//
    
   
    std::string fileFine = basemeshpath + "/Case1_Cilamce/case1_fine5400.msh";
    std::string fileCoarse = basemeshpath + "/Case1_Cilamce/case1_coarse5400.msh";
    
//    std::string fileFine = basemeshpath + "/Case1_Cilamce/case1_fine9000.msh";
//    std::string fileCoarse = basemeshpath + "/Case1_Cilamce/case1_coarse.msh";
    
//    std::string fileFine = basemeshpath + "/Case1_Cilamce/case1_fine7000.msh";
//    std::string fileCoarse = basemeshpath + "/Case1_Cilamce/case1_coarse7000.msh";
    
//    std::string fileCoarse = basemeshpath +"/Case1_Cilamce/case1_coarse.msh";
//    std::string fileFine =basemeshpath + "/Case1_Cilamce/case1_fine5400.msh";
    

//    std::string fileFine = basemeshpath + "/jose6_fine.msh";
//    std::string fileCoarse = basemeshpath +"/jose6_coarse.msh";
////

   
//    std::string fileFine = basemeshpath + "/embedFrac_subWithFrac.msh";
//    std::string fileCoarse = basemeshpath + "/embedFrac_coarse.msh";

    
    
//    std::string fileCoarse("../../FracMeshes/flem_case1_Coarse_BC.msh");
//    std::string fileFine("../../FracMeshes/flem_case1_Fine_BC.msh");
//
    
    
//    std::string fileCoarse("../../FracMeshes/jose10_coarse.msh");
//    std::string fileFine("../../FracMeshes/jose10_fine.msh");

    
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4); // From 0D to 3D
    /*
     2 4 "inlet"
     2 5 "outlet"
     2 6 "noflux"
     3 3 "k33"
     3 10 "k31"
     */
    dim_name_and_physical_tagCoarse[3]["k33"] = 1;
    dim_name_and_physical_tagCoarse[3]["k31"] = 2;
    dim_name_and_physical_tagCoarse[2]["inlet"] = -2;
    dim_name_and_physical_tagCoarse[2]["outlet"] = -4;
    dim_name_and_physical_tagCoarse[2]["noflux"] = -1;
    /*
     2 2 "Fractures"
     3 1 "c1"
     */
    dim_name_and_physical_tagFine[2]["Fractures"] = 10;
    dim_name_and_physical_tagFine[2]["Fracture"] = 10;
    dim_name_and_physical_tagFine[2]["Fracture2"] = 10;
    dim_name_and_physical_tagFine[2]["Fracture10"] = 10;
    dim_name_and_physical_tagFine[2]["frac0"] = 10;
    dim_name_and_physical_tagFine[1]["BCfrac0"] = -11;
    
    
    dim_name_and_physical_tagFine[1]["nofluxFrac"] = -11;
    dim_name_and_physical_tagFine[1]["outletFrac"] = -12;
    dim_name_and_physical_tagFine[1]["inletFrac"] = -13;
    for(int i=1; i<=200; i++)
    {
        std::stringstream sout;
        sout << "c" << i;
        dim_name_and_physical_tagFine[3][sout.str()] = i+9;
    }
    TPZGmshReader GeometryCoarse, GeometryFine;
    TPZGeoMesh *gmeshCoarse, *gmeshFine;
    {
        REAL l = 1.0;
        GeometryCoarse.SetCharacteristiclength(l);
        GeometryCoarse.SetDimNamePhysical(dim_name_and_physical_tagCoarse);
        gmeshCoarse = GeometryCoarse.GeometricGmshMesh(fileCoarse);
        GeometryCoarse.PrintPartitionSummary(std::cout);
    }
    {
        REAL l = 1.0;
        GeometryFine.SetCharacteristiclength(l);
        GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
        gmeshFine = GeometryFine.GeometricGmshMesh(fileFine);
        GeometryCoarse.PrintPartitionSummary(std::cout);
    }
    {
        std::ofstream fileCoarse("mesh3dCoarse.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshCoarse, fileCoarse);
    }
    {
        std::ofstream fileFine("mesh3dFine.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFine);
    }
    int nels =gmeshFine->NElements();
    for (int iel =0; iel<nels; iel++) {
        TPZGeoEl *gel = gmeshFine->Element(iel);
        TPZVec<REAL> qsi(3,0.0);
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        REAL detjac;
        TPZFMatrix<REAL> jacinv;
        gel->Jacobian(qsi, jac, axes, detjac, jacinv);
        if(detjac<0.0){
            DebugStop();
        }
        
        
    }
    
    
    MergeMeshes(gmeshFine, gmeshCoarse, subdomain);
    {
        std::ofstream fileFine("mesh3dFineMerge.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFine,subdomain);
    }
    delete gmeshCoarse;
    return gmeshFine;
}

TPZGeoMesh *ReadFractureMesh(){
//    std::string fileFine("../../FracMeshes/case1_Tetra.msh");
//    std::string fileFine("../../FracMeshes/simple_tetra.msh");
    
//    std::string fileFine("../../FracMeshes/CaseCube.msh");
    
    
//    std::string fileFine("../../FracMeshes/case1_simple_matchnew.msh");
//    std::string fileFine("../../FracMeshes/flem_case1_Submesh_Fractures.msh");
//    std::string fileFine("../../FracMeshes/case_1bas.msh");
//    std::string fileFine("../../FracMeshes/case1_withoutFrac.msh");
// std::string fileFine("../../FracMeshes/case1_withFrac.msh");
    
    
    
    std::string fileFine("../../FracMeshes/jose_simple0.msh");
    
//    std::string fileFine("../../FracMeshes/case1_simple_5elem.msh");
    
    
    
//    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4); // From 0D to 3D
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
//
//    dim_name_and_physical_tagFine[1]["BCFracInlet"] = -12;
//    dim_name_and_physical_tagFine[1]["BCFracOutlet"] = -13;
//    dim_name_and_physical_tagFine[1]["FracIntersection"] = -14;
    
 
    /*
     2 2 "Fractures"
     3 1 "c1"
     */

//    for(int i=1; i<=100; i++)
//    {
//        std::stringstream sout;
//        sout << "c" << i;
//        dim_name_and_physical_tagFine[3][sout.str()] = i+9;
//    }
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;

        REAL l = 1.0;
        GeometryFine.SetCharacteristiclength(l);        

        GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
        gmeshFine = GeometryFine.GeometricGmshMesh(fileFine);
        GeometryFine.PrintPartitionSummary(std::cout);
        std::ofstream fileFinevtk("mesh3dFine.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshFine, fileFinevtk);
    
   
   
    return gmeshFine;
}

void MergeMeshes(TPZGeoMesh *finemesh, TPZGeoMesh *coarsemesh, TPZVec<int64_t> &subdomain)
{
    int fine_skeleton_matid = 18;
    int coarse_skeleton_matid = 19;
    std::map<int,int64_t> MatFinetoCoarseElIndex;
    std::map<int64_t,int64_t> NodeCoarseToNodeFine;
    std::map<int64_t,int64_t> ElCoarseToElFine;
    int temp_bc_mat = -10;
    // create boundary elements for elements without neighbour
    {
        std::map<int, int> num_created;
        int64_t nel_fine = finemesh->NElements();
        for (int64_t el = 0; el<nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            int dim = gel->Dimension();
            int nsides = gel->NSides();
            int firstside = nsides-gel->NSides(dim-1)-1;
            for (int side = firstside; side<nsides-1; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour.Element()->Dimension() != dim)
                {
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == gelside)
                {
                    TPZGeoElBC gelbc(gelside,temp_bc_mat);
                    num_created[gel->MaterialId()]++;
                }
            }
        }
#ifdef PZDEBUG
        for (auto it : num_created) {
            std::cout << "For matid " << it.first << " number of elements created " << it.second << std::endl;
        }
#endif
    }
    // find the correspondence between coarse nodes and fine nodes (MOST EXPENSIVE OPERATION)
    {
        int64_t nnode_coarse = coarsemesh->NNodes();
        for (int64_t n = 0; n<nnode_coarse; n++) {
            TPZGeoNode &no = coarsemesh->NodeVec()[n];
            if(no.Id() == -1) continue;
            TPZManVector<REAL,3> co(3);
            no.GetCoordinates(co);
            int64_t fineindex;
            TPZGeoNode *finenode = finemesh->FindNode(co,fineindex);
            NodeCoarseToNodeFine[n] = fineindex;
        }
    }
    // identify the correspondence between the material id of the fine mesh and the coarse element index
    // of the coarse mesh
    // this also defines the subdomain of the elements
    {
        int64_t first3DCoarse = 0;
        int dim = coarsemesh->Dimension();
        {
            int64_t nelcoarse = coarsemesh->NElements();
            for (int64_t el=0; el<nelcoarse; el++) {
                TPZGeoEl *gel = coarsemesh->Element(el);
                if(gel->Dimension() == dim)
                {
                    first3DCoarse = el;
                    break;
                }
            }
        }
        int64_t nel_fine = finemesh->NElements();
        subdomain.Resize(nel_fine);
        subdomain.Fill(-1);
        for (int64_t el = 0; el<nel_fine; el++) {
            auto *gel = finemesh->Element(el);
            if(gel->Dimension() != dim) continue;
            int matid = gel->MaterialId();
            subdomain[el] = matid-10+first3DCoarse;
#ifdef PZDEBUG
            if(MatFinetoCoarseElIndex.find(matid) == MatFinetoCoarseElIndex.end())
            {
                TPZManVector<REAL,3> xcenter(3);
                TPZGeoElSide gelside(gel);
                gelside.CenterX(xcenter);
                TPZManVector<REAL,3> qsi(dim,0.);
                int64_t coarse_index = 0;
                
                TPZGeoEl *gelcoarse = coarsemesh->FindElementCaju(xcenter, qsi, coarse_index, dim);
                if(coarse_index-first3DCoarse != matid -10) DebugStop();
                MatFinetoCoarseElIndex[matid] = coarse_index;
            }
#else
            MatFinetoCoarseElIndex[matid] = matid-10+first3DCoarse;
#endif
        }
#ifdef PZDEBUG
        for(auto it : MatFinetoCoarseElIndex)
        {
            std::cout << "Fine mat id " << it.first << " coarse element index " << it.second << std::endl;
        }
#endif
    }
    // modify the material id of the boundary elements of the fine mesh (EXPENSIVE OPERATION)
    {
        int64_t nel_fine = finemesh->NElements();
        int meshdim = finemesh->Dimension();
        std::map<int,int> created_by_mat;
        for (int64_t el = 0; el < nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            if(gel->MaterialId() == temp_bc_mat)
            {
                int dim = gel->Dimension();
                if(dim != meshdim-1) continue;
                TPZGeoElSide gelside(gel);
                TPZManVector<REAL,3> xcenter(3);
                gelside.CenterX(xcenter);
                int64_t elindex3D = 0;
                TPZManVector<REAL, 3> qsi3D(3,0.);
                coarsemesh->FindElementCaju(xcenter, qsi3D, elindex3D, meshdim);
                TPZGeoEl *coarsegel3D = coarsemesh->Element(elindex3D);
                int coarseside3D = coarsegel3D->WhichSide(qsi3D);
                if(coarsegel3D->SideDimension(coarseside3D) != dim) DebugStop();
                TPZGeoElSide BCSide(coarsegel3D,coarseside3D);
                TPZGeoElSide neighbour = BCSide.Neighbour();
                while(neighbour != BCSide)
                {
                    if(neighbour.Element()->Dimension() == dim) break;
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == BCSide) DebugStop();
                int bc_id = neighbour.Element()->MaterialId();
                created_by_mat[bc_id]++;
                gel->SetMaterialId(bc_id);
            }
        }
        for (auto it : created_by_mat) {
            std::cout << "For matid " << it.first << " number of elements created " << it.second << std::endl;
        }
        auto subsize = subdomain.size();
        auto finesize = finemesh->NElements();
        subdomain.Resize(finesize, -1);
    }
    // create a Skeleton element between the large elements of the coarse mesh
    std::map<std::pair<int64_t,int64_t>, int64_t> CoarseFaceEl;
    {
        int64_t nel = coarsemesh->NElements();
        int dim = coarsemesh->Dimension();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZGeoEl *gel = coarsemesh->Element(el);
            int geldim = gel->Dimension();
            if(geldim != 3) continue;
            int firstside = gel->NSides()-gel->NSides(dim-1)-1;
            for (int side = firstside; side < gel->NSides()-1; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension() == dim) break;
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == gelside) continue;
                int64_t neighindex = neighbour.Element()->Index();
                TPZGeoElBC gelbc(gelside,coarse_skeleton_matid);
                std::pair<int64_t, int64_t> leftright(el,neighindex);
                if(neighindex < el) leftright = std::pair<int64_t, int64_t>(neighindex,el);
                CoarseFaceEl[leftright] = gelbc.CreatedElement()->Index();
            }
        }
    }
    
    // duplicate the skeleton elements of the coarse mesh within the fine mesh
    {
        int64_t nelcoarse = coarsemesh->NElements();
        int meshdim = coarsemesh->Dimension();
        for (int64_t el = 0; el<nelcoarse; el++) {
            auto gel = coarsemesh->Element(el);
            int matid = gel->MaterialId();
            if(matid != coarse_skeleton_matid) continue;
            int nnode = gel->NNodes();
            TPZManVector<int64_t, 8> nodeindices(nnode);
            for(int n=0; n<nnode; n++)
            {
                int64_t node_index_coarse = gel->NodeIndex(n);
                int64_t node_index_fine = NodeCoarseToNodeFine[node_index_coarse];
                nodeindices[n] = node_index_fine;
            }
            auto eltype = gel->Type();
            int64_t fine_index;
            finemesh->CreateGeoElement(eltype, nodeindices, matid, fine_index);
            ElCoarseToElFine[el] = fine_index;
        }
    }
    // the pair represents the subdomain indices of the elements
    // the integer is the element index of the skeleton element in the fine mesh
    std::map<std::pair<int64_t,int64_t>, int64_t> FineFaceEl;
    
    for(auto it : CoarseFaceEl)
    {
        auto coarsepair = it.first;
        int64_t coarseface = it.second;
        //        if(ElCoarseToElFine.find(coarsepair.first) == ElCoarseToElFine.end()) DebugStop();
        //        if(ElCoarseToElFine.find(coarsepair.second) == ElCoarseToElFine.end()) DebugStop();
        if(ElCoarseToElFine.find(coarseface) == ElCoarseToElFine.end()) DebugStop();
        FineFaceEl[coarsepair] = ElCoarseToElFine[coarseface];
    }
    finemesh->BuildConnectivity();
    // modify the material id of the volumetric elements of the fine mesh
    {
        int64_t nel_fine = finemesh->NElements();
        int dim = finemesh->Dimension();
        for (int64_t el = 0; el < nel_fine; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            if(gel->Dimension() != dim) continue;
            int matid = gel->MaterialId();
            if(MatFinetoCoarseElIndex.find(matid) == MatFinetoCoarseElIndex.end())
            {
                continue;
            }
            int64_t coarse_index = MatFinetoCoarseElIndex[matid];
            int64_t fine_index = ElCoarseToElFine[coarse_index];
            TPZGeoEl *father = coarsemesh->Element(coarse_index);
            int fathermatid = father->MaterialId();
            gel->SetMaterialId(fathermatid);
        }
    }
    // create face elements along the small elements as sons of macroscopic faces
    {
        // identify lists of element/sides that connect two subdomains (in the fine mesh)
        std::map<std::pair<int64_t,int64_t>, std::list<int64_t>> facelist;
        {
            int64_t nel = finemesh->NElements();
            int dim = finemesh->Dimension();
            for(int64_t el = 0; el<nel; el++)
            {
                TPZGeoEl *gel = finemesh->Element(el);
                if(gel->Dimension() != dim) continue;
                int64_t domain = subdomain[el];
                if(domain == -1) continue;
                int firstside = gel->NSides()-gel->NSides(dim-1)-1;
                for (int side = firstside; side < gel->NSides()-1; side++) {
                    TPZGeoElSide gelside(gel,side);
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    while(neighbour != gelside)
                    {
                        if(neighbour.Element()->Dimension() == dim) break;
                        neighbour = neighbour.Neighbour();
                    }
                    if(neighbour == gelside) continue;
                    int64_t neighdomain = subdomain[neighbour.Element()->Index()];
                    if(neighdomain == -1) DebugStop();
                    if(neighdomain < domain)
                    {
                        TPZGeoElBC gbc(neighbour,fine_skeleton_matid);
                        std::pair<int64_t,int64_t> leftright(neighdomain,domain);
                        facelist[leftright].push_back(gbc.CreatedElement()->Index());
                    }
                }
            }
        }
        // create the refinement patterns between small element/side and skeleton elements
        // facelist : key : left/right domain
        // second : list of geometric element indexes of (dim-1) face elements
        for (auto it : facelist) {
            if(FineFaceEl.find(it.first) == FineFaceEl.end()) DebugStop();
            int64_t fine_skel = FineFaceEl[it.first];
            int nelmesh = it.second.size()+1;
            TPZVec<TPZGeoEl *> gelvec(nelmesh);
            gelvec[0] = finemesh->Element(fine_skel);
            int64_t count = 1;
            for(auto itel : it.second) gelvec[count++] = finemesh->Element(itel);
#ifdef PZDEBUG2
            REAL Area = gelvec[0]->Volume();
            REAL Sum = 0.;
            for(int i=1; i<gelvec.size(); i++) Sum += gelvec[i]->Volume();
            REAL diff = Area-Sum;
            std::cout << "Skeleton area of el " << fine_skel << " area " << Area << " sum of small " << Sum << std::endl;
#endif
            TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::GetRefPatternBasedOnRealMeshElements(gelvec);
            TPZGeoEl *gelcoarse = finemesh->Element(fine_skel);
            gelcoarse->SetRefPattern(refpat);
            for(int i=1; i<gelvec.size(); i++){
                gelvec[i]->SetFather(gelvec[0]);
                gelvec[0]->SetSubElement(i-1, gelvec[i]);
            }
        }
        auto subsize = subdomain.size();
        auto finesize = finemesh->NElements();
        subdomain.Resize(finesize, -1);
    }
    // complement the domain of the lower dimensional elements. If all volumetric neighbours share
    // the same subdomain, the element belongs to "that" domain
    {
        int64_t nel = finemesh->NElements();
        int dim = finemesh->Dimension();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            int geldim = gel->Dimension();
            int64_t domain = subdomain[el];
            if(geldim == dim && domain == -1) DebugStop();
            if(geldim < dim && domain != -1) continue;
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            std::set<int64_t> neighdomains;
            while(neighbour != gelside)
            {
                int64_t locdomain = subdomain[neighbour.Element()->Index()];
                if(locdomain != -1) neighdomains.insert(locdomain);
                neighbour = neighbour.Neighbour();
            }
            if(neighdomains.size() == 1)
            {
                subdomain[el] = *neighdomains.begin();
            }
        }
    }
    // set the boundary of the fractures to no flow **** WATCH OUT FOR THIS **** TO BE ADJUSTED
    {
        int64_t nel = finemesh->NElements();
        int dim = finemesh->Dimension();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            int matid = gel->MaterialId();
            if(matid == temp_bc_mat && gel->Dimension() != dim-2)
            {
                std::cout << "gel index " << gel->Index() << " dim " << gel->Dimension() << " matid " << matid << std::endl;
            }
            if(gel->Dimension() != dim-2) continue;
            if(matid == temp_bc_mat) matid = -11;
            gel->SetMaterialId(matid);
        }

    }
#ifdef PZDEBUG
    {
        int64_t nel = finemesh->NElements();
        std::map<int,int> numels;
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = finemesh->Element(el);
            int matid = gel->MaterialId();
            numels[matid]++;
        }
        for(auto it: numels)
        {
            std::cout << "For matid " << it.first << " number of elements " << it.second << std::endl;
        }
    }
#endif
}

TMRSDataTransfer SettingBenchmarkCase1(){
    
    TMRSDataTransfer sim_data;
    
    //    dim_name_and_physical_tagCoarse[3]["k33"] = 1;
    //    dim_name_and_physical_tagCoarse[3]["k31"] = 2;
    //    dim_name_and_physical_tagCoarse[2]["inlet"] = -2;
    //    dim_name_and_physical_tagCoarse[2]["outlet"] = -4;
    //    dim_name_and_physical_tagCoarse[2]["noflux"] = -1;
    
    sim_data.mTGeometry.mDomainNameAndMatId["k33"] = 1;
    sim_data.mTGeometry.mDomainNameAndMatId["k31"] = 2;
    sim_data.mTGeometry.mDomainFracNameAndMatId["Fractures"] = 10;
    sim_data.mTGeometry.mInterface_material_id = 100;
    sim_data.mTGeometry.mInterface_material_idFracInf = 101;
    sim_data.mTGeometry.mInterface_material_idFracSup = 102;
    sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
    sim_data.mTGeometry.mInterface_material_idFracBound = 104;
    
    int D_Type = 0;
    int N_Type = 1;
    int zero_flux=0.0;
    REAL pressure_in = 4.0 ;
    REAL pressure_out = 1.0 ;
    REAL flux_int = -1.5;
    
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[-1] = std::make_pair(N_Type, zero_flux);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[-2] = std::make_pair(D_Type,pressure_in);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[-4] = std::make_pair(D_Type,pressure_out);
	sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue[-5] = std::make_pair(N_Type,zero_flux);
	
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[-11] = std::make_pair(N_Type, zero_flux);
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[-13] = std::make_pair(D_Type, pressure_in);
	sim_data.mTBoundaryConditions.mBCFlowFracMatIdToTypeValue[-12] = std::make_pair(D_Type, pressure_out);
    
    //Transport boundary Conditions
    int bc_inlet = 0;
    int bc_outlet = 1;
    REAL sat_in = 0.01;	
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[-1] = std::make_pair(bc_outlet, 0.0);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[-2] = std::make_pair(bc_inlet, sat_in);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[-4] = std::make_pair(bc_outlet, 0.0);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[-11] = std::make_pair(bc_outlet, sat_in);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[-13] = std::make_pair(bc_inlet, sat_in);
	sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue[-12] = std::make_pair(bc_outlet, 0.0);
    
    //Fluid Properties
    sim_data.mTFluidProperties.mWaterViscosity = 0.1;
    sim_data.mTFluidProperties.mOilViscosity = 0.1;
    sim_data.mTFluidProperties.mWaterDensityRef = 1000.0;
    sim_data.mTFluidProperties.mOilDensityRef = 1000.0;
    
    // Numerical controls
    sim_data.mTNumerics.m_max_iter_mixed = 5;
    sim_data.mTNumerics.m_max_iter_transport = 1;
    sim_data.mTNumerics.m_max_iter_sfi = 1;
    //BorderElementOrder
    sim_data.mTNumerics.m_MortarBorderElementPresOrder=1;
    sim_data.mTNumerics.m_MortarBorderElementFluxOrder=1;
    
    sim_data.mTGeometry.mSkeletonDiv = 0;
    sim_data.mTNumerics.m_sfi_tol = 0.0001;
    sim_data.mTNumerics.m_res_tol_transport = 0.0001;
    sim_data.mTNumerics.m_corr_tol_transport = 0.0001;
    sim_data.mTNumerics.m_n_steps = 100;
    REAL day = 86400.0;
    sim_data.mTNumerics.m_dt      = 1.0e7;//*day;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q          = true;
    std::vector<REAL> grav(3,0.0);
    grav[1] = 0.0;//-9.8*(1.0e-6); // hor
    sim_data.mTNumerics.m_gravity = grav;
    sim_data.mTNumerics.m_ISLinearKrModelQ = true;
    sim_data.mTNumerics.m_nThreadsMixedProblem = 0;
    
    //FracProperties
    //FracAndReservoirProperties
    sim_data.mTFracProperties.m_Permeability = 1.0e-3;
    REAL kappa1=1.0e-5;
    REAL kappa2=1.0e-6;

    
    int  id1=1;
    int  id2=2;
    std::map<int, REAL>idPerm;
    idPerm[id1]= kappa1;
    idPerm[id2]= kappa2;
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
    sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
    sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
    
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    TPZStack<REAL,100> reporting_times;
    REAL time = sim_data.mTPostProcess.m_file_time_step;
    int n_reporting_times =(n_steps)/(time/dt) + 1;
    REAL r_time =0.0;
    for (int i =1; i<= n_reporting_times; i++) {
        r_time += dt*(time/dt);
        reporting_times.push_back(r_time);
    }
    sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
    return sim_data;
}


void IdentifySubDomain()
{
    
}




void BenchmarkCase1()
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
    TPZGeoMesh *gmesh = ReadFractureMesh(subdomain);
//    TPZGeoMesh *gmesh = ReadFractureMesh();
    
    TMRSApproxSpaceGenerator aspace;
    TMRSDataTransfer sim_data  = SettingBenchmarkCase1();
    sim_data.mTFracProperties.m_matid = 10;
    sim_data.mTGeometry.mSkeletonDiv =0;
    sim_data.mTGeometry.m_skeletonMatId = 19;
    sim_data.mTNumerics.m_four_approx_spaces_Q = true;
    sim_data.mTNumerics.m_mhm_mixed_Q = false;
    sim_data.mTNumerics.m_UseSubstructures_Q = false ;
    sim_data.mTNumerics.m_SpaceType = TMRSDataTransfer::TNumerics::E4SpaceMortar;
    //mSimData.mTGeometry.mDomainNameAndMatId
    aspace.SetGeometry(gmesh);
  
    aspace.SetSubdomainIndexes(subdomain);
//    std::ofstream name("fractureTest.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, name, subdomain);
   
    
    aspace.SetDataTransfer(sim_data);
    

    
  
    int order = 1;
    aspace.BuildMixedMultiPhysicsCompMesh(order);
    TPZMultiphysicsCompMesh * mixed_operator = aspace.GetMixedOperator();
    
#ifdef USING_BOOST
    boost::posix_time::ptime tsim2 = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration deltat = tsim2-tsim1;
    std::cout << "Mixed:: OverHead " << deltat << std::endl;
#endif
    
    
    std::ofstream fileinform("Infomation.txt");
    int neq = mixed_operator->NEquations();
    int nels = mixed_operator->NElements();
    fileinform<<"NeqMixed: "<<neq<<std::endl;
    fileinform<<"NeElements: "<<nels<<std::endl;
    
    bool must_opt_band_width_Q = true;
    int n_threads = 0;
    
    //This parameter should be always "true"
    bool UsePardiso_Q = true;
    
//    mixedAnal->Configure(n_threads, UsePardiso_Q, UsingPzSparse);
//    mixedAnal->SetDataTransfer(&sim_data);
//    aspace.BuildTransportMultiPhysicsCompMesh();
    aspace.BuildAuxTransportCmesh();
    TPZCompMesh * transport_operator = aspace.GetTransportOperator();
    std::ofstream file("transportyMult.vtk");
    std::ofstream file2("transportyMult.txt");
    transport_operator->Print(file2);
    TPZVTKGeoMesh::PrintCMeshVTK(transport_operator, file);
    
    TMRSSFIAnalysis * sfi_analysis = new TMRSSFIAnalysis(mixed_operator,transport_operator,must_opt_band_width_Q);
    sfi_analysis->SetDataTransferAndBuildAlgDatStruct(&sim_data);
    
    bool usingpzSparse = false;
    sfi_analysis->Configure(n_threads, UsePardiso_Q, usingpzSparse);

    //If the parameter "UsingPzSparse" is true, it uses the pz sparse matrix, otherwise it uses eigen sparse matrix
//    bool usingpzSparse = false;
    
    //The parallelism is just implemented for the "UsingPzSparse=True" case, with eigen for now is running in serial (the next task to do)
    sfi_analysis->Configure(n_threads, UsePardiso_Q, usingpzSparse);
    int n_steps = sim_data.mTNumerics.m_n_steps;
    REAL dt = sim_data.mTNumerics.m_dt;
    
    
    TPZStack<REAL,100> reporting_times;
    reporting_times = sim_data.mTPostProcess.m_vec_reporting_times;
    
    REAL sim_time = 0.0;
    int pos =0;
    REAL current_report_time = reporting_times[pos];
    int npos = reporting_times.size();
    
    sfi_analysis->m_transport_module->UpdateInitialSolutionFromCellsData();
    REAL initial_mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
    
    std::cout << "Mass report at time : " << 0.0 << std::endl;
    std::cout << "Mass integral :  " << initial_mass << std::endl;
    std::ofstream fileCilamce("IntegratedSat.txt");
    TPZFastCondensedElement::fSkipLoadSolution = false;
    bool first=true;
    for (int it = 1; it <= n_steps; it++) {
        sim_time = it*dt;
        if (sim_time >=  current_report_time) {
            TPZFastCondensedElement::fSkipLoadSolution = false;
        }
        
        sfi_analysis->m_transport_module->SetCurrentTime(dt);
        sfi_analysis->RunTimeStep();
        mixed_operator->LoadSolution(mixed_operator->Solution());
        if (sim_time >=  current_report_time) {
            std::cout << "Time step number:  " << it << std::endl;
            std::cout << "PostProcess over the reporting time:  " << sim_time << std::endl;
            if(first==true){
                sfi_analysis->PostProcessTimeStep(1);
                first=false;
            }
            if(pos >=npos-2){
                sfi_analysis->PostProcessTimeStep(2);
            }
//            sfi_analysis->PostProcessTimeStep(2);
            pos++;
            current_report_time =reporting_times[pos];
             REAL InntMassFrac=sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMassById(10);
            fileCilamce<<current_report_time/(86400*365)<<", "<<InntMassFrac<<std::endl;
           
            REAL mass = sfi_analysis->m_transport_module->fAlgebraicTransport.CalculateMass();
            std::cout << "Mass report at time : " << sim_time << std::endl;
            std::cout << "Mass integral :  " << mass << std::endl;
            TPZFastCondensedElement::fSkipLoadSolution = true;
        }
    }
    
//
//    TPZAlgebraicDataTransfer transfer;
//    transfer.SetMeshes(*mixed_operator, *transport_operator);
//    TPZAlgebraicTransport transport;
//    transfer.BuildTransportDataStructure(transport);
//
//    TMRSTransportAnalysis *anal = new TMRSTransportAnalysis(transport_operator, true);
    
    delete gmesh;
}

void TransferLamdasToCondensedCompel(TPZMultiphysicsCompMesh *mixed_operator){
   TPZCompMesh *fluxCmesh = mixed_operator->MeshVector()[0];
 
    
    int nels = mixed_operator->NElements();
    
    for (int iel = 0; iel<nels; iel++) {
        TPZCompEl *cel = mixed_operator->Element(iel);
        TPZStack<TPZFastCondensedElement *> condensedelements;
        TPZReservoirTools::FindCondensed(cel, condensedelements);
        for (auto condensed:condensedelements) {
            TPZGeoEl *gel = condensed->Reference();
            int dim = gel->Dimension();
            int side = gel->NSides()-1;
            TPZVec<REAL> coordFast(3), xifast(dim,0);
            gel->CenterPoint(side, xifast);
            gel->X(xifast, coordFast);
            double lambda = (coordFast[0])/10.0;
            condensed->SetLambda(lambda);
         
        }
        
    }
}
