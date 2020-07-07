
//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"
#include <TPZRefPattern.h>
#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"
#include "TPZGenGrid2D.h"
#include <time.h>
#include <stdio.h>

#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pznoderep.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"

#include "TMRSApproxSpaceGenerator.h"
#include "TMRSDataTransfer.h"
#include "TMRSSFIAnalysis.h"
#include "TMRSMixedAnalysis.h"
#include "TMRSTransportAnalysis.h"
#include "TPZExtendGridDimension.h"
#include "TPZFastCondensedElement.h"
#include "TPZReservoirTools.h"
#include "pzcondensedcompel.h"

#include "TPZAlgebraicDataTransfer.h"

#include "TMRSPropertiesFunctions.h"
#include "TRMSpatialPropertiesMap.h"
#include <libInterpolate/Interpolate.hpp>
#include <libInterpolate/AnyInterpolator.hpp>

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "TPZFileStream.h"
#include "TPZBFileStream.h"


TPZGeoMesh * CreateGeoMeshWithTopeAndBase(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ);
TPZGeoMesh * CreateGeoMeshMHM3DTest(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ);

void ModifyTopeAndBase(TPZGeoMesh * gmesh, std::string filename);
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers);
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);

void SimpleTest3D();
void UNISIMGeoMesh();



void PostProcessResProps(TPZMultiphysicsCompMesh *cmesh, TPZAlgebraicTransport *alg);

TMRSDataTransfer SettingSimple2D();
void SimpleTest2D();
//
int main(){
    InitializePZLOG();
    UNISIMGeoMesh();
    return 0;
}
void SimpleTest2D(){
    TMRSApproxSpaceGenerator aspace;
    aspace.CreateUniformMesh(10, 10, 10, 10);
    std::string name = "2D_geo";
    aspace.PrintGeometry(name);
    aspace.ApplyUniformRefinement(1);
    std::string name_ref = "2D_ref_geo";
    aspace.PrintGeometry(name_ref);

}


void UNISIMGeoMesh(){
    
    TMRSApproxSpaceGenerator aspace;
    std::string geometry_file2D ="gmsh/UNISIMT4R8P2p5.msh";
    int nLayers = 5;
    bool is3DQ = true;
    bool print3DMesh = true;
    gRefDBase.InitializeAllUniformRefPatterns();
    TPZGeoMesh *gmesh = CreateGeoMeshWithTopeAndBase( geometry_file2D,  nLayers, print3DMesh, is3DQ);
    aspace.SetGeometry(gmesh);
    aspace.ApplyUniformRefinement(0);
    std::ofstream file("geo.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    TPZFileStream buf;
    TPZBFileStream buf2;
    std::string fileName("test.txt");
    buf.OpenWrite(fileName);
    gmesh->Write(buf, 1);
    
//    TPZGeoMesh *newmesh = new TPZGeoMesh;
//    newmesh->CleanUp();
//    std::string fileName2("test.txt");
//    buf2.OpenRead(fileName2);
//    buf2.OpenWrite(fileName2);
//    void *val =0;
//    gmesh->Read(buf2,0);
//    std::ofstream file("testprint.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(newmesh, file);
 
    
    
}

TPZGeoMesh * CreateGeoMeshWithTopeAndBase(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["RockMatrix"] = 1;
    dim_name_and_physical_tag[2]["RockMatrix2"] = 1;
    dim_name_and_physical_tag[1]["Injectors"] = -2;
    dim_name_and_physical_tag[1]["Productors"] = -4;
    dim_name_and_physical_tag[1]["ZeroFlux"] = -1;
    
    TPZGmshReader Geometry;
    TPZGeoMesh * gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);
    double w = 200.0;
    
    std::string name2D("mesh2d.vtk");
    
    
    int topID= -3;
    int baseID = -1;
    TPZGeoMesh * returnedMesh = nullptr;
    if (Is3DQ) {
        TPZExtendGridDimension extend(gmesh2d, w);
        extend.SetElType(1);
        returnedMesh = extend.ExtendedMesh(nLayers,topID,baseID);
        ModifyTopeAndBase2(returnedMesh ,nLayers);
        
    }
    if (!Is3DQ){
        std::string filename1 = "Reservoir/tope_unisim2.txt";
        ModifyTopeAndBase(gmesh2d,filename1 );
        returnedMesh =gmesh2d;
        
    }
    return returnedMesh;
    
}

TPZGeoMesh * CreateGeoMeshMHM3DTest(std::string geometry_file2D, int nLayers, bool print3DMesh, bool Is3DQ){
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[2]["Reservoir"] = 1;
    dim_name_and_physical_tag[2]["Wellbore_p"] = 1;
    dim_name_and_physical_tag[2]["Wellbore_i"] = 1;
    dim_name_and_physical_tag[1]["producers"] = -4;
    dim_name_and_physical_tag[1]["injectors"] = -2;
    dim_name_and_physical_tag[1]["Reservoir_south"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_West"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_north"] = -1;
    dim_name_and_physical_tag[1]["Reservoir_East"] = -1;
    
    TPZGmshReader Geometry;
    TPZGeoMesh * gmesh2d;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh2d = Geometry.GeometricGmshMesh(geometry_file2D);
    Geometry.PrintPartitionSummary(std::cout);
    double w = 25.0;
    
    std::string name2D("mesh2d.vtk");
    
    int topID= -3;
    int baseID = -1;
    TPZGeoMesh * returnedMesh = nullptr;
    TPZExtendGridDimension extend(gmesh2d, w);
    extend.SetElType(1);
    returnedMesh = extend.ExtendedMesh(nLayers,topID,baseID);
    return returnedMesh;
    
}
void ModifyTopeAndBase(TPZGeoMesh * gmesh, std::string filename){
    
    
    std::vector<double> x, y, z;
    ReadData(filename, true, x, y, z);
    
    _2D::ThinPlateSplineInterpolator <double> interp;
    interp.setData(x,y,z);
    
    int nCoordinates = gmesh->NodeVec().NElements();
    double sum=0.0;
    for (auto val:z) {
        sum += val;
    }
    double val_storage= sum / z.size();
    
    for (int icoord=0; icoord< nCoordinates; icoord++) {
        TPZGeoNode node = gmesh->NodeVec()[icoord];
        TPZVec<REAL> co(3);
        node.GetCoordinates(co);
        double val_interp =interp(co[0],co[1]);
        
        if (val_interp==0.0) {
            co[2] =val_storage;
        }
        if (val_interp>1.0) {
            co[2] =val_interp;
        }
        
        gmesh->NodeVec()[icoord].SetCoord(co);
    }
    
    
}
void ModifyTopeAndBase2(TPZGeoMesh * gmesh ,int nlayers){
//    std::string filename2 = "Reservoir/base_unisimMOD.txt";
    
     std::string filename1 = "Reservoir/topeMOD.txt";
    std::string filename2 = "Reservoir/baseMOD.txt";
    std::vector<double> x, y, z, x1,y1,z1;
    ReadData(filename1, true, x, y, z);
    ReadData(filename2, true, x1, y1, z1);

    _2D::ThinPlateSplineInterpolator <double> interpTope;
    _2D::ThinPlateSplineInterpolator <double> interpBase;

    interpTope.setData(x,y,z);
    interpBase.setData(x1,y1,z1);

    int nCoordinates = gmesh->NodeVec().NElements();
    double sum=0.0;
    for (auto val:z) {
        sum += val;
    }
    double val_tope= sum / z.size();
    sum=0.0;
    for (auto val:z1) {
        sum += val;
    }
    double val_base= sum / z1.size();
//    val_base = 1000;
//    val_tope = 5000;
//
//    val_tope = 3000;
//    val_base = 3000;
    int npointsPerLayer = nCoordinates/(nlayers+1);
    double valinter=0.0;
    for (int ilay = 1; ilay <= nlayers+1; ilay++) {
        for (int ipoint = (ilay-1)*npointsPerLayer; ipoint<(ilay)*npointsPerLayer; ipoint++) {
            TPZGeoNode node = gmesh->NodeVec()[ipoint];
            TPZVec<REAL> co(3);
            node.GetCoordinates(co);
            double topeinterpol =interpTope(co[0],co[1]);
            double baseinterpol = interpBase(co[0],co[1]);
            if (topeinterpol==0) {
                topeinterpol = val_tope;
                if (co[0]>1000.00) {
                    topeinterpol -= 120;
                }

//                if (ipoint==0) {
//                    bool find = 0;
//                    int i =1;
//                    while (find==0) {
//                        TPZGeoNode node = gmesh->NodeVec()[ipoint+i];
//                        TPZVec<REAL> coaux(3);
//                        node.GetCoordinates(coaux);
//                        topeinterpol =interpTope(coaux[0],coaux[1]);
//                        if (topeinterpol!=0) {
//                            find=1;
//                        }
//                        i=i+1;
////                        topeinterpol =coaux[2];
//                    }
//                }
//                if (ipoint != 0) {
//                    TPZGeoNode node = gmesh->NodeVec()[ipoint-1];
//                    TPZVec<REAL> coaux(3);
//                    node.GetCoordinates(coaux);
//                    REAL valz = coaux[2];
//                    topeinterpol =coaux[2];
//
//                }
                
                
            }
            if (baseinterpol==0) {
                
                baseinterpol = val_base;
                if (co[0]>1000.00) {
                   baseinterpol = val_base-80;
                }
//                std::cout<<"{"<<co[0]<<","<<co[1]<<"};"<<std::endl;
//                    if (ipoint==npointsPerLayer) {
//                        bool find = 0;
//                        int i =1;
//                        while (find==0) {
//                            TPZGeoNode node = gmesh->NodeVec()[ipoint+i];
//                            TPZVec<REAL> coaux(3);
//                            node.GetCoordinates(coaux);
//                            baseinterpol =interpBase(coaux[0],coaux[1]);
//                            if (baseinterpol!=0) {
//                                find=1;
//                            }
//                            i=i+1;
//    //                        topeinterpol =coaux[2];
//                        }
//                    }
//                    if (ipoint != npointsPerLayer && ilay!=1) {
//                        TPZGeoNode node = gmesh->NodeVec()[ipoint-1];
//                        TPZVec<REAL> coaux(3);
//                        node.GetCoordinates(coaux);
//                        REAL valz = coaux[2];
//                        baseinterpol =coaux[2];
//
//                    }


            }

            if (ilay==1) {
                valinter=topeinterpol;
//                valinter = 3500;
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
            if (ilay==nlayers+1) {
                valinter = baseinterpol;
//                valinter = 2850;
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
            if (ilay>1   && ilay < nlayers+1) {
                valinter = topeinterpol + (ilay-1)*(baseinterpol - topeinterpol)/(nlayers);
                co[2]=valinter;
                gmesh->NodeVec()[ipoint].SetCoord(co);
            }
//            if(co[2]==0.0){
//                DebugStop();
//            }
//

        }
    }
}
void ReadData(std::string name, bool print_table_Q, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z){
    
    bool modpoints = true;
    std::ifstream file;
    file.open(name);
    int i=1;
    
    
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::istringstream issText(line);
        char l = line[0];
        if(l != '/'){
            i=i+1;
            int val = i%15;
            if(val ==0){
                double a, b, c;
                if(iss >> a >> b >> c) ;
                if (modpoints) {
                    x.push_back(a - 350808.47);
                    y.push_back(b - 7.51376238e6);
                    z.push_back(c);
                }
                else{
                x.push_back(a);
                y.push_back(b);
                z.push_back(c);
                }
            };
        };
    };
    
    if(x.size() == 0){
        std::cout<<"No data read."<<std::endl;
        
        DebugStop();
    }
    if(print_table_Q){
        std::cout<<"*************************"<<std::endl;
        std::cout<<"Reading file... ok!"<<std::endl;
        std::cout<<"*************************"<<std::endl;
        std::cout<<x.size()<<std::endl;
        std::cout<<y.size()<<std::endl;
        std::cout<<z.size()<<std::endl;
    }
    
}

void PostProcessResProps(TPZMultiphysicsCompMesh *cmesh, TPZAlgebraicTransport *alTransport){
    
    
     TPZAnalysis *an = new TPZAnalysis(cmesh,false);
     int ncells = alTransport->fCellsData.fVolume.size();
     TPZFMatrix<STATE> phi(ncells,1,0.1);
     TPZFMatrix<STATE> Kx(ncells,1,1.0);
     TPZFMatrix<STATE> Ky(ncells,1,2.0);
     TPZFMatrix<STATE> Kz(ncells,1,3.0);
    
    for (int ivol = 0; ivol<ncells; ivol++) {
        int eq = alTransport->fCellsData.fEqNumber[ivol];
        REAL phival = alTransport->fCellsData.fporosity[ivol];
        REAL Kxval = alTransport->fCellsData.fKx[ivol];
        REAL Kyval = alTransport->fCellsData.fKy[ivol];
        REAL Kzval = alTransport->fCellsData.fKz[ivol];
        phi(ivol,0) =phival;
        Kx(ivol,0) =Kxval;
        Ky(ivol,0) =Kyval;
        Kz(ivol,0) =Kzval;
    }
//    cmesh->MeshVector()[0]->Solution().Resize(ncells, 1);
//    cmesh->MeshVector()[1]->Solution().Resize(ncells, 1);
//    cmesh->MeshVector()[2]->Solution().Resize(ncells, 1);
//    cmesh->MeshVector()[3]->Solution().Resize(ncells, 1);
     cmesh->MeshVector()[0]->Solution() = phi;
     cmesh->MeshVector()[1]->Solution() = Kx;
     cmesh->MeshVector()[2]->Solution() = Ky;
     cmesh->MeshVector()[3]->Solution() = Kz;
     int dim = cmesh->Reference()->Dimension();
     TPZStack<std::string,10> scalnames, vecnames;
//     vecnames.Push("q");
     scalnames.Push("Porosity");
     scalnames.Push("Permeability_x");
     scalnames.Push("Permeability_y");
     scalnames.Push("Permeability_z");
     std::string file("Props.vtk");

     an->DefineGraphMesh(dim,scalnames,vecnames,file);
     an->PostProcess(0,dim);

    
}
