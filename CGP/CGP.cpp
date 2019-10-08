//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

//#include "pzgeoel.h"
//#include "pzgnode.h"
//#include "pzgmesh.h"
//#include "pzbiharmonic.h"
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
//#include "pzelasmat.h"
//#include "pzlog.h"
#include "pzgengrid.h"

//#include <time.h>
//#include <stdio.h>

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

//#include "TPZMatLaplacian.h"
#include "pzpoisson3d.h"
#include "TPZNullMaterial.h"
#include "pzbndcond.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzsolve.h"
#include "TPZPersistenceManager.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#include "../TPZMixedDarcyWithFourSpaces.h"
#include "pzcmesh.h"
#include "pzbdstrmatrix.h"

#include "TPZHdivTransfer.h"
#include "pzseqsolver.h"
#include "pzmgsolver.h"
#include "TPZTimer.h"

#ifdef _AUTODIFF
//#include "tfad.h"
//#include "fad.h"
//#include "pzextractval.h"
#endif


//Creating geometric 1D, 2D and 3D mesh
TPZGeoMesh * GenerateGmesh1D(int nx, double l);
TPZGeoMesh * GenerateGmesh2D(int nx, int ny, double l, double h);
TPZGeoMesh * GenerateGmesh3D(int nx, int ny, int nz, double l, double h, double w);

//Creating computational meshes
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order_internal);
TPZCompMesh * GenerateConstantCmesh(TPZGeoMesh *Gmesh, bool third_LM);
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *Gmesh, int order_internal, int order_border);
TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int two_d_Q);

//Stablish de force fuction for 1D and 2D mesh
void Ladoderecho_1D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Ladoderecho_2D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void Ladoderecho_3D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//Creates index vector
void IndexVectorCoFi(TPZMultiphysicsCompMesh *Coarse_sol, TPZMultiphysicsCompMesh *Fine_sol, TPZVec<int64_t> & indexvec);

//Hdiv Test
void HDiv(int nx, int order_small, int order_high, bool condense_equations_Q, int dim);

//Transfer DOF from coarse mesh to fine mesh
void TransferDegreeOfFreedom(TPZFMatrix<STATE> & CoarseDoF, TPZFMatrix<STATE> & FineDoF, TPZVec<int64_t> & DoFIndexes);

//Analysis configuration
void ConfigurateAnalyses(TPZCompMesh * cmesh_c, TPZCompMesh * cmesh_f, bool must_opt_band_width_Q, int number_threads, TPZAnalysis *an_c,TPZAnalysis *an_f, bool UsePardiso_Q);

//Stablish an exact solution
void SolExact(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

//Split connects from a certain mesh
void SplitConnects(TPZCompMesh *fluxmesh, TPZGeoEl *gel ,int j);

//Shows shape functions for a certain element
void ShowShape(TPZCompMesh * cmesh, int element, int funcion,std::string plotname);


std::ofstream log_file("Results.txt");
using namespace std;

int main(int argc, char **argv){
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int fine_order_max = 2;
    int coarse_order = 1;
    int max_nx = 4;
    
    log_file<<"Order    NEls    qC   pC  qCav    pCav    mixC    qCCd    pCCd    qavCCd  pavCCd  mixCCd  qF   pF  qFav    pFav    mixF    qFCd    pFCd    qavFCd  pavFCd  mixFCd    NIter"<<std::endl;

    for (int fineorder=1; fineorder<= fine_order_max; fineorder++){
        for (int nx=1; nx<=max_nx; nx= nx+1) {
            std::cout<<"-------------------------------------------"<<endl;
            std::cout<<"*******************************************"<<endl;
            std::cout<<"Simulation: "<<(fineorder-1)*(max_nx) + nx <<" / "<<fine_order_max*max_nx<<std::endl;
            std::cout<<"*******************************************"<<endl;
            std::cout<<"-------------------------------------------"<<endl;
            log_file<<fineorder;
            HDiv(nx, coarse_order, fineorder, true, 2);
            log_file<<endl;
        }
    }
}


/**
 * @brief Runs a HDiv problem with 4 spaces for 1D 2D or 3D case
 * @param order_small: Low order for border elements
 * @param order_high: High order for internal elements
 * @param condense_equations_Q: Bool wether the problem is condensed or not
 * @param dim: Dimension of the problem
 */
void HDiv(int nx, int order_small, int order_high, bool condense_equations_Q, int dim){
    
    //Show Shapes functions
//    TPZGeoMesh *gmesh_1D = GenerateGmesh2D(nx, nx, 8, 8);      // Generates a 2D geo mesh
//    TPZCompMesh *flux = GenerateFluxCmesh(gmesh_1D, 1, 1);
//    int el_index=1;
//    int nfun=flux->Element(el_index)->NEquations();
//    for (int i=0; i<nfun; i++) {
//
//        std::string filename("elementFunc.vtk");
//        std::string file(filename+std::to_string(i)+".vtk");
//        ShowShape(flux,el_index,i,file);
//    };
//
//End Shapes functions

    bool KeepOneLagrangian = false;
    bool KeepMatrix = false;
    bool must_opt_band_width_Q = true;
    int number_threads = 0;
    
    TPZGeoMesh *gmesh;
    
    //Creates a geometric mesh with a given dimension
    switch (dim) {
        case 1:
            gmesh = GenerateGmesh1D(nx, 1);                 // 1D
            break;
        case 2:
            gmesh = GenerateGmesh2D(nx, nx, 1, 1);          // 2D
            break;
        case 3:
            gmesh = GenerateGmesh3D(nx, nx, nx, 1, 1, 1);   // 3D
            break;
    }
    int num=gmesh->NElements()-(6*nx);
   log_file<<" "<<num;
    
    TPZMultiphysicsCompMesh *MixedMesh_c = 0;
    TPZManVector<TPZCompMesh *> vecmesh_c(4);      //Vector for coarse mesh case (4 spaces)
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_small, order_small);
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_small);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh_c[0] = q_cmesh;              //Flux
        vecmesh_c[1] = p_cmesh;              //Pressure
        vecmesh_c[2] = gavg_cmesh;           //Average distribute flux
        vecmesh_c[3] = pavg_cmesh;           //Average pressure
        MixedMesh_c = GenerateMixedCmesh(vecmesh_c, dim);
        
        log_file<<" "<<q_cmesh->NEquations();
        log_file<<" "<<p_cmesh->NEquations();
        log_file<<" "<<gavg_cmesh->NEquations();
        log_file<<" "<<pavg_cmesh->NEquations();
        log_file<<" "<<MixedMesh_c->NEquations();
    
    if (condense_equations_Q) {             //Asks if you want to condesate the problem
        MixedMesh_c->ComputeNodElCon();
        int dim = MixedMesh_c->Dimension();
        int64_t nel = MixedMesh_c->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_c->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
        

        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_c, KeepOneLagrangian, KeepMatrix);
        int nconnects = MixedMesh_c->NConnects();
        int eqflux = 0;
        int eqpress = 0;
        int eqqav = 0;
        int eqpav = 0;
        for (int icon = 0; icon< nconnects; icon++) {
            TPZConnect &conect = MixedMesh_c->ConnectVec()[icon];
            if (conect.LagrangeMultiplier() == 0 && conect.IsCondensed()==0) {
                eqflux += conect.NShape();
            }
            if (conect.LagrangeMultiplier() == 1 && conect.IsCondensed()==0) {
                eqpress += conect.NShape();;
            }
            if (conect.LagrangeMultiplier() == 2 && conect.IsCondensed()==0) {
                eqqav +=conect.NShape();
            }
            if (conect.LagrangeMultiplier() == 3 && conect.IsCondensed()==0) {
                eqpav += conect.NShape();
            }
        }
        log_file<<" "<<eqflux;
        log_file<<" "<<eqpress;
        log_file<<" "<<eqqav;
        log_file<<" "<<eqpav;
        log_file<<" "<<MixedMesh_c->NEquations();
        }
    }
    
    TPZMultiphysicsCompMesh * MixedMesh_f = 0;
    TPZManVector<TPZCompMesh *> vecmesh_f(4);      //Vector for fine mesh case (4 spaces)
    {
        TPZCompMesh *q_cmesh = GenerateFluxCmesh(gmesh, order_high, order_small);
        
        TPZCompMesh *p_cmesh = GeneratePressureCmesh(gmesh, order_high);
        TPZCompMesh *gavg_cmesh = GenerateConstantCmesh(gmesh,false);
        TPZCompMesh *pavg_cmesh = GenerateConstantCmesh(gmesh,true);
        vecmesh_f[0] = q_cmesh;              //Flux
        vecmesh_f[1] = p_cmesh;              //Pressure
        vecmesh_f[2] = gavg_cmesh;           //Average distribute flux
        vecmesh_f[3] = pavg_cmesh;           //Average pressure
        
        MixedMesh_f = GenerateMixedCmesh(vecmesh_f, dim);
        log_file<<" "<<q_cmesh->NEquations();
        log_file<<" "<<p_cmesh->NEquations();
        log_file<<" "<<gavg_cmesh->NEquations();
        log_file<<" "<<pavg_cmesh->NEquations();
        log_file<<" "<<MixedMesh_f->NEquations();
        
    
    
    //Asks if you want to condesate the problem
    if (condense_equations_Q) {
        
        MixedMesh_f->ComputeNodElCon();
        int dim = MixedMesh_f->Dimension();
        int64_t nel = MixedMesh_f->NElements();
        for (int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = MixedMesh_f->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if(!gel) continue;
            if(gel->Dimension() != dim) continue;
            int nc = cel->NConnects();
            cel->Connect(nc-1).IncrementElConnected();
        }
        
        // Created condensed elements for the elements that have internal nodes
        TPZCompMeshTools::CreatedCondensedElements(MixedMesh_f, KeepOneLagrangian, KeepMatrix);
        
        int nconnects = MixedMesh_c->NConnects();
        int eqflux = 0;
        int eqpress = 0;
        int eqqav = 0;
        int eqpav = 0;
        for (int icon = 0; icon< nconnects; icon++) {
            TPZConnect &conect = MixedMesh_c->ConnectVec()[icon];
            if (conect.LagrangeMultiplier() == 0 && conect.IsCondensed()==0) {
                eqflux += conect.NShape();
            }
            
            if (conect.LagrangeMultiplier() == 1 && conect.IsCondensed()==0) {
                eqpress += conect.NShape();;
            }
            if (conect.LagrangeMultiplier() == 2 && conect.IsCondensed()==0) {
                eqqav +=conect.NShape();
            }
            if (conect.LagrangeMultiplier() == 3 && conect.IsCondensed()==0) {
                eqpav += conect.NShape();
            }
        }
        log_file<<" "<<eqflux;
        log_file<<" "<<eqpress;
        log_file<<" "<<eqqav;
        log_file<<" "<<eqpav;
        log_file<<" "<<MixedMesh_f->NEquations();
    }
    }
    
    //Solving the system:
    MixedMesh_c->InitializeBlock();    //Resequence the block object, remove unconnected connect objects
    MixedMesh_f->InitializeBlock();    //and reset the dimension of the solution vector
    TPZAnalysis *an_c = new TPZAnalysis;
    TPZAnalysis *an_f = new TPZAnalysis;
    //True to use pardiso
    ConfigurateAnalyses(MixedMesh_c, MixedMesh_f, must_opt_band_width_Q, number_threads, an_c, an_f, true);

//    std::cout<<"------------------------------"<<std::endl;
//    std::cout<<"Analysis configuration done"<<std::endl;
//    std::cout<<"------------------------------"<<std::endl;
    
    // Assembly fine operator
    an_f->Assemble();
    
//    std::cout<<"------------------------------"<<std::endl;
//    std::cout<<"Assembly Fine"<<std::endl;
//    std::cout<<"------------------------------"<<std::endl;
    
    // Assembly for coarse operator
    an_c->Assemble();

//    std::cout<<"------------------------------"<<std::endl;
//    std::cout<<"Assembly Coarse"<<std::endl;
//    std::cout<<"------------------------------"<<std::endl;
    
    // An iterative solution
    {

        // Constructing block diagonal
        if(1){
            TPZBlockDiagonalStructMatrix bdstr(MixedMesh_f);     //Give the fine mesh
            TPZBlockDiagonal<STATE> * sp = new TPZBlockDiagonal<STATE>();
            bdstr.AssembleBlockDiagonal(*sp);

            TPZAutoPointer<TPZMatrix<STATE> > sp_auto(sp);
            int64_t n_con = MixedMesh_f->NConnects();
            for (int ic = 0; ic < n_con; ic++) {
                TPZConnect & con = MixedMesh_f->ConnectVec()[ic];
                bool check = con.IsCondensed() || con.HasDependency() || con.LagrangeMultiplier() == 0;
                if (check) {
                    continue;
                }

                int64_t seqnum = con.SequenceNumber();
                int block_size = MixedMesh_f->Block().Size(seqnum);
                if (block_size != 1) {
                    continue;
                }

                int64_t pos = MixedMesh_f->Block().Position(seqnum);
                (*sp).PutVal(pos, pos, 1.0);
            }
            
//            std::cout<<"------------------------------"<<std::endl;
//            std::cout<<"Diagonal block constructed"<<std::endl;
//            std::cout<<"------------------------------"<<std::endl;
            
            TPZVec<int64_t> Indexes;
            IndexVectorCoFi(MixedMesh_c, MixedMesh_f, Indexes);
            int64_t neq_coarse = MixedMesh_c->NEquations();
            int64_t neq_fine = MixedMesh_f->NEquations();
            TPZHdivTransfer<STATE> *transfer = new TPZHdivTransfer<STATE>(neq_coarse, neq_fine, Indexes);
            TPZFMatrix<STATE> coarsesol(neq_coarse,1,1.), finesol(neq_fine,1,1.);
            transfer->Multiply(finesol, coarsesol,0); //It mutiplies itself by TPZMatrix<TVar>A adding the result in res
            transfer->Multiply(coarsesol, finesol,1); //z = beta * y(coarse) + alpha * opt(this)*x (fine)
            
            //Transfers the solution from coarse to fine mesh
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            an_c->Solver().Solve(coarsesol,coarsesol); //Force vector, solution
            TPZMGSolver<STATE> mgsolve(transfer,an_c->Solver(),1);
            mgsolve.SetMatrix(an_f->Solver().Matrix());
            mgsolve.Solve(finesol, finesol);
            
            //End of transferation
            //            finesol.Print(std::cout);
            //            coarsesol.Print(std::cout);
            
//            std::cout<<"------------------------------"<<std::endl;
//            std::cout<<"Transferation done"<<std::endl;
//            std::cout<<"------------------------------"<<std::endl;

            //Iterative method process
            
            TPZFMatrix<STATE> rhscoarse = an_c->Rhs();
            TPZFMatrix<STATE> rhsfine = an_f->Rhs();
            TPZStepSolver<STATE> BDSolve(sp_auto);
            BDSolve.SetDirect(ELU);
            TPZSequenceSolver<STATE> seqsolver;
            seqsolver.SetMatrix(an_f->Solver().Matrix());
            seqsolver.AppendSolver(mgsolve); //Updates the values of the preconditioner based on the values of the matrix
            seqsolver.AppendSolver(BDSolve);
            seqsolver.AppendSolver(mgsolve);
            seqsolver.Solve(rhsfine, finesol);
            TPZStepSolver<STATE> cg_solve(an_f->Solver().Matrix());
            cg_solve.SetCG(200, seqsolver, 1.e-10, 0);
            finesol.Zero();
            cg_solve.Solve(rhsfine, finesol);
            log_file<<" "<<cg_solve.NumIterations();
//            std::cout<<"------------------------------"<<std::endl;
//            std::cout<<"Iterative method done"<<std::endl;
//            std::cout<<"------------------------------"<<std::endl;
        }
        

     if(1){

            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_c, MixedMesh_c);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(vecmesh_f, MixedMesh_f);

            //PostProcess
            TPZStack<std::string> scalar, vectors;
            TPZManVector<std::string,10> scalnames(4), vecnames(1);
            vecnames[0]  = "q";
            scalnames[0] = "p";
            scalnames[1] = "kappa";
            scalnames[1] = "div_q";
            scalnames[2] = "g_average";
            scalnames[3] = "u_average";
            
//            std::ofstream filePrint_coarse("MixedHdiv_coarse.txt");
//            MixedMesh_c->Print(filePrint_coarse);
//            std::string name_coarse = "MixedHdiv_coarse.vtk";
//
         
//            std::ofstream filePrint_fine("MixedHdiv_fine.txt");
//            MixedMesh_f->Print(filePrint_fine);
//            std::string name_fine = "MixedHdiv_fine.vtk";
       

//            an_c->DefineGraphMesh(dim, scalnames, vecnames, name_coarse);
            an_c->PostProcess(0,dim);
           
//            an_f->DefineGraphMesh(dim, scalnames, vecnames, name_fine);
            an_f->PostProcess(0,dim);

         
//         std::cout<<"------------------------------"<<std::endl;
//         std::cout<<"Postprocess done"<<std::endl;
//         std::cout<<"------------------------------"<<std::endl;
        }
    }
}

/**
 * @brief Generates an 1D geometric mesh
 * @param nx: Partitions number on x
 * @param l: Mesh lenght
 * @return Geometric 1D mesh
 */

TPZGeoMesh * GenerateGmesh1D(int nx, double l){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    double h = l/nx;
    int Domain_Mat_Id = 1;
    int Inlet_bc_Id = -1;
    int Outletbc_Id = -2;
    TPZVec<REAL> xp(3,0.0);   //Creates vector nodes

    gmesh->NodeVec().Resize(nx+1);
    for (int64_t i=0; i<nx+1; i++) {
        xp[0] =(i)*h;
        gmesh->NodeVec()[i]= TPZGeoNode(i, xp, *gmesh);
    }
    
    TPZVec<int64_t> cornerindexes(2);   //Creates elements
    for (int64_t iel=0; iel<nx; iel++) {
        cornerindexes[0]=iel;
        cornerindexes[1]=iel+1;
        gmesh->CreateGeoElement(EOned, cornerindexes, Domain_Mat_Id, iel);
    }
    gmesh->Element(0)->CreateBCGeoEl(0, Inlet_bc_Id);     //Sets BC
    gmesh->Element(nx-1)->CreateBCGeoEl(1, Outletbc_Id);
    gmesh->SetDimension(1);
    gmesh->BuildConnectivity();
    return gmesh;
}

/**
 * @brief Generates a 2D geometric mesh
 * @param nx: number of partions on x
 * @param ny: number of partions on y
 * @param l: lenght
 * @param h: height
 * @return Geometric mesh
 */

TPZGeoMesh * GenerateGmesh2D(int nx, int ny, double l, double h){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y
    
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,l);
    x1[1]=h;
    x1[2]=0;
    
    //Setting boundary conditions (negative numbers to recognize them)
    TPZGenGrid gen(nels,x0,x1);
    gen.SetElementType(EQuadrilateral);
    gen.Read(gmesh);
    gen.SetBC(gmesh, 4, -1);
    gen.SetBC(gmesh, 5, -2);
    gen.SetBC(gmesh, 6, -3);
    gen.SetBC(gmesh, 7, -4);
    return gmesh;
}

/**
 * @brief Generates a 3D geometric mesh
 * @param nx: Partitions number on x
 * @param ny: Partitions number on y
 * @param nz: Partitions number on z
 * @param l: Mesh lenght
 * @param h: Mesh heigh
 * @param w: Mesh width
 * @return Geometric 3D mesh
 */

TPZGeoMesh * GenerateGmesh3D(int nx, int ny, int nz, double l, double h, double w){
    
    TPZVec<int> nels(3,0);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y
    
    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,l);
    x1[1]=h;
    x1[2]=0;
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    TPZGenGrid gen(nels,x0,x1);
    gen.SetElementType(EQuadrilateral);
    gen.Read(gmesh);
    double var = w/nz;
    TPZExtendGridDimension extend(gmesh,var);
    extend.SetElType(1);
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(nz);
    
    for (auto gel:gmesh3d->ElementVec()) {
        TPZFMatrix<REAL> coordinates;
        gel->NodesCoordinates(coordinates);
        if(coordinates(2,0)==0){
            gel->CreateBCGeoEl(20, -1);
        }
        if(coordinates(2,4)==w){
            gel->CreateBCGeoEl(25, -2);
        }
        
        if(coordinates(0,0)==0.0 ){
            gel->CreateBCGeoEl(24, -3);
        }
        if(coordinates(1,0)==0.0 ){
            gel->CreateBCGeoEl(21, -4);
        }
        
        if(coordinates(0,1)== l ){
            gel->CreateBCGeoEl(22, -5);
        }
        if(coordinates(1,3)==h){
            gel->CreateBCGeoEl(23, -6);
        }
    };
   
    gmesh3d->SetDimension(3);
    gmesh3d->BuildConnectivity();

    return gmesh3d;
}

/**
 * @brief Generates the pressure computational mesh
 * @param Geometric mesh
 * @param Order
 * @return Pressure computational mesh
 */
TPZCompMesh * GeneratePressureCmesh(TPZGeoMesh *Gmesh, int order){
    
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(order);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    Cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    STATE Permeability=1;
    
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(MaterialId, dimen);
    
    //No convection
    REAL conv=0;
    TPZVec<REAL> convdir(dimen, 0);
    mat->SetParameters(Permeability, conv, convdir);
    
    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    //Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    
    //Set Lagrange multiplier
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    return Cmesh;
}


/**
 * @brief Generates the constant computational mesh
 * @param Gmesh: Geometric mesh
 * @param third_LM: Bool Third Lagrange multiplier
 * @return Constant computational mesh
 */
TPZCompMesh * GenerateConstantCmesh(TPZGeoMesh *Gmesh, bool third_LM)
{
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(0);
    Cmesh->SetAllCreateFunctionsDiscontinuous();

    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    
    TPZNullMaterial *mat =new TPZNullMaterial(MaterialId);
    mat->SetDimension(dimen);
    mat->SetNStateVariables(1);

    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    //Autobuild
    Cmesh->AutoBuild();
    
    int ncon = Cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = Cmesh->ConnectVec()[i];
        if (third_LM) {
            newnod.SetLagrangeMultiplier(3);
        }else{
            newnod.SetLagrangeMultiplier(2);
        }
    }
    return Cmesh;
}


/**
 * @brief Generates the flux computational mesh
 * @param mesh: Geometric mesh
 * @param order_internal: Order used for internal elements
 * @param order_border: Order used for border elements
 * @return Flux computational mesh
 */
TPZCompMesh * GenerateFluxCmesh(TPZGeoMesh *mesh, int order_internal, int order_border){
    
    int dimen = mesh->Dimension();
    TPZCompMesh *Cmesh = new TPZCompMesh(mesh);
    Cmesh->SetDimModel(dimen);
    Cmesh->SetDefaultOrder(order_border);
    
    //Definition of the approximation space
    int perm=1;
    REAL conv=0;
    REAL perme=1;
    TPZVec<REAL> convdir(dimen , 0.0);
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(perm , dimen);
    mat->SetParameters(perme, conv, convdir);
    
    //Inserting volumetric materials objects
    Cmesh->InsertMaterialObject(mat);
    
    //Create H(div) functions
 
    Cmesh->SetAllCreateFunctionsHDiv();
    
    //Insert boundary conditions
    int D=0;
    int N=1;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
     val2(0,0)=0;
    TPZMaterial *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    Cmesh->InsertMaterialObject(bc1);

    int BC2=-2;
    val2(0,0)=0;
    TPZMaterial *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    Cmesh->InsertMaterialObject(bc2);

    int BC3=-3;
    val2(0,0)=0;
    TPZMaterial *bc3 = mat->CreateBC(mat, BC3,D, val1, val2);
    Cmesh->InsertMaterialObject(bc3);

    int BC4=-4;
    val2(0,0)=0;
    TPZMaterial *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    Cmesh->InsertMaterialObject(bc4);
    
    int BC5=-5;
    val2(0,0)=0;
    TPZMaterial *bc5 = mat->CreateBC(mat, BC5, D, val1, val2);
    Cmesh->InsertMaterialObject(bc5);
    
    int BC6=-6;
    val2(0,0)=0;
    TPZMaterial *bc6 = mat->CreateBC(mat, BC6, D, val1, val2);
    Cmesh->InsertMaterialObject(bc6);
    

    Cmesh->AutoBuild();
    
    int64_t nel = Cmesh->NElements();
    for (int el=0; el<nel; el++) {
        TPZCompEl *cel = Cmesh->Element(el);
        if(!cel) continue;
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) DebugStop();
        TPZGeoEl *gel = intel->Reference();
        intel->SetSideOrder(gel->NSides()-1, order_internal);
    }
    Cmesh->ExpandSolution();
    
    return Cmesh;
}

/**
 * @brief Generates the mixed computational mesh
 * @param fvecmesh: Vector thats contains flux and pressure computational mesh
 * @param order: Ordertwo_d_Q
 * @param two_d_Q: Bool wether is 1D (false) or 2D (true)
 * @return Mixed computational mesh
 */

TPZMultiphysicsCompMesh * GenerateMixedCmesh(TPZVec<TPZCompMesh *> fvecmesh, int dim){
    TPZGeoMesh *gmesh = fvecmesh[1]->Reference();
    TPZMultiphysicsCompMesh *MixedMesh = new TPZMultiphysicsCompMesh(gmesh);
    
    //Definition of the approximation space
    int dimen= gmesh->Dimension();
    int matnum=1;
    REAL perm=1;
    
    //Inserting material
    TPZMixedDarcyWithFourSpaces * mat = new TPZMixedDarcyWithFourSpaces(matnum, dimen);
    mat->SetPermeability(perm);
    TPZAutoPointer<TPZFunction<STATE> > sourceterm ;
    switch (dim) {
        case 1:
           sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_1D, 10);
            mat->SetForcingFunction(sourceterm);
            break;
        case 2:
            sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_2D, 10);
            mat->SetForcingFunction(sourceterm);
            break;
        case 3:
             sourceterm = new TPZDummyFunction<STATE>(Ladoderecho_3D, 10);
            mat->SetForcingFunction(sourceterm);
        break;
    }
    
    
    //Inserting volumetric materials objects
    MixedMesh->InsertMaterialObject(mat);
    
    //Boundary conditions
    int D=0;
    int N=1;
    int BC1=-1;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZFMatrix<STATE> val2(1,1,0.0);
    
    val2(0,0)=0.0;
    TPZMaterial *bc1 = mat->CreateBC(mat, BC1, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc1);
    
    int BC2=-2;
    val2(0,0)=0;
    TPZMaterial *bc2 = mat->CreateBC(mat, BC2, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc2);
    
    int BC3=-3;
    val2(0,0)=0;
    TPZMaterial *bc3 = mat->CreateBC(mat, BC3, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc3);

    int BC4=-4;
    val2(0,0)=0;
    TPZMaterial *bc4 = mat->CreateBC(mat, BC4, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc4);
    
    int BC5=-5;
    val2(0,0)=0;
    TPZMaterial *bc5 = mat->CreateBC(mat, BC5, D, val1, val2);
    MixedMesh->InsertMaterialObject(bc5);
    
    int BC6=-6;
    val2(0,0)=0;
    TPZMaterial *bc6 = mat->CreateBC(mat, BC6,D, val1, val2);
    MixedMesh->InsertMaterialObject(bc6);
    
    
    MixedMesh->SetAllCreateFunctionsMultiphysicElem();
    MixedMesh->SetDimModel(dimen);
    
    //Autobuild
    TPZManVector<int,5> active_approx_spaces(4); /// 4 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    active_approx_spaces[3] = 1;
    MixedMesh->BuildMultiphysicsSpace(active_approx_spaces,fvecmesh);
    
    TPZBuildMultiphysicsMesh::AddElements(fvecmesh, MixedMesh);
    TPZBuildMultiphysicsMesh::AddConnects(fvecmesh,MixedMesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fvecmesh, MixedMesh);
    
    return MixedMesh;
};

/**
 * @brief Generates a force function for the 1D case
 * @param pt: Points values
 * @param disp: Vector
 */
void Ladoderecho_1D (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    double fx= 4*M_PI*M_PI*sin(2*M_PI*x);        //Force function definition
    disp[0]=fx;
    
}

/**
 * @brief Generates a force function for the 2D case
 * @param pt: Points values
 * @param disp: Vector
 */

void Ladoderecho_2D (const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    STATE x = pt[0];
    STATE y = pt[1];
    
//    double fx= 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);        //Force function definition
    
//    double fx= 2*(-1 + x)*x*(-1 + y)*y + 2*(-1 + x)*x + 2*(-1 + y)*y;
    
//    double fx =-4144.653167389283*pow(10,
//                                      -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
//    (-2*M_PI + 15.*x) + 4771.70829943056*
//    pow(10,-pow(-2*M_PI + 15.*x,2) -
//        pow(-2*M_PI + 15.*y,2))*pow(-2*M_PI + 15.*x,3) +
//    4771.70829943056*pow(10,
//                         -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
//    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);

    double fx = -4144.653167389282*pow(2,2 - pow(-2*M_PI + 15.*x,2) -
                                       pow(-2*M_PI + 15.*y,2))*
    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    (-2*M_PI + 15.*x) + 4771.708299430558*
    pow(2,2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    pow(-2*M_PI + 15.*x,3) +
    4771.708299430558*pow(2,
                          2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
    
    disp[0]=fx;
}

/**
 * @brief Generates a force function for the 3D case
 * @param pt: Points values
 * @param disp: Vector
 */

void Ladoderecho_3D(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    STATE x = pt[0];
    STATE y = pt[1];
    STATE z = pt[2];
    
    double fx= 2*(-1 + x)*x*(-1 + y)*y + 2*(-1 + x)*x*(-1 + z)*z + 2*(-1 + y)*y*(-1 + z)*z;      //Force function definition
    
    disp[0]=fx;
}

/**
 * @brief Generates a index vector which relates coarse and fine mesh
 * @param Coarse_sol: Coarse mesh
 * @param Fine_sol: Fine mesh
 * @param indexvec:
 */
void IndexVectorCoFi(TPZMultiphysicsCompMesh *Coarse_sol, TPZMultiphysicsCompMesh *Fine_sol, TPZVec<int64_t> & indexvec)
{
    int64_t maxcone = Coarse_sol->NConnects();
    
    int64_t indexvecsize = 0;
    for (int j=0; j<maxcone; j++) {
        bool is_condensed = Coarse_sol->ConnectVec()[j].IsCondensed();
        if (is_condensed == true ) continue;

        int64_t sequence_coarse = Coarse_sol->ConnectVec()[j].SequenceNumber();
        int blocksize_coarse = Coarse_sol->Block().Size(sequence_coarse);

        indexvecsize += blocksize_coarse;
    }
    indexvec.Resize(indexvecsize,-1);
    
    for (int j=0; j<maxcone; j++) {
        bool is_condensed = Coarse_sol->ConnectVec()[j].IsCondensed();
        if (is_condensed == true ) continue;

        int64_t sequence_coarse = Coarse_sol->ConnectVec()[j].SequenceNumber();
        int64_t sequence_fine = Fine_sol->ConnectVec()[j].SequenceNumber();
        int blocksize_coarse = Coarse_sol->Block().Size(sequence_coarse);
        int blocksize_fine = Fine_sol->Block().Size(sequence_fine);

        if (blocksize_coarse > blocksize_fine) {
            DebugStop();

        }
    
        for(int i=0; i<blocksize_coarse; i++){
            int64_t pos_coarse = Coarse_sol->Block().Position(sequence_coarse);
            int64_t pos_fine = Fine_sol->Block().Position(sequence_fine);
            indexvec[pos_coarse+i] = pos_fine+i;
        }
    }
}

/**
 * @brief Transfer DOF from coarse to fine mesh
 * @param CoarseDoF: Solution coarse matrix
 * @param FineDoF: Solution fine matrix
 * @param DoFIndexes: DOF index vector
 */
void TransferDegreeOfFreedom(TPZFMatrix<STATE> & CoarseDoF, TPZFMatrix<STATE> & FineDoF, TPZVec<int64_t> & DoFIndexes){
    
    int64_t n_data = DoFIndexes.size();
    for (int64_t i = 0 ; i < n_data; i++) {
        FineDoF(DoFIndexes[i],0) = CoarseDoF(i,0);
    }
}

/**
 * @brief Configurate the matrix for analysis
 * @param cmesh_c: Computational coarse mesh
 * @param cmesh_f: Computational fine mesh
 * @param must_opt_band_width_Q: Wether the band width is optimized or not (it is neccesaty to rearrange the matrix in order to be less sparse)
 * @param number_threads: 
 * @param an_c: Coarse analysis mesh
 * @param an_f: Fine analysis mesh
 * @param UsePardiso_Q: Wether using Pardiso or not
 */
void ConfigurateAnalyses(TPZCompMesh * cmesh_c, TPZCompMesh * cmesh_f, bool must_opt_band_width_Q, int number_threads, TPZAnalysis *an_c,TPZAnalysis *an_f, bool UsePardiso_Q){
    
        an_c->SetCompMesh(cmesh_c,must_opt_band_width_Q);
        an_f->SetCompMesh(cmesh_f,must_opt_band_width_Q);
        TPZStepSolver<STATE> step;
        if (UsePardiso_Q) {
    
            TPZSymetricSpStructMatrix sparse_matrix_coarse(cmesh_c);
            TPZSymetricSpStructMatrix sparse_matrix_fine(cmesh_f);
            sparse_matrix_coarse.SetNumThreads(number_threads);
            sparse_matrix_fine.SetNumThreads(number_threads);
            an_c->SetStructuralMatrix(sparse_matrix_coarse);
            an_f->SetStructuralMatrix(sparse_matrix_fine);
    
        }else{
    
            TPZSkylineStructMatrix sparse_matrix_coarse(cmesh_c);
            TPZSkylineStructMatrix sparse_matrix_fine(cmesh_f);
            sparse_matrix_coarse.SetNumThreads(number_threads);
            sparse_matrix_fine.SetNumThreads(number_threads);
            an_c->SetStructuralMatrix(sparse_matrix_coarse);
            an_f->SetStructuralMatrix(sparse_matrix_fine);
    
        }
        step.SetDirect(ELDLt);
        an_c->SetSolver(step);
        an_f->SetSolver(step);
    }

/**
 * @brief Shows shape functions of a certain element
 * @param cmesh: Computational mesh
 * @param element: Element number to show it shapes functions
 * @param funcion: Shape function number to show
 * @return plotfile: Name for the plot to print
 */
void ShowShape(TPZCompMesh * cmesh, int element, int funcion, std::string plotfile){
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    
    int nels = cmesh->NElements();
    for (int iel =0 ; iel<nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        TPZGeoEl *gel = cel->Reference();
        if (gel->Index() == element) {
            gel->SetMaterialId(100);
        }
    }
    
    TPZMaterial  * mat_2(cmesh->MaterialVec()[1]);
    std::map<int, TPZMaterial *> matvec;
    mat_2->Clone(matvec);
    cmesh->CleanUp();
 
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(100 , 2);
    TPZMaterial *perf_mat( matvec[1]);
    TPZMatPoisson3d *aux_mat = dynamic_cast<TPZMatPoisson3d *>(perf_mat);
    
    cmesh->InsertMaterialObject(mat);
    cmesh->AutoBuild();
  
    int gdl = cmesh->Solution().Rows();
    int nfuncols = cmesh->Solution().Cols();
    int nfunrows = cmesh->Solution().Rows();
    TPZFMatrix<STATE> solu(nfunrows,nfuncols,1.0);
    TPZFMatrix<STATE> sol(gdl,1,0.0);
  
    sol.Resize(gdl, 1);
    int nel = cmesh->NElements();
    int val =1;
    for (int i=0; i<nel; i++) {
        TPZCompEl *cel = cmesh->Element(i);
        TPZGeoEl *gel =cel->Reference();
        int index = gel->Index();
        if (index != element) {
            continue;
        }

        int nconnects = cel->NConnects();
        int acum =0;
       
        for(int j=0; j<nconnects; j++){
            int sec_number = cel->Connect(j).fSequenceNumber;
            if(sec_number > -1)
            {
                acum = acum + cmesh->Block().Size(sec_number);
                if (acum > funcion) {

//                    if (val) {
                    
                        //SplitConnects(cmesh, gel ,j);
//                        cmesh->CleanUp();
                        val=0;
                        
                        int delta = funcion + cmesh->Block().Size(sec_number) - acum;
                        int64_t pos = cmesh->Block().Position(sec_number);
                        pos = pos +  delta;
                        sol.Zero();
                        sol(pos,0)=-1.0;
                    if (funcion == 3 || funcion == 6 ){sol(pos,0)=1.0;};
                        cmesh->LoadSolution(sol);

                        TPZAnalysis *an = new TPZAnalysis(cmesh,false);
                        {
                            const int dim = an->Mesh()->Dimension();
                            int div = 3;
                            // std::string plotfile = "SHAPES.vtk";
                            TPZStack<std::string> scalar_names, vec_names;
                            
                            scalar_names.push_back("Solution");
                            an->DefineGraphMesh(dim,scalar_names,vec_names,plotfile);
                            an->PostProcess(div,dim);
                            std::cout << "The function :"<<funcion<<" of element "<<element<<" has been printed on .vtk file" << std::endl;
                            
                        }
                        int nelss = gmesh->NElements();
                        for (int iel =0 ; iel<nelss; iel++) {
                            TPZGeoEl *gel = gmesh->Element(iel);
                            if (gel->Index() == element) {
                                gel->SetMaterialId(1);
                            }
                        cmesh->CleanUp();
                        cmesh->InsertMaterialObject(aux_mat);
                        cmesh->AutoBuild();
                        
                        return;

                    }
                }
            }
        }
    }
    
    
}
void SplitConnects(TPZCompMesh *fluxmesh, TPZGeoEl *gel , int j){
    TPZCompEl *cel = gel->Reference();
    int nelconnect = cel->Connect(j).NElConnected();
    if (nelconnect>1) {
        
    }
   
    int iside = j + 4;
    if (iside >= gel->NSides()-1) {
        return;
    }
        TPZGeoElSide neighb(gel, iside);
        TPZStack<TPZGeoElSide> allneigh;
        TPZGeoElSide neig = neighb.Neighbour();
        if (neig.Element()->Dimension() !=2) {
            return;
        }
    
        TPZCompElSide left(cel,iside);
        TPZGeoElSide gelr = gel->Neighbour(iside);
        TPZCompElSide right(gelr.Element()->Reference(), gelr.Side());
    
        TPZGeoElSide gleft(gel,iside);
        TPZGeoElSide gright(gleft.Neighbour());
        TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
        TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *> (right.Element());
        intelleft->SetSideOrient(left.Side(), 1);
        intelright->SetSideOrient(right.Side(), 1);
        TPZStack<TPZCompElSide> equalright;
        TPZConnect &cleft = intelleft->SideConnect(0, left.Side());
        int conr = right.Element()->Connect(right.Side()-4).SequenceNumber();
       int conl = left.Element()->Connect(left.Side()-4).SequenceNumber();
    
    if (conr != conl) {
        return;
    }
    
        int64_t index = fluxmesh->AllocateNewConnect(cleft);
        TPZConnect &newcon = fluxmesh->ConnectVec()[index];
        cleft.DecrementElConnected();
        newcon.ResetElConnected();
        newcon.IncrementElConnected();
        newcon.SetSequenceNumber(fluxmesh->NConnects() - 1);
        
        int rightlocindex = intelright->SideConnectLocId(0, right.Side());
        intelright->SetConnectIndex(rightlocindex, index);
        
        int sideorder = cleft.Order();
        fluxmesh->SetDefaultOrder(sideorder);
        
        int gdl = fluxmesh->Solution().Rows();
        fluxmesh->Solution().Resize(gdl+2, 1);
        TPZFMatrix<STATE> sol(gdl+2,1,0.0);
        fluxmesh->LoadSolution(sol);

    return;
}
