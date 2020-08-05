//
//  TPZFastCondensedElement.cpp
//  Fast Mixed Finite Elements
//
//  Created by PHILIPPE DEVLOO on 28/4/2020.
//  Copyright © 2020 PHILIPPE DEVLOO. All rights reserved.
//


#include "TPZFastCondensedElement.h"
#include "pzcmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMixedDarcyWithFourSpaces.h"
#include "TPZDarcyFlowWithMem.h"
#include "TPZSSpMatrixEigen.h"
bool TPZFastCondensedElement::fSkipLoadSolution = true;

/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
void TPZFastCondensedElement::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef)
{
    
    if(this->fMatrixComputed == false)
    {
     
//        TPZMixedDarcyWithFourSpaces *matDarcy = dynamic_cast<TPZMixedDarcyWithFourSpaces *>(Material());
        TPZDarcyFlowWithMem *matDarcy = dynamic_cast<TPZDarcyFlowWithMem *>(Material());
        if (!matDarcy) {
            DebugStop();
        }
        
//        matDarcy->SetPermeability(fPermeabilityTensor);
        TPZCondensedCompEl::CalcStiff(ek, ef);
        ShrinkElementMatrix(ek, fEK);
        ShrinkElementMatrix(ef, fEF);
        this->fMatrixComputed = true;
    }
    
    ek = fEK;
    ef = fEF;
    
    int nrows = ek.fMat.Rows();
    int ncols = ek.fMat.Rows();
    

    REAL Glambda = fMixedDensity;

    ek.fMat *= (1./fLambda);
    for (int icol=0; icol<ncols; icol++) {
        ek.fMat(nrows-1,icol) *= fLambda;
    }
    for (int irow=0; irow<nrows; irow++) {
        ek.fMat(irow,ncols-1) *= fLambda;
    }
    ek.fMat(nrows-1,ncols-1) *=fLambda;
    
    TPZFMatrix<STATE> solvec(fEK.fMat.Rows(),1,0.);
    GetSolutionVector(solvec);
    

    ef.fMat *= -1.0*0.0*Glambda;
    REAL norm =Norm(ef.fMat);
    REAL normek = Norm(ek.fMat);
    if (norm>1000) {
        DebugStop();
    }
    /** @brief Computes z = alpha * opt(this)*x + beta * y */
    /** @note z and x cannot overlap in memory */
    //    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
    //                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
    STATE alpha = 1.;
//    TPZSYsmpMatrixEigen<REAL>* mat = dynamic_cast<TPZSYsmpMatrixEigen<REAL> >(ek.fMat) ;

    ek.fMat.MultAdd(solvec, ef.fMat, ef.fMat, alpha, 1);
    REAL norm3 =Norm(ef.fMat);
    REAL normek3 = Norm(ek.fMat);
    if (norm3>1000) {
        ef.fMat *=0.0;
        ek.fMat.Print("Ek=", std::cout, EMathematicaInput);
        solvec.Print("Solvec=", std::cout, EMathematicaInput);
        solvec.Print("Ef=", std::cout, EMathematicaInput);
        ek.fMat.MultAdd(solvec, ef.fMat, ef.fMat, alpha, 1);
//        DebugStop();
    }
    int ol=0;
}

// extract the solution vector of the condensed element
void TPZFastCondensedElement::GetSolutionVector(TPZFMatrix<STATE> &solvec)
{
    int nc = fEK.fConnect.size();
    TPZCompMesh *cmesh = Mesh();
    int64_t vecsize = fEK.fMat.Rows();
    int count = 0;
    for(int ic=0; ic<nc; ic++)
    {
        int64_t cindex = fEK.fConnect[ic];
        TPZConnect &c = cmesh->ConnectVec()[cindex];
        int64_t seqnum = c.SequenceNumber();
        int blsize = c.NShape()*c.NState();
        for(int dof=0; dof<blsize; dof++)
        {
            solvec(count+dof,0) = cmesh->Block()(seqnum,0,dof,0);
        }
        count += blsize;
    }
    if(count != vecsize) DebugStop();
}

/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
void TPZFastCondensedElement::CalcResidual(TPZElementMatrix &ef)
{
      TPZElementMatrix ek;
      CalcStiff(ek, ef);

}

void TPZFastCondensedElement::ShrinkElementMatrix(TPZElementMatrix &input, TPZElementMatrix &output)
{
    output.fType = input.fType;
    output.fMesh = input.fMesh;
    int nindep = 0;
    int64_t condense_size = 0;
    int ncon = input.fConnect.size();
    output.fConnect.Resize(ncon);
    int firstnon = -1;
    for (int ic = 0; ic<ncon; ic++) {
        int64_t cindex = input.fConnect[ic];
        TPZConnect &c = input.fMesh->ConnectVec()[cindex];
        if(c.IsCondensed() && firstnon != -1)
        {
            DebugStop();
        }
        if(c.IsCondensed()) continue;
        if(firstnon == -1) firstnon = ic;
        output.fConnect[nindep] = cindex;
        condense_size += c.NShape()*c.NState();
        nindep++;
    }
    output.fConnect.Resize(nindep);
    if(input.fType == TPZElementMatrix::EK)
    {
        output.fMat.Redim(condense_size,condense_size);
    } else
    {
        output.fMat.Redim(condense_size,input.fMat.Cols());
    }
    int row_orig = input.fMat.Rows();
    if(input.fType == TPZElementMatrix::EK)
    {
        input.fMat.GetSub(row_orig-condense_size, row_orig-condense_size, condense_size, condense_size, output.fMat);
    } else
    {
        input.fMat.GetSub(row_orig-condense_size, 0, condense_size,
                          input.fMat.Cols(), output.fMat);
    }
    output.fBlock.SetNBlocks(nindep);
    
    for (int ic = ncon-nindep; ic < ncon; ic++) {
        output.fBlock.Set(ic-ncon+nindep, input.fBlock.Size(ic));
    }
    output.fBlock.Resequence();
    
    
    
}

void TPZFastCondensedElement::SetLambda(REAL lambda){
    fLambda = lambda;
}
REAL TPZFastCondensedElement::GetLambda(){
    return fLambda;
}
void TPZFastCondensedElement::SetMixedDensity(REAL mdensity){
    fMixedDensity = mdensity;
}
REAL TPZFastCondensedElement::GetMixedDensity(){
    return fMixedDensity;
}
void TPZFastCondensedElement::SetPermTensorAndInv(TPZFNMatrix<9, REAL> &PermeabilityTensor, TPZFNMatrix<9, REAL> &InvPerm){
    fPermeabilityTensor = PermeabilityTensor;
    fInvPerm = InvPerm;
}
TPZFMatrix<REAL> &TPZFastCondensedElement::GetPermTensor(){
    return  fPermeabilityTensor;
}

/**
 * @brief Calculates the solution - sol - for the variable var
 * at point qsi, where qsi is expressed in terms of the
 * master element coordinates
 * @param qsi master element coordinate
 * @param var variable name
 * @param sol vetor for the solution
 */
void TPZFastCondensedElement::Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    switch (var) {
        case 7:
            sol[0] = fPermeabilityTensor(0,0);
            break;
        case 8:
            sol[0] = fPermeabilityTensor(1,1);
            break;
        case 9:
            sol[0] = fPermeabilityTensor(2,2);
            break;
        case 10:
            sol[0] = fLambda;
            break;
        default:
            TPZCondensedCompEl::Solution(qsi, var, sol);
            break;
    }
}

/** @brief Loads the solution within the internal data structure of the element */
/**
 * Is used to initialize the solution of connect objects with dependency \n
 * Is also used to load the solution within SuperElements
 */
void TPZFastCondensedElement::LoadSolution()
{
    if(fSkipLoadSolution)
    {
        TPZCompEl::LoadSolution();
    }
    else
    {
        AdjustPressureCoefficients();
    }
}

// global indices of the pressure equations
void TPZFastCondensedElement::PressureEquations(TPZVec<int64_t> &eqs)
{
    TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(fReferenceCompEl);
    if(!mphys) DebugStop();
    int nel = mphys->NMeshes();
    if(nel != 4 && nel != 6) DebugStop();
    int firstpressure_connect = mphys->Element(0)->NConnects();
    int numpressure_connects = mphys->Element(1)->NConnects();
    int numpressure_equations = 0;
    for(int ic = firstpressure_connect; ic < firstpressure_connect+numpressure_connects; ic++)
    {
        TPZConnect &c = fReferenceCompEl->Connect(ic);
        int neq = c.NShape()*c.NState();
        numpressure_equations += neq;
    }
    TPZCompMesh *cmesh = Mesh();
    TPZBlock<STATE> &block = cmesh->Block();
    eqs.Resize(numpressure_equations, 0);
    int count = 0;
    for(int ic = firstpressure_connect; ic < firstpressure_connect+numpressure_connects; ic++)
    {
        TPZConnect &c = fReferenceCompEl->Connect(ic);
        int neq = c.NShape()*c.NState();
        int64_t seqnum = c.SequenceNumber();
        int64_t firsteq = block.Position(seqnum);
        for(int i = 0; i<neq; i++)
        {
            eqs[count] = firsteq+i;
            count++;
        }
    }
}

// global index of the average pressure equation
int64_t TPZFastCondensedElement::AveragePressureEquation()
{
    TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(fReferenceCompEl);
    if(!mphys) DebugStop();
    int nel = mphys->NMeshes();
    if(nel != 4 && nel != 6) DebugStop();
    int firstpressure_connect = mphys->Element(0)->NConnects()+mphys->Element(1)->NConnects()
        +mphys->Element(2)->NConnects();
    int numpressure_connects = mphys->Element(3)->NConnects();
    if(numpressure_connects != 1) DebugStop();
    TPZConnect &c = mphys->Connect(firstpressure_connect);
    int64_t seq_num = c.SequenceNumber();
    int64_t globeq = Mesh()->Block().Position(seq_num);
    return globeq;
}

// global indices of the boundary flux equations
void TPZFastCondensedElement::BoundaryFluxEquations(TPZVec<int64_t> &eqs)
{
    TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(fReferenceCompEl);
    if(!mphys) DebugStop();
    int nel = mphys->NMeshes();
    if(nel != 4 && nel != 6) DebugStop();
    int firstflux_connect = 0;
    int numflux_connects = mphys->Element(0)->NConnects();
    int numflux_equations = 0;
    TPZCompMesh *cmesh = Mesh();

    // skip the internal flux connect
    for(int ic = firstflux_connect; ic < firstflux_connect+numflux_connects-1; ic++)
    {
        TPZConnect &c = fReferenceCompEl->Connect(ic);
        int neq = c.NShape()*c.NState();
        numflux_equations += neq;
        if(c.HasDependency())
        {
            int64_t cindex = c.FirstDepend()->fDepConnectIndex;
            TPZConnect &cdep = cmesh->ConnectVec()[cindex];
            neq = cdep.NShape()*cdep.NState();
            numflux_equations += neq;
        }
    }
    TPZBlock<STATE> &block = cmesh->Block();
    eqs.Resize(numflux_equations, 0);
    int count = 0;
    for(int ic = firstflux_connect; ic < firstflux_connect+numflux_connects-1; ic++)
    {
        TPZConnect &c = fReferenceCompEl->Connect(ic);
        int neq = c.NShape()*c.NState();
        int64_t seqnum = c.SequenceNumber();
        int64_t firsteq = block.Position(seqnum);
        for(int i = 0; i<neq; i++)
        {
            eqs[count] = firsteq+i;
            count++;
        }
        if(c.HasDependency())
        {
            int64_t cindex = c.FirstDepend()->fDepConnectIndex;
            TPZConnect &cdep = cmesh->ConnectVec()[cindex];
            seqnum = cdep.SequenceNumber();
            firsteq = block.Position(seqnum);
            for(int i = 0; i<neq; i++)
            {
                eqs[count] = firsteq+i;
                count++;
            }
        }
    }
}

// adjust the multiplying coeficients of the pressure equations
void TPZFastCondensedElement::AdjustPressureCoefficients()
{
    TPZManVector<int64_t,20> fluxeqs, pressureqs;
    PressureEquations(pressureqs);
    BoundaryFluxEquations(fluxeqs);
    int64_t averagepressureq = AveragePressureEquation();
    TPZMatrix<STATE> &solution = *(Mesh()->Block().Matrix());
    int npres = pressureqs.size();
    int nflux = fluxeqs.size();
    TPZManVector<STATE> gravity_pressure(npres,0.);
    TPZManVector<STATE> average_pressure(npres,0.);
    TPZManVector<STATE> flux_pressure(npres,0.);
    TPZManVector<STATE> boundary_fluxes(nflux,0.);
    STATE average = solution(averagepressureq,0);
    for (int ifl = 0; ifl<nflux; ifl++) {
        boundary_fluxes[ifl] = solution(fluxeqs[ifl],0);
    }
    // build the gravity pressure coefs
    solution(averagepressureq,0) = 0.;
    // zero de boundary fluxes
    for(int ifl=0; ifl < nflux; ifl++) solution(fluxeqs[ifl],0) = 0.;
    TPZCondensedCompEl::LoadSolution();
    for(int ipr=0; ipr < npres; ipr++) gravity_pressure[ipr] = solution(pressureqs[ipr],0);
    
    // build the pressure for constant pressure
    // zero the boundary fluxes
    for(int ifl=0; ifl < nflux; ifl++) solution(fluxeqs[ifl],0) = 0.;
    solution(averagepressureq,0) = average;
    TPZCondensedCompEl::LoadSolution();
    for(int ipr=0; ipr < npres; ipr++) average_pressure[ipr] = solution(pressureqs[ipr],0)-gravity_pressure[ipr];

    // build the pressure solution due to boundary fluxes
    // restore the boundary fluxes
    for(int ifl=0; ifl < nflux; ifl++) solution(fluxeqs[ifl],0) = boundary_fluxes[ifl];
    solution(averagepressureq,0) = 0.;
    TPZCondensedCompEl::LoadSolution();
    for(int ipr=0; ipr < npres; ipr++) flux_pressure[ipr] = solution(pressureqs[ipr],0)-gravity_pressure[ipr];
    // compose the final pressure coeficients
    solution(averagepressureq,0) = average;
    for (int ipr = 0; ipr < npres; ipr++) {
        solution(pressureqs[ipr],0) = average_pressure[ipr]-
            fMixedDensity*gravity_pressure[ipr]+
            flux_pressure[ipr]/fLambda;
    }
}
