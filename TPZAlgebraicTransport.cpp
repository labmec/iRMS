//
//  AlgebraicTransport.cpp
//  ALL_BUILD
//
//  Created by Jose on 6/1/20.
//

#include "TPZAlgebraicTransport.h"
#include "pzelmat.h"
/// Default constructor
TPZAlgebraicTransport::TPZAlgebraicTransport(){
    
}

/// Copy constructor
TPZAlgebraicTransport::TPZAlgebraicTransport(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
}

/// Assignement constructor
const TPZAlgebraicTransport & TPZAlgebraicTransport::operator=(const TPZAlgebraicTransport & other){
    fNFluxCoefficients = other.fNFluxCoefficients;
    fNVolumesTransport = other.fNVolumesTransport;
    fCellsData = other.fCellsData;
    fInterfaceData = other.fInterfaceData;
    return *this;
}
TPZAlgebraicTransport::~TPZAlgebraicTransport(){
    
}



void TPZAlgebraicTransport::CalcLambdas(){
//    int nlamdas = fPressure.size();
//    flambdas.resize(nlamdas);
//    for (int ilambda =0; ilambda<nlamdas; ilambda++) {
//        flambdas[ilambda] = CalcLambda(fSaturation[ilambda], fPressure[ilambda], fDensityOil[ilambda], fDensityWater[ilambda] );
//    }
}
void TPZAlgebraicTransport::CalcDensities(){
//    int ndensities = fPressure.size();
//    fDensityOil.resize(ndensities);
//    fDensityWater.resize(ndensities);
//    double rhoo_ref = 760.00;
//    double rhow_ref = 1000.00;
//    double p_ref = 1.013E5;
//    for (int idensity = 0; idensity< ndensities; idensity++) {
//        fDensityWater[idensity] = CalcDensity(fPressure[idensity], fCompressibility[0],rhow_ref, p_ref);
//        fDensityOil[idensity] = CalcDensity(fPressure[idensity],fCompressibility[1], rhoo_ref, p_ref);
//    }
}

double TPZAlgebraicTransport::CalcLambda(double sw, double paverage, double densityo, double densityw){
    double kro = (1 - sw)*(1 - sw);
    double krw = sw*sw;
    double muo = 0.05;
    double muw = 0.01;
    double lambda;
    lambda = ((densityo*kro)/muo) + ((densityw*krw)/muw);
    
    return lambda;
}

double TPZAlgebraicTransport::CalcDensity(double paverage,double compress, double reff, double pref){
    double density = reff*(1+ compress*(paverage-pref));
    return density;
}



void TPZAlgebraicTransport::BuildDataStructures(TPZMultiphysicsCompMesh &transportmesh)
{
//    TPZMultiphysicsCompMesh *transp = &transportmesh;
//    int neltransport = transp->NElements();
//    for (int iel = 0; iel< neltransport; iel++) {
//        TPZCompEl * cel = transp->Element(iel);
//        if (!cel) {
//            continue;
//        }
//        if (cel->Dimension() != transp->Dimension()) {
//            continue;
//        }
//        TPZElementMatrix ek;
//        TPZElementMatrix ef;
//        cel->CalcStiff(ek, ef);
//        double vol_val = ek.fMat.Get(0, 0);
//        fVolumes.push_back(vol_val);
//    }
//    int nvols = fVolumes.size();
//
//    //these could be default values
//    fSaturation.resize(nvols,0.0);
//    fDensityOil.resize(nvols,800.00);
//    fDensityWater.resize(nvols,1000.00);
//    fPressure.resize(nvols, 1.03E5);
//    flambdas.resize(nvols,1.0);
//    fCompressibility.resize(2);
//    fCompressibility[0] = 0.0; //Water compressibility
//    fCompressibility[1] = 0.0; //Oil compressibility
}
void TPZAlgebraicTransport::Assamble(TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef ){
    if (fNVolumesTransport <= 0) {
        DebugStop();
    }
    ek.Resize(fNVolumesTransport, fNVolumesTransport);
    ek.Zero();
    std::cout<<"NRows: "<<ek.Rows()<<"NCols: "<<ek.Cols();
    this->Contribute(ek, ef);
    ek.Print(std::cout);
}
void TPZAlgebraicTransport::AssambleResidual(TPZFNMatrix<100, REAL> &ef){
    
}
void TPZAlgebraicTransport::Contribute(TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef){
    
    
    for (auto celDatabyID: fCellsData) {
        for (auto CelData : celDatabyID.second) {
            CelData.Print(std::cout);
            int index = CelData.index;
            if (CelData.fVolume<0) {
                continue;
            }
            REAL vol = CelData.fVolume;
            ek(index, index) = vol;
        }
    }
}
void TPZAlgebraicTransport::ContributeInterface(TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef){
    
}
void TPZAlgebraicTransport::ContributeBCInterface(TPZFNMatrix<100, REAL> &ek,TPZFNMatrix<100, REAL> &ef){
    
}

void TPZAlgebraicTransport::TCellData::Print(std::ostream &out){
    
    out<<"Element index: "<<this->index<<" : "<<std::endl;
    out<<"Volume = "<<this->fVolume<<std::endl;
    out<<"Pressure = "<<this->fPressure<<std::endl;
    out<<"DensityWater = "<<this->fDensityWater<<std::endl;
    out<<"DensityOil = "<<this->fDensityOil<<std::endl;
    out<<"Saturation = "<<this->fSaturation<<std::endl;
    out<<"Lambda = "<<this->flambda<<std::endl;
    out <<std::endl;
}
void TPZAlgebraicTransport::TInterfaceDataTransport::Print(std::ostream &out){
    
    out<<"Gel index: "<<this->gelIndex<<" : "<<std::endl;
    out << "fCoefficientsFlux :";
    for (auto const& value : this->fCoefficientsFlux) out << value << ' ';
    out << std::endl;
    out << "fIntegralFluxFunctions :";
    for (auto const& value : this->fIntegralFluxFunctions) out << value << ' ';
    out << std::endl;
    out << "fIntegralFlux :";
    for (auto const& value : this->fIntegralFlux) out << value << ' ';
    out << std::endl;
    out<<"fFluxSing = "<<this->fFluxSing<<std::endl;
    out << "fNormalFaceDirection :";
    for (auto const& value : this->fNormalFaceDirection) out << value << ' ';
    out << std::endl;
    out <<" " <<std::endl;
}
