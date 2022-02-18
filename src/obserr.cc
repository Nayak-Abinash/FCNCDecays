#include "obserr.h"

//Bd->Kstr,ll:
double BdtoKstrll_obserr::J1s(double qsq, double ml, double unv[]){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml,unv)*AaRIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*AaRRe(qsq,ml,unv) + ApLIm(qsq,ml,unv)*ApRIm(qsq,ml,unv) +
                             ApLRe(qsq,ml,unv)*ApRRe(qsq,ml,unv)))/qsq +
    ((pow(AaLIm(qsq,ml,unv),2.0) + pow(AaLRe(qsq,ml,unv),2.0) + pow(AaRIm(qsq,ml,unv),2.0) + pow(AaRRe(qsq,ml,unv),2.0) +
      pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double BdtoKstrll_obserr::J1c(double qsq, double ml, double unv[]){
    return pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml,unv),2.0) + pow(AtRe(qsq,ml,unv),2.0) +
                      2.0*(AzLIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) + AzLRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))))/qsq +
    (pow(ASIm(qsq,ml,unv),2.0) + pow(ASRe(qsq,ml,unv),2.0))*pow(betal(qsq,ml),2.0) ;}

double BdtoKstrll_obserr::J2s(double qsq, double ml, double unv[]){
    return ((pow(AaLIm(qsq,ml,unv),2.0) + pow(AaLRe(qsq,ml,unv),2.0) + pow(AaRIm(qsq,ml,unv),2.0) + pow(AaRRe(qsq,ml,unv),2.0) +
             pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double BdtoKstrll_obserr::J2c(double qsq, double ml, double unv[]){
    return -((pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double BdtoKstrll_obserr::J3(double qsq, double ml, double unv[]){
    return ((-pow(AaLIm(qsq,ml,unv),2.0) - pow(AaLRe(qsq,ml,unv),2.0) - pow(AaRIm(qsq,ml,unv),2.0) - pow(AaRRe(qsq,ml,unv),2.0) +
             pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double BdtoKstrll_obserr::J4(double qsq, double ml, double unv[]){
    return ((AaLIm(qsq,ml,unv)*AzLIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*AzLRe(qsq,ml,unv) + AaRIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) +
             AaRRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obserr::J5(double qsq, double ml, double unv[]){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml,unv)*ASIm(qsq,ml,unv) + AaRIm(qsq,ml,unv)*ASIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*ASRe(qsq,ml,unv)
                             + AaRRe(qsq,ml,unv)*ASRe(qsq,ml,unv)))/sqrt(qsq)) + ApLIm(qsq,ml,unv)*AzLIm(qsq,ml,unv) +
                      ApLRe(qsq,ml,unv)*AzLRe(qsq,ml,unv) - ApRIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) - ApRRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))*betal(qsq,ml) ;}

double BdtoKstrll_obserr::J6s(double qsq, double ml, double unv[]){
    return 2.0*(AaLIm(qsq,ml,unv)*ApLIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*ApLRe(qsq,ml,unv) - AaRIm(qsq,ml,unv)*ApRIm(qsq,ml,unv) -
                AaRRe(qsq,ml,unv)*ApRRe(qsq,ml,unv))*betal(qsq,ml) ;}

double BdtoKstrll_obserr::J6c(double qsq, double ml, double unv[]){
    return (4.0*ml*(ASIm(qsq,ml,unv)*AzLIm(qsq,ml,unv) + ASRe(qsq,ml,unv)*AzLRe(qsq,ml,unv) + ASIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) +
                    ASRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))*betal(qsq,ml))/sqrt(qsq) ; }

double BdtoKstrll_obserr::J7(double qsq, double ml, double unv[]){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml,unv)*ASIm(qsq,ml,unv)) - ApRRe(qsq,ml,unv)*ASIm(qsq,ml,unv) + ApLIm(qsq,ml,unv)*ASRe(qsq,ml,unv)
                           + ApRIm(qsq,ml,unv)*ASRe(qsq,ml,unv)))/sqrt(qsq) + AaLRe(qsq,ml,unv)*AzLIm(qsq,ml,unv) -
                      AaLIm(qsq,ml,unv)*AzLRe(qsq,ml,unv) - AaRRe(qsq,ml,unv)*AzRIm(qsq,ml,unv) + AaRIm(qsq,ml,unv)*AzRRe(qsq,ml,unv))*betal(qsq,ml) ;}

double BdtoKstrll_obserr::J8(double qsq, double ml, double unv[]){
    return ((ApLRe(qsq,ml,unv)*AzLIm(qsq,ml,unv) - ApLIm(qsq,ml,unv)*AzLRe(qsq,ml,unv) + ApRRe(qsq,ml,unv)*AzRIm(qsq,ml,unv) -
             ApRIm(qsq,ml,unv)*AzRRe(qsq,ml,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obserr::J9(double qsq, double ml, double unv[]){
    return (AaLRe(qsq,ml,unv)*ApLIm(qsq,ml,unv) - AaLIm(qsq,ml,unv)*ApLRe(qsq,ml,unv) + AaRRe(qsq,ml,unv)*ApRIm(qsq,ml,unv) -
            AaRIm(qsq,ml,unv)*ApRRe(qsq,ml,unv))*pow(betal(qsq,ml),2.0) ;}

double BdtoKstrll_obserr::diffWidth(double qsq, double ml, double unv[]){
    return (3.0*(J1c(qsq,ml,unv) + 2.0*J1s(qsq,ml,unv) + (-J2c(qsq,ml,unv) - 2.0*J2s(qsq,ml,unv))/3.0))/4.0;}

double BdtoKstrll_obserr::FL(double qsq, double ml, double unv[]){
    return (pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0))/
    (pow(AaLIm(qsq,ml,unv),2.0) + pow(AaLRe(qsq,ml,unv),2.0) + pow(AaRIm(qsq,ml,unv),2.0) + pow(AaRRe(qsq,ml,unv),2.0) +
     pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0) +
     pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0));}

double BdtoKstrll_obserr::AFB(double qsq, double ml, double unv[]){
    return (3.0*(J6c(qsq,ml,unv) + 2.0*J6s(qsq,ml,unv)))/(8.0*diffWidth(qsq,ml,unv));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll:
double Bstophill_obserr::J1s(double qsq, double ml, double unv[]){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml,unv)*AaRIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*AaRRe(qsq,ml,unv) + ApLIm(qsq,ml,unv)*ApRIm(qsq,ml,unv)
                             + ApLRe(qsq,ml,unv)*ApRRe(qsq,ml,unv)))/qsq +
    ((pow(AaLIm(qsq,ml,unv),2.0) + pow(AaLRe(qsq,ml,unv),2.0) + pow(AaRIm(qsq,ml,unv),2.0) + pow(AaRRe(qsq,ml,unv),2.0) +
      pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double Bstophill_obserr::J1c(double qsq, double ml, double unv[]){
    return pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml,unv),2.0) + pow(AtRe(qsq,ml,unv),2.0) +
                      2.0*(AzLIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) + AzLRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))))/qsq +
    (pow(ASIm(qsq,ml,unv),2.0) + pow(ASRe(qsq,ml,unv),2.0))*pow(betal(qsq,ml),2.0) ;}

double Bstophill_obserr::J2s(double qsq, double ml, double unv[]){
    return ((pow(AaLIm(qsq,ml,unv),2.0) + pow(AaLRe(qsq,ml,unv),2.0) + pow(AaRIm(qsq,ml,unv),2.0) + pow(AaRRe(qsq,ml,unv),2.0) +
             pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double Bstophill_obserr::J2c(double qsq, double ml, double unv[]){
    return -((pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double Bstophill_obserr::J3(double qsq, double ml, double unv[]){
    return ((-pow(AaLIm(qsq,ml,unv),2.0) - pow(AaLRe(qsq,ml,unv),2.0) - pow(AaRIm(qsq,ml,unv),2.0) - pow(AaRRe(qsq,ml,unv),2.0) +
             pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double Bstophill_obserr::J4(double qsq, double ml, double unv[]){
    return ((AaLIm(qsq,ml,unv)*AzLIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*AzLRe(qsq,ml,unv) + AaRIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) +
             AaRRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obserr::J5(double qsq, double ml, double unv[]){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml,unv)*ASIm(qsq,ml,unv) + AaRIm(qsq,ml,unv)*ASIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*ASRe(qsq,ml,unv)
                             + AaRRe(qsq,ml,unv)*ASRe(qsq,ml,unv)))/sqrt(qsq)) + ApLIm(qsq,ml,unv)*AzLIm(qsq,ml,unv) +
        ApLRe(qsq,ml,unv)*AzLRe(qsq,ml,unv) - ApRIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) - ApRRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))*betal(qsq,ml) ;}

double Bstophill_obserr::J6s(double qsq, double ml, double unv[]){
    return 2.0*(AaLIm(qsq,ml,unv)*ApLIm(qsq,ml,unv) + AaLRe(qsq,ml,unv)*ApLRe(qsq,ml,unv) - AaRIm(qsq,ml,unv)*ApRIm(qsq,ml,unv) -
                AaRRe(qsq,ml,unv)*ApRRe(qsq,ml,unv))*betal(qsq,ml) ;}

double Bstophill_obserr::J6c(double qsq, double ml, double unv[]){
    return (4.0*ml*(ASIm(qsq,ml,unv)*AzLIm(qsq,ml,unv) + ASRe(qsq,ml,unv)*AzLRe(qsq,ml,unv) + ASIm(qsq,ml,unv)*AzRIm(qsq,ml,unv) +
                    ASRe(qsq,ml,unv)*AzRRe(qsq,ml,unv))*betal(qsq,ml))/sqrt(qsq) ; }

double Bstophill_obserr::J7(double qsq, double ml, double unv[]){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml,unv)*ASIm(qsq,ml,unv)) - ApRRe(qsq,ml,unv)*ASIm(qsq,ml,unv) + ApLIm(qsq,ml,unv)*ASRe(qsq,ml,unv)
                           + ApRIm(qsq,ml,unv)*ASRe(qsq,ml,unv)))/sqrt(qsq) + AaLRe(qsq,ml,unv)*AzLIm(qsq,ml,unv) -
        AaLIm(qsq,ml,unv)*AzLRe(qsq,ml,unv) - AaRRe(qsq,ml,unv)*AzRIm(qsq,ml,unv) + AaRIm(qsq,ml,unv)*AzRRe(qsq,ml,unv))*betal(qsq,ml) ;}

double Bstophill_obserr::J8(double qsq, double ml, double unv[]){
    return ((ApLRe(qsq,ml,unv)*AzLIm(qsq,ml,unv) - ApLIm(qsq,ml,unv)*AzLRe(qsq,ml,unv) + ApRRe(qsq,ml,unv)*AzRIm(qsq,ml,unv) -
             ApRIm(qsq,ml,unv)*AzRRe(qsq,ml,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obserr::J9(double qsq, double ml, double unv[]){
    return (AaLRe(qsq,ml,unv)*ApLIm(qsq,ml,unv) - AaLIm(qsq,ml,unv)*ApLRe(qsq,ml,unv) + AaRRe(qsq,ml,unv)*ApRIm(qsq,ml,unv) -
            AaRIm(qsq,ml,unv)*ApRRe(qsq,ml,unv))*pow(betal(qsq,ml),2.0) ;}

double Bstophill_obserr::diffWidth(double qsq, double ml, double unv[]){
    return (3.0*(J1c(qsq,ml,unv) + 2.0*J1s(qsq,ml,unv) + (-J2c(qsq,ml,unv) - 2.0*J2s(qsq,ml,unv))/3.0))/4.0;}

double Bstophill_obserr::FL(double qsq, double ml, double unv[]){
    return (pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0))/
    (pow(AaLIm(qsq,ml,unv),2.0) + pow(AaLRe(qsq,ml,unv),2.0) + pow(AaRIm(qsq,ml,unv),2.0) + pow(AaRRe(qsq,ml,unv),2.0) +
     pow(ApLIm(qsq,ml,unv),2.0) + pow(ApLRe(qsq,ml,unv),2.0) + pow(ApRIm(qsq,ml,unv),2.0) + pow(ApRRe(qsq,ml,unv),2.0) +
     pow(AzLIm(qsq,ml,unv),2.0) + pow(AzLRe(qsq,ml,unv),2.0) + pow(AzRIm(qsq,ml,unv),2.0) + pow(AzRRe(qsq,ml,unv),2.0));}

double Bstophill_obserr::AFB(double qsq, double ml, double unv[]){
    return (3.0*(J6c(qsq,ml,unv) + 2.0*J6s(qsq,ml,unv)))/(8.0*diffWidth(qsq,ml,unv));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_obserr::alNP(double qsq, double ml, double unv[]){
    return (4.0*pow(mBd(),2.0)*pow(ml,2.0)*(pow(FAIm(qsq,ml,unv),2.0) + pow(FARe(qsq,ml,unv),2.0)) +
            2.0*ml*(pow(mBd(),2.0) - pow(mK(),2.0) + qsq)*(FAIm(qsq,ml,unv)*FPIm(qsq,ml,unv) + FARe(qsq,ml,unv)*FPRe(qsq,ml,unv)) +
            qsq*(pow(FPIm(qsq,ml,unv),2.0) + pow(FPRe(qsq,ml,unv),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FSIm(qsq,ml,unv),2.0) + pow(FSRe(qsq,ml,unv),2.0))) +
            ((pow(FAIm(qsq,ml,unv),2.0) + pow(FARe(qsq,ml,unv),2.0) + pow(FVIm(qsq,ml,unv),2.0) + pow(FVRe(qsq,ml,unv),2.0))*
             lambda(qsq))/4.0)*nf(qsq,ml)*pow(xiP(qsq,unv),2.0);}

double BdtoKll_obserr::blNP(double qsq, double ml, double unv[]){
    return 2.0*(qsq*(FPIm(qsq,ml,unv)*FT5Im(qsq,ml,unv) + FPRe(qsq,ml,unv)*FT5Re(qsq,ml,unv) +
                     pow(betal(qsq,ml),2.0)*(FSIm(qsq,ml,unv)*FTIm(qsq,ml,unv) + FSRe(qsq,ml,unv)*FTRe(qsq,ml,unv))) +
                ml*((pow(mBd(),2.0) - pow(mK(),2.0) + qsq)*(FAIm(qsq,ml,unv)*FT5Im(qsq,ml,unv) + FARe(qsq,ml,unv)*FT5Re(qsq,ml,unv)) +
                    betal(qsq,ml)*(FSIm(qsq,ml,unv)*FVIm(qsq,ml,unv) + FSRe(qsq,ml,unv)*FVRe(qsq,ml,unv))*sqrt(lambda(qsq))))*
    nf(qsq,ml)*pow(xiP(qsq,unv),2.0);}

double BdtoKll_obserr::clNP(double qsq, double ml, double unv[]){
    return (qsq*(pow(FT5Im(qsq,ml,unv),2.0) + pow(FT5Re(qsq,ml,unv),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FTIm(qsq,ml,unv),2.0) + pow(FTRe(qsq,ml,unv),2.0))) +
            2.0*ml*betal(qsq,ml)*(FTIm(qsq,ml,unv)*FVIm(qsq,ml,unv) + FTRe(qsq,ml,unv)*FVRe(qsq,ml,unv))*sqrt(lambda(qsq)) -
            (pow(betal(qsq,ml),2.0)*(pow(FAIm(qsq,ml,unv),2.0) + pow(FARe(qsq,ml,unv),2.0) + pow(FVIm(qsq,ml,unv),2.0) +
                                     pow(FVRe(qsq,ml,unv),2.0))*lambda(qsq))/4.0)*nf(qsq,ml)*pow(xiP(qsq,unv),2.0);}

double BdtoKll_obserr::diffWidth(double qsq, double ml, double unv[]){
    return 2.0*(alNP(qsq,ml,unv) + 1.0/3.0*clNP(qsq,ml,unv)); }

double BdtoKll_obserr::diffBrnch(double qsq, double ml, double unv[]){
    return tauBd()/2.0*diffWidth(qsq,ml,unv); }

double BdtoKll_obserr::diffAFB(double qsq, double ml, double unv[]){
    return blNP(qsq,ml,unv)/diffWidth(qsq,ml,unv); }

double BdtoKll_obserr::diffFH(double qsq, double ml, double unv[]){
    return 2.0*(alNP(qsq,ml,unv)+clNP(qsq,ml,unv))/diffWidth(qsq,ml,unv); }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
//B->ll
double Btoll_obserr :: Btoll_nf(double mBq){
    if (mBq==mBd()){
        return pow(GF()*alpha_e()*fBd()*absVtbVtdStr(),2.0)*tauBd()/(64.0*pow(mBd()*pi,3.0)); }
    if (mBq==mBs()){
        return pow(GF()*alpha_e()*fBs()*absVtbVtsStr(),2.0)*tauBs()/(64.0*pow(mBs()*pi,3.0)); }
    return 0.0;}

double Btoll_obserr :: ADeltaGammaf(double mBq, double mq, double ml1, double ml2){
    return (-pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0) + pow(ampSIm(mBq,mq,ml1,ml2),2.0) -
            pow(ampSRe(mBq,mq,ml1,ml2),2.0))/
    (pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0) + pow(ampSIm(mBq,mq,ml1,ml2),2.0) +
     pow(ampSRe(mBq,mq,ml1,ml2),2.0));}

double Btoll_obserr :: CorrctnFctr(double mBq, double mq, double ml1, double ml2){
    return (1.0 - pow(y(mBq),2.0))/(1.0 + ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq));}

double Btoll_obserr :: BrInst(double mBq, double mq, double ml1, double ml2){
    return ((pow(mBq,2.0) - pow(ml1 - ml2,2.0))*(pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0)) + (pow(mBq,2.0) - pow(ml1 + ml2,2.0))*(pow(ampSIm(mBq,mq,ml1,ml2),2.0) +
        pow(ampSRe(mBq,mq,ml1,ml2),2.0)))*sqrt(Btoll_lambda(mBq,ml1,ml2))*Btoll_nf(mBq);}

double Btoll_obserr :: BrTimeIntgratd(double mBq, double mq, double ml1, double ml2){
    return BrInst(mBq,mq,ml1,ml2)/CorrctnFctr(mBq,mq,ml1,ml2);}

double Btoll_obserr :: efftau(double mBq, double mq, double ml1, double ml2){
    return (tauB(mBq)*(1.0 + 2.0*ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq) + pow(y(mBq),2.0)))/
    ((1.0 + ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq))*(1.0 - pow(y(mBq),2.0)));}
*/











