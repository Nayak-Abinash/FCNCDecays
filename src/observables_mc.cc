#include "observables_mc.h"

//Bd->Kstr,ll:
double BdtoKstrll_obserr::J1s(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml,smwc,npwc,unv)*AaRIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*AaRRe(qsq,ml,smwc,npwc,unv) + ApLIm(qsq,ml,smwc,npwc,unv)*ApRIm(qsq,ml,smwc,npwc,unv) +
                             ApLRe(qsq,ml,smwc,npwc,unv)*ApRRe(qsq,ml,smwc,npwc,unv)))/qsq +
    ((pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
      pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double BdtoKstrll_obserr::J1c(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AtRe(qsq,ml,smwc,npwc,unv),2.0) +
                      2.0*(AzLIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) + AzLRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))))/qsq +
    (pow(ASIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ASRe(qsq,ml,smwc,npwc,unv),2.0))*pow(betal(qsq,ml),2.0) ;}

double BdtoKstrll_obserr::J2s(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
             pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double BdtoKstrll_obserr::J2c(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -((pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double BdtoKstrll_obserr::J3(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((-pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) - pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) - pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) - pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
             pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double BdtoKstrll_obserr::J4(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((AaLIm(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) + AaRIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) +
             AaRRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obserr::J5(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv) + AaRIm(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)
                             + AaRRe(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)))/sqrt(qsq)) + ApLIm(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) +
                      ApLRe(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) - ApRIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) - ApRRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml) ;}

double BdtoKstrll_obserr::J6s(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*(AaLIm(qsq,ml,smwc,npwc,unv)*ApLIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*ApLRe(qsq,ml,smwc,npwc,unv) - AaRIm(qsq,ml,smwc,npwc,unv)*ApRIm(qsq,ml,smwc,npwc,unv) -
                AaRRe(qsq,ml,smwc,npwc,unv)*ApRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml) ;}

double BdtoKstrll_obserr::J6c(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (4.0*ml*(ASIm(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) + ASRe(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) + ASIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) +
                    ASRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml))/sqrt(qsq) ; }

double BdtoKstrll_obserr::J7(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv)) - ApRRe(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv) + ApLIm(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)
                           + ApRIm(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)))/sqrt(qsq) + AaLRe(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) -
                      AaLIm(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) - AaRRe(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) + AaRIm(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml) ;}

double BdtoKstrll_obserr::J8(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((ApLRe(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) - ApLIm(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) + ApRRe(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) -
             ApRIm(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obserr::J9(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (AaLRe(qsq,ml,smwc,npwc,unv)*ApLIm(qsq,ml,smwc,npwc,unv) - AaLIm(qsq,ml,smwc,npwc,unv)*ApLRe(qsq,ml,smwc,npwc,unv) + AaRRe(qsq,ml,smwc,npwc,unv)*ApRIm(qsq,ml,smwc,npwc,unv) -
            AaRIm(qsq,ml,smwc,npwc,unv)*ApRRe(qsq,ml,smwc,npwc,unv))*pow(betal(qsq,ml),2.0) ;}

double BdtoKstrll_obserr::diffWidth(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (3.0*(J1c(qsq,ml,smwc,npwc,unv) + 2.0*J1s(qsq,ml,smwc,npwc,unv) + (-J2c(qsq,ml,smwc,npwc,unv) - 2.0*J2s(qsq,ml,smwc,npwc,unv))/3.0))/4.0;}

double BdtoKstrll_obserr::FL(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0))/
    (pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
     pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0) +
     pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0));}

double BdtoKstrll_obserr::AFB(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (3.0*(J6c(qsq,ml,smwc,npwc,unv) + 2.0*J6s(qsq,ml,smwc,npwc,unv)))/(8.0*diffWidth(qsq,ml,smwc,npwc,unv));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll:
double Bstophill_obserr::J1s(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml,smwc,npwc,unv)*AaRIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*AaRRe(qsq,ml,smwc,npwc,unv) + ApLIm(qsq,ml,smwc,npwc,unv)*ApRIm(qsq,ml,smwc,npwc,unv)
                             + ApLRe(qsq,ml,smwc,npwc,unv)*ApRRe(qsq,ml,smwc,npwc,unv)))/qsq +
    ((pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
      pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double Bstophill_obserr::J1c(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AtRe(qsq,ml,smwc,npwc,unv),2.0) +
                      2.0*(AzLIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) + AzLRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))))/qsq +
    (pow(ASIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ASRe(qsq,ml,smwc,npwc,unv),2.0))*pow(betal(qsq,ml),2.0) ;}

double Bstophill_obserr::J2s(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
             pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double Bstophill_obserr::J2c(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -((pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double Bstophill_obserr::J3(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((-pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) - pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) - pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) - pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
             pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double Bstophill_obserr::J4(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((AaLIm(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) + AaRIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) +
             AaRRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obserr::J5(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv) + AaRIm(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)
                             + AaRRe(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)))/sqrt(qsq)) + ApLIm(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) +
        ApLRe(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) - ApRIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) - ApRRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml) ;}

double Bstophill_obserr::J6s(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*(AaLIm(qsq,ml,smwc,npwc,unv)*ApLIm(qsq,ml,smwc,npwc,unv) + AaLRe(qsq,ml,smwc,npwc,unv)*ApLRe(qsq,ml,smwc,npwc,unv) - AaRIm(qsq,ml,smwc,npwc,unv)*ApRIm(qsq,ml,smwc,npwc,unv) -
                AaRRe(qsq,ml,smwc,npwc,unv)*ApRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml) ;}

double Bstophill_obserr::J6c(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (4.0*ml*(ASIm(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) + ASRe(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) + ASIm(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) +
                    ASRe(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml))/sqrt(qsq) ; }

double Bstophill_obserr::J7(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv)) - ApRRe(qsq,ml,smwc,npwc,unv)*ASIm(qsq,ml,smwc,npwc,unv) + ApLIm(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)
                           + ApRIm(qsq,ml,smwc,npwc,unv)*ASRe(qsq,ml,smwc,npwc,unv)))/sqrt(qsq) + AaLRe(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) -
        AaLIm(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) - AaRRe(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) + AaRIm(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*betal(qsq,ml) ;}

double Bstophill_obserr::J8(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return ((ApLRe(qsq,ml,smwc,npwc,unv)*AzLIm(qsq,ml,smwc,npwc,unv) - ApLIm(qsq,ml,smwc,npwc,unv)*AzLRe(qsq,ml,smwc,npwc,unv) + ApRRe(qsq,ml,smwc,npwc,unv)*AzRIm(qsq,ml,smwc,npwc,unv) -
             ApRIm(qsq,ml,smwc,npwc,unv)*AzRRe(qsq,ml,smwc,npwc,unv))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obserr::J9(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (AaLRe(qsq,ml,smwc,npwc,unv)*ApLIm(qsq,ml,smwc,npwc,unv) - AaLIm(qsq,ml,smwc,npwc,unv)*ApLRe(qsq,ml,smwc,npwc,unv) + AaRRe(qsq,ml,smwc,npwc,unv)*ApRIm(qsq,ml,smwc,npwc,unv) -
            AaRIm(qsq,ml,smwc,npwc,unv)*ApRRe(qsq,ml,smwc,npwc,unv))*pow(betal(qsq,ml),2.0) ;}

double Bstophill_obserr::diffWidth(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (3.0*(J1c(qsq,ml,smwc,npwc,unv) + 2.0*J1s(qsq,ml,smwc,npwc,unv) + (-J2c(qsq,ml,smwc,npwc,unv) - 2.0*J2s(qsq,ml,smwc,npwc,unv))/3.0))/4.0;}

double Bstophill_obserr::FL(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0))/
    (pow(AaLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AaRRe(qsq,ml,smwc,npwc,unv),2.0) +
     pow(ApLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(ApRRe(qsq,ml,smwc,npwc,unv),2.0) +
     pow(AzLIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzLRe(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRIm(qsq,ml,smwc,npwc,unv),2.0) + pow(AzRRe(qsq,ml,smwc,npwc,unv),2.0));}

double Bstophill_obserr::AFB(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (3.0*(J6c(qsq,ml,smwc,npwc,unv) + 2.0*J6s(qsq,ml,smwc,npwc,unv)))/(8.0*diffWidth(qsq,ml,smwc,npwc,unv));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_obserr::alNP(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (4.0*pow(mBd(unv),2.0)*pow(ml,2.0)*(pow(FAIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FARe(qsq,ml,smwc,npwc,unv),2.0)) +
            2.0*ml*(pow(mBd(unv),2.0) - pow(mK(unv),2.0) + qsq)*(FAIm(qsq,ml,smwc,npwc,unv)*FPIm(qsq,ml,smwc,npwc,unv) + FARe(qsq,ml,smwc,npwc,unv)*FPRe(qsq,ml,smwc,npwc,unv)) +
            qsq*(pow(FPIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FPRe(qsq,ml,smwc,npwc,unv),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FSIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FSRe(qsq,ml,smwc,npwc,unv),2.0))) +
            ((pow(FAIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FARe(qsq,ml,smwc,npwc,unv),2.0) + pow(FVIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FVRe(qsq,ml,smwc,npwc,unv),2.0))*
             lambda(qsq,unv))/4.0)*nf(qsq,ml,unv)*pow(xiP(qsq,smwc,npwc,unv),2.0);}

double BdtoKll_obserr::blNP(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*(qsq*(FPIm(qsq,ml,smwc,npwc,unv)*FT5Im(qsq,ml,smwc,npwc,unv) + FPRe(qsq,ml,smwc,npwc,unv)*FT5Re(qsq,ml,smwc,npwc,unv) +
                     pow(betal(qsq,ml),2.0)*(FSIm(qsq,ml,smwc,npwc,unv)*FTIm(qsq,ml,smwc,npwc,unv) + FSRe(qsq,ml,smwc,npwc,unv)*FTRe(qsq,ml,smwc,npwc,unv))) +
                ml*((pow(mBd(unv),2.0) - pow(mK(unv),2.0) + qsq)*(FAIm(qsq,ml,smwc,npwc,unv)*FT5Im(qsq,ml,smwc,npwc,unv) + FARe(qsq,ml,smwc,npwc,unv)*FT5Re(qsq,ml,smwc,npwc,unv)) +
                    betal(qsq,ml)*(FSIm(qsq,ml,smwc,npwc,unv)*FVIm(qsq,ml,smwc,npwc,unv) + FSRe(qsq,ml,smwc,npwc,unv)*FVRe(qsq,ml,smwc,npwc,unv))*sqrt(lambda(qsq,unv))))*
    nf(qsq,ml,unv)*pow(xiP(qsq,smwc,npwc,unv),2.0);}

double BdtoKll_obserr::clNP(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (qsq*(pow(FT5Im(qsq,ml,smwc,npwc,unv),2.0) + pow(FT5Re(qsq,ml,smwc,npwc,unv),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FTIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FTRe(qsq,ml,smwc,npwc,unv),2.0))) +
            2.0*ml*betal(qsq,ml)*(FTIm(qsq,ml,smwc,npwc,unv)*FVIm(qsq,ml,smwc,npwc,unv) + FTRe(qsq,ml,smwc,npwc,unv)*FVRe(qsq,ml,smwc,npwc,unv))*sqrt(lambda(qsq,unv)) -
            (pow(betal(qsq,ml),2.0)*(pow(FAIm(qsq,ml,smwc,npwc,unv),2.0) + pow(FARe(qsq,ml,smwc,npwc,unv),2.0) + pow(FVIm(qsq,ml,smwc,npwc,unv),2.0) +
                                     pow(FVRe(qsq,ml,smwc,npwc,unv),2.0))*lambda(qsq,unv))/4.0)*nf(qsq,ml,unv)*pow(xiP(qsq,smwc,npwc,unv),2.0);}

double BdtoKll_obserr::diffWidth(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*(alNP(qsq,ml,smwc,npwc,unv) + 1.0/3.0*clNP(qsq,ml,smwc,npwc,unv)); }

double BdtoKll_obserr::diffBrnch(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return tauBd(unv)/2.0*diffWidth(qsq,ml,smwc,npwc,unv); }

double BdtoKll_obserr::diffAFB(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return blNP(qsq,ml,smwc,npwc,unv)/diffWidth(qsq,ml,smwc,npwc,unv); }

double BdtoKll_obserr::diffFH(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*(alNP(qsq,ml,smwc,npwc,unv)+clNP(qsq,ml,smwc,npwc,unv))/diffWidth(qsq,ml,smwc,npwc,unv); }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//B->ll
double Btoll_obserr :: Btoll_nf(double mBq, double unv[]){
    if (mBq==mBd(unv)){
        return pow(GF()*alpha_e()*fBd(unv)*absVtbVtdStr(unv),2.0)*tauBd(unv)/(64.0*pow(mBd(unv)*pi,3.0)); }
    if (mBq==mBs(unv)){
        return pow(GF()*alpha_e()*fBs(unv)*absVtbVtsStr(unv),2.0)*tauBs(unv)/(64.0*pow(mBs(unv)*pi,3.0)); }
    return 0.0;}

double Btoll_obserr :: ADeltaGammaf(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return (-pow(ampPIm(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) + pow(ampPRe(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) + pow(ampSIm(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) -
            pow(ampSRe(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0))/
    (pow(ampPIm(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) + pow(ampPRe(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) + pow(ampSIm(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) +
     pow(ampSRe(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0));}

double Btoll_obserr :: CorrctnFctr(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return (1.0 - pow(y(mBq,smwc,npwc,unv),2.0))/(1.0 + ADeltaGammaf(mBq,mq,ml1,ml2,smwc,npwc,unv)*y(mBq,smwc,npwc,unv));}

double Btoll_obserr :: BrInst(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return ((pow(mBq,2.0) - pow(ml1 - ml2,2.0))*(pow(ampPIm(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0) + pow(ampPRe(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0)) + (pow(mBq,2.0)
            - pow(ml1 + ml2,2.0))*(pow(ampSIm(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0)
            + pow(ampSRe(mBq,mq,ml1,ml2,smwc,npwc,unv),2.0)))*sqrt(Btoll_lambda(mBq,ml1,ml2,smwc,npwc,unv))*Btoll_nf(mBq,unv);}

double Btoll_obserr :: BrTimeIntgratd(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return BrInst(mBq,mq,ml1,ml2,smwc,npwc,unv)/CorrctnFctr(mBq,mq,ml1,ml2,smwc,npwc,unv);}

double Btoll_obserr :: efftau(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return (tauB(mBq,smwc,npwc,unv)*(1.0 + 2.0*ADeltaGammaf(mBq,mq,ml1,ml2,smwc,npwc,unv)*y(mBq,smwc,npwc,unv) + pow(y(mBq,smwc,npwc,unv),2.0)))/
    ((1.0 + ADeltaGammaf(mBq,mq,ml1,ml2,smwc,npwc,unv)*y(mBq,smwc,npwc,unv))*(1.0 - pow(y(mBq,smwc,npwc,unv),2.0)));}












