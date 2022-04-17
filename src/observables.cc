#include "observables.h"

//Bd->Kstr,ll:
double BdtoKstrll_obs::J1s(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml)*AaRIm(qsq,ml) + AaLRe(qsq,ml)*AaRRe(qsq,ml) + ApLIm(qsq,ml)*ApRIm(qsq,ml) + ApLRe(qsq,ml)*ApRRe(qsq,ml)))/qsq +
    ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
      pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double BdtoKstrll_obs::J1c(double qsq, double ml){
    return pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml),2.0) + pow(AtRe(qsq,ml),2.0) +
                      2.0*(AzLIm(qsq,ml)*AzRIm(qsq,ml) + AzLRe(qsq,ml)*AzRRe(qsq,ml))))/qsq +
    (pow(ASIm(qsq,ml),2.0) + pow(ASRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ;}

double BdtoKstrll_obs::J1(double qsq, double ml){
    return 2.0*J1s(qsq,ml) + J1c(qsq,ml);}

double BdtoKstrll_obs::J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double BdtoKstrll_obs::J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double BdtoKstrll_obs::J2(double qsq, double ml){
    return 2.0*J2s(qsq,ml) + J2c(qsq,ml);}

double BdtoKstrll_obs::J3(double qsq, double ml){
    return ((-pow(AaLIm(qsq,ml),2.0) - pow(AaLRe(qsq,ml),2.0) - pow(AaRIm(qsq,ml),2.0) - pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double BdtoKstrll_obs::J4(double qsq, double ml){
    return ((AaLIm(qsq,ml)*AzLIm(qsq,ml) + AaLRe(qsq,ml)*AzLRe(qsq,ml) + AaRIm(qsq,ml)*AzRIm(qsq,ml) +
             AaRRe(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obs::J5(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml)*ASIm(qsq,ml) + AaRIm(qsq,ml)*ASIm(qsq,ml) + AaLRe(qsq,ml)*ASRe(qsq,ml) + AaRRe(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq)) + ApLIm(qsq,ml)*AzLIm(qsq,ml) +
                      ApLRe(qsq,ml)*AzLRe(qsq,ml) - ApRIm(qsq,ml)*AzRIm(qsq,ml) - ApRRe(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double BdtoKstrll_obs::J6s(double qsq, double ml){
    return 2.0*(AaLIm(qsq,ml)*ApLIm(qsq,ml) + AaLRe(qsq,ml)*ApLRe(qsq,ml) - AaRIm(qsq,ml)*ApRIm(qsq,ml) -
                AaRRe(qsq,ml)*ApRRe(qsq,ml))*betal(qsq,ml) ;}

double BdtoKstrll_obs::J6c(double qsq, double ml){
    return (4.0*ml*(ASIm(qsq,ml)*AzLIm(qsq,ml) + ASRe(qsq,ml)*AzLRe(qsq,ml) + ASIm(qsq,ml)*AzRIm(qsq,ml) +
                    ASRe(qsq,ml)*AzRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double BdtoKstrll_obs::J6(double qsq, double ml){
    return 2.0*J6s(qsq,ml) + J6c(qsq,ml);}

double BdtoKstrll_obs::J7(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml)*ASIm(qsq,ml)) - ApRRe(qsq,ml)*ASIm(qsq,ml) + ApLIm(qsq,ml)*ASRe(qsq,ml) + ApRIm(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq) + AaLRe(qsq,ml)*AzLIm(qsq,ml) -
                      AaLIm(qsq,ml)*AzLRe(qsq,ml) - AaRRe(qsq,ml)*AzRIm(qsq,ml) + AaRIm(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double BdtoKstrll_obs::J8(double qsq, double ml){
    return ((ApLRe(qsq,ml)*AzLIm(qsq,ml) - ApLIm(qsq,ml)*AzLRe(qsq,ml) + ApRRe(qsq,ml)*AzRIm(qsq,ml) -
             ApRIm(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obs::J9(double qsq, double ml){
    return (AaLRe(qsq,ml)*ApLIm(qsq,ml) - AaLIm(qsq,ml)*ApLRe(qsq,ml) + AaRRe(qsq,ml)*ApRIm(qsq,ml) -
            AaRIm(qsq,ml)*ApRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ;}

//Conjugate J's
double BdtoKstrll_obs::J1bs(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AabLIm(qsq,ml)*AabRIm(qsq,ml) + AabLRe(qsq,ml)*AabRRe(qsq,ml) + ApbLIm(qsq,ml)*ApbRIm(qsq,ml) + ApbLRe(qsq,ml)*ApbRRe(qsq,ml)))/qsq +
    ((pow(AabLIm(qsq,ml),2.0) + pow(AabLRe(qsq,ml),2.0) + pow(AabRIm(qsq,ml),2.0) + pow(AabRRe(qsq,ml),2.0) +
      pow(ApbLIm(qsq,ml),2.0) + pow(ApbLRe(qsq,ml),2.0) + pow(ApbRIm(qsq,ml),2.0) + pow(ApbRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double BdtoKstrll_obs::J1bc(double qsq, double ml){
    return pow(AzbLIm(qsq,ml),2.0) + pow(AzbLRe(qsq,ml),2.0) + pow(AzbRIm(qsq,ml),2.0) + pow(AzbRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtbIm(qsq,ml),2.0) + pow(AtbRe(qsq,ml),2.0) +
                      2.0*(AzbLIm(qsq,ml)*AzbRIm(qsq,ml) + AzbLRe(qsq,ml)*AzbRRe(qsq,ml))))/qsq +
    (pow(ASbIm(qsq,ml),2.0) + pow(ASbRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ;}

double BdtoKstrll_obs::J1b(double qsq, double ml){
    return 2.0*J1bs(qsq,ml) + J1bc(qsq,ml);}

double BdtoKstrll_obs::J2bs(double qsq, double ml){
    return ((pow(AabLIm(qsq,ml),2.0) + pow(AabLRe(qsq,ml),2.0) + pow(AabRIm(qsq,ml),2.0) + pow(AabRRe(qsq,ml),2.0) +
             pow(ApbLIm(qsq,ml),2.0) + pow(ApbLRe(qsq,ml),2.0) + pow(ApbRIm(qsq,ml),2.0) + pow(ApbRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double BdtoKstrll_obs::J2bc(double qsq, double ml){
    return -((pow(AzbLIm(qsq,ml),2.0) + pow(AzbLRe(qsq,ml),2.0) + pow(AzbRIm(qsq,ml),2.0) + pow(AzbRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double BdtoKstrll_obs::J2b(double qsq, double ml){
    return 2.0*J2bs(qsq,ml) + J2bc(qsq,ml);}

double BdtoKstrll_obs::J3b(double qsq, double ml){
    return ((-pow(AabLIm(qsq,ml),2.0) - pow(AabLRe(qsq,ml),2.0) - pow(AabRIm(qsq,ml),2.0) - pow(AabRRe(qsq,ml),2.0) +
             pow(ApbLIm(qsq,ml),2.0) + pow(ApbLRe(qsq,ml),2.0) + pow(ApbRIm(qsq,ml),2.0) + pow(ApbRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double BdtoKstrll_obs::J4b(double qsq, double ml){
    return ((AabLIm(qsq,ml)*AzbLIm(qsq,ml) + AabLRe(qsq,ml)*AzbLRe(qsq,ml) + AabRIm(qsq,ml)*AzbRIm(qsq,ml) +
             AabRRe(qsq,ml)*AzbRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obs::J5b(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AabLIm(qsq,ml)*ASbIm(qsq,ml) + AabRIm(qsq,ml)*ASbIm(qsq,ml) + AabLRe(qsq,ml)*ASbRe(qsq,ml)
                             + AabRRe(qsq,ml)*ASbRe(qsq,ml)))/sqrt(qsq)) + ApbLIm(qsq,ml)*AzbLIm(qsq,ml) +
                      ApbLRe(qsq,ml)*AzbLRe(qsq,ml) - ApbRIm(qsq,ml)*AzbRIm(qsq,ml) - ApbRRe(qsq,ml)*AzbRRe(qsq,ml))*betal(qsq,ml) ;}

double BdtoKstrll_obs::J6bs(double qsq, double ml){
    return 2.0*(AabLIm(qsq,ml)*ApbLIm(qsq,ml) + AabLRe(qsq,ml)*ApbLRe(qsq,ml) - AabRIm(qsq,ml)*ApbRIm(qsq,ml) -
                AabRRe(qsq,ml)*ApbRRe(qsq,ml))*betal(qsq,ml) ;}

double BdtoKstrll_obs::J6bc(double qsq, double ml){
    return (4.0*ml*(ASbIm(qsq,ml)*AzbLIm(qsq,ml) + ASbRe(qsq,ml)*AzbLRe(qsq,ml) + ASbIm(qsq,ml)*AzbRIm(qsq,ml) +
                    ASbRe(qsq,ml)*AzbRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double BdtoKstrll_obs::J6b(double qsq, double ml){
    return 2.0*J6bs(qsq,ml) + J6bc(qsq,ml);}

double BdtoKstrll_obs::J7b(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApbLRe(qsq,ml)*ASbIm(qsq,ml)) - ApbRRe(qsq,ml)*ASbIm(qsq,ml) + ApbLIm(qsq,ml)*ASbRe(qsq,ml)
                           + ApbRIm(qsq,ml)*ASbRe(qsq,ml)))/sqrt(qsq) + AabLRe(qsq,ml)*AzbLIm(qsq,ml) -
                      AabLIm(qsq,ml)*AzbLRe(qsq,ml) - AabRRe(qsq,ml)*AzbRIm(qsq,ml) + AabRIm(qsq,ml)*AzbRRe(qsq,ml))*betal(qsq,ml) ;}

double BdtoKstrll_obs::J8b(double qsq, double ml){
    return ((ApbLRe(qsq,ml)*AzbLIm(qsq,ml) - ApbLIm(qsq,ml)*AzbLRe(qsq,ml) + ApbRRe(qsq,ml)*AzbRIm(qsq,ml) -
             ApbRIm(qsq,ml)*AzbRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double BdtoKstrll_obs::J9b(double qsq, double ml){
    return (AabLRe(qsq,ml)*ApbLIm(qsq,ml) - AabLIm(qsq,ml)*ApbLRe(qsq,ml) + AabRRe(qsq,ml)*ApbRIm(qsq,ml) -
            AabRIm(qsq,ml)*ApbRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ;}


double BdtoKstrll_obs::diffWidth(double qsq, double ml){
    return 3.0/4.0*(J1(qsq,ml) - J2(qsq,ml)/3.0);}

double BdtoKstrll_obs::diffWidthConj(double qsq, double ml){
    return 3.0/4.0*(J1b(qsq,ml) - J2b(qsq,ml)/3.0);}

double BdtoKstrll_obs::FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0));}

double BdtoKstrll_obs::AFB(double qsq, double ml){
    return 3.0*J6(qsq,ml)/(8.0*diffWidth(qsq,ml));}

double BdtoKstrll_obs::P1(double qsq, double ml){return 1.0/2.0*(J3(qsq,ml)+J3b(qsq,ml))/(J2s(qsq,ml)+J2bs(qsq,ml));}

double BdtoKstrll_obs::P2(double qsq, double ml){return 1.0/8.0*(J6s(qsq,ml)+J6bs(qsq,ml))/(J2s(qsq,ml)+J2bs(qsq,ml));}

double BdtoKstrll_obs::P4p(double qsq, double ml){return 1.0/sqrt(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml)))*(J4(qsq,ml)+J4b(qsq,ml));}

double BdtoKstrll_obs::P5p(double qsq, double ml){return 1.0/(2.0*sqrt(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml))))*(J5(qsq,ml)+J5b(qsq,ml));}

double BdtoKstrll_obs::P6p(double qsq, double ml){return -1.0/(2.0*sqrt(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml))))*(J7(qsq,ml)+J7b(qsq,ml));}

double BdtoKstrll_obs::P8p(double qsq, double ml){return -1.0/sqrt(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml)))*(J8(qsq,ml)+J8b(qsq,ml));}

double BdtoKstrll_obs::S1(double qsq, double ml){return (J1s(qsq,ml) + J1c(qsq,ml) + J1bs(qsq,ml) + J1bc(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S2(double qsq, double ml){return (J2s(qsq,ml) + J2c(qsq,ml) + J2bs(qsq,ml) + J2bc(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S3(double qsq, double ml){return (J3(qsq,ml) + J3b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S4(double qsq, double ml){return (J4(qsq,ml) + J4b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S5(double qsq, double ml){return (J5(qsq,ml) + J5b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S6(double qsq, double ml){return (J6s(qsq,ml) + J6c(qsq,ml) + J6bs(qsq,ml) + J6bc(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S7(double qsq, double ml){return (J7(qsq,ml) + J7b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S8(double qsq, double ml){return (J8(qsq,ml) + J8b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double BdtoKstrll_obs::S9(double qsq, double ml){return (J9(qsq,ml) + J9b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll:
double Bstophill_obs::J1s(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml)*AaRIm(qsq,ml) + AaLRe(qsq,ml)*AaRRe(qsq,ml) + ApLIm(qsq,ml)*ApRIm(qsq,ml) + ApLRe(qsq,ml)*ApRRe(qsq,ml)))/qsq +
    ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
      pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double Bstophill_obs::J1c(double qsq, double ml){
    return pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml),2.0) + pow(AtRe(qsq,ml),2.0) +
                      2.0*(AzLIm(qsq,ml)*AzRIm(qsq,ml) + AzLRe(qsq,ml)*AzRRe(qsq,ml))))/qsq +
    (pow(ASIm(qsq,ml),2.0) + pow(ASRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ;}

double Bstophill_obs::J1(double qsq, double ml){
    return 2.0*J1s(qsq,ml) + J1c(qsq,ml);}

double Bstophill_obs::J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double Bstophill_obs::J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double Bstophill_obs::J2(double qsq, double ml){
    return 2.0*J2s(qsq,ml) + J2c(qsq,ml);}

double Bstophill_obs::J3(double qsq, double ml){
    return ((-pow(AaLIm(qsq,ml),2.0) - pow(AaLRe(qsq,ml),2.0) - pow(AaRIm(qsq,ml),2.0) - pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double Bstophill_obs::J4(double qsq, double ml){
    return ((AaLIm(qsq,ml)*AzLIm(qsq,ml) + AaLRe(qsq,ml)*AzLRe(qsq,ml) + AaRIm(qsq,ml)*AzRIm(qsq,ml) +
             AaRRe(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obs::J5(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml)*ASIm(qsq,ml) + AaRIm(qsq,ml)*ASIm(qsq,ml) + AaLRe(qsq,ml)*ASRe(qsq,ml) + AaRRe(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq)) + ApLIm(qsq,ml)*AzLIm(qsq,ml) +
        ApLRe(qsq,ml)*AzLRe(qsq,ml) - ApRIm(qsq,ml)*AzRIm(qsq,ml) - ApRRe(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double Bstophill_obs::J6s(double qsq, double ml){
    return 2.0*(AaLIm(qsq,ml)*ApLIm(qsq,ml) + AaLRe(qsq,ml)*ApLRe(qsq,ml) - AaRIm(qsq,ml)*ApRIm(qsq,ml) -
                AaRRe(qsq,ml)*ApRRe(qsq,ml))*betal(qsq,ml) ;}

double Bstophill_obs::J6c(double qsq, double ml){
    return (4.0*ml*(ASIm(qsq,ml)*AzLIm(qsq,ml) + ASRe(qsq,ml)*AzLRe(qsq,ml) + ASIm(qsq,ml)*AzRIm(qsq,ml) +
                    ASRe(qsq,ml)*AzRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double Bstophill_obs::J6(double qsq, double ml){
    return 2.0*J6s(qsq,ml) + J6c(qsq,ml);}

double Bstophill_obs::J7(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml)*ASIm(qsq,ml)) - ApRRe(qsq,ml)*ASIm(qsq,ml) + ApLIm(qsq,ml)*ASRe(qsq,ml) + ApRIm(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq) + AaLRe(qsq,ml)*AzLIm(qsq,ml) -
        AaLIm(qsq,ml)*AzLRe(qsq,ml) - AaRRe(qsq,ml)*AzRIm(qsq,ml) + AaRIm(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double Bstophill_obs::J8(double qsq, double ml){
    return ((ApLRe(qsq,ml)*AzLIm(qsq,ml) - ApLIm(qsq,ml)*AzLRe(qsq,ml) + ApRRe(qsq,ml)*AzRIm(qsq,ml) -
             ApRIm(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obs::J9(double qsq, double ml){
    return (AaLRe(qsq,ml)*ApLIm(qsq,ml) - AaLIm(qsq,ml)*ApLRe(qsq,ml) + AaRRe(qsq,ml)*ApRIm(qsq,ml) -
            AaRIm(qsq,ml)*ApRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ;}

//Conjugate J's
double Bstophill_obs::J1bs(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AabLIm(qsq,ml)*AabRIm(qsq,ml) + AabLRe(qsq,ml)*AabRRe(qsq,ml) + ApbLIm(qsq,ml)*ApbRIm(qsq,ml) + ApbLRe(qsq,ml)*ApbRRe(qsq,ml)))/qsq +
    ((pow(AabLIm(qsq,ml),2.0) + pow(AabLRe(qsq,ml),2.0) + pow(AabRIm(qsq,ml),2.0) + pow(AabRRe(qsq,ml),2.0) +
      pow(ApbLIm(qsq,ml),2.0) + pow(ApbLRe(qsq,ml),2.0) + pow(ApbRIm(qsq,ml),2.0) + pow(ApbRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double Bstophill_obs::J1bc(double qsq, double ml){
    return pow(AzbLIm(qsq,ml),2.0) + pow(AzbLRe(qsq,ml),2.0) + pow(AzbRIm(qsq,ml),2.0) + pow(AzbRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtbIm(qsq,ml),2.0) + pow(AtbRe(qsq,ml),2.0) +
                      2.0*(AzbLIm(qsq,ml)*AzbRIm(qsq,ml) + AzbLRe(qsq,ml)*AzbRRe(qsq,ml))))/qsq +
    (pow(ASbIm(qsq,ml),2.0) + pow(ASbRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ;}

double Bstophill_obs::J1b(double qsq, double ml){
    return 2.0*J1bs(qsq,ml) + J1bc(qsq,ml);}

double Bstophill_obs::J2bs(double qsq, double ml){
    return ((pow(AabLIm(qsq,ml),2.0) + pow(AabLRe(qsq,ml),2.0) + pow(AabRIm(qsq,ml),2.0) + pow(AabRRe(qsq,ml),2.0) +
             pow(ApbLIm(qsq,ml),2.0) + pow(ApbLRe(qsq,ml),2.0) + pow(ApbRIm(qsq,ml),2.0) + pow(ApbRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double Bstophill_obs::J2bc(double qsq, double ml){
    return -((pow(AzbLIm(qsq,ml),2.0) + pow(AzbLRe(qsq,ml),2.0) + pow(AzbRIm(qsq,ml),2.0) + pow(AzbRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double Bstophill_obs::J2b(double qsq, double ml){
    return 2.0*J2bs(qsq,ml) + J2bc(qsq,ml);}

double Bstophill_obs::J3b(double qsq, double ml){
    return ((-pow(AabLIm(qsq,ml),2.0) - pow(AabLRe(qsq,ml),2.0) - pow(AabRIm(qsq,ml),2.0) - pow(AabRRe(qsq,ml),2.0) +
             pow(ApbLIm(qsq,ml),2.0) + pow(ApbLRe(qsq,ml),2.0) + pow(ApbRIm(qsq,ml),2.0) + pow(ApbRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double Bstophill_obs::J4b(double qsq, double ml){
    return ((AabLIm(qsq,ml)*AzbLIm(qsq,ml) + AabLRe(qsq,ml)*AzbLRe(qsq,ml) + AabRIm(qsq,ml)*AzbRIm(qsq,ml) +
             AabRRe(qsq,ml)*AzbRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obs::J5b(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AabLIm(qsq,ml)*ASbIm(qsq,ml) + AabRIm(qsq,ml)*ASbIm(qsq,ml) + AabLRe(qsq,ml)*ASbRe(qsq,ml)
                             + AabRRe(qsq,ml)*ASbRe(qsq,ml)))/sqrt(qsq)) + ApbLIm(qsq,ml)*AzbLIm(qsq,ml) +
                      ApbLRe(qsq,ml)*AzbLRe(qsq,ml) - ApbRIm(qsq,ml)*AzbRIm(qsq,ml) - ApbRRe(qsq,ml)*AzbRRe(qsq,ml))*betal(qsq,ml) ;}

double Bstophill_obs::J6bs(double qsq, double ml){
    return 2.0*(AabLIm(qsq,ml)*ApbLIm(qsq,ml) + AabLRe(qsq,ml)*ApbLRe(qsq,ml) - AabRIm(qsq,ml)*ApbRIm(qsq,ml) -
                AabRRe(qsq,ml)*ApbRRe(qsq,ml))*betal(qsq,ml) ;}

double Bstophill_obs::J6bc(double qsq, double ml){
    return (4.0*ml*(ASbIm(qsq,ml)*AzbLIm(qsq,ml) + ASbRe(qsq,ml)*AzbLRe(qsq,ml) + ASbIm(qsq,ml)*AzbRIm(qsq,ml) +
                    ASbRe(qsq,ml)*AzbRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double Bstophill_obs::J6b(double qsq, double ml){
    return 2.0*J6bs(qsq,ml) + J6bc(qsq,ml);}

double Bstophill_obs::J7b(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApbLRe(qsq,ml)*ASbIm(qsq,ml)) - ApbRRe(qsq,ml)*ASbIm(qsq,ml) + ApbLIm(qsq,ml)*ASbRe(qsq,ml)
                           + ApbRIm(qsq,ml)*ASbRe(qsq,ml)))/sqrt(qsq) + AabLRe(qsq,ml)*AzbLIm(qsq,ml) -
                      AabLIm(qsq,ml)*AzbLRe(qsq,ml) - AabRRe(qsq,ml)*AzbRIm(qsq,ml) + AabRIm(qsq,ml)*AzbRRe(qsq,ml))*betal(qsq,ml) ;}

double Bstophill_obs::J8b(double qsq, double ml){
    return ((ApbLRe(qsq,ml)*AzbLIm(qsq,ml) - ApbLIm(qsq,ml)*AzbLRe(qsq,ml) + ApbRRe(qsq,ml)*AzbRIm(qsq,ml) -
             ApbRIm(qsq,ml)*AzbRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double Bstophill_obs::J9b(double qsq, double ml){
    return (AabLRe(qsq,ml)*ApbLIm(qsq,ml) - AabLIm(qsq,ml)*ApbLRe(qsq,ml) + AabRRe(qsq,ml)*ApbRIm(qsq,ml) -
            AabRIm(qsq,ml)*ApbRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ;}


double Bstophill_obs::diffWidth(double qsq, double ml){
    return 3.0/4.0*(J1(qsq,ml) - J2(qsq,ml)/3.0);}

double Bstophill_obs::diffWidthConj(double qsq, double ml){
    return 3.0/4.0*(J1b(qsq,ml) - J2b(qsq,ml)/3.0);}

double Bstophill_obs::FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0));}

double Bstophill_obs::AFB(double qsq, double ml){
    return (3.0*(J6c(qsq,ml) + 2.0*J6s(qsq,ml)))/(8.0*diffWidth(qsq,ml));}

double Bstophill_obs::P1(double qsq, double ml){return 1.0/2.0*(J3(qsq,ml)+J3b(qsq,ml))/(J2s(qsq,ml)+J2bs(qsq,ml));}

double Bstophill_obs::P2(double qsq, double ml){return 1.0/8.0*(J6s(qsq,ml)+J6bs(qsq,ml))/(J2s(qsq,ml)+J2bs(qsq,ml));}

double Bstophill_obs::P4p(double qsq, double ml){return 1.0/(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml)))*(J4(qsq,ml)+J4b(qsq,ml));}

double Bstophill_obs::P5p(double qsq, double ml){return 1.0/(2.0*(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml))))*(J5(qsq,ml)+J5b(qsq,ml));}

double Bstophill_obs::P6p(double qsq, double ml){return -1.0/(2.0*(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml))))*(J7(qsq,ml)+J7b(qsq,ml));}

double Bstophill_obs::P8p(double qsq, double ml){return -1.0/(-(J2s(qsq,ml)+J2bs(qsq,ml))*(J2c(qsq,ml)+J2bc(qsq,ml)))*(J8(qsq,ml)+J8b(qsq,ml));}

double Bstophill_obs::S1(double qsq, double ml){return (J1s(qsq,ml) + J1c(qsq,ml) + J1bs(qsq,ml) + J1bc(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S2(double qsq, double ml){return (J2s(qsq,ml) + J2c(qsq,ml) + J2bs(qsq,ml) + J2bc(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S3(double qsq, double ml){return (J3(qsq,ml) + J3b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S4(double qsq, double ml){return (J4(qsq,ml) + J4b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S5(double qsq, double ml){return (J5(qsq,ml) + J5b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S6(double qsq, double ml){return (J6s(qsq,ml) + J6c(qsq,ml) + J6bs(qsq,ml) + J6bc(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S7(double qsq, double ml){return (J7(qsq,ml) + J7b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S8(double qsq, double ml){return (J8(qsq,ml) + J8b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}

double Bstophill_obs::S9(double qsq, double ml){return (J9(qsq,ml) + J9b(qsq,ml))/(diffWidth(qsq,ml)+diffWidthConj(qsq,ml));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_obs::alNP(double qsq, double ml){
    return (4.0*pow(mBd(),2.0)*pow(ml,2.0)*(pow(FAIm(qsq,ml),2.0) + pow(FARe(qsq,ml),2.0)) +
            2.0*ml*(pow(mBd(),2.0) - pow(mK(),2.0) + qsq)*(FAIm(qsq,ml)*FPIm(qsq,ml) + FARe(qsq,ml)*FPRe(qsq,ml)) +
            qsq*(pow(FPIm(qsq,ml),2.0) + pow(FPRe(qsq,ml),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FSIm(qsq,ml),2.0) + pow(FSRe(qsq,ml),2.0))) +
            ((pow(FAIm(qsq,ml),2.0) + pow(FARe(qsq,ml),2.0) + pow(FVIm(qsq,ml),2.0) + pow(FVRe(qsq,ml),2.0))*
             lambda(qsq))/4.0)*nf(qsq,ml)*pow(xiP(qsq),2.0);}

double BdtoKll_obs::blNP(double qsq, double ml){
    return 2.0*(qsq*(FPIm(qsq,ml)*FT5Im(qsq,ml) + FPRe(qsq,ml)*FT5Re(qsq,ml) +
                     pow(betal(qsq,ml),2.0)*(FSIm(qsq,ml)*FTIm(qsq,ml) + FSRe(qsq,ml)*FTRe(qsq,ml))) +
                ml*((pow(mBd(),2.0) - pow(mK(),2.0) + qsq)*(FAIm(qsq,ml)*FT5Im(qsq,ml) + FARe(qsq,ml)*FT5Re(qsq,ml)) +
                    betal(qsq,ml)*(FSIm(qsq,ml)*FVIm(qsq,ml) + FSRe(qsq,ml)*FVRe(qsq,ml))*sqrt(lambda(qsq))))*
    nf(qsq,ml)*pow(xiP(qsq),2.0);}

double BdtoKll_obs::clNP(double qsq, double ml){
    return (qsq*(pow(FT5Im(qsq,ml),2.0) + pow(FT5Re(qsq,ml),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FTIm(qsq,ml),2.0) + pow(FTRe(qsq,ml),2.0))) +
            2.0*ml*betal(qsq,ml)*(FTIm(qsq,ml)*FVIm(qsq,ml) + FTRe(qsq,ml)*FVRe(qsq,ml))*sqrt(lambda(qsq)) -
            (pow(betal(qsq,ml),2.0)*(pow(FAIm(qsq,ml),2.0) + pow(FARe(qsq,ml),2.0) + pow(FVIm(qsq,ml),2.0) +
                                     pow(FVRe(qsq,ml),2.0))*lambda(qsq))/4.0)*nf(qsq,ml)*pow(xiP(qsq),2.0);}

double BdtoKll_obs::diffWidth(double qsq, double ml){
    return 2.0*(alNP(qsq,ml) + 1.0/3.0*clNP(qsq,ml)); }

double BdtoKll_obs::diffBrnch(double qsq, double ml){
    return tauBd()/2.0*diffWidth(qsq,ml); }

double BdtoKll_obs::diffAFB(double qsq, double ml){
    return blNP(qsq,ml)/diffWidth(qsq,ml); }

double BdtoKll_obs::diffFH(double qsq, double ml){
    return 2.0*(alNP(qsq,ml)+clNP(qsq,ml))/diffWidth(qsq,ml); }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//B->ll
double Btoll_obs :: Btoll_nf(double mBq){
    if (mBq==mBd()){
        return pow(GF()*alpha_e()*fBd()*absVtbVtdStr(),2.0)*tauBd()/(64.0*pow(mBd()*pi,3.0)); }
    if (mBq==mBs()){
        return pow(GF()*alpha_e()*fBs()*absVtbVtsStr(),2.0)*tauBs()/(64.0*pow(mBs()*pi,3.0)); }
    return 0.0;}

double Btoll_obs :: ADeltaGammaf(double mBq, double mq, double ml1, double ml2){
    return (-pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0) + pow(ampSIm(mBq,mq,ml1,ml2),2.0) -
            pow(ampSRe(mBq,mq,ml1,ml2),2.0))/
    (pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0) + pow(ampSIm(mBq,mq,ml1,ml2),2.0) +
     pow(ampSRe(mBq,mq,ml1,ml2),2.0));}

double Btoll_obs :: CorrctnFctr(double mBq, double mq, double ml1, double ml2){
    return (1.0 - pow(y(mBq),2.0))/(1.0 + ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq));}

double Btoll_obs :: BrInst(double mBq, double mq, double ml1, double ml2){
    return ((pow(mBq,2.0) - pow(ml1 - ml2,2.0))*(pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0)) + (pow(mBq,2.0) - pow(ml1 + ml2,2.0))*(pow(ampSIm(mBq,mq,ml1,ml2),2.0) +
        pow(ampSRe(mBq,mq,ml1,ml2),2.0)))*sqrt(Btoll_lambda(mBq,ml1,ml2))*Btoll_nf(mBq);}

double Btoll_obs :: BrTimeIntgratd(double mBq, double mq, double ml1, double ml2){
    return BrInst(mBq,mq,ml1,ml2)/CorrctnFctr(mBq,mq,ml1,ml2);}

double Btoll_obs :: efftau(double mBq, double mq, double ml1, double ml2){
    return (tauB(mBq)*(1.0 + 2.0*ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq) + pow(y(mBq),2.0)))/
    ((1.0 + ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq))*(1.0 - pow(y(mBq),2.0)));}











