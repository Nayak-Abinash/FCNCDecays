#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include "amp.h"
#include "obs.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bs->phill:
double obs::J1s(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml)*AaRIm(qsq,ml) + AaLRe(qsq,ml)*AaRRe(qsq,ml) + ApLIm(qsq,ml)*ApRIm(qsq,ml) +
                             ApLRe(qsq,ml)*ApRRe(qsq,ml)))/qsq +
    ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
      pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ;}

double obs::J1c(double qsq, double ml){
    return pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml),2.0) + pow(AtRe(qsq,ml),2.0) +
                      2.0*(AzLIm(qsq,ml)*AzRIm(qsq,ml) + AzLRe(qsq,ml)*AzRRe(qsq,ml))))/qsq +
    (pow(ASIm(qsq,ml),2.0) + pow(ASRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ;}

double obs::J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double obs::J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

double obs::J3(double qsq, double ml){
    return ((-pow(AaLIm(qsq,ml),2.0) - pow(AaLRe(qsq,ml),2.0) - pow(AaRIm(qsq,ml),2.0) - pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ;}

double obs::J4(double qsq, double ml){
    return ((AaLIm(qsq,ml)*AzLIm(qsq,ml) + AaLRe(qsq,ml)*AzLRe(qsq,ml) + AaRIm(qsq,ml)*AzRIm(qsq,ml) +
             AaRRe(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double obs::J5(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml)*ASIm(qsq,ml) + AaRIm(qsq,ml)*ASIm(qsq,ml) + AaLRe(qsq,ml)*ASRe(qsq,ml) +
                             AaRRe(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq)) + ApLIm(qsq,ml)*AzLIm(qsq,ml) +
                      ApLRe(qsq,ml)*AzLRe(qsq,ml) - ApRIm(qsq,ml)*AzRIm(qsq,ml) - ApRRe(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double obs::J6s(double qsq, double ml){
    return 2.0*(AaLIm(qsq,ml)*ApLIm(qsq,ml) + AaLRe(qsq,ml)*ApLRe(qsq,ml) - AaRIm(qsq,ml)*ApRIm(qsq,ml) -
                AaRRe(qsq,ml)*ApRRe(qsq,ml))*betal(qsq,ml) ;}

double obs::J6c(double qsq, double ml){
    return (4.0*ml*(ASIm(qsq,ml)*AzLIm(qsq,ml) + ASRe(qsq,ml)*AzLRe(qsq,ml) + ASIm(qsq,ml)*AzRIm(qsq,ml) +
                    ASRe(qsq,ml)*AzRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double obs::J7(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml)*ASIm(qsq,ml)) - ApRRe(qsq,ml)*ASIm(qsq,ml) + ApLIm(qsq,ml)*ASRe(qsq,ml) +
                           ApRIm(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq) + AaLRe(qsq,ml)*AzLIm(qsq,ml) -
                      AaLIm(qsq,ml)*AzLRe(qsq,ml) - AaRRe(qsq,ml)*AzRIm(qsq,ml) + AaRIm(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ;}

double obs::J8(double qsq, double ml){
    return ((ApLRe(qsq,ml)*AzLIm(qsq,ml) - ApLIm(qsq,ml)*AzLRe(qsq,ml) + ApRRe(qsq,ml)*AzRIm(qsq,ml) -
             ApRIm(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ;}

double obs::J9(double qsq, double ml){
    return (AaLRe(qsq,ml)*ApLIm(qsq,ml) - AaLIm(qsq,ml)*ApLRe(qsq,ml) + AaRRe(qsq,ml)*ApRIm(qsq,ml) -
            AaRIm(qsq,ml)*ApRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ;}

double obs::diffWidth(double qsq, double ml){
    return (3.0*(J1c(qsq,ml) + 2.0*J1s(qsq,ml) + (-J2c(qsq,ml) - 2.0*J2s(qsq,ml))/3.0))/4.0;}

double obs::FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0));}

double obs::AFB(double qsq, double ml){
    return (3.0*(J6c(qsq,ml) + 2.0*J6s(qsq,ml)))/(8.0*diffWidth(qsq,ml));}











