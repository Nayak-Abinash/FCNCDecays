#include "amplitudes_mc.h"

//Bd->Kstr,ll:
double BdtoKstrll_amperr::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq)
                -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq)
                -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double BdtoKstrll_amperr::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double BdtoKstrll_amperr::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0)); }

double BdtoKstrll_amperr::hImZ(double qsq){
    return 4.0/9.0*pi; }

double BdtoKstrll_amperr::YRe(double qsq, double smwc[], double npwc[], double unv[]){
    return hRe(qsq,mc(unv))*(4.0/3.0*C1(smwc,unv) + C2(smwc,unv) + 6.0*C3(smwc,unv) + 60.0*C5(smwc,unv))
            - 1.0/2.0*hRe(qsq,mb(unv))*(7.0*C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 76.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            - 1.0/2.0*hReZ(qsq)*(C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 16.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            + 4.0/3.0*C3(smwc,unv) + 64.0/9.0*C5(smwc,unv) + 64.0/27.0*C6(smwc,unv); }

double BdtoKstrll_amperr::YIm(double qsq, double smwc[], double npwc[], double unv[]){
    return hIm(qsq,mc(unv))*(4.0/3.0*C1(smwc,unv) + C2(smwc,unv) + 6.0*C3(smwc,unv) + 60.0*C5(smwc,unv))
            - 1.0/2.0*hIm(qsq,mb(unv))*(7.0*C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 76.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            - 1.0/2.0*hImZ(qsq)*(C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 16.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv)); }

double BdtoKstrll_amperr::C9effRe(double qsq, double smwc[], double npwc[], double unv[]){return C9Re(smwc,unv) + YRe(qsq, smwc, npwc, unv);}
double BdtoKstrll_amperr::C9effIm(double qsq, double smwc[], double npwc[], double unv[]){return C9Im(smwc,unv) + YIm(qsq, smwc, npwc, unv);}

double BdtoKstrll_amperr::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double BdtoKstrll_amperr::lambda(double qsq, double unv[]){
    return pow(qsq,2.0)+pow(mBd(unv),4.0)+pow(mKst(unv),4.0)-2.0*(pow(mBd(unv)*mKst(unv),2.0)+qsq*pow(mBd(unv),2.0)+qsq*pow(mKst(unv),2.0));}

double BdtoKstrll_amperr:: Renf(double qsq, double ml, double unv[]){
    return ReVtbVtsStr(unv)*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq, unv))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBd(unv),3.0)) );}

double BdtoKstrll_amperr:: Imnf(double qsq, double ml, double unv[]){
    return ImVtbVtsStr(unv)*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq, unv))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBd(unv),3.0)) );}

double BdtoKstrll_amperr::ApLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double BdtoKstrll_amperr::ApLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double BdtoKstrll_amperr::ApRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double BdtoKstrll_amperr::ApRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double BdtoKstrll_amperr::AaLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double BdtoKstrll_amperr::AaLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double BdtoKstrll_amperr::AaRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double BdtoKstrll_amperr::AaRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double BdtoKstrll_amperr::AzLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}
double BdtoKstrll_amperr::AzLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}

double BdtoKstrll_amperr::AzRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}
double BdtoKstrll_amperr::AzRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}

double BdtoKstrll_amperr::AtRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv)
    - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv);}
double BdtoKstrll_amperr::AtIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv)
    + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv);}

double BdtoKstrll_amperr::ASRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv)
                 - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv));}
double BdtoKstrll_amperr::ASIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv)
                 + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv));}

//Conjugate Amplitudes
double BdtoKstrll_amperr::ApbLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double BdtoKstrll_amperr::ApbLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double BdtoKstrll_amperr::ApbRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double BdtoKstrll_amperr::ApbRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBd(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double BdtoKstrll_amperr::AabLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double BdtoKstrll_amperr::AabLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double BdtoKstrll_amperr::AabRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double BdtoKstrll_amperr::AabRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBd(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double BdtoKstrll_amperr::AzbLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}
double BdtoKstrll_amperr::AzbLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}

double BdtoKstrll_amperr::AzbRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}
double BdtoKstrll_amperr::AzbRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBd(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBd(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBd(unv)*pow(mKst(unv),2.0)*(mBd(unv)-mKst(unv))*T23(qsq,unv)) ));}

double BdtoKstrll_amperr::AtbRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv)
    + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv);}
double BdtoKstrll_amperr::AtbIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv)
    - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv);}

double BdtoKstrll_amperr::ASbRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv)
                 + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv));}
double BdtoKstrll_amperr::ASbIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv)
                 - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv));}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll:
double Bstophill_amperr::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq)
                -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq)
                -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double Bstophill_amperr::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double Bstophill_amperr::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0)); }

double Bstophill_amperr::hImZ(double qsq){
    return 4.0/9.0*pi; }

double Bstophill_amperr::YRe(double qsq, double smwc[], double npwc[], double unv[]){
    return hRe(qsq,mc(unv))*(4.0/3.0*C1(smwc,unv) + C2(smwc,unv) + 6.0*C3(smwc,unv) + 60.0*C5(smwc,unv))
            - 1.0/2.0*hRe(qsq,mb(unv))*(7.0*C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 76.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            - 1.0/2.0*hReZ(qsq)*(C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 16.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            + 4.0/3.0*C3(smwc,unv) + 64.0/9.0*C5(smwc,unv) + 64.0/27.0*C6(smwc,unv); }

double Bstophill_amperr::YIm(double qsq, double smwc[], double npwc[], double unv[]){
    return hIm(qsq,mc(unv))*(4.0/3.0*C1(smwc,unv) + C2(smwc,unv) + 6.0*C3(smwc,unv) + 60.0*C5(smwc,unv))
            - 1.0/2.0*hIm(qsq,mb(unv))*(7.0*C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 76.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            - 1.0/2.0*hImZ(qsq)*(C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 16.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv)); }

double Bstophill_amperr::C9effRe(double qsq, double smwc[], double npwc[], double unv[]){return C9Re(smwc,unv) + YRe(qsq,smwc,npwc,unv);}
double Bstophill_amperr::C9effIm(double qsq, double smwc[], double npwc[], double unv[]){return C9Im(smwc,unv) + YIm(qsq,smwc,npwc,unv);}

double Bstophill_amperr::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double Bstophill_amperr::lambda(double qsq, double unv[]){
    return pow(qsq,2.0)+pow(mBs(unv),4.0)+pow(mKst(unv),4.0)-2.0*(pow(mBs(unv)*mKst(unv),2.0)+qsq*pow(mBs(unv),2.0)+qsq*pow(mKst(unv),2.0));}

double Bstophill_amperr:: Renf(double qsq, double ml, double unv[]){
    return ReVtbVtsStr(unv)*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq, unv))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(unv),3.0)) );}

double Bstophill_amperr:: Imnf(double qsq, double ml, double unv[]){
    return ImVtbVtsStr(unv)*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq, unv))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(unv),3.0)) );}

double Bstophill_amperr::ApLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double Bstophill_amperr::ApLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double Bstophill_amperr::ApRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double Bstophill_amperr::ApRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double Bstophill_amperr::AaLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double Bstophill_amperr::AaLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double Bstophill_amperr::AaRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double Bstophill_amperr::AaRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double Bstophill_amperr::AzLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}
double Bstophill_amperr::AzLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}

double Bstophill_amperr::AzRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}
double Bstophill_amperr::AzRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}

double Bstophill_amperr::AtRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv)
    - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv);}
double Bstophill_amperr::AtIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv)
    + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv);}

double Bstophill_amperr::ASRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv)
                 - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv));}
double Bstophill_amperr::ASIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv)
                 + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv));}

//Conjugate Amplitudes
double Bstophill_amperr::ApbLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double Bstophill_amperr::ApbLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))-(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))-(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double Bstophill_amperr::ApbRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv))
    + Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv));}
double Bstophill_amperr::ApbRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effIm(qsq,smwc,npwc,unv)+C9RHIm(npwc,unv))+(C10Im(smwc,unv)+C10RHIm(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)+C7RHIm(npwc,unv))*T1(qsq,unv))
    - Imnf(qsq,ml,unv)*sqrt(2.0*lambda(qsq,unv))*(((C9effRe(qsq,smwc,npwc,unv)+C9RHRe(npwc,unv))+(C10Re(smwc,unv)+C10RHRe(npwc,unv)))*V(qsq,unv)/(mBs(unv)+mKst(unv))
                                                     + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)+C7RHRe(npwc,unv))*T1(qsq,unv));}

double Bstophill_amperr::AabLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double Bstophill_amperr::AabLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        -(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        -(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double Bstophill_amperr::AabRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv))
                 + Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv)));}
double Bstophill_amperr::AabRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))
                                        +(C10Im(smwc,unv)-C10RHIm(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))*T2(qsq,unv))
                 - Imnf(qsq,ml,unv)*sqrt(2.0)*(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))
                                        +(C10Re(smwc,unv)-C10RHRe(npwc,unv)))*A1(qsq,unv)/(mBs(unv)-mKst(unv))
                                        + 2.0*mb(unv)/qsq*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))*T2(qsq,unv)));}

double Bstophill_amperr::AzbLRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}
double Bstophill_amperr::AzbLIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))-(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))-(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}

double Bstophill_amperr::AzbRRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 + Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}
double Bstophill_amperr::AzbRIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -1.0*(Renf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effIm(qsq,smwc,npwc,unv)-C9RHIm(npwc,unv))+(C10Im(smwc,unv)-C10RHIm(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effIm(smwc,unv)-C7RHIm(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) )
                 - Imnf(qsq,ml,unv)/(2.0*mKst(unv)*sqrt(qsq))*( ((C9effRe(qsq,smwc,npwc,unv)-C9RHRe(npwc,unv))+(C10Re(smwc,unv)-C10RHRe(npwc,unv)))
                                                          *16.0*mBs(unv)*mKst(unv)*mKst(unv)*A12(qsq,unv)
                    + 2.0*mb(unv)/(pow(mBs(unv),2.0)-pow(mKst(unv),2.0))*(C7effRe(smwc,unv)-C7RHRe(npwc,unv))
                                                          *(8.0*mBs(unv)*pow(mKst(unv),2.0)*(mBs(unv)-mKst(unv))*T23(qsq,unv)) ));}

double Bstophill_amperr::AtbRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv)
    + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv);}
double Bstophill_amperr::AtbIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Im(smwc,unv)-C10RHIm(npwc,unv)) + qsq/ml*(CPIm(npwc,unv)-CPRHIm(npwc,unv)))*A0(qsq,unv)
    - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv)/qsq)*(2.0*(C10Re(smwc,unv)-C10RHRe(npwc,unv)) + qsq/ml*(CPRe(npwc,unv)-CPRHRe(npwc,unv)))*A0(qsq,unv);}

double Bstophill_amperr::ASbRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv)
                 + Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv));}
double Bstophill_amperr::ASbIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return -2.0*(Renf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSIm(npwc,unv)-CSRHIm(npwc,unv))*A0(qsq,unv)
                 - Imnf(qsq,ml,unv)*sqrt(lambda(qsq,unv))*(CSRe(npwc,unv)-CSRHRe(npwc,unv))*A0(qsq,unv));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll
double BdtoKll_amperr::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq)
                -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq)
                -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double BdtoKll_amperr::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double BdtoKll_amperr::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0));}

double BdtoKll_amperr::hImZ(double qsq){
    return 4.0/9.0*pi;}

double BdtoKll_amperr::YRe(double qsq, double smwc[], double npwc[], double unv[]){
    return hRe(qsq,mc(unv))*(4.0/3.0*C1(smwc,unv) + C2(smwc,unv) + 6.0*C3(smwc,unv) + 60.0*C5(smwc,unv))
            - 1.0/2.0*hRe(qsq,mb(unv))*(7.0*C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 76.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            - 1.0/2.0*hReZ(qsq)*(C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 16.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            + 4.0/3.0*C3(smwc,unv) + 64.0/9.0*C5(smwc,unv) + 64.0/27.0*C6(smwc,unv);}

double BdtoKll_amperr::YIm(double qsq, double smwc[], double npwc[], double unv[]){
    return hIm(qsq,mc(unv))*(4.0/3.0*C1(smwc,unv) + C2(smwc,unv) + 6.0*C3(smwc,unv) + 60.0*C5(smwc,unv))
            - 1.0/2.0*hIm(qsq,mb(unv))*(7.0*C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 76.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv))
            - 1.0/2.0*hImZ(qsq)*(C3(smwc,unv) + 4.0/3.0*C4(smwc,unv) + 16.0*C5(smwc,unv) + 64.0/3.0*C6(smwc,unv));}

double BdtoKll_amperr::C9effRe(double qsq, double smwc[], double npwc[], double unv[]){
    return C9Re(smwc,unv) + YRe(qsq,smwc,npwc,unv);}
double BdtoKll_amperr::C9effIm(double qsq, double smwc[], double npwc[], double unv[]){
    return C9Im(smwc,unv) + YIm(qsq,smwc,npwc,unv);}

double BdtoKll_amperr::betal(double qsq, double ml){
    return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double BdtoKll_amperr::lambda(double qsq, double unv[]){
    return pow(qsq,2.0)+pow(mBd(unv),4.0)+pow(mK(unv),4.0)-2.0*(pow(mBd(unv)*mK(unv),2.0)+qsq*pow(mBd(unv),2.0)+qsq*pow(mK(unv),2.0));}

double BdtoKll_amperr::nf(double qsq, double ml, double unv[]){
    return pow(GF()*alpha_e()*absVtbVtsStr(unv),2.0)/(512.0*pow(pi,5.0)*pow(mBd(unv),3.0))*betal(qsq,ml)*sqrt(lambda(qsq,unv));}

double BdtoKll_amperr::xiP(double qsq, double smwc[], double npwc[], double unv[]){return fp(qsq,unv);}

double BdtoKll_amperr::TauPRe(double qsq, double smwc[], double npwc[], double unv[]){
    return xiP(qsq,smwc,npwc,unv)*(C7effRe(smwc,unv) + mBd(unv)/(2.0*mb(unv))*YRe(qsq,smwc,npwc,unv));}

double BdtoKll_amperr::TauPIm(double qsq, double smwc[], double npwc[], double unv[]){
    return xiP(qsq,smwc,npwc,unv)*(C7effIm(smwc,unv) + mBd(unv)/(2.0*mb(unv))*YIm(qsq,smwc,npwc,unv));}

double BdtoKll_amperr::FVRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (C9Re(smwc,unv) + C9RHRe(npwc,unv))+ 2.0*mb(unv)/mBd(unv)*TauPRe(qsq,smwc,npwc,unv)/xiP(qsq,smwc,npwc,unv) + 8.0*ml/(mBd(unv)+mK(unv))*fT(qsq,unv)/fp(qsq,unv)*CTRe(npwc,unv);}
double BdtoKll_amperr::FVIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (C9Im(smwc,unv) + C9RHIm(npwc,unv))+ 2.0*mb(unv)/mBd(unv)*TauPIm(qsq,smwc,npwc,unv)/xiP(qsq,smwc,npwc,unv) + 8.0*ml/(mBd(unv)+mK(unv))*fT(qsq,unv)/fp(qsq,unv)*CTIm(npwc,unv);}

double BdtoKll_amperr::FARe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (C10Re(smwc,unv) + C10RHRe(npwc,unv));}
double BdtoKll_amperr::FAIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return (C10Im(smwc,unv) + C10RHIm(npwc,unv));}

double BdtoKll_amperr::FSRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 1.0/2.0*(pow(mBd(unv),2.0)-pow(mK(unv),2.0))/(mb(unv)-ms(unv))*fz(qsq,unv)/fp(qsq,unv)*(CSRe(npwc,unv) + CSRHRe(npwc,unv));}
double BdtoKll_amperr::FSIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 1.0/2.0*(pow(mBd(unv),2.0)-pow(mK(unv),2.0))/(mb(unv)-ms(unv))*fz(qsq,unv)/fp(qsq,unv)*(CSIm(npwc,unv) + CSRHIm(npwc,unv));}

double BdtoKll_amperr::FPRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 1.0/2.0*(pow(mBd(unv),2.0)-pow(mK(unv),2.0))/(mb(unv)-ms(unv))*fz(qsq,unv)/fp(qsq,unv)*(CPRe(npwc,unv) + CPRHRe(npwc,unv)) +
    ml*C10Re(smwc,unv)*((pow(mBd(unv),2.0)-pow(mK(unv),2.0))/qsq*(fz(qsq,unv)/fp(qsq,unv)-1.0)-1.0);}
double BdtoKll_amperr::FPIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 1.0/2.0*(pow(mBd(unv),2.0)-pow(mK(unv),2.0))/(mb(unv)-ms(unv))*fz(qsq,unv)/fp(qsq,unv)*(CPIm(npwc,unv) + CPRHIm(npwc,unv)) +
    ml*C10Im(smwc,unv)*((pow(mBd(unv),2.0)-pow(mK(unv),2.0))/qsq*(fz(qsq,unv)/fp(qsq,unv)-1.0)-1.0);}


double BdtoKll_amperr::FTRe(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*sqrt(lambda(qsq,unv))*betal(qsq,ml)/(mBd(unv)+mK(unv))*fT(qsq,unv)/fp(qsq,unv)*CTRe(npwc,unv);}
double BdtoKll_amperr::FTIm(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*sqrt(lambda(qsq,unv))*betal(qsq,ml)/(mBd(unv)+mK(unv))*fT(qsq,unv)/fp(qsq,unv)*CTIm(npwc,unv);}

double BdtoKll_amperr::FT5Re(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*sqrt(lambda(qsq,unv))*betal(qsq,ml)/(mBd(unv)+mK(unv))*fT(qsq,unv)/fp(qsq,unv)*CT5Re(npwc,unv);}
double BdtoKll_amperr::FT5Im(double qsq, double ml, double smwc[], double npwc[], double unv[]){
    return 2.0*sqrt(lambda(qsq,unv))*betal(qsq,ml)/(mBd(unv)+mK(unv))*fT(qsq,unv)/fp(qsq,unv)*CT5Im(npwc,unv);}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//B->ll
double Btoll_amperr :: ampSRe(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return ((C9Re(smwc,unv)-C9RHRe(npwc,unv))*(ml1-ml2) + (CSRe(npwc,unv)-CSRHRe(npwc,unv))*pow(mBq,2.0)/(mb(unv)+mq));}
double Btoll_amperr :: ampSIm(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return ((C9Im(smwc,unv)-C9RHIm(npwc,unv))*(ml1-ml2) + (CSIm(npwc,unv)-CSRHIm(npwc,unv))*pow(mBq,2.0)/(mb(unv)+mq));}

double Btoll_amperr :: ampPRe(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return ((C10Re(smwc,unv)-C10RHRe(npwc,unv))*(ml1+ml2) + (CPRe(npwc,unv)-CPRHRe(npwc,unv))*pow(mBq,2.0)/(mb(unv)+mq));}
double Btoll_amperr :: ampPIm(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return ((C10Im(smwc,unv)-C10RHIm(npwc,unv))*(ml1+ml2) + (CPIm(npwc,unv)-CPRHIm(npwc,unv))*pow(mBq,2.0)/(mb(unv)+mq));}

double Btoll_amperr :: Btoll_lambda(double mBq, double ml1, double ml2, double smwc[], double npwc[], double unv[]){
    return (pow(mBq,2.0)-pow(ml1-ml2,2.0))*(pow(mBq,2.0)-pow(ml1+ml2,2.0));}

double Btoll_amperr :: y(double mBq, double smwc[], double npwc[], double unv[]){
    if (mBq==mBs(unv)){
        return DGamma_sbar(unv)/(2.0); }
    if (mBq==mBd(unv)){
        return DGamma_dbar(unv)/(2.0); }
    return 0.0;}

double Btoll_amperr :: tauB(double mBq, double smwc[], double npwc[], double unv[]){
    if (mBq==mBs(unv)){
        return tauBs(unv); }
    if (mBq==mBd(unv)){
        return tauBd(unv); }
    return 0.0;}


















