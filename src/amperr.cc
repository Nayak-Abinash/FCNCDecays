#include "amperr.h"

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

double BdtoKstrll_amperr::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6())
            - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6(); }

double BdtoKstrll_amperr::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6())
            - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()); }

double BdtoKstrll_amperr::C9effRe(double qsq){return 4.211 + YRe(qsq);}
double BdtoKstrll_amperr::C9effIm(double qsq){return YIm(qsq);}

double BdtoKstrll_amperr::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double BdtoKstrll_amperr::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBd(),4.0)+pow(mKst(),4.0)-2.0*(pow(mBd()*mKst(),2.0)+qsq*pow(mBd(),2.0)+qsq*pow(mKst(),2.0));}

double BdtoKstrll_amperr:: nf(double qsq, double ml){
    return absVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBd(),3.0)) );}

double BdtoKstrll_amperr::ApLRe(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq,unv)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq,unv));}
double BdtoKstrll_amperr::ApLIm(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq,unv)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq,unv));}

double BdtoKstrll_amperr::ApRRe(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq,unv)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq,unv));}
double BdtoKstrll_amperr::ApRIm(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq,unv)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq,unv));}

double BdtoKstrll_amperr::AaLRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq,unv)/(mBd()-mKst())
                                        + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq,unv));}
double BdtoKstrll_amperr::AaLIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq,unv)/(mBd()-mKst())
                                        + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq,unv));}

double BdtoKstrll_amperr::AaRRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq,unv)/(mBd()-mKst())
                                        + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq,unv));}
double BdtoKstrll_amperr::AaRIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq,unv)/(mBd()-mKst())
                                        + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq,unv));}

double BdtoKstrll_amperr::AzLRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq,unv)) );}
double BdtoKstrll_amperr::AzLIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq,unv)) );}

double BdtoKstrll_amperr::AzRRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq,unv)) );}
double BdtoKstrll_amperr::AzRIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq,unv)) );}

double BdtoKstrll_amperr::AtRe(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq,unv);}
double BdtoKstrll_amperr::AtIm(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq,unv);}

double BdtoKstrll_amperr::ASRe(double qsq, double ml, double unv[]){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq,unv);}
double BdtoKstrll_amperr::ASIm(double qsq, double ml, double unv[]){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq,unv);}
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

double Bstophill_amperr::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6())
            - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6(); }

double Bstophill_amperr::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6())
            - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()); }

double Bstophill_amperr::C9effRe(double qsq){return 4.211 + YRe(qsq);}
double Bstophill_amperr::C9effIm(double qsq){return YIm(qsq);}

double Bstophill_amperr::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double Bstophill_amperr::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBs(),4.0)+pow(mKst(),4.0)-2.0*(pow(mBs()*mKst(),2.0)+qsq*pow(mBs(),2.0)+qsq*pow(mKst(),2.0));}

double Bstophill_amperr:: nf(double qsq, double ml){
    return absVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(),3.0)) );}

double Bstophill_amperr::ApLRe(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq,unv)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq,unv));}
double Bstophill_amperr::ApLIm(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq,unv)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq,unv));}

double Bstophill_amperr::ApRRe(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq,unv)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq,unv));}
double Bstophill_amperr::ApRIm(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq,unv)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq,unv));}

double Bstophill_amperr::AaLRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq,unv)/(mBs()-mKst())
                                        + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq,unv));}
double Bstophill_amperr::AaLIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq,unv)/(mBs()-mKst())
                                        + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq,unv));}

double Bstophill_amperr::AaRRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq,unv)/(mBs()-mKst())
                                        + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq,unv));}
double Bstophill_amperr::AaRIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq,unv)/(mBs()-mKst())
                                        + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq,unv));}

double Bstophill_amperr::AzLRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq,unv)) );}
double Bstophill_amperr::AzLIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq,unv)) );}

double Bstophill_amperr::AzRRe(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq,unv)) );}
double Bstophill_amperr::AzRIm(double qsq, double ml, double unv[]){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq,unv)
                                        + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq,unv)) );}

double Bstophill_amperr::AtRe(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq,unv);}
double Bstophill_amperr::AtIm(double qsq, double ml, double unv[]){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq,unv);}

double Bstophill_amperr::ASRe(double qsq, double ml, double unv[]){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq,unv);}
double Bstophill_amperr::ASIm(double qsq, double ml, double unv[]){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq,unv);}
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

double BdtoKll_amperr::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6())
            - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6();}

double BdtoKll_amperr::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6())
            - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6());}

double BdtoKll_amperr::C9effRe(double qsq){
    return 4.211 + YRe(qsq);}
double BdtoKll_amperr::C9effIm(double qsq){
    return YIm(qsq);}

double BdtoKll_amperr::betal(double qsq, double ml){
    return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double BdtoKll_amperr::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBd(),4.0)+pow(mK(),4.0)-2.0*(pow(mBd()*mK(),2.0)+qsq*pow(mBd(),2.0)+qsq*pow(mK(),2.0));}

double BdtoKll_amperr::nf(double qsq, double ml){
    return pow(GF()*alpha_e()*absVtbVtsStr(),2.0)/(512.0*pow(pi,5.0)*pow(mBd(),3.0))*betal(qsq,ml)*sqrt(lambda(qsq));}

double BdtoKll_amperr::xiP(double qsq, double unv[]){return fp(qsq,unv);}

double BdtoKll_amperr::TauPRe(double qsq, double unv[]){
    return xiP(qsq,unv)*(C7effRe() + mBd()/(2.0*mb())*YRe(qsq));}

double BdtoKll_amperr::TauPIm(double qsq, double unv[]){
    return xiP(qsq,unv)*(C7effIm() + mBd()/(2.0*mb())*YIm(qsq));}

double BdtoKll_amperr::FVRe(double qsq, double ml, double unv[]){
    return (C9Re() + C9RHRe())+ 2.0*mb()/mBd()*TauPRe(qsq,unv)/xiP(qsq,unv) + 8.0*ml/(mBd()+mK())*fT(qsq,unv)/fp(qsq,unv)*CTRe();}
double BdtoKll_amperr::FVIm(double qsq, double ml, double unv[]){
    return (C9Im() + C9RHIm())+ 2.0*mb()/mBd()*TauPIm(qsq,unv)/xiP(qsq,unv) + 8.0*ml/(mBd()+mK())*fT(qsq,unv)/fp(qsq,unv)*CTIm();}

double BdtoKll_amperr::FARe(double qsq, double ml, double unv[]){
    return (C10Re() + C10RHRe());}
double BdtoKll_amperr::FAIm(double qsq, double ml, double unv[]){
    return (C10Im() + C10RHIm());}

double BdtoKll_amperr::FSRe(double qsq, double ml, double unv[]){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq,unv)/fp(qsq,unv)*(CSRe() + CSRHRe());}
double BdtoKll_amperr::FSIm(double qsq, double ml, double unv[]){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq,unv)/fp(qsq,unv)*(CSIm() + CSRHIm());}

double BdtoKll_amperr::FPRe(double qsq, double ml, double unv[]){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq,unv)/fp(qsq,unv)*(CPRe() + CPRHRe()) +
    ml*C10Re()*((pow(mBd(),2.0)-pow(mK(),2.0))/qsq*(fz(qsq,unv)/fp(qsq,unv)-1.0)-1.0);}
double BdtoKll_amperr::FPIm(double qsq, double ml, double unv[]){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq,unv)/fp(qsq,unv)*(CPIm() + CPRHIm()) +
    ml*C10Im()*((pow(mBd(),2.0)-pow(mK(),2.0))/qsq*(fz(qsq,unv)/fp(qsq,unv)-1.0)-1.0);}


double BdtoKll_amperr::FTRe(double qsq, double ml, double unv[]){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq,unv)/fp(qsq,unv)*CTRe();}
double BdtoKll_amperr::FTIm(double qsq, double ml, double unv[]){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq,unv)/fp(qsq,unv)*CTIm();}

double BdtoKll_amperr::FT5Re(double qsq, double ml, double unv[]){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq,unv)/fp(qsq,unv)*CT5Re();}
double BdtoKll_amperr::FT5Im(double qsq, double ml, double unv[]){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq,unv)/fp(qsq,unv)*CT5Im();}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
//B->ll
double Btoll_amperr :: ampSRe(double mBq, double mq, double ml1, double ml2){
    return ((C9Re()-C9RHRe())*(ml1-ml2) + (CSRe()-CSRHRe())*pow(mBq,2.0)/(mb()+mq));}
double Btoll_amperr :: ampSIm(double mBq, double mq, double ml1, double ml2){
    return ((C9Im()-C9RHIm())*(ml1-ml2) + (CSIm()-CSRHIm())*pow(mBq,2.0)/(mb()+mq));}

double Btoll_amperr :: ampPRe(double mBq, double mq, double ml1, double ml2){
    return ((C10Re()-C10RHRe())*(ml1+ml2) + (CPRe()-CPRHRe())*pow(mBq,2.0)/(mb()+mq));}
double Btoll_amperr :: ampPIm(double mBq, double mq, double ml1, double ml2){
    return ((C10Im()-C10RHIm())*(ml1+ml2) + (CPIm()-CPRHIm())*pow(mBq,2.0)/(mb()+mq));}

double Btoll_amperr :: Btoll_lambda(double mBq, double ml1, double ml2){
    return (pow(mBq,2.0)-pow(ml1-ml2,2.0))*(pow(mBq,2.0)-pow(ml1+ml2,2.0));}

double Btoll_amperr :: y(double mBq){
    if (mBq==mBs()){
        return DGamma_sbar()/(2.0); }
    if (mBq==mBd()){
        return DGamma_dbar()/(2.0); }
    return 0.0;}

double Btoll_amperr :: tauB(double mBq){
    if (mBq==mBs()){
        return tauBs(); }
    if (mBq==mBd()){
        return tauBd(); }
    return 0.0;}
*/

















