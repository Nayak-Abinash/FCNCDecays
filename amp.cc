#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include "amp.h"
#include <iostream>
#include <cmath>
using namespace std;

//Bs->phill:
double amp::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double amp::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double amp::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0)); }

double amp::hImZ(double qsq){
    return 4.0/9.0*pi; }

double amp::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6(); }

double amp::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()); }

double amp::C9effRe(double qsq){return 4.211 + YRe(qsq);}
double amp::C9effIm(double qsq){return YIm(qsq);}

double amp::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double amp::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBs(),4.0)+pow(mKst(),4.0)-2.0*(pow(mBs()*mKst(),2.0)+qsq*pow(mBs(),2.0)+qsq*pow(mKst(),2.0));}

double amp:: nf(double qsq, double ml){
    return absVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(),3.0)) );}

double amp::ApLRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}
double amp::ApLIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}

double amp::ApRRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}
double amp::ApRIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}

double amp::AaLRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq));}
double amp::AaLIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq));}

double amp::AaRRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq));}
double amp::AaRIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*(((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq));}

double amp::AzLRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}
double amp::AzLIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}

double amp::AzRRe(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}
double amp::AzRIm(double qsq, double ml){
    return -1.0*nf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq) + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) );}

double amp::AtRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq);}
double amp::AtIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq);}

double amp::ASRe(double qsq, double ml){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq);}
double amp::ASIm(double qsq, double ml){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq);}



















