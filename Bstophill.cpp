#include<iostream>
#include<string>
#include<cmath>
using namespace std;

//parameters:
double pi=M_PI;
double mB=5.27953,mKst=0.896,mu=4.8,mc=1.27,mb=4.8;
double C1=-0.257,C2=1.009,C3=-0.005,C4=-0.078,C5=0.000,C6=0.001,C7effRe=-0.304,/*C8eff=-0.167,*/C10effRe=-4.103,
        tauB=1.519*1.52*pow(10.0,12.0),alpha_e=1.0/133.0,GF=1.16637*pow(10.0,-5.0),absVtbVtsStr=0.0409,
        me=0.511*pow(10.0,-3.0),mmu=105.658*pow(10.0,-3.0);
//SMsub:
double C10effIm=0.0,C7effIm=0.0,C7RHRe=-0.006,C7RHIm=0,C9RHRe=0.0,C9RHIm=0.0,C10RHRe=0.0,C10RHIm=0.0,CPRe=0.0,CPIm=0.0,CPRHRe=0.0,CPRHIm=0.0,CSRe=0.0,CSIm=0.0,CSRHRe=0.0,CSRHIm=0.0,CTRe=0.0,CTIm=0.0,CT5Re=0.0,CT5Im=0.0;

double tp = pow(mB+mKst,2.0),tm = pow(mB-mKst,2.0),tz = tp*(1.0-sqrt(1.0-tm/tp));
double z(double qsq){
    return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));   }

//extrapolatedLCSR-Latticeformfactors:(see,arXiv:1503.05534)
double mPSbs=5.366,mVbs=5.415,mAbs=5.829;
double Va0=0.376313,Va1=-1.16597,Va2=2.42443,
        A0a0=0.369196,A0a1=-1.36584,A0a2=0.128191,
        A1a0=0.29725,A1a1=0.392378,A1a2=1.18916,
        A12a0=0.265375,A12a1=0.533638,A12a2=0.483166,
        T1a0=0.312055,T1a1=-1.00893,T1a2=1.5272,
        T2a0=0.312055,T2a1=0.496846,T2a2=1.61431,
        T23a0=0.667412,T23a1=1.31812,T23a2=3.82334;

double V(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs,2.0))*( Va0*pow(z(qsq)-z(0.0),0.0) + Va1*pow(z(qsq)-z(0.0),1.0) + Va2*pow(z(qsq)-z(0.0),2.0) ); }

double A0(double qsq){
    return 1.0/(1.0-qsq/pow(mPSbs,2.0))*( A0a0*pow(z(qsq)-z(0.0),0.0) + A0a1*pow(z(qsq)-z(0.0),1.0) + A0a2*pow(z(qsq)-z(0.0),2.0) ); }

double A1(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs,2.0))*( A1a0*pow(z(qsq)-z(0.0),0.0) + A1a1*pow(z(qsq)-z(0.0),1.0) + A1a2*pow(z(qsq)-z(0.0),2.0) ); }

double A12(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs,2.0))*( A12a0*pow(z(qsq)-z(0.0),0.0) + A12a1*pow(z(qsq)-z(0.0),1.0) + A12a2*pow(z(qsq)-z(0.0),2.0) ); }

double T1(double qsq){
    return 1.0/(1.0-qsq/pow(mVbs,2.0))*( T1a0*pow(z(qsq)-z(0.0),0.0) + T1a1*pow(z(qsq)-z(0.0),1.0) + T1a2*pow(z(qsq)-z(0.0),2.0) ); }

double T2(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs,2.0))*( T2a0*pow(z(qsq)-z(0.0),0.0) + T2a1*pow(z(qsq)-z(0.0),1.0) + T2a2*pow(z(qsq)-z(0.0),2.0) ); }

double T23(double qsq){
    return 1.0/(1.0-qsq/pow(mAbs,2.0))*( T23a0*pow(z(qsq)-z(0.0),0.0) + T23a1*pow(z(qsq)-z(0.0),1.0) + T23a2*pow(z(qsq)-z(0.0),2.0) ); }


//C9eff:(see,arXiv:0811.1214)
double hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu,2.0)); }

double hImZ(double qsq){
    return 4.0/9.0*pi; }

double YRe(double qsq){
    return hRe(qsq,mc)*(4.0/3.0*C1 + C2 + 6.0*C3 + 60.0*C5) - 1.0/2.0*hRe(qsq,mb)*(7.0*C3 + 4.0/3.0*C4 + 76.0*C5 + 64.0/3.0*C6) - 1.0/2.0*hReZ(qsq)*(C3 + 4.0/3.0*C4 + 16.0*C5 + 64.0/3.0*C6) + 4.0/3.0*C3 + 64.0/9.0*C5 + 64.0/27.0*C6; }

double YIm(double qsq){
    return hIm(qsq,mc)*(4.0/3.0*C1 + C2 + 6.0*C3 + 60.0*C5) - 1.0/2.0*hIm(qsq,mb)*(7.0*C3 + 4.0/3.0*C4 + 76.0*C5 + 64.0/3.0*C6) - 1.0/2.0*hImZ(qsq)*(C3 + 4.0/3.0*C4 + 16.0*C5 + 64.0/3.0*C6); }

double C9effRe(double qsq){
    return 4.211 + YRe(qsq); }
double C9effIm(double qsq){
    return YIm(qsq); }

//amplitudes:(see,arXiv:0811.1214)
double betal(double qsq, double ml){
    return sqrt(1.0-4.0*pow(ml,2.0)/qsq); }
double lambda(double qsq){
    return pow(qsq,2.0)+pow(mB,4.0)+pow(mKst,4.0)-2.0*(pow(mB*mKst,2.0)+qsq*pow(mB,2.0)+qsq*pow(mKst,2.0)); }
double nf(double qsq, double ml){
    return absVtbVtsStr*sqrt( pow(GF*alpha_e,2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mB,3.0)) ); }

//TransversityAmplitudes_SM+NP(S(CiS)+PS(CiPS)+RHC(CiRH)):
double ApLRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe)-(C10effRe+C10RHRe))*V(qsq)/(mB+mKst) + 2.0*mb/qsq*(C7effRe+C7RHRe)*T1(qsq)); }
double ApLIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm)-(C10effIm+C10RHIm))*V(qsq)/(mB+mKst) + 2.0*mb/qsq*(C7effIm+C7RHIm)*T1(qsq)); }

double ApRRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe)+(C10effRe+C10RHRe))*V(qsq)/(mB+mKst) + 2.0*mb/qsq*(C7effRe+C7RHRe)*T1(qsq)); }
double ApRIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm)+(C10effIm+C10RHIm))*V(qsq)/(mB+mKst) + 2.0*mb/qsq*(C7effIm+C7RHIm)*T1(qsq)); }

double AaLRe(double qsq, double ml){
    return -nf(qsq,ml)*sqrt(2.0)*(pow(mB,2.0)-pow(mKst,2.0))*(((C9effRe(qsq)-C9RHRe)-(C10effRe-C10RHRe))*A1(qsq)/(mB-mKst) + 2.0*mb/qsq*(C7effRe-C7RHRe)*T2(qsq)); }
double AaLIm(double qsq, double ml){
    return -nf(qsq,ml)*sqrt(2.0)*(pow(mB,2.0)-pow(mKst,2.0))*(((C9effIm(qsq)-C9RHIm)-(C10effIm-C10RHIm))*A1(qsq)/(mB-mKst) + 2.0*mb/qsq*(C7effIm-C7RHIm)*T2(qsq)); }

double AaRRe(double qsq, double ml){
    return -nf(qsq,ml)*sqrt(2.0)*(pow(mB,2.0)-pow(mKst,2.0))*(((C9effRe(qsq)-C9RHRe)+(C10effRe-C10RHRe))*A1(qsq)/(mB-mKst) + 2.0*mb/qsq*(C7effRe-C7RHRe)*T2(qsq)); }
double AaRIm(double qsq, double ml){
    return -nf(qsq,ml)*sqrt(2.0)*(pow(mB,2.0)-pow(mKst,2.0))*(((C9effIm(qsq)-C9RHIm)+(C10effIm-C10RHIm))*A1(qsq)/(mB-mKst) + 2.0*mb/qsq*(C7effIm-C7RHIm)*T2(qsq)); }

double AzLRe(double qsq, double ml){
    return -nf(qsq,ml)/(2.0*mKst*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe)-(C10effRe-C10RHRe))*16.0*mB*mKst*mKst*A12(qsq) + 2.0*mb/(pow(mB,2.0)-pow(mKst,2.0))*(C7effRe-C7RHRe)*(8.0*mB*pow(mKst,2.0)*(mB-mKst)*T23(qsq)) ); }
double AzLIm(double qsq, double ml){
    return -nf(qsq,ml)/(2.0*mKst*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm)-(C10effIm-C10RHIm))*16.0*mB*mKst*mKst*A12(qsq) + 2.0*mb/(pow(mB,2.0)-pow(mKst,2.0))*(C7effIm-C7RHIm)*(8.0*mB*pow(mKst,2.0)*(mB-mKst)*T23(qsq)) ); }

double AzRRe(double qsq, double ml){
    return -nf(qsq,ml)/(2.0*mKst*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe)+(C10effRe-C10RHRe))*16.0*mB*mKst*mKst*A12(qsq) + 2.0*mb/(pow(mB,2.0)-pow(mKst,2.0))*(C7effRe-C7RHRe)*(8.0*mB*pow(mKst,2.0)*(mB-mKst)*T23(qsq)) ); }
double AzRIm(double qsq, double ml){
    return -nf(qsq,ml)/(2.0*mKst*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm)+(C10effIm-C10RHIm))*16.0*mB*mKst*mKst*A12(qsq) + 2.0*mb/(pow(mB,2.0)-pow(mKst,2.0))*(C7effIm-C7RHIm)*(8.0*mB*pow(mKst,2.0)*(mB-mKst)*T23(qsq)) ); }

double AtRe(double qsq, double ml){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10effRe-C10RHRe) + qsq/ml*(CPRe-CPRHRe))*A0(qsq); }
double AtIm(double qsq, double ml){
    return nf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10effIm-C10RHIm) + qsq/ml*(CPIm-CPRHIm))*A0(qsq); }

double ASRe(double qsq, double ml){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSRe-CSRHRe)*A0(qsq); }
double ASIm(double qsq, double ml){
    return -2.0*nf(qsq,ml)*sqrt(lambda(qsq))*(CSIm-CSRHIm)*A0(qsq); }

//AngularCoefficients:(see,complexexpand.nb)
double J1s(double qsq, double ml){
    return (4.0*pow(ml,2.0)*(AaLIm(qsq,ml)*AaRIm(qsq,ml) + AaLRe(qsq,ml)*AaRRe(qsq,ml) + ApLIm(qsq,ml)*ApRIm(qsq,ml) +
                           ApLRe(qsq,ml)*ApRRe(qsq,ml)))/qsq +
    ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
      pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))
     *(2.0 + pow(betal(qsq,ml),2.0)))/4.0 ; }

double J1c(double qsq, double ml){
    return pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) +
    (4.0*pow(ml,2.0)*(pow(AtIm(qsq,ml),2.0) + pow(AtRe(qsq,ml),2.0) +
                    2.0*(AzLIm(qsq,ml)*AzRIm(qsq,ml) + AzLRe(qsq,ml)*AzRRe(qsq,ml))))/qsq +
    (pow(ASIm(qsq,ml),2.0) + pow(ASRe(qsq,ml),2.0))*pow(betal(qsq,ml),2.0) ; }

double J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ; }

double J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ; }

double J3(double qsq, double ml){
    return ((-pow(AaLIm(qsq,ml),2.0) - pow(AaLRe(qsq,ml),2.0) - pow(AaRIm(qsq,ml),2.0) - pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/2.0 ; }

double J4(double qsq, double ml){
    return ((AaLIm(qsq,ml)*AzLIm(qsq,ml) + AaLRe(qsq,ml)*AzLRe(qsq,ml) + AaRIm(qsq,ml)*AzRIm(qsq,ml) +
             AaRRe(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ; }

double J5(double qsq, double ml){
    return sqrt(2.0)*(-((ml*(AaLIm(qsq,ml)*ASIm(qsq,ml) + AaRIm(qsq,ml)*ASIm(qsq,ml) + AaLRe(qsq,ml)*ASRe(qsq,ml) +
                           AaRRe(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq)) + ApLIm(qsq,ml)*AzLIm(qsq,ml) +
                    ApLRe(qsq,ml)*AzLRe(qsq,ml) - ApRIm(qsq,ml)*AzRIm(qsq,ml) - ApRRe(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ; }

double J6s(double qsq, double ml){
    return 2.0*(AaLIm(qsq,ml)*ApLIm(qsq,ml) + AaLRe(qsq,ml)*ApLRe(qsq,ml) - AaRIm(qsq,ml)*ApRIm(qsq,ml) -
              AaRRe(qsq,ml)*ApRRe(qsq,ml))*betal(qsq,ml) ; }

double J6c(double qsq, double ml){
    return (4.0*ml*(ASIm(qsq,ml)*AzLIm(qsq,ml) + ASRe(qsq,ml)*AzLRe(qsq,ml) + ASIm(qsq,ml)*AzRIm(qsq,ml) +
                  ASRe(qsq,ml)*AzRRe(qsq,ml))*betal(qsq,ml))/sqrt(qsq) ; }

double J7(double qsq, double ml){
    return sqrt(2.0)*((ml*(-(ApLRe(qsq,ml)*ASIm(qsq,ml)) - ApRRe(qsq,ml)*ASIm(qsq,ml) + ApLIm(qsq,ml)*ASRe(qsq,ml) +
                         ApRIm(qsq,ml)*ASRe(qsq,ml)))/sqrt(qsq) + AaLRe(qsq,ml)*AzLIm(qsq,ml) -
                    AaLIm(qsq,ml)*AzLRe(qsq,ml) - AaRRe(qsq,ml)*AzRIm(qsq,ml) + AaRIm(qsq,ml)*AzRRe(qsq,ml))*
    betal(qsq,ml) ; }

double J8(double qsq, double ml){
    return ((ApLRe(qsq,ml)*AzLIm(qsq,ml) - ApLIm(qsq,ml)*AzLRe(qsq,ml) + ApRRe(qsq,ml)*AzRIm(qsq,ml) -
             ApRIm(qsq,ml)*AzRRe(qsq,ml))*pow(betal(qsq,ml),2.0))/sqrt(2.0) ; }

double J9(double qsq, double ml){
    return (AaLRe(qsq,ml)*ApLIm(qsq,ml) - AaLIm(qsq,ml)*ApLRe(qsq,ml) + AaRRe(qsq,ml)*ApRIm(qsq,ml) -
            AaRIm(qsq,ml)*ApRRe(qsq,ml))*pow(betal(qsq,ml),2.0) ; }

/*
//AngularObservables:(see,arXiv:1409.3088)
double diffWidth(double qsq, double ml){
    return pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
    pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
    pow(ASIm(qsq,ml),2.0) + pow(ASRe(qsq,ml),2.0) + pow(AtIm(qsq,ml),2.0) + pow(AtRe(qsq,ml),2.0) +
    pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) ; }

double diffBrnch(double qsq, double ml){
    return tauB/2.0*diffWidth(qsq,ml); }

double FL(double qsq, double ml){
    return 1.0/diffWidth(qsq,ml)*( pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0) ) ; }

double FP(double qsq, double ml){
    return 1.0/diffWidth(qsq,ml)*( pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) ) ; }

double AFB(double qsq, double ml){
    return (3.0*(AaLIm(qsq,ml)*ApLIm(qsq,ml) + AaLRe(qsq,ml)*ApLRe(qsq,ml) - AaRIm(qsq,ml)*ApRIm(qsq,ml) -
               AaRRe(qsq,ml)*ApRRe(qsq,ml)))/(2.0*diffWidth(qsq,ml)) ; }

double A4(double qsq, double ml){
    return (sqrt(2.0)*(AaLIm(qsq,ml)*AzLIm(qsq,ml) + AaLRe(qsq,ml)*AzLRe(qsq,ml) + AaRIm(qsq,ml)*AzRIm(qsq,ml) +
                     AaRRe(qsq,ml)*AzRRe(qsq,ml)))/(pi*diffWidth(qsq,ml)) ; }

double A5(double qsq, double ml){
    return (3.0*(ApLIm(qsq,ml)*AzLIm(qsq,ml) + ApLRe(qsq,ml)*AzLRe(qsq,ml) - ApRIm(qsq,ml)*AzRIm(qsq,ml) -
               ApRRe(qsq,ml)*AzRRe(qsq,ml)))/(2.0*sqrt(2.0)*diffWidth(qsq,ml)) ; }

double A7(double qsq, double ml){
    return (3.0*(AaLRe(qsq,ml)*AzLIm(qsq,ml) - AaLIm(qsq,ml)*AzLRe(qsq,ml) - AaRRe(qsq,ml)*AzRIm(qsq,ml) +
               AaRIm(qsq,ml)*AzRRe(qsq,ml)))/(2.0*sqrt(2.0)*diffWidth(qsq,ml)) ; }

double A8(double qsq, double ml){
    return (sqrt(2.0)*(ApLRe(qsq,ml)*AzLIm(qsq,ml) - ApLIm(qsq,ml)*AzLRe(qsq,ml) + ApRRe(qsq,ml)*AzRIm(qsq,ml) -
                     ApRIm(qsq,ml)*AzRRe(qsq,ml)))/(pi*diffWidth(qsq,ml)) ; }

double A9(double qsq, double ml){
    return (3.0*(AaLRe(qsq,ml)*ApLIm(qsq,ml) - AaLIm(qsq,ml)*ApLRe(qsq,ml) + AaRRe(qsq,ml)*ApRIm(qsq,ml) -
               AaRIm(qsq,ml)*ApRRe(qsq,ml)))/(2.0*pi*diffWidth(qsq,ml)) ; }

//LHCbObservables:
double S3(double qsq, double ml){
    return (2.0*FP(qsq,ml) + FL(qsq,ml) - 1.0)/2.0; }

double S4(double qsq, double ml){
    return -pi/2.0*A4(qsq,ml); }

double S5(double qsq, double ml){
    return 4.0/3.0*A5(qsq,ml); }

double S7(double qsq, double ml){
    return 4.0/3.0*A7(qsq,ml); }

double S8(double qsq, double ml){
    return -pi/2.0*A8(qsq,ml); }

double S9(double qsq, double ml){
    return 2.0*pi/3.0*A9(qsq,ml); }
*/

//observables:(see,arXiv:2006.03489)

double diffWidth(double qsq, double ml){
    return (3.0*(J1c(qsq,ml) + 2.0*J1s(qsq,ml) + (-J2c(qsq,ml) - 2.0*J2s(qsq,ml))/3.0))/4.0; }

double FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0)); }

double AFB(double qsq, double ml){
    return (3.0*(J6c(qsq,ml) + 2.0*J6s(qsq,ml)))/(8.0*diffWidth(qsq,ml)); }




int main(){
    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;
    cout << "Observable values at this point are:"
    << diffWidth(qsq,mmu) << "," << FL(qsq,mmu) << "," << AFB(qsq,mmu) << endl;
    return 0;
}










    









