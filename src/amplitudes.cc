#include "amplitudes.h"

//Bd->Kstr,ll:
double BdtoKstrll_amp::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double BdtoKstrll_amp::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double BdtoKstrll_amp::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0)); }

double BdtoKstrll_amp::hImZ(double qsq){
    return 4.0/9.0*pi; }

double BdtoKstrll_amp::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6(); }

double BdtoKstrll_amp::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()); }

double BdtoKstrll_amp::C9effRe(double qsq){return C9Re() + YRe(qsq);}
double BdtoKstrll_amp::C9effIm(double qsq){return C9Im() + YIm(qsq);}

double BdtoKstrll_amp::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double BdtoKstrll_amp::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBd(),4.0)+pow(mKst(),4.0)-2.0*(pow(mBd()*mKst(),2.0)+qsq*pow(mBd(),2.0)+qsq*pow(mKst(),2.0));}

double BdtoKstrll_amp:: Renf(double qsq, double ml){
    return ReVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBd(),3.0)) );}

double BdtoKstrll_amp:: Imnf(double qsq, double ml){
    return ImVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBd(),3.0)) );}

double BdtoKstrll_amp::ApLRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double BdtoKstrll_amp::ApLIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double BdtoKstrll_amp::ApRRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double BdtoKstrll_amp::ApRIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double BdtoKstrll_amp::AaLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double BdtoKstrll_amp::AaLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 + Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double BdtoKstrll_amp::AaRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                 - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double BdtoKstrll_amp::AaRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 + Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double BdtoKstrll_amp::AzLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}
double BdtoKstrll_amp::AzLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}

double BdtoKstrll_amp::AzRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}
double BdtoKstrll_amp::AzRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}

double BdtoKstrll_amp::AtRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq)
            - Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq);}
double BdtoKstrll_amp::AtIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq)
    + Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq);}

double BdtoKstrll_amp::ASRe(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq)
                 - Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq));}
double BdtoKstrll_amp::ASIm(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq)
                 + Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq));}


/////////Conjugate mode
double BdtoKstrll_amp::ApbLRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double BdtoKstrll_amp::ApbLIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double BdtoKstrll_amp::ApbRRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double BdtoKstrll_amp::ApbRIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBd()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double BdtoKstrll_amp::AabLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                + Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double BdtoKstrll_amp::AabLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double BdtoKstrll_amp::AabRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                 + Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double BdtoKstrll_amp::AabRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBd(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBd()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double BdtoKstrll_amp::AzbLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}
double BdtoKstrll_amp::AzbLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}

double BdtoKstrll_amp::AzbRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}
double BdtoKstrll_amp::AzbRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBd()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBd(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBd()*pow(mKst(),2.0)*(mBd()-mKst())*T23(qsq)) ));}

double BdtoKstrll_amp::AtbRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq)
            + Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq);}
double BdtoKstrll_amp::AtbIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq)
    - Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq);}

double BdtoKstrll_amp::ASbRe(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq)
                 + Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq));}
double BdtoKstrll_amp::ASbIm(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq)
                 - Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bs->phi,ll:
double Bstophill_amp::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double Bstophill_amp::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double Bstophill_amp::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0)); }

double Bstophill_amp::hImZ(double qsq){
    return 4.0/9.0*pi; }

double Bstophill_amp::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6(); }

double Bstophill_amp::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()); }

double Bstophill_amp::C9effRe(double qsq){return C9Re() + YRe(qsq);}
double Bstophill_amp::C9effIm(double qsq){return C9Im() + YIm(qsq);}

double Bstophill_amp::betal(double qsq, double ml){return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double Bstophill_amp::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBs(),4.0)+pow(mKst(),4.0)-2.0*(pow(mBs()*mKst(),2.0)+qsq*pow(mBs(),2.0)+qsq*pow(mKst(),2.0));}

double Bstophill_amp:: Renf(double qsq, double ml){
    return ReVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(),3.0)) );}

double Bstophill_amp:: Imnf(double qsq, double ml){
    return ImVtbVtsStr()*sqrt( pow(GF()*alpha_e(),2.0)*qsq*sqrt(lambda(qsq))*betal(qsq,ml)/(3.0*pow(2.0,10.0)*pow(pi,5.0)*pow(mBs(),3.0)) );}

double Bstophill_amp::ApLRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double Bstophill_amp::ApLIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double Bstophill_amp::ApRRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double Bstophill_amp::ApRIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double Bstophill_amp::AaLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double Bstophill_amp::AaLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 + Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double Bstophill_amp::AaRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                 - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double Bstophill_amp::AaRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 + Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double Bstophill_amp::AzLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}
double Bstophill_amp::AzLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}

double Bstophill_amp::AzRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}
double Bstophill_amp::AzRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}

double Bstophill_amp::AtRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq)
            - Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq);}
double Bstophill_amp::AtIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq)
            + Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq);}

double Bstophill_amp::ASRe(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq)
                 - Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq));}
double Bstophill_amp::ASIm(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq)
                 + Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq));}


/////////Conjugate mode
double Bstophill_amp::ApbLRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double Bstophill_amp::ApbLIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())-(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())-(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double Bstophill_amp::ApbRRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq))
    + Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq));}
double Bstophill_amp::ApbRIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effIm(qsq)+C9RHIm())+(C10Im()+C10RHIm()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effIm()+C7RHIm())*T1(qsq))
    - Imnf(qsq,ml)*sqrt(2.0*lambda(qsq))*(((C9effRe(qsq)+C9RHRe())+(C10Re()+C10RHRe()))*V(qsq)/(mBs()+mKst()) + 2.0*mb()/qsq*(C7effRe()+C7RHRe())*T1(qsq));}

double Bstophill_amp::AabLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                + Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double Bstophill_amp::AabLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double Bstophill_amp::AabRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq))
                 + Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq)));}
double Bstophill_amp::AabRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effIm()-C7RHIm())*T2(qsq))
                 - Imnf(qsq,ml)*sqrt(2.0)*(pow(mBs(),2.0)-pow(mKst(),2.0))*
                    (((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*A1(qsq)/(mBs()-mKst()) + 2.0*mb()/qsq*(C7effRe()-C7RHRe())*T2(qsq)));}

double Bstophill_amp::AzbLRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}
double Bstophill_amp::AzbLIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())-(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())-(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}

double Bstophill_amp::AzbRRe(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 + Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}
double Bstophill_amp::AzbRIm(double qsq, double ml){
    return -1.0*(Renf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effIm(qsq)-C9RHIm())+(C10Im()-C10RHIm()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effIm()-C7RHIm())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) )
                 - Imnf(qsq,ml)/(2.0*mKst()*sqrt(qsq))*( ((C9effRe(qsq)-C9RHRe())+(C10Re()-C10RHRe()))*16.0*mBs()*mKst()*mKst()*A12(qsq)
                    + 2.0*mb()/(pow(mBs(),2.0)-pow(mKst(),2.0))*(C7effRe()-C7RHRe())*(8.0*mBs()*pow(mKst(),2.0)*(mBs()-mKst())*T23(qsq)) ));}

double Bstophill_amp::AtbRe(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq)
            + Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq);}
double Bstophill_amp::AtbIm(double qsq, double ml){
    return Renf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Im()-C10RHIm()) + qsq/ml*(CPIm()-CPRHIm()))*A0(qsq)
            - Imnf(qsq,ml)*sqrt(lambda(qsq)/qsq)*(2.0*(C10Re()-C10RHRe()) + qsq/ml*(CPRe()-CPRHRe()))*A0(qsq);}

double Bstophill_amp::ASbRe(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq)
                 + Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq));}
double Bstophill_amp::ASbIm(double qsq, double ml){
    return -2.0*(Renf(qsq,ml)*sqrt(lambda(qsq))*(CSIm()-CSRHIm())*A0(qsq)
                 - Imnf(qsq,ml)*sqrt(lambda(qsq))*(CSRe()-CSRHRe())*A0(qsq));}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Bd->K,ll//
double BdtoKll_amp::hRe(double qsq, double mq){
    if (qsq < 4.0*pow(mq,2.0)){
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(4.0*pow(mq,2.0)/qsq-1.0)* atan(1.0/sqrt(4.0*pow(mq,2.0)/qsq-1.0)); }
    else {
        return -4.0/9.0*( 2.0*log(mq/mu()) -2.0/3.0 -4.0*pow(mq,2.0)/qsq) -4.0/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq)* log((1.0+sqrt(1.0-4.0*pow(mq,2.0)/qsq))/sqrt(4.0*pow(mq,2.0)/qsq)); } }

double BdtoKll_amp::hIm(double qsq, double mq){
    if (qsq >= 4.0*pow(mq,2.0)){
        return 2.0*pi/9.0*( 2.0+ 4.0*pow(mq,2.0)/qsq)* sqrt(1.0-4.0*pow(mq,2.0)/qsq); }
    else {
        return 0.0; } }

double BdtoKll_amp::hReZ(double qsq){
    return 8.0/27.0 - 4.0/9.0*log(qsq/pow(mu(),2.0));}

double BdtoKll_amp::hImZ(double qsq){
    return 4.0/9.0*pi;}

double BdtoKll_amp::YRe(double qsq){
    return hRe(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hRe(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hReZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6()) + 4.0/3.0*C3() + 64.0/9.0*C5() + 64.0/27.0*C6();}

double BdtoKll_amp::YIm(double qsq){
    return hIm(qsq,mc())*(4.0/3.0*C1() + C2() + 6.0*C3() + 60.0*C5()) - 1.0/2.0*hIm(qsq,mb())*(7.0*C3() + 4.0/3.0*C4() + 76.0*C5() + 64.0/3.0*C6()) - 1.0/2.0*hImZ(qsq)*(C3() + 4.0/3.0*C4() + 16.0*C5() + 64.0/3.0*C6());}

double BdtoKll_amp::C9effRe(double qsq){
    return C9Re() + YRe(qsq);}
double BdtoKll_amp::C9effIm(double qsq){
    return C9Im() + YIm(qsq);}

double BdtoKll_amp::betal(double qsq, double ml){
    return sqrt(1.0-4.0*pow(ml,2.0)/qsq);}

double BdtoKll_amp::lambda(double qsq){
    return pow(qsq,2.0)+pow(mBd(),4.0)+pow(mK(),4.0)-2.0*(pow(mBd()*mK(),2.0)+qsq*pow(mBd(),2.0)+qsq*pow(mK(),2.0));}

double BdtoKll_amp::nf(double qsq, double ml){
    return pow(GF()*alpha_e()*absVtbVtsStr(),2.0)/(512.0*pow(pi,5.0)*pow(mBd(),3.0))*betal(qsq,ml)*sqrt(lambda(qsq));}

double BdtoKll_amp::xiP(double qsq){return fp(qsq);}

double BdtoKll_amp::TauPRe(double qsq){
    return xiP(qsq)*(C7effRe() + mBd()/(2.0*mb())*YRe(qsq));}

double BdtoKll_amp::TauPIm(double qsq){
    return xiP(qsq)*(C7effIm() + mBd()/(2.0*mb())*YIm(qsq));}

double BdtoKll_amp::FVRe(double qsq, double ml){
    return (C9Re() + C9RHRe())+ 2.0*mb()/mBd()*TauPRe(qsq)/xiP(qsq) + 8.0*ml/(mBd()+mK())*fT(qsq)/fp(qsq)*CTRe();}
double BdtoKll_amp::FVIm(double qsq, double ml){
    return (C9Im() + C9RHIm())+ 2.0*mb()/mBd()*TauPIm(qsq)/xiP(qsq) + 8.0*ml/(mBd()+mK())*fT(qsq)/fp(qsq)*CTIm();}

double BdtoKll_amp::FARe(double qsq, double ml){
    return (C10Re() + C10RHRe());}
double BdtoKll_amp::FAIm(double qsq, double ml){
    return (C10Im() + C10RHIm());}

double BdtoKll_amp::FSRe(double qsq, double ml){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq)/fp(qsq)*(CSRe() + CSRHRe());}
double BdtoKll_amp::FSIm(double qsq, double ml){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq)/fp(qsq)*(CSIm() + CSRHIm());}

double BdtoKll_amp::FPRe(double qsq, double ml){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq)/fp(qsq)*(CPRe() + CPRHRe()) +
    ml*C10Re()*((pow(mBd(),2.0)-pow(mK(),2.0))/qsq*(fz(qsq)/fp(qsq)-1.0)-1.0);}
double BdtoKll_amp::FPIm(double qsq, double ml){
    return 1.0/2.0*(pow(mBd(),2.0)-pow(mK(),2.0))/(mb()-ms())*fz(qsq)/fp(qsq)*(CPIm() + CPRHIm()) +
    ml*C10Im()*((pow(mBd(),2.0)-pow(mK(),2.0))/qsq*(fz(qsq)/fp(qsq)-1.0)-1.0);}


double BdtoKll_amp::FTRe(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq)/fp(qsq)*CTRe();}
double BdtoKll_amp::FTIm(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq)/fp(qsq)*CTIm();}

double BdtoKll_amp::FT5Re(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq)/fp(qsq)*CT5Re();}
double BdtoKll_amp::FT5Im(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mBd()+mK())*fT(qsq)/fp(qsq)*CT5Im();}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//B->ll
double Btoll_amp :: ampSRe(double mBq, double mq, double ml1, double ml2){
    return ((C9Re()-C9RHRe())*(ml1-ml2) + (CSRe()-CSRHRe())*pow(mBq,2.0)/(mb()+mq));}
double Btoll_amp :: ampSIm(double mBq, double mq, double ml1, double ml2){
    return ((C9Im()-C9RHIm())*(ml1-ml2) + (CSIm()-CSRHIm())*pow(mBq,2.0)/(mb()+mq));}

double Btoll_amp :: ampPRe(double mBq, double mq, double ml1, double ml2){
    return ((C10Re()-C10RHRe())*(ml1+ml2) + (CPRe()-CPRHRe())*pow(mBq,2.0)/(mb()+mq));}
double Btoll_amp :: ampPIm(double mBq, double mq, double ml1, double ml2){
    return ((C10Im()-C10RHIm())*(ml1+ml2) + (CPIm()-CPRHIm())*pow(mBq,2.0)/(mb()+mq));}

double Btoll_amp :: Btoll_lambda(double mBq, double ml1, double ml2){
    return (pow(mBq,2.0)-pow(ml1-ml2,2.0))*(pow(mBq,2.0)-pow(ml1+ml2,2.0));}

double Btoll_amp :: y(double mBq){
    if (mBq==mBs()){
        return DGamma_sbar()/(2.0); }
    if (mBq==mBd()){
        return DGamma_dbar()/(2.0); }
    return 0.0;}

double Btoll_amp :: tauB(double mBq){
    if (mBq==mBs()){
        return tauBs(); }
    if (mBq==mBd()){
        return tauBd(); }
    return 0.0;}

















