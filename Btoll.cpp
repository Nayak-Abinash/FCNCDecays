#include <iostream>
#include <cmath>
using namespace std;

////header(#include):////

class par {
public:
    double
    pi,mb,md,ms,mBd,mBs,tauBd,tauBs,DGamma_dbar,DGamma_sbar,fBd,fBs,absVtbVtdStr,absVtbVtsStr,GF,alpha_e,mmu,me,
    mtau,C9Re,C10Re,C9Im,C9RHRe,C9RHIm,C10Im,C10RHRe,C10RHIm,CSRe,CSIm,CSRHRe,CSRHIm,CPRe,CPIm,CPRHRe,CPRHIm;
    par();
    };

class Btoll_amp : public par {
public:
    double ampSRe(double mBq, double mq, double ml1, double ml2);
    double ampSIm(double mBq, double mq, double ml1, double ml2);
    double ampPRe(double mBq, double mq, double ml1, double ml2);
    double ampPIm(double mBq, double mq, double ml1, double ml2);
    double lambda(double mBq, double ml1, double ml2);
    double y(double mBq);
    double tauB(double mBq);
    };
class Btoll_obs : public Btoll_amp {
public:
    double nf(double mBq);
    double ADeltaGammaf(double mBq, double mq, double ml1, double ml2);
    double CorrctnFctr(double mBq, double mq, double ml1, double ml2);
    double BrInst(double mBq, double mq, double ml1, double ml2);
    double BrTimeIntgratd(double mBq, double mq, double ml1, double ml2);
    double efftau(double mBq, double mq, double ml1, double ml2);
};

////source(#src):////

//parameter values:
par :: par() {
    pi=M_PI;
    mBd=5.27965,mBs=5.36688,mb=4.18,md=0.00467,ms=0.093,tauBd=1.519*1.52*pow(10.0,12.0),
    tauBs=1.527*1.52*pow(10.0,12.0),DGamma_dbar=0.001,DGamma_sbar=0.124,fBd=0.1905,fBs=0.2303,
    absVtbVtdStr=0.00853871,absVtbVtsStr=0.0408726;
    
    GF=1.1663787*pow(10.0,-5.0),alpha_e=1.0/127.944,mmu=0.1056583745,me=0.0005109989461,mtau=1.77686;
    
    //SMsubstitutions:
    C9Re=4.211,C10Re=-4.103;
    
    C9Im=0.0,C9RHRe=0.0,C9RHIm=0.0,C10Im=0.0,C10RHRe=0.0,C10RHIm=0.0,CSRe=0.0,CSIm=0.0,CSRHRe=0.0,CSRHIm=0.0,
    CPRe=0.0,CPIm=0.0,CPRHRe=0.0,CPRHIm=0.0;
    }
//amplitudes:(see,arXiv:1204.1735; arXiv:1204.1737; arXiv:1602.00881;)
double Btoll_amp :: ampSRe(double mBq, double mq, double ml1, double ml2){
        return ((C9Re-C9RHRe)*(ml1-ml2) + (CSRe-CSRHRe)*pow(mBq,2.0)/(mb+mq));
        }
double Btoll_amp :: ampSIm(double mBq, double mq, double ml1, double ml2){
        return ((C9Im-C9RHIm)*(ml1-ml2) + (CSIm-CSRHIm)*pow(mBq,2.0)/(mb+mq));
        }
double Btoll_amp :: ampPRe(double mBq, double mq, double ml1, double ml2){
        return ((C10Re-C10RHRe)*(ml1+ml2) + (CPRe-CPRHRe)*pow(mBq,2.0)/(mb+mq));
        }
double Btoll_amp :: ampPIm(double mBq, double mq, double ml1, double ml2){
        return ((C10Im-C10RHIm)*(ml1+ml2) + (CPIm-CPRHIm)*pow(mBq,2.0)/(mb+mq));
        }
double Btoll_amp :: lambda(double mBq, double ml1, double ml2){
        return (pow(mBq,2.0)-pow(ml1-ml2,2.0))*(pow(mBq,2.0)-pow(ml1+ml2,2.0));
        }
double Btoll_amp :: y(double mBq){
        if (mBq==mBs){
            return DGamma_sbar/(2.0); }
        if (mBq==mBd){
            return DGamma_dbar/(2.0); }
        return 0.0;
        }
double Btoll_amp :: tauB(double mBq){
        if (mBq==mBs){
            return tauBs; }
        if (mBq==mBd){
            return tauBd; }
        return 0.0;
        }
//observables:
double Btoll_obs :: nf(double mBq){
        if (mBq==mBd){
            return pow(GF*alpha_e*fBd*absVtbVtdStr,2.0)*tauBd/(64.0*pow(mBd*pi,3.0)); }
        if (mBq==mBs){
            return pow(GF*alpha_e*fBs*absVtbVtsStr,2.0)*tauBs/(64.0*pow(mBs*pi,3.0)); }
        return 0.0;
        }
double Btoll_obs :: ADeltaGammaf(double mBq, double mq, double ml1, double ml2){
    return (-pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0) + pow(ampSIm(mBq,mq,ml1,ml2),2.0) -
            pow(ampSRe(mBq,mq,ml1,ml2),2.0))/
    (pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0) + pow(ampSIm(mBq,mq,ml1,ml2),2.0) +
     pow(ampSRe(mBq,mq,ml1,ml2),2.0));
    }
double Btoll_obs :: CorrctnFctr(double mBq, double mq, double ml1, double ml2){
    return (1.0 - pow(y(mBq),2.0))/(1.0 + ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq));
    }
double Btoll_obs :: BrInst(double mBq, double mq, double ml1, double ml2){
    return ((pow(mBq,2.0) - pow(ml1 - ml2,2.0))*(pow(ampPIm(mBq,mq,ml1,ml2),2.0) + pow(ampPRe(mBq,mq,ml1,ml2),2.0))             + (pow(mBq,2.0) - pow(ml1 + ml2,2.0))*(pow(ampSIm(mBq,mq,ml1,ml2),2.0) +
        pow(ampSRe(mBq,mq,ml1,ml2),2.0)))*sqrt(lambda(mBq,ml1,ml2))*nf(mBq);
    }
double Btoll_obs :: BrTimeIntgratd(double mBq, double mq, double ml1, double ml2){
    return BrInst(mBq,mq,ml1,ml2)/CorrctnFctr(mBq,mq,ml1,ml2);
    }
double Btoll_obs :: efftau(double mBq, double mq, double ml1, double ml2){
    return (tauB(mBq)*(1.0 + 2.0*ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq) + pow(y(mBq),2.0)))/
    ((1.0 + ADeltaGammaf(mBq,mq,ml1,ml2)*y(mBq))*(1.0 - pow(y(mBq),2.0)));
    }

////main(#main):////


int main(){
    par p;
    Btoll_obs obs;
    cout << obs.BrTimeIntgratd(p.mBs,p.ms,p.mmu,p.mmu) << endl;
    cout << obs.efftau(p.mBs,p.ms,p.mmu,p.mmu) << endl;
    cout << obs.ADeltaGammaf(p.mBs,p.ms,p.mmu,p.mmu) << endl;
    return 0;
}











