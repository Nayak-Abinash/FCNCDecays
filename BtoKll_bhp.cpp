#include<iostream>
#include<string>
#include<cmath>
using namespace std;

//parameters:
double pi=M_PI;
double mB=5.27958,mK=0.497614,mu=4.18,mc=1.27,mb=4.18,ms=0.093;
double C1=-0.257,C2=1.009,C3=-0.005,C4=-0.078,C5=0.000,C6=0.001,C7effRe=-0.304,C9Re=4.211,C10effRe=-4.103,C10Re=-4.103,tauB=1.519*1.52*pow(10.0,12.0),alpha_e=1.0/133.0,GF=1.16637*pow(10.0,-5.0),absVtbVtsStr=0.0409,
    me=0.511*pow(10.0,-3.0),mmu=105.658*pow(10.0,-3.0);
//SMsub:
double C7effIm=0.0,C7RHRe=0.0,C7RHIm=0.0,C9Im=0.0,C9RHRe=0.0,C9RHIm=0.0,C10Im=0.0,C10effIm=0.0,C10RHRe=0.0,C10RHIm=0.0,CPRe=0.0,CPIm=0.0,CPRHRe=0.0,CPRHIm=0.0,CSRe=0.0,CSIm=0.0,CSRHRe=0.0,CSRHIm=0.0,CTRe=0.0,CTIm=0.0,CT5Re=0.0,CT5Im=0.0;

double tp = pow(mB+mK,2.0),tm = pow(mB-mK,2.0),tz = tp*(1.0-sqrt(1.0-tm/tp));
double z(double qsq){
    return (sqrt(tp-qsq)-sqrt(tp-tz))/(sqrt(tp-qsq)+sqrt(tp-tz));   }

//extrapolatedLCSRformfactors:(see, arXiv:1811.00983)
double lcs_mVBs=5.412,lcs_mSBs=5.630;
double cfp0=0.32909,cfp1=-0.866947,cfp2=0.00609567,cfz0=0.0,cfz1=0.195117,cfz2=-0.446126,cft0=0.299383,cft1=-0.773546,cft2=0.00955438;

double lcsfp(double qsq){
    return 1.0/(1.0-qsq/pow(lcs_mVBs,2.0))*( cfp0*pow(z(qsq)-z(0.0),0.0) + cfp1*pow(z(qsq)-z(0.0),1.0) + cfp2*pow(z(qsq)-z(0.0),2.0) ); }

double lcsfz(double qsq){
    return 1.0/(1.0-qsq/pow(lcs_mSBs,2.0))*( cfz0*pow(z(qsq)-z(0.0),0.0) + cfz1*pow(z(qsq)-z(0.0),1.0) + cfz2*pow(z(qsq)-z(0.0),2.0) ); }

double lcsft(double qsq){
    return 1.0/(1.0-qsq/pow(lcs_mVBs,2.0))*( cft0*pow(z(qsq)-z(0.0),0.0) + cft1*pow(z(qsq)-z(0.0),1.0) + cft2*pow(z(qsq)-z(0.0),2.0) ); }

//extrapolatedlatticeformfactors:(see, arXiv:1509.06235)
double lat_mVBs=5.4154,lat_mSBs=5.711;
double b0p=0.466,b1p=-0.885,b2p=-0.213,b0z=0.292,b1z=0.281,b2z=0.150,b0t=0.460,b1t=-1.089,b2t=-1.114;

double latfp(double qsq){
    return 1.0/(1.0-qsq/pow(lat_mVBs,2.0))*(b0p + b1p*(z(qsq)-1.0/3.0*pow(z(qsq),3.0)) + b2p*(pow(z(qsq),2.0) + 2.0/3.0*pow(z(qsq),3.0))); }

double latfz(double qsq){
    return 1.0/(1.0-qsq/pow(lat_mSBs,2.0))*(b0z + b1z*z(qsq) + b2z*pow(z(qsq),2.0)); }

double latft(double qsq){
    return 1.0/(1.0-qsq/pow(lat_mVBs,2.0))*(b0t + b1t*(z(qsq)-1.0/3.0*pow(z(qsq),3.0)) + b2t*(pow(z(qsq),2.0) + 2.0/3.0*pow(z(qsq),3.0))); }

//formfactorfunction:
double fp(double qsq){
    if (qsq<4.0*pow(mc,2.0)){
        return lcsfp(qsq); }
    else {
        return latfp(qsq); } }

double fz(double qsq){
    if (qsq<4.0*pow(mc,2.0)){
        return lcsfz(qsq); }
    else {
        return latfz(qsq); } }

double ft(double qsq){
    if (qsq<4.0*pow(mc,2.0)){
        return lcsft(qsq); }
    else {
        return latft(qsq); } }

//C9eff:(see,arXiv:0709.4174)
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

//amplitudes:(see,arXiv:0709.4174)
double betal(double qsq, double ml){
    return sqrt(1.0-4.0*pow(ml,2.0)/qsq); }
double lambda(double qsq){
    return pow(qsq,2.0)+pow(mB,4.0)+pow(mK,4.0)-2.0*(pow(mB*mK,2.0)+qsq*pow(mB,2.0)+qsq*pow(mK,2.0)); }
double nf(double qsq, double ml){
    return pow(GF*alpha_e*absVtbVtsStr,2.0)/(512.0*pow(pi,5.0)*pow(mB,3.0))*betal(qsq,ml)*sqrt(lambda(qsq)); }

//SM+NP(S(CiS)+PS(CiPS)+T(CiT)+RHC(CiRH))(complex C_i's):
double xiP(double qsq){
    return fp(qsq); }
double TauPRe(double qsq){
    return xiP(qsq)*(C7effRe + mB/(2.0*mb)*YRe(qsq)); }
double TauPIm(double qsq){
    return xiP(qsq)*(C7effIm + mB/(2.0*mb)*YIm(qsq)); }
    
double FVRe(double qsq, double ml){
    return (C9Re + C9RHRe)+ 2.0*mb/mB*TauPRe(qsq)/xiP(qsq) + 8.0*ml/(mB+mK)*ft(qsq)/fp(qsq)*CTRe; }
double FVIm(double qsq, double ml){
    return (C9Im + C9RHIm)+ 2.0*mb/mB*TauPIm(qsq)/xiP(qsq) + 8.0*ml/(mB+mK)*ft(qsq)/fp(qsq)*CTIm; }


double FARe(double qsq, double ml){
    return (C10Re + C10RHRe); }
double FAIm(double qsq, double ml){
    return (C10Im + C10RHIm); }

double FSRe(double qsq, double ml){
    return 1.0/2.0*(pow(mB,2.0)-pow(mK,2.0))/(mb-ms)*fz(qsq)/fp(qsq)*(CSRe + CSRHRe); }
double FSIm(double qsq, double ml){
    return 1.0/2.0*(pow(mB,2.0)-pow(mK,2.0))/(mb-ms)*fz(qsq)/fp(qsq)*(CSIm + CSRHIm); }

double FPRe(double qsq, double ml){
    return 1.0/2.0*(pow(mB,2.0)-pow(mK,2.0))/(mb-ms)*fz(qsq)/fp(qsq)*(CPRe + CPRHRe) +
        ml*C10Re*((pow(mB,2.0)-pow(mK,2.0))/qsq*(fz(qsq)/fp(qsq)-1.0)-1.0); }
double FPIm(double qsq, double ml){
    return 1.0/2.0*(pow(mB,2.0)-pow(mK,2.0))/(mb-ms)*fz(qsq)/fp(qsq)*(CPIm + CPRHIm) +
        ml*C10Im*((pow(mB,2.0)-pow(mK,2.0))/qsq*(fz(qsq)/fp(qsq)-1.0)-1.0); }


double FTRe(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mB+mK)*ft(qsq)/fp(qsq)*CTRe; }
double FTIm(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mB+mK)*ft(qsq)/fp(qsq)*CTIm; }

double FT5Re(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mB+mK)*ft(qsq)/fp(qsq)*CT5Re; }
double FT5Im(double qsq, double ml){
    return 2.0*sqrt(lambda(qsq))*betal(qsq,ml)/(mB+mK)*ft(qsq)/fp(qsq)*CT5Im; }

//angular coefficients:
double alNP(double qsq, double ml){
    return (4.0*pow(mB,2.0)*pow(ml,2.0)*(pow(FAIm(qsq,ml),2.0) + pow(FARe(qsq,ml),2.0)) +
            2.0*ml*(pow(mB,2.0) - pow(mK,2.0) + qsq)*(FAIm(qsq,ml)*FPIm(qsq,ml) + FARe(qsq,ml)*FPRe(qsq,ml)) +
            qsq*(pow(FPIm(qsq,ml),2.0) + pow(FPRe(qsq,ml),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FSIm(qsq,ml),2.0) + pow(FSRe(qsq,ml),2.0))) +
            ((pow(FAIm(qsq,ml),2.0) + pow(FARe(qsq,ml),2.0) + pow(FVIm(qsq,ml),2.0) + pow(FVRe(qsq,ml),2.0))*
             lambda(qsq))/4.0)*nf(qsq,ml)*pow(xiP(qsq),2.0); }

double blNP(double qsq, double ml){
    return 2.0*(qsq*(FPIm(qsq,ml)*FT5Im(qsq,ml) + FPRe(qsq,ml)*FT5Re(qsq,ml) +
                   pow(betal(qsq,ml),2.0)*(FSIm(qsq,ml)*FTIm(qsq,ml) + FSRe(qsq,ml)*FTRe(qsq,ml))) +
              ml*((pow(mB,2.0) - pow(mK,2.0) + qsq)*(FAIm(qsq,ml)*FT5Im(qsq,ml) + FARe(qsq,ml)*FT5Re(qsq,ml)) +
                  betal(qsq,ml)*(FSIm(qsq,ml)*FVIm(qsq,ml) + FSRe(qsq,ml)*FVRe(qsq,ml))*sqrt(lambda(qsq))))*
    nf(qsq,ml)*pow(xiP(qsq),2.0); }

double clNP(double qsq, double ml){
    return (qsq*(pow(FT5Im(qsq,ml),2.0) + pow(FT5Re(qsq,ml),2.0) +
                 pow(betal(qsq,ml),2.0)*(pow(FTIm(qsq,ml),2.0) + pow(FTRe(qsq,ml),2.0))) +
            2.0*ml*betal(qsq,ml)*(FTIm(qsq,ml)*FVIm(qsq,ml) + FTRe(qsq,ml)*FVRe(qsq,ml))*sqrt(lambda(qsq)) -
            (pow(betal(qsq,ml),2.0)*(pow(FAIm(qsq,ml),2.0) + pow(FARe(qsq,ml),2.0) + pow(FVIm(qsq,ml),2.0) +
                                     pow(FVRe(qsq,ml),2.0))*lambda(qsq))/4.0)*nf(qsq,ml)*pow(xiP(qsq),2.0); }

//observables:
double diffWidth(double qsq, double ml){
    return 2.0*(alNP(qsq,ml) + 1.0/3.0*clNP(qsq,ml)); }

double diffBrnch(double qsq, double ml){
    return tauB/2.0*diffWidth(qsq,ml); }

double diffAFB(double qsq, double ml){
    return blNP(qsq,ml)/diffWidth(qsq,ml); }

double diffFH(double qsq, double ml){
    return 2.0*(alNP(qsq,ml)+clNP(qsq,ml))/diffWidth(qsq,ml); }


int main(){
    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;
    cout << "Observable values at this point are:"
    << diffBrnch(qsq,mmu) << "," << diffAFB(qsq,mmu) << "," << diffFH(qsq,mmu) << endl;
    return 0;
}










    









