#include "obs.h"

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

double BdtoKstrll_obs::J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double BdtoKstrll_obs::J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

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

double BdtoKstrll_obs::diffWidth(double qsq, double ml){
    return (3.0*(J1c(qsq,ml) + 2.0*J1s(qsq,ml) + (-J2c(qsq,ml) - 2.0*J2s(qsq,ml))/3.0))/4.0;}

double BdtoKstrll_obs::FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0));}

double BdtoKstrll_obs::AFB(double qsq, double ml){
    return (3.0*(J6c(qsq,ml) + 2.0*J6s(qsq,ml)))/(8.0*diffWidth(qsq,ml));}
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

double Bstophill_obs::J2s(double qsq, double ml){
    return ((pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
             pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0))*
            pow(betal(qsq,ml),2.0))/4.0 ;}

double Bstophill_obs::J2c(double qsq, double ml){
    return -((pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))*
             pow(betal(qsq,ml),2.0)) ;}

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

double Bstophill_obs::diffWidth(double qsq, double ml){
    return (3.0*(J1c(qsq,ml) + 2.0*J1s(qsq,ml) + (-J2c(qsq,ml) - 2.0*J2s(qsq,ml))/3.0))/4.0;}

double Bstophill_obs::FL(double qsq, double ml){
    return (pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0))/
    (pow(AaLIm(qsq,ml),2.0) + pow(AaLRe(qsq,ml),2.0) + pow(AaRIm(qsq,ml),2.0) + pow(AaRRe(qsq,ml),2.0) +
     pow(ApLIm(qsq,ml),2.0) + pow(ApLRe(qsq,ml),2.0) + pow(ApRIm(qsq,ml),2.0) + pow(ApRRe(qsq,ml),2.0) +
     pow(AzLIm(qsq,ml),2.0) + pow(AzLRe(qsq,ml),2.0) + pow(AzRIm(qsq,ml),2.0) + pow(AzRRe(qsq,ml),2.0));}

double Bstophill_obs::AFB(double qsq, double ml){
    return (3.0*(J6c(qsq,ml) + 2.0*J6s(qsq,ml)))/(8.0*diffWidth(qsq,ml));}
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











