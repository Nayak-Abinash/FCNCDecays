#ifndef observables_h
#define observables_h

#include "amplitudes.h"

//Bd->Kstr,ll:
class BdtoKstrll_obs : public BdtoKstrll_amp {
public:
    double J1s(double qsq, double ml);
    double J1c(double qsq, double ml);
    double J1(double qsq, double ml);
    double J2s(double qsq, double ml);
    double J2c(double qsq, double ml);
    double J2(double qsq, double ml);
    double J3(double qsq, double ml);
    double J4(double qsq, double ml);
    double J5(double qsq, double ml);
    double J6s(double qsq, double ml);
    double J6c(double qsq, double ml);
    double J6(double qsq, double ml);
    double J7(double qsq, double ml);
    double J8(double qsq, double ml);
    double J9(double qsq, double ml);
    //Conjugate J's
    double J1bs(double qsq, double ml);
    double J1bc(double qsq, double ml);
    double J1b(double qsq, double ml);
    double J2bs(double qsq, double ml);
    double J2bc(double qsq, double ml);
    double J2b(double qsq, double ml);
    double J3b(double qsq, double ml);
    double J4b(double qsq, double ml);
    double J5b(double qsq, double ml);
    double J6bs(double qsq, double ml);
    double J6bc(double qsq, double ml);
    double J6b(double qsq, double ml);
    double J7b(double qsq, double ml);
    double J8b(double qsq, double ml);
    double J9b(double qsq, double ml);
    //Observables
    double diffWidth(double qsq, double ml);
    double diffWidthConj(double qsq, double ml);
    double FL(double qsq, double ml);
    double AFB(double qsq, double ml);
    double P1(double qsq, double ml);
    double P2(double qsq, double ml);
    double P4p(double qsq, double ml);
    double P5p(double qsq, double ml);
    double P6p(double qsq, double ml);
    double P8p(double qsq, double ml);
    double S1(double qsq, double ml);
    double S2(double qsq, double ml);
    double S3(double qsq, double ml);
    double S4(double qsq, double ml);
    double S5(double qsq, double ml);
    double S6(double qsq, double ml);
    double S7(double qsq, double ml);
    double S8(double qsq, double ml);
    double S9(double qsq, double ml);
};

//Bs->phi,ll:
class Bstophill_obs : public Bstophill_amp {
public:
    double J1s(double qsq, double ml);
    double J1c(double qsq, double ml);
    double J1(double qsq, double ml);
    double J2s(double qsq, double ml);
    double J2c(double qsq, double ml);
    double J2(double qsq, double ml);
    double J3(double qsq, double ml);
    double J4(double qsq, double ml);
    double J5(double qsq, double ml);
    double J6s(double qsq, double ml);
    double J6c(double qsq, double ml);
    double J6(double qsq, double ml);
    double J7(double qsq, double ml);
    double J8(double qsq, double ml);
    double J9(double qsq, double ml);
    //Conjugate J's
    double J1bs(double qsq, double ml);
    double J1bc(double qsq, double ml);
    double J1b(double qsq, double ml);
    double J2bs(double qsq, double ml);
    double J2bc(double qsq, double ml);
    double J2b(double qsq, double ml);
    double J3b(double qsq, double ml);
    double J4b(double qsq, double ml);
    double J5b(double qsq, double ml);
    double J6bs(double qsq, double ml);
    double J6bc(double qsq, double ml);
    double J6b(double qsq, double ml);
    double J7b(double qsq, double ml);
    double J8b(double qsq, double ml);
    double J9b(double qsq, double ml);
    //Observables
    double diffWidth(double qsq, double ml);
    double diffWidthConj(double qsq, double ml);
    double FL(double qsq, double ml);
    double AFB(double qsq, double ml);
    double P1(double qsq, double ml);
    double P2(double qsq, double ml);
    double P4p(double qsq, double ml);
    double P5p(double qsq, double ml);
    double P6p(double qsq, double ml);
    double P8p(double qsq, double ml);
    double S1(double qsq, double ml);
    double S2(double qsq, double ml);
    double S3(double qsq, double ml);
    double S4(double qsq, double ml);
    double S5(double qsq, double ml);
    double S6(double qsq, double ml);
    double S7(double qsq, double ml);
    double S8(double qsq, double ml);
    double S9(double qsq, double ml);
};

//Bd->K,ll
class BdtoKll_obs : public BdtoKll_amp {
public:
    double alNP(double qsq, double ml);
    double blNP(double qsq, double ml);
    double clNP(double qsq, double ml);
    double diffWidth(double qsq, double ml);
    double diffBrnch(double qsq, double ml);
    double diffAFB(double qsq, double ml);
    double diffFH(double qsq, double ml);
};

//B->ll
class Btoll_obs : public Btoll_amp {
public:
    double Btoll_nf(double mBq);
    double ADeltaGammaf(double mBq, double mq, double ml1, double ml2);
    double CorrctnFctr(double mBq, double mq, double ml1, double ml2);
    double BrInst(double mBq, double mq, double ml1, double ml2);
    double BrTimeIntgratd(double mBq, double mq, double ml1, double ml2);
    double efftau(double mBq, double mq, double ml1, double ml2);
};

#endif

