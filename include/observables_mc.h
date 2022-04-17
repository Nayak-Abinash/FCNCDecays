#ifndef observables_mc_h
#define observables_mc_h

#include "amplitudes_mc.h"


//Bd->Kstr,ll:
class BdtoKstrll_obserr : public BdtoKstrll_amperr {
public:
    double J1s(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1c(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2s(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2c(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J3(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J4(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J5(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6s(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6c(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J7(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J8(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J9(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    //Conjugate J's
    double J1bs(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1bc(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2bs(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2bc(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J3b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J4b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J5b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6bs(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6bc(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J7b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J8b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J9b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    //Observables
    double diffWidth(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double diffWidthConj(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FL(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AFB(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P1(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P2(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P4p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P5p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P6p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P8p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S1(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S2(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S3(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S4(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S5(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S6(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S7(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S8(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S9(double qsq, double ml, double smwc[], double npwc[], double unv[]);
};

//Bs->phi,ll:
class Bstophill_obserr : public Bstophill_amperr {
public:
    double J1s(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1c(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2s(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2c(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J3(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J4(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J5(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6s(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6c(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J7(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J8(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J9(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    //Conjugate J's
    double J1bs(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1bc(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J1b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2bs(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2bc(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J2b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J3b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J4b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J5b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6bs(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6bc(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J6b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J7b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J8b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double J9b(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    //Observables
    double diffWidth(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double diffWidthConj(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double FL(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double AFB(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P1(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P2(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P4p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P5p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P6p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double P8p(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S1(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S2(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S3(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S4(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S5(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S6(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S7(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S8(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double S9(double qsq, double ml, double smwc[], double npwc[], double unv[]);
};

//Bd->K,ll
class BdtoKll_obserr : public BdtoKll_amperr {
public:
    double alNP(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double blNP(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double clNP(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double diffWidth(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double diffBrnch(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double diffAFB(double qsq, double ml, double smwc[], double npwc[], double unv[]);
    double diffFH(double qsq, double ml, double smwc[], double npwc[], double unv[]);
};


//B->ll
class Btoll_obserr : public Btoll_amperr {
public:
    double Btoll_nf(double mBq, double unv[]);
    double ADeltaGammaf(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double CorrctnFctr(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double BrInst(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double BrTimeIntgratd(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
    double efftau(double mBq, double mq, double ml1, double ml2, double smwc[], double npwc[], double unv[]);
};


#endif

