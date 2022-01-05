#ifndef fferrpar_h
#define fferrpar_h

#include "smpar.h"
#include "myfun.h"
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
#include <iostream>
using namespace std;

//Bd->Kstr,ll
class BdtoKstrll_fferrpar : public smpar, public myfun {
public:
    BdtoKstrll_fferrpar();
    double tpd,tmd,tzd;
    double zd(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),
    T2a2(),T23a0(),T23a1(),T23a2();
private:
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,
    T2_a2,T23_a0,T23_a1,T23_a2;
    double  eV_a0,eV_a1,eV_a2,eA0_a0,eA0_a1,eA0_a2,eA1_a0,eA1_a1,eA1_a2,eA12_a0,eA12_a1,eA12_a2,eT1_a0,eT1_a1,eT1_a2,eT2_a0,eT2_a1,
    eT2_a2,eT23_a0,eT23_a1,eT23_a2;
};

//Bs->phi,ll
class Bstophill_fferrpar : public smpar, public myfun {
public:
    Bstophill_fferrpar();
    double tps,tms,tzs;
    double zs(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),
    T2a2(),T23a0(),T23a1(),T23a2();
private:
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,
    T2_a2,T23_a0,T23_a1,T23_a2;
    double  eV_a0,eV_a1,eV_a2,eA0_a0,eA0_a1,eA0_a2,eA1_a0,eA1_a1,eA1_a2,eA12_a0,eA12_a1,eA12_a2,eT1_a0,eT1_a1,eT1_a2,eT2_a0,eT2_a1,
    eT2_a2,eT23_a0,eT23_a1,eT23_a2;
};

//Bd->K,ll
class BdtoKll_fferrpar : public smpar, public myfun {
public:
    BdtoKll_fferrpar();
    double tp,tm,tz;
    double z(double qsq);
    //LCSR
    double lcsr_mVbs(),lcsr_mSbs();
    double cfp0(),cfp1(),cfp2(),cfz0(),cfz1(),cfz2(),cft0(),cft1(),cft2();
    //Lattice
    double lat_mVbs(),lat_mSbs();
    double b0p(),b1p(),b2p(),b0z(),b1z(),b2z(),b0t(),b1t(),b2t();
private:
    //LCSR
    double lcsr_m_Vbs,lcsr_m_Sbs;
    double c_fp0,c_fp1,c_fp2,c_fz0,c_fz1,c_fz2,c_ft0,c_ft1,c_ft2;
    double ec_fp0,ec_fp1,ec_fp2,ec_fz0,ec_fz1,ec_fz2,ec_ft0,ec_ft1,ec_ft2;
    //Lattice
    double lat_m_Vbs,lat_m_Sbs;
    double b0_p,b1_p,b2_p,b0_z,b1_z,b2_z,b0_t,b1_t,b2_t;
    double eb0_p,eb1_p,eb2_p,eb0_z,eb1_z,eb2_z,eb0_t,eb1_t,eb2_t;
};

#endif
