#ifndef ffpar_h
#define ffpar_h

#include "smpar.h"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

//Bd->Kstr,ll
class BdtoKstrll_ffpar : public smpar {
public:
    BdtoKstrll_ffpar();
    double tpd,tmd,tzd;
    double zd(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),T2a2(),T23a0(),T23a1(),T23a2();
private:
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,T2_a2,T23_a0,T23_a1,T23_a2;
};

//Bs->phi,ll
class Bstophill_ffpar : public smpar {
public:
    Bstophill_ffpar();
    double tps,tms,tzs;
    double zs(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),T2a2(),T23a0(),T23a1(),T23a2();
private:
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,T2_a2,T23_a0,T23_a1,T23_a2;
};

//Bd->K,ll
class BdtoKll_ffpar : public smpar {
public:
    BdtoKll_ffpar();
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
    //Lattice
    double lat_m_Vbs,lat_m_Sbs;
    double b0_p,b1_p,b2_p,b0_z,b1_z,b2_z,b0_t,b1_t,b2_t;
};

#endif