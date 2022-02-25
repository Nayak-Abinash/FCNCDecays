#ifndef ffpar_h
#define ffpar_h

#include "smpar.h"

//Bd->Kstr,ll
class BdtoKstrll_ffpar : public smpar, public ref_fun {
public:
    BdtoKstrll_ffpar();
    double tpd,tmd,tzd;
    double zd(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),
    T2a2(),T23a0(),T23a1(),T23a2();
    static double cov_V[mxdm][mxdm]; static double cov_A0[mxdm][mxdm]; static double cov_A1[mxdm][mxdm]; static double cov_A12[mxdm][mxdm];
    static double cov_T1[mxdm][mxdm]; static double cov_T2[mxdm][mxdm]; static double cov_T23[mxdm][mxdm];
private:
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,T2_a2,T23_a0,T23_a1,T23_a2;
};

//Bs->phi,ll
class Bstophill_ffpar : public smpar, public ref_fun {
public:
    Bstophill_ffpar();
    double tps,tms,tzs;
    double zs(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),
    T2a2(),T23a0(),T23a1(),T23a2();
    static double cov_V[mxdm][mxdm]; static double cov_A0[mxdm][mxdm]; static double cov_A1[mxdm][mxdm]; static double cov_A12[mxdm][mxdm];
    static double cov_T1[mxdm][mxdm]; static double cov_T2[mxdm][mxdm]; static double cov_T23[mxdm][mxdm];
private:
    double m_PSbs,m_Vbs,m_Abs;
    double  V_a0,V_a1,V_a2,A0_a0,A0_a1,A0_a2,A1_a0,A1_a1,A1_a2,A12_a0,A12_a1,A12_a2,T1_a0,T1_a1,T1_a2,T2_a0,T2_a1,T2_a2,T23_a0,T23_a1,T23_a2;
};

//Bd->K,ll
class BdtoKll_ffpar : public smpar, public ref_fun {
public:
    BdtoKll_ffpar();
    double tp,tm,tz;
    double z(double qsq);
    //LCSR+Lattice
    double mVbs(),mSbs();
    double fza1(),fza2(),fTa0(),fTa1(),fTa2(),fpa0(),fpa1(),fpa2();
    static double cov_fz[mxdm][mxdm]; static double cov_fT[mxdm][mxdm]; static double cov_fp[mxdm][mxdm];
private:
    //LCSR+Lattice
    double m_Vbs,m_Sbs;
    double fz_a1,fz_a2,fT_a0,fT_a1,fT_a2,fp_a0,fp_a1,fp_a2;
};

#endif
