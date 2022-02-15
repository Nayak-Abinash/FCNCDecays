#ifndef fferrpar_h
#define fferrpar_h

#include "myfun.h"
#include "smpar.h"


//Bd->Kstr,ll
class BdtoKstrll_fferrpar : public smpar, public myfun {
public:
    BdtoKstrll_fferrpar();
    double tpd,tmd,tzd;
    double zd(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(),Va1(),Va2(),A0a0(),A0a1(),A0a2(),A1a0(),A1a1(),A1a2(),A12a0(),A12a1(),A12a2(),T1a0(),T1a1(),T1a2(),T2a0(),T2a1(),
    T2a2(),T23a0(),T23a1(),T23a2();
    static double cen_FF[MXdm]; static double unc_FF[MXdm]; static double chd_cov[MXdm][MXdm];
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
    static double cen_FF[MXdm]; static double unc_FF[MXdm]; static double chd_cov[MXdm][MXdm];
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
    //LCSR+Lattice
    double mVbs(),mSbs();
    double fza1(),fza2(),fTa0(),fTa1(),fTa2(),fpa0(),fpa1(),fpa2();
    static double cen_FF[MXdm]; static double unc_FF[MXdm]; static double chd_cov[MXdm][MXdm];
private:
    //LCSR+Lattice
    double m_Vbs,m_Sbs;
    double fz_a1,fz_a2,fT_a0,fT_a1,fT_a2,fp_a0,fp_a1,fp_a2;
    double efz_a1,efz_a2,efT_a0,efT_a1,efT_a2,efp_a0,efp_a1,efp_a2;
};

#endif
