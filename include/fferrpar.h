#ifndef fferrpar_h
#define fferrpar_h

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
};

#endif
