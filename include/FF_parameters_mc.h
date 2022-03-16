#ifndef FF_parameters_mc_h
#define FF_parameters_mc_h

#include "SM_parameters_mc.h"

//Bd->Kstr,ll
class BdtoKstrll_fferrpar : public smpar_mc {
public:
    BdtoKstrll_fferrpar();
    double tpd,tmd,tzd;
    double zd(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(double unv[]), Va1(double unv[]), Va2(double unv[]), A0a0(double unv[]), A0a1(double unv[]), A0a2(double unv[]),
        A1a0(double unv[]), A1a1(double unv[]), A1a2(double unv[]), A12a0(double unv[]), A12a1(double unv[]), A12a2(double unv[]),
        T1a0(double unv[]), T1a1(double unv[]), T1a2(double unv[]), T2a0(double unv[]), T2a1(double unv[]), T2a2(double unv[]),
        T23a0(double unv[]), T23a1(double unv[]), T23a2(double unv[]);
    double static cen_FF[MXdm]; double static unc_FF[MXdm]; double static chd_cov[MXdm][MXdm];

private:
    double m_PSbs,m_Vbs,m_Abs;
};

//Bs->phi,ll
class Bstophill_fferrpar : public smpar_mc {
public:
    Bstophill_fferrpar();
    double tps,tms,tzs;
    double zs(double qsq);
    double mPSbs(),mVbs(),mAbs();
    double Va0(double unv[]), Va1(double unv[]), Va2(double unv[]), A0a0(double unv[]), A0a1(double unv[]), A0a2(double unv[]),
        A1a0(double unv[]), A1a1(double unv[]), A1a2(double unv[]), A12a0(double unv[]), A12a1(double unv[]), A12a2(double unv[]),
        T1a0(double unv[]), T1a1(double unv[]), T1a2(double unv[]), T2a0(double unv[]), T2a1(double unv[]), T2a2(double unv[]),
        T23a0(double unv[]), T23a1(double unv[]), T23a2(double unv[]);
    double static cen_FF[MXdm]; double static unc_FF[MXdm]; double static chd_cov[MXdm][MXdm];

private:
    double m_PSbs,m_Vbs,m_Abs;
};

//Bd->K,ll
class BdtoKll_fferrpar : public smpar_mc {
public:
    BdtoKll_fferrpar();
    double tp,tm,tz;
    double z(double qsq);
    //LCSR+Lattice
    double mVbs(),mSbs();
    double fza1(double unv[]), fza2(double unv[]), fTa0(double unv[]), fTa1(double unv[]), fTa2(double unv[]),
        fpa0(double unv[]), fpa1(double unv[]), fpa2(double unv[]);
    double static cen_FF[MXdm]; double static unc_FF[MXdm]; double static chd_cov[MXdm][MXdm];

private:
    //LCSR+Lattice
    double m_Vbs,m_Sbs;
};

#endif
