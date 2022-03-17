#ifndef SM_parameters_mc_h
#define SM_parameters_mc_h

#include "ckm_elements_mc.h"

class smpar_mc : public ckm_elements_mc
{
public:
    smpar_mc ();
    double pi;
    double mu(), GF(),alpha_e();
    double md(double unv[]),mc(double unv[]),ms(double unv[]),mb(double unv[]);
    double static cen_smpar[MXdm]; double static unc_smpar[MXdm]; double static chd_smpar[MXdm][MXdm];
    double static cen_wilsonSM[MXdm]; double static unc_wilsonSM[MXdm]; double static chd_wilsonSM[MXdm][MXdm];
    double static cen_wilsonNP[MXdm]; double static unc_wilsonNP[MXdm]; double static chd_wilsonNP[MXdm][MXdm];
    double me(double unv[]),mmu(double unv[]),mtau(double unv[]);
    double mBd(double unv[]),mBs(double unv[]),mK(double unv[]),mKst(double unv[]),mphi(double unv[]);
    double tauBd(double unv[]),tauBs(double unv[]),DGamma_dbar(double unv[]),DGamma_sbar(double unv[]),fBd(double unv[]),fBs(double unv[]);
    double C1(double unv[]),C2(double unv[]),C3(double unv[]),C4(double unv[]),C5(double unv[]),C6(double unv[]),C7effRe(double unv[]),
        C7effIm(double unv[]),C8eff(double unv[]),C9Re(double unv[]),C9Im(double unv[]),C10Re(double unv[]),C10Im(double unv[]);
    double C7RHRe(double unv[]),C7RHIm(double unv[]),C9RHRe(double unv[]),C9RHIm(double unv[]),C10RHRe(double unv[]),C10RHIm(double unv[]),
        CSRe(double unv[]),CSIm(double unv[]),CSRHRe(double unv[]),CSRHIm(double unv[]),CPRe(double unv[]),CPIm(double unv[]),CPRHRe(double unv[]),
        CPRHIm(double unv[]),CTRe(double unv[]),CTIm(double unv[]),CT5Re(double unv[]),CT5Im(double unv[]);

private:
    double scale_mu;
    double G_F,alphae;
};

#endif

