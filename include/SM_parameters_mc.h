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
    //double static cen_wilsonSM[MXdm]; double static unc_wilsonSM[MXdm];
    double static chd_wilsonSM[MXdm][MXdm];
    //double static cen_wilsonNP[MXdm]; double static unc_wilsonNP[MXdm];
    double static chd_wilsonNP[MXdm][MXdm];
    double me(double unv[]),mmu(double unv[]),mtau(double unv[]);
    double mBd(double unv[]),mBs(double unv[]),mK(double unv[]),mKst(double unv[]),mphi(double unv[]);
    double tauBd(double unv[]),tauBs(double unv[]),DGamma_dbar(double unv[]),DGamma_sbar(double unv[]),fBd(double unv[]),fBs(double unv[]);
    double C1(double smwc[], double unv[]),C2(double smwc[], double unv[]),C3(double smwc[], double unv[]),C4(double smwc[], double unv[]),
        C5(double smwc[], double unv[]),C6(double smwc[], double unv[]),C7effRe(double smwc[], double unv[]), C7effIm(double smwc[], double unv[]),
        C8eff(double smwc[], double unv[]),C9Re(double smwc[], double unv[]),C9Im(double smwc[], double unv[]),C10Re(double smwc[], double unv[]),
        C10Im(double smwc[], double unv[]);
    double C7RHRe(double npwc[], double unv[]),C7RHIm(double npwc[], double unv[]),C9RHRe(double npwc[], double unv[]),C9RHIm(double npwc[], double unv[]),
        C10RHRe(double npwc[], double unv[]),C10RHIm(double npwc[], double unv[]),CSRe(double npwc[], double unv[]),CSIm(double npwc[], double unv[]),
        CSRHRe(double npwc[], double unv[]),CSRHIm(double npwc[], double unv[]),CPRe(double npwc[], double unv[]),CPIm(double npwc[], double unv[]),
        CPRHRe(double npwc[], double unv[]),CPRHIm(double npwc[], double unv[]),CTRe(double npwc[], double unv[]),CTIm(double npwc[], double unv[]),
        CT5Re(double npwc[], double unv[]),CT5Im(double npwc[], double unv[]);

private:
    double scale_mu;
    double G_F,alphae;
};

#endif

