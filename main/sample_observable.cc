//ProjectFiles
#include "observables.h"
#include "observables_mc.h"

int main()
{
    Btoll_obs o1; Btoll_obserr eo1;
    double unv[] = {/*A0_a0*/eo1.mnd_default(), /*A0_a1*/eo1.mnd_default(), /*A0_a2*/eo1.mnd_default(), /*A1_a0*/eo1.mnd_default(),
                            /*A1_a1*/eo1.mnd_default(), /*A1_a2*/eo1.mnd_default(), /*A12_a1*/eo1.mnd_default(), /*A12_a2*/eo1.mnd_default(),
                            /*V_a0*/eo1.mnd_default(), /*V_a1*/eo1.mnd_default(), /*V_a2*/eo1.mnd_default(), /*T1_a0*/eo1.mnd_default(),
                            /*T1_a1*/eo1.mnd_default(), /*T1_a2*/eo1.mnd_default(), /*T2_a1*/eo1.mnd_default(), /*T2_a2*/eo1.mnd_default(),
                            /*T23_a0*/eo1.mnd_default(), /*T23_a1*/eo1.mnd_default(), /*T23_a2*/eo1.mnd_default(), /*CKM_A*/eo1.mnd_default(),
                            /*CKM_lambda*/eo1.mnd_default(), /*CKM_rhobr*/eo1.mnd_default(), /*CKM_etabr*/eo1.mnd_default(), /*C1*/eo1.mnd_default(),
                            /*C2*/eo1.mnd_default(), /*C3*/eo1.mnd_default(), /*C4*/eo1.mnd_default(), /*C5*/eo1.mnd_default(),
                            /*C6*/eo1.mnd_default(), /*C7effRe*/eo1.mnd_default(), /*C7effIm*/eo1.mnd_default(), /*C8eff*/eo1.mnd_default(),
                            /*C9Re*/eo1.mnd_default(), /*C9Im*/eo1.mnd_default(), /*C10Re*/eo1.mnd_default(), /*C10Im*/eo1.mnd_default(),
                            /*C7RHRe*/eo1.mnd_default(), /*C7RHIm*/eo1.mnd_default(), /*C9RHRe*/eo1.mnd_default(), /*C9RHIm*/eo1.mnd_default(),
                            /*C10RHRe*/eo1.mnd_default(), /*C10RHIm*/eo1.mnd_default(), /*CSRe*/eo1.mnd_default(), /*CSIm*/eo1.mnd_default(),
                            /*CSRHRe*/eo1.mnd_default(), /*CSRHIm*/eo1.mnd_default(), /*CPRe*/eo1.mnd_default(), /*CPIm*/eo1.mnd_default(),
                            /*CPRHRe*/eo1.mnd_default(), /*CPRHIm*/eo1.mnd_default(), /*CTRe*/eo1.mnd_default(), /*CTIm*/eo1.mnd_default(),
                            /*CT5Re*/eo1.mnd_default(), /*CT5Im*/eo1.mnd_default(), /*md*/eo1.mnd_default(), /*mc*/eo1.mnd_default(),
                            /*ms*/eo1.mnd_default(), /*mb*/eo1.mnd_default(), /*me*/eo1.mnd_default(), /*mmu*/eo1.mnd_default(),
                            /*mtau*/eo1.mnd_default(), /*mBd*/eo1.mnd_default(), /*mBs*/eo1.mnd_default(), /*mK*/eo1.mnd_default(),
                            /*mKst*/eo1.mnd_default(), /*mphi*/eo1.mnd_default(), /*tauBd*/eo1.mnd_default(), /*tauBs*/eo1.mnd_default(),
                            /*DGamma_dbar*/eo1.mnd_default(), /*DGamma_sbar*/eo1.mnd_default(), /*fBd*/eo1.mnd_default(), /*fBs*/eo1.mnd_default()};

    cout << o1.BrTimeIntgratd(o1.mBs(),o1.ms(),o1.mmu(),o1.mmu()) << endl;
    cout << eo1.BrTimeIntgratd(eo1.mBs(unv),eo1.ms(unv),eo1.mmu(unv),eo1.mmu(unv),unv) << endl;
    return 0;
}

