#include "observables_binavg.h"

//Bd->Kstr,ll:
void obs_binavg::BtoKstrll_obsavg(double smwc[], double npwc[]){
    //B->Vll
    BdtoKstrll_obserr eo1;
    double bin_a[] = {0.10, 1.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 11.00, 11.75, 15.0, 16.0, 17.0, 18.0};
    double bin_b[] = {0.98, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 11.75, 12.50, 16.0, 17.0, 18.0, 19.0};

    double unv[72]; int bin_num = 3; int bin_prolf = 5;

    int n_sample = 20;
    double q2, dq2;
    dq2 = (bin_b[bin_num] - bin_a[bin_num])/bin_prolf;

    double simdata_afb[n_sample], simdata_fl[n_sample];
    double bindata_afbmean[bin_prolf], bindata_afbsd[bin_prolf], bindata_flmean[bin_prolf], bindata_flsd[bin_prolf];
    double bindatacov_afb_fl[bin_prolf];

    for(int i=0; i<bin_prolf; i++){
        q2 = bin_a[bin_num] + dq2*(2*i + 1)/2;
        for(int j=0; j<n_sample; j++){
            for(int k=0; k<72; k++){unv[k] = eo1.mnd_default();};
            simdata_afb[j] = eo1.AFB(q2,eo1.mmu(unv),smwc,npwc,unv);
            simdata_fl[j] = eo1.FL(q2,eo1.mmu(unv),smwc,npwc,unv);
            };
        bindata_afbmean[i] = eo1.mean_model(simdata_afb,n_sample);
        bindata_afbsd[i] = eo1.sd_model(simdata_afb,n_sample);
        bindata_flmean[i] = eo1.mean_model(simdata_fl,n_sample);
        bindata_flsd[i] = eo1.sd_model(simdata_fl,n_sample);
        bindatacov_afb_fl[i] = eo1.cov_model(simdata_afb,simdata_fl,n_sample);
    }
    double binavg_afbmean = eo1.mean_model(bindata_afbmean,bin_prolf);
    double binavg_afbsd = eo1.sd_model(bindata_afbsd,bin_prolf);
    double binavg_flmean = eo1.mean_model(bindata_flmean,bin_prolf);
    double binavg_flsd = eo1.sd_model(bindata_flsd,bin_prolf);
    double binavgcov_afb_fl = eo1.mean_model(bindatacov_afb_fl,bin_prolf);

    cout << binavg_afbmean << endl;
    cout << binavg_afbsd << endl;
    cout << binavg_flmean << endl;
    cout << binavg_flsd << endl;
    cout << binavgcov_afb_fl << endl;
}
