#ifndef observables_binavg_h
#define observables_binavg_h

#include "observables_mc.h"

//#define nbin 20

//Bd->Kstr,ll:
class obs_binavg : public BdtoKstrll_obserr, public Bstophill_obserr, public BdtoKll_obserr
{
public:
    //Observables
    void BtoKstrll_obsavg(double smwc[], double npwc[]);
    //void BtoPll_obsavg(string s, double smwc[], double npwc[]);
};
#endif // observables_binavg_h


