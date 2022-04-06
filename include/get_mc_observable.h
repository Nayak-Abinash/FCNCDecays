#ifndef get_mc_observable_h
#define get_mc_observable_h

#include "observables_mc.h"

class get_obserr : public BdtoKstrll_obserr, public Bstophill_obserr, public BdtoKll_obserr, public Btoll_obserr
{
public:
    void obsval(string s, double qsq, double smwc[], double npwc[]);
    //double static smwc[MXdm]; double static npwc[MXdm];
};

#endif // get_mc_observable_h

