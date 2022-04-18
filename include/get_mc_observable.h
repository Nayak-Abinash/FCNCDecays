#ifndef get_mc_observable_h
#define get_mc_observable_h

#include "observables_mc.h"

class get_obserr : public BdtoKstrll_obserr, public Bstophill_obserr, public BdtoKll_obserr, public Btoll_obserr
{
public:
    void BtoVll_obsval(string s, double qsq, double smwc[], double npwc[]);
    void BtoPll_obsval(string s, double qsq, double smwc[], double npwc[]);
    void Btoll_obsval(string s, double smwc[], double npwc[]);
};

#endif // get_mc_observable_h

