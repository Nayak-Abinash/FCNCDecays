#ifndef get_mcplot_h
#define get_mcplot_h

#include "observables_mc.h"

class get_obsplot : public BdtoKstrll_obserr, public Bstophill_obserr, public BdtoKll_obserr, public Btoll_obserr
{
public:
    void BtoVll_obsplot(string s, double smwc[], double npwc[]);
    void BtoPll_obsplot(string s, double smwc[], double npwc[]);
    //void Btoll_obsplot(string s, double smwc[], double npwc[]);
};

#endif // get_mc_observable_h

