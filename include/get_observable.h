#ifndef get_observable_h
#define get_observable_h

#include "observables.h"

class get_obs : public BdtoKstrll_obs, public Bstophill_obs, public BdtoKll_obs, public Btoll_obs
{
public:
    //void BtoVll_obsval(string s, double qsq);
    void BtoPll_obsval(string s, double qsq);
    void Btoll_obsval(string s);
};

#endif // get_observable_h
