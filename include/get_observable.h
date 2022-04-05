#ifndef get_observable_h
#define get_observable_h

#include "observables.h"

class get_obs : public BdtoKstrll_obs, public Bstophill_obs, public BdtoKll_obs, public Btoll_obs
{
public:
    //get_obs();
    void obsval(string s, double qsq);
};

#endif // get_observable_h
