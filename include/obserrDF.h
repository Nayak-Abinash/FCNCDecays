#ifndef obserrDF_h
#define obserrDF_h

#include "myfun.h"
#include "smpar.h"
#include "ffpar.h"
#include "fffun.h"
#include "amp.h"

//Bd->K,ll
class BdtoKll_obserrDF : public BdtoKll_amp {
public:
    double Er_diffWidth(double qsq, double ml);
    //double diffBrnch(double qsq, double ml);
    //double diffAFB(double qsq, double ml);
    //double diffFH(double qsq, double ml);
};

#endif // obserrDF_h
