#include "smpar.h"
#include "myfun.h"
#include "ffpar.h"
#include "fferrpar.h"
#include "fffun.h"
#include "fferrfun.h"
#include "amp.h"
#include "amperr.h"
#include "obs.h"
#include "obserr.h"


int main(){
    BdtoKstrll_obs o1;
    BdtoKstrll_obserr eo1;

    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;

    cout << "(" << o1.V(qsq) << "," << eo1.ErV(qsq) << ")" << endl;

///////////Error///////////////
////68.2%ConfidenceLevel (point counting from boundary values)////////////////
    double cval = o1.V(qsq);
    int iter=10000;
    int sdcl = int(iter*15.9/100);
    double data[iter];
    for(int i=0; i<iter; i++)
        {
            data[i]= eo1.V(qsq);
        }
    sort(data,data+iter);
    cout << "Observable at this qsq: " << cval << "(+" << data[iter-sdcl-1]- cval << ",-" << cval-data[sdcl] << ")" << endl;

    return 0;
}

