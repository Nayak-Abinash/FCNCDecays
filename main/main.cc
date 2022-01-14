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

////Histogram/////////////////////
    /*double minv = *min_element(data,data+iter);
    double maxv = *max_element(data,data+iter);
    double nbin = 20;
    double init_val = minv;
    for (int i=0; i<nbin; ++i)
        {
            int npnt = mylength(data,iter,minv+(i)*(maxv-minv)/nbin,minv+(i+1)*(maxv-minv)/nbin);
            if(i<10)
                {
                    cout << "-" << "bin-0" << i << ": ";
                    cout << string(npnt/(iter/100),'*') << "\t \t" << npnt <<endl;
                }
            else
                {
                    cout << "-" << "bin-" << i << ": ";
                    cout << string(npnt/(iter/100),'*') << "\t \t" << npnt <<endl;
                }
        }

    cout << "(" << maxv << "," << minv << ")" << endl;
////OutPut///////////////////////
    cout << "The Observable value at this qsq:" << o1.AFB(qsq,p.mmu()) << endl;
    cout << "The Mean and Standard Deviation at this qsq:" << "(" << meanv << "," << stdv <<")" << endl;*/

    return 0;
}

