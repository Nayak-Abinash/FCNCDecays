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
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <string>
using namespace std;

/*int mylength(double lst[], int iter, double a, double b){
    double num=0;
    for(int i=0; i<iter ; ++i)
    {
        if(a <= lst[i] && lst[i] < b)
            num = num+1;
    }
    return num;}
*/

int main(){
    smpar p;
    BdtoKstrll_obs o1;
    Bstophill_obs o2;
    Btoll_obs o3;
    BdtoKll_obs o4;
    BdtoKstrll_obserr eo1;
    double qsq;
    cout << "Type a value for qsq:";
    cin >> qsq;
    double obsval = o1.AFB(qsq,p.mmu());

    double simAFB = eo1.AFB(qsq,p.mmu());
    cout << simAFB << endl;

    cout << obsval << endl;

///////////Error///////////////
////68.2%ConfidenceLevel(via point counting from boundary values)////////////////
    int iter=10000;
    int sdcl = int(iter*15.9/100);
    double data[iter];
    for(int i=0; i<iter; i++)
        {
            data[i]= simAFB;
        }
    sort(data,data+iter);
    cout << "Observable at this qsq: " << obsval << "(+" << data[iter-sdcl-1]- obsval << "," << obsval-data[sdcl] << ")" << endl;

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

