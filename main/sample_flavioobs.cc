//ProjectFiles
#include "get_observable.h"

int main()
{
    get_obs obj;
    string s; double qsq;
    cout << "Choose one of the following: B->Vll, B->Pll, B->ll" << endl;
    cin >> s;
    if(s=="B->Vll")
        obj.BtoPll_obsval(s,qsq);
    else if(s=="B->Pll")
        obj.BtoPll_obsval(s,qsq);
    else if(s=="B->ll")
        obj.Btoll_obsval(s);
    return 0;
}
