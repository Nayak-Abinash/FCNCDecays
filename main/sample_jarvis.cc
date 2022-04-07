#include "get_observable.h"
#include "get_mc_observable.h"
#include "get_mcplot.h"

double smwc[] = {/*C_1*/-0.257, /*C_2*/1.009, /*C_3*/-0.005, /*C_4*/-0.078, /*C_5*/0.000, /*C_6*/0.001, /*C_7effRe*/-0.304, /*C_7effIm*/0.0,
                                    /*C_8eff*/-0.167, /*C_9Re*/4.211, /*C_9Im*/0.0, /*C_10Re*/-4.103, /*C_10Im*/0.0};

double npwc[] = {/*C_7RHRe*/-0.006, /*C_7RHIm*/0.0, /*C_9RHRe*/0.0, /*C_9RHIm*/0.0, /*C_10RHRe*/0.0, /*C_10RHIm*/0.0, /*C_SRe*/0.0,
                                    /*C_SIm*/0.0, /*C_SRHRe*/0.0, /*C_SRHIm*/0.0, /*C_PRe*/0.0, /*C_PIm*/0.0, /*C_PRHRe*/0.0, /*C_PRHIm*/0.0,
                                    /*C_TRe*/0.0, /*C_TIm*/0.0, /*C_T5Re*/0.0, /*C_T5Im*/0.0};

int main()
{
    get_obs obj1; get_obserr obj2; get_obsplot obj3;
    string s; double qsq;
    cout << "Welcome to FlavBit!" << endl;
    cout << "Choose one of the following decay modes: B->Vll, B->Pll, B->ll" << endl;
    cin >> s;
    if(s=="B->Pll")
    {
        cout << "Specify an action: obsval, obserr or obsplot" << endl;
        cin >> s;
        if(s=="obsval")
        {
            obj1.obsval(s,qsq);
        }
        else if(s=="obserr")
        {
            obj2.obsval(s,qsq,smwc,npwc);
        }
        else if(s=="obsplot")
        {
            obj3.obsplot(s,smwc,npwc);
        }
        else
            cout << "Error!! Try again." << endl;
    }
    else if(s=="B->ll")
    {
        cout << "Specify an action: obsval or obserr" << endl;
        cin >> s;
        if(s=="obsval")
        {
            obj1.obsval(s,qsq);
        }
        else if(s=="obserr")
        {
            obj2.obsval(s,qsq,smwc,npwc);
        }
        else
            cout << "Error!! Try again." << endl;
    }
    else
        cout << "Error!! Try again." << endl;
    return 0;
}


