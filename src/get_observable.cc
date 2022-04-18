#include "get_observable.h"

void get_obs::BtoVll_obsval(string s, double qsq){
    //B->Vll
    BdtoKstrll_obs o1; Bstophill_obs o2;
    //string input
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->K*ll #######//
    //B0->K*mumu
    if(s=="AFB(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.AFB(qsq,o1.mmu()) << endl;}
    else if(s=="FL(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.FL(qsq,o1.mmu()) << endl;}
    else if(s=="P1(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P1(qsq,o1.mmu()) << endl;}
    else if(s=="P2(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P2(qsq,o1.mmu()) << endl;}
    else if(s=="P4p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P4p(qsq,o1.mmu()) << endl;}
    else if(s=="P5p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P5p(qsq,o1.mmu()) << endl;}
    else if(s=="P6p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P6p(qsq,o1.mmu()) << endl;}
    else if(s=="P8p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P8p(qsq,o1.mmu()) << endl;}
    else if(s=="Rmue(B0->K*ll)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.diffWidth(qsq,o1.mmu())/o1.diffWidth(qsq,o1.me()) << endl;}
    else if(s=="Rtaumu(B0->K*ll)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.diffWidth(qsq,o1.mtau())/o1.diffWidth(qsq,o1.mmu()) << endl;}
    else if(s=="S1(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S1(qsq,o1.mmu()) << endl;}
    else if(s=="S2(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S2(qsq,o1.mmu()) << endl;}
    else if(s=="S3(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S3(qsq,o1.mmu()) << endl;}
    else if(s=="S4(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S4(qsq,o1.mmu()) << endl;}
    else if(s=="S5(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S5(qsq,o1.mmu()) << endl;}
    else if(s=="S6(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S6(qsq,o1.mmu()) << endl;}
    else if(s=="S7(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S7(qsq,o1.mmu()) << endl;}
    else if(s=="S8(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S8(qsq,o1.mmu()) << endl;}
    else if(s=="S9(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S9(qsq,o1.mmu()) << endl;}
    else if(s=="dBR/dq2(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.diffWidth(qsq,o1.mmu()) << endl;}
    //B0->K*tautau
    else if(s=="AFB(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.AFB(qsq,o1.mtau()) << endl;}
    else if(s=="FL(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.FL(qsq,o1.mtau()) << endl;}
    else if(s=="P1(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P1(qsq,o1.mtau()) << endl;}
    else if(s=="P2(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P2(qsq,o1.mtau()) << endl;}
    else if(s=="P4p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P4p(qsq,o1.mtau()) << endl;}
    else if(s=="P5p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P5p(qsq,o1.mtau()) << endl;}
    else if(s=="P6p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P6p(qsq,o1.mtau()) << endl;}
    else if(s=="P8p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P8p(qsq,o1.mtau()) << endl;}
    else if(s=="S1(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S1(qsq,o1.mtau()) << endl;}
    else if(s=="S2(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S2(qsq,o1.mtau()) << endl;}
    else if(s=="S3(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S3(qsq,o1.mtau()) << endl;}
    else if(s=="S4(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S4(qsq,o1.mtau()) << endl;}
    else if(s=="S5(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S5(qsq,o1.mtau()) << endl;}
    else if(s=="S6(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S6(qsq,o1.mtau()) << endl;}
    else if(s=="S7(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S7(qsq,o1.mtau()) << endl;}
    else if(s=="S8(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S8(qsq,o1.mtau()) << endl;}
    else if(s=="S9(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S9(qsq,o1.mtau()) << endl;}
    else if(s=="dBR/dq2(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.diffWidth(qsq,o1.mtau()) << endl;}
    //B0->K*ee
    else if(s=="AFB(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.AFB(qsq,o1.me()) << endl;}
    else if(s=="FL(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.FL(qsq,o1.me()) << endl;}
    else if(s=="P1(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P1(qsq,o1.me()) << endl;}
    else if(s=="P2(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P2(qsq,o1.me()) << endl;}
    else if(s=="P4p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P4p(qsq,o1.me()) << endl;}
    else if(s=="P5p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P5p(qsq,o1.me()) << endl;}
    else if(s=="P6p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P6p(qsq,o1.me()) << endl;}
    else if(s=="P8p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.P8p(qsq,o1.me()) << endl;}
    else if(s=="S1(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S1(qsq,o1.me()) << endl;}
    else if(s=="S2(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S2(qsq,o1.me()) << endl;}
    else if(s=="S3(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S3(qsq,o1.me()) << endl;}
    else if(s=="S4(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S4(qsq,o1.me()) << endl;}
    else if(s=="S5(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S5(qsq,o1.me()) << endl;}
    else if(s=="S6(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S6(qsq,o1.me()) << endl;}
    else if(s=="S7(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S7(qsq,o1.me()) << endl;}
    else if(s=="S8(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S8(qsq,o1.me()) << endl;}
    else if(s=="S9(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.S9(qsq,o1.me()) << endl;}
    else if(s=="dBR/dq2(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o1.diffWidth(qsq,o1.me()) << endl;}

//####### Bs->phill #######//
    //Bs->phimumu
    else if(s=="AFB(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.AFB(qsq,o2.mmu()) << endl;}
    else if(s=="FL(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.FL(qsq,o2.mmu()) << endl;}
    else if(s=="P1(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P1(qsq,o2.mmu()) << endl;}
    else if(s=="P2(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P2(qsq,o2.mmu()) << endl;}
    else if(s=="P4p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P4p(qsq,o2.mmu()) << endl;}
    else if(s=="P5p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P5p(qsq,o2.mmu()) << endl;}
    else if(s=="P6p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P6p(qsq,o2.mmu()) << endl;}
    else if(s=="P8p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P8p(qsq,o2.mmu()) << endl;}
    else if(s=="Rmue(Bs->phill)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.diffWidth(qsq,o2.mmu())/o2.diffWidth(qsq,o2.me()) << endl;}
    else if(s=="Rtaumu(Bs->phill)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.diffWidth(qsq,o2.mtau())/o2.diffWidth(qsq,o2.mmu()) << endl;}
    else if(s=="S1(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S1(qsq,o2.mmu()) << endl;}
    else if(s=="S2(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S2(qsq,o2.mmu()) << endl;}
    else if(s=="S3(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S3(qsq,o2.mmu()) << endl;}
    else if(s=="S4(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S4(qsq,o2.mmu()) << endl;}
    else if(s=="S5(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S5(qsq,o2.mmu()) << endl;}
    else if(s=="S6(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S6(qsq,o2.mmu()) << endl;}
    else if(s=="S7(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S7(qsq,o2.mmu()) << endl;}
    else if(s=="S8(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S8(qsq,o2.mmu()) << endl;}
    else if(s=="S9(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S9(qsq,o2.mmu()) << endl;}
    else if(s=="dBR/dq2(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.diffWidth(qsq,o2.mmu()) << endl;}
    //Bs->phitautau
    else if(s=="AFB(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.AFB(qsq,o2.mtau()) << endl;}
    else if(s=="FL(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.FL(qsq,o2.mtau()) << endl;}
    else if(s=="P1(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P1(qsq,o2.mtau()) << endl;}
    else if(s=="P2(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P2(qsq,o2.mtau()) << endl;}
    else if(s=="P4p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P4p(qsq,o2.mtau()) << endl;}
    else if(s=="P5p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P5p(qsq,o2.mtau()) << endl;}
    else if(s=="P6p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P6p(qsq,o2.mtau()) << endl;}
    else if(s=="P8p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P8p(qsq,o2.mtau()) << endl;}
    else if(s=="S1(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S1(qsq,o2.mtau()) << endl;}
    else if(s=="S2(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S2(qsq,o2.mtau()) << endl;}
    else if(s=="S3(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S3(qsq,o2.mtau()) << endl;}
    else if(s=="S4(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S4(qsq,o2.mtau()) << endl;}
    else if(s=="S5(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S5(qsq,o2.mtau()) << endl;}
    else if(s=="S6(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S6(qsq,o2.mtau()) << endl;}
    else if(s=="S7(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S7(qsq,o2.mtau()) << endl;}
    else if(s=="S8(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S8(qsq,o2.mtau()) << endl;}
    else if(s=="S9(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S9(qsq,o2.mtau()) << endl;}
    else if(s=="dBR/dq2(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.diffWidth(qsq,o2.mtau()) << endl;}
    //Bs->phiee
    else if(s=="AFB(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.AFB(qsq,o2.me()) << endl;}
    else if(s=="FL(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.FL(qsq,o2.me()) << endl;}
    else if(s=="P1(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P1(qsq,o2.me()) << endl;}
    else if(s=="P2(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P2(qsq,o2.me()) << endl;}
    else if(s=="P4p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P4p(qsq,o2.me()) << endl;}
    else if(s=="P5p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P5p(qsq,o2.me()) << endl;}
    else if(s=="P6p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P6p(qsq,o2.me()) << endl;}
    else if(s=="P8p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.P8p(qsq,o2.me()) << endl;}
    else if(s=="S1(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S1(qsq,o2.me()) << endl;}
    else if(s=="S2(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S2(qsq,o2.me()) << endl;}
    else if(s=="S3(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S3(qsq,o2.me()) << endl;}
    else if(s=="S4(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S4(qsq,o2.me()) << endl;}
    else if(s=="S5(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S5(qsq,o2.me()) << endl;}
    else if(s=="S6(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S6(qsq,o2.me()) << endl;}
    else if(s=="S7(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S7(qsq,o2.me()) << endl;}
    else if(s=="S8(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S8(qsq,o2.me()) << endl;}
    else if(s=="S9(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.S9(qsq,o2.me()) << endl;}
    else if(s=="dBR/dq2(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o2.diffWidth(qsq,o2.me()) << endl;}

//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}

void get_obs::BtoPll_obsval(string s, double qsq){
    //B->Pll
    BdtoKll_obs o3;
    //string input
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->Kll(SM) #######//
    //B0->Kmumu
    if(s=="AFB(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffAFB(qsq,o3.mmu()) << endl;}
    else if(s=="FH(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffFH(qsq,o3.mmu()) << endl;}
    else if(s=="Rmue(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mmu())/o3.diffBrnch(qsq,o3.me()) << endl;}
    else if(s=="Rtaumu(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mtau())/o3.diffBrnch(qsq,o3.mmu()) << endl;}
    else if(s=="dBR/dq2(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mmu()) << endl;}
    //B0->Ktautau
    else if(s=="AFB(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffAFB(qsq,o3.mtau()) << endl;}
    else if(s=="FH(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffFH(qsq,o3.mtau()) << endl;}
    else if(s=="dBR/dq2(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.mtau()) << endl;}
    //B0->Kee
    else if(s=="AFB(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffAFB(qsq,o3.me()) << endl;}
    else if(s=="FH(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffFH(qsq,o3.me()) << endl;}
    else if(s=="dBR/dq2(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq; cout << o3.diffBrnch(qsq,o3.me()) << endl;}

//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}

void get_obs::Btoll_obsval(string s){
    //B->ll
    Btoll_obs o4;
    //string input
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->ll(SM) #######//
    //B0->mumu
    if(s=="BR(B0->mumu)")
        cout << o4.BrTimeIntgratd(o4.mBd(),o4.md(),o4.mmu(),o4.mmu()) << endl;
    //B0->tautau
    else if(s=="BR(B0->tautau)")
        cout << o4.BrTimeIntgratd(o4.mBd(),o4.md(),o4.mtau(),o4.mtau()) << endl;
    //B0->ee
    else if(s=="BR(B0->ee)")
        cout << o4.BrTimeIntgratd(o4.mBd(),o4.md(),o4.me(),o4.me()) << endl;
//####### Bs->ll(SM) #######//
    //Bs->mumu
    else if(s=="ADeltaGamma(Bs->mumu)")
        cout << o4.ADeltaGammaf(o4.mBs(),o4.ms(),o4.mmu(),o4.mmu()) << endl;
    else if(s=="BR(Bs->mumu)")
        cout << o4.BrTimeIntgratd(o4.mBs(),o4.ms(),o4.mmu(),o4.mmu()) << endl;
    else if(s=="tau_mumu")
        cout << o4.efftau(o4.mBs(),o4.ms(),o4.mmu(),o4.mmu()) << endl;
    //Bs->tautau
    else if(s=="ADeltaGamma(Bs->tautau)")
        cout << o4.ADeltaGammaf(o4.mBs(),o4.ms(),o4.mtau(),o4.mtau()) << endl;
    else if(s=="BR(Bs->tautau)")
        cout << o4.BrTimeIntgratd(o4.mBs(),o4.ms(),o4.mtau(),o4.mtau()) << endl;
    else if(s=="tau_tautau")
        cout << o4.efftau(o4.mBs(),o4.ms(),o4.mtau(),o4.mtau()) << endl;
    //Bs->ee
    else if(s=="ADeltaGamma(Bs->ee)")
        cout << o4.ADeltaGammaf(o4.mBs(),o4.ms(),o4.me(),o4.me()) << endl;
    else if(s=="BR(Bs->ee)")
        cout << o4.BrTimeIntgratd(o4.mBs(),o4.ms(),o4.me(),o4.me()) << endl;
    else if(s=="tau_ee")
        cout << o4.efftau(o4.mBs(),o4.ms(),o4.me(),o4.me()) << endl;

//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}
