#include "get_mc_observable.h"

void get_obserr::BtoVll_obsval(string s, double qsq, double smwc[], double npwc[]){
    //B->Vll
    BdtoKstrll_obserr eo1; Bstophill_obserr eo2;

    int iter=100;
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->K*ll #######//
    //B0->K*mumu
    if(s=="AFB(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.AFB(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="FL(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.FL(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P1(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P1(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P2(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P2(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P4p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P4p(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P5p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P5p(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P6p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P6p(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P8p(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P8p(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="Rmue(B0->K*ll)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.diffWidth(qsq,eo1.mmu(unv),smwc,npwc,unv)/eo1.diffWidth(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="Rtaumu(B0->K*ll)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.diffWidth(qsq,eo1.mtau(unv),smwc,npwc,unv)/eo1.diffWidth(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S1(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S1(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S2(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S2(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S3(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S3(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S4(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S4(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S5(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S5(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S6(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S6(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S7(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S7(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S8(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S8(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S9(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S9(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="dBR/dq2(B0->K*mumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.diffWidth(qsq,eo1.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}

    //B0->K*tautau
    else if(s=="AFB(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.AFB(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="FL(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.FL(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P1(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P1(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P2(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P2(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P4p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P4p(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P5p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P5p(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P6p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P6p(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P8p(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P8p(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S1(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S1(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S2(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S2(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S3(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S3(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S4(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S4(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S5(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S5(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S6(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S6(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S7(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S7(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S8(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S8(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S9(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S9(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="dBR/dq2(B0->K*tautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.diffWidth(qsq,eo1.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}

    //B0->K*ee
    else if(s=="AFB(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.AFB(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="FL(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.FL(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P1(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P1(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P2(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P2(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P4p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P4p(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P5p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P5p(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P6p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P6p(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="P8p(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.P8p(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S1(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S1(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S2(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S2(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S3(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S3(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S4(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S4(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S5(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S5(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S6(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S6(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S7(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S7(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S8(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S8(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="S9(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.S9(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}
    else if(s=="dBR/dq2(B0->K*ee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo1.mnd_default();}
                data[i]=eo1.diffWidth(qsq,eo1.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo1.mean_model(data,iter) << endl;
        cout << "SD: " << eo1.sd_model(data,iter) << endl;}

//####### Bs->phill #######//
    //Bs->phimumu
    else if(s=="AFB(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.AFB(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="FL(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.FL(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P1(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P1(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P2(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P2(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P4p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P4p(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P5p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P5p(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P6p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P6p(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P8p(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P8p(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="Rmue(Bs->phill)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.diffWidth(qsq,eo2.mmu(unv),smwc,npwc,unv)/eo2.diffWidth(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="Rtaumu(Bs->phill)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.diffWidth(qsq,eo2.mtau(unv),smwc,npwc,unv)/eo2.diffWidth(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S1(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S1(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S2(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S2(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S3(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S3(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S4(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S4(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S5(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S5(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S6(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S6(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S7(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S7(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S8(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S8(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S9(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S9(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="dBR/dq2(Bs->phimumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.diffWidth(qsq,eo2.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}

    //Bs->phitautau
    else if(s=="AFB(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.AFB(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="FL(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.FL(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P1(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P1(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P2(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P2(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P4p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P4p(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P5p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P5p(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P6p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P6p(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P8p(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P8p(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S1(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S1(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S2(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S2(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S3(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S3(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S4(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S4(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S5(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S5(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S6(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S6(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S7(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S7(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S8(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S8(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S9(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S9(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="dBR/dq2(Bs->phitautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.diffWidth(qsq,eo2.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}

    //Bs->phiee
    else if(s=="AFB(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.AFB(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="FL(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.FL(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P1(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P1(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P2(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P2(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P4p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P4p(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P5p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P5p(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P6p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P6p(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="P8p(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.P8p(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S1(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S1(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S2(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S2(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S3(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S3(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S4(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S4(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S5(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S5(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S6(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S6(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S7(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S7(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S8(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S8(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="S9(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.S9(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}
    else if(s=="dBR/dq2(Bs->phiee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo2.mnd_default();}
                data[i]=eo2.diffWidth(qsq,eo2.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo2.mean_model(data,iter) << endl;
        cout << "SD: " << eo2.sd_model(data,iter) << endl;}

//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}

void get_obserr::BtoPll_obsval(string s, double qsq, double smwc[], double npwc[]){
    //B->Pll
    BdtoKll_obserr eo3;

    int iter=100;
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->Kll #######//
    //B0->Kmumu
    if(s=="AFB(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffAFB(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="FH(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffFH(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="Rmue(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mmu(unv),smwc,npwc,unv)/eo3.diffBrnch(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="Rtaumu(B0->Kll)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mtau(unv),smwc,npwc,unv)/eo3.diffBrnch(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="dBR/dq2(B0->Kmumu)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    //B0->Ktautau
    else if(s=="AFB(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffAFB(qsq,eo3.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="FH(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffFH(qsq,eo3.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="dBR/dq2(B0->Ktautau)")
        {cout << "Enter a value of qsq:"; cin >> qsq; double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    //B0->Kee
    else if(s=="AFB(B0->Kee)")
        {cout << "Enter a value of qsq:";
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffAFB(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="FH(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffFH(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}

    else if(s=="dBR/dq2(B0->Kee)")
        {cout << "Enter a value of qsq:"; cin >> qsq;
        double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo3.mnd_default();}
                data[i]=eo3.diffBrnch(qsq,eo3.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo3.mean_model(data,iter) << endl;
        cout << "SD: " << eo3.sd_model(data,iter) << endl;}


//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}





void get_obserr::Btoll_obsval(string s, double smwc[], double npwc[]){
    //B->ll
    Btoll_obserr eo4;

    int iter=100;
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->ll #######//
    //B0->mumu
    if(s=="BR(B0->mumu)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBd(unv),eo4.md(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //B0->tautau
    else if(s=="BR(B0->tautau)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBd(unv),eo4.md(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //B0->ee
    else if(s=="BR(B0->ee)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBd(unv),eo4.md(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

//####### Bs->ll #######//
    //Bs->mumu
    else if(s=="ADeltaGamma(Bs->mumu)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.ADeltaGammaf(eo4.mBs(unv),eo4.ms(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="BR(Bs->mumu)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBs(unv),eo4.ms(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="tau_mumu")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.efftau(eo4.mBs(unv),eo4.ms(unv),eo4.mmu(unv),eo4.mmu(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //Bs->tautau
    else if(s=="ADeltaGamma(Bs->tautau)")
         {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.ADeltaGammaf(eo4.mBs(unv),eo4.ms(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="BR(Bs->tautau)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBs(unv),eo4.ms(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="tau_tautau")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.efftau(eo4.mBs(unv),eo4.ms(unv),eo4.mtau(unv),eo4.mtau(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    //Bs->ee
    else if(s=="ADeltaGamma(Bs->ee)")
         {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.ADeltaGammaf(eo4.mBs(unv),eo4.ms(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="BR(Bs->ee)")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.BrTimeIntgratd(eo4.mBs(unv),eo4.ms(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

    else if(s=="tau_ee")
        {double data[iter]; double unv[72];
        for(int i=0;i<iter;i++){
                for(int j=0;j<72;j++){unv[j] = eo4.mnd_default();}
                data[i]=eo4.efftau(eo4.mBs(unv),eo4.ms(unv),eo4.me(unv),eo4.me(unv),smwc,npwc,unv);}
        cout << "Mean: " << eo4.mean_model(data,iter) << endl;
        cout << "SD: " << eo4.sd_model(data,iter) << endl;}

//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}
