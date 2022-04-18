#include "get_mcplot.h"


void get_obsplot::BtoVll_obsplot(string s, double smwc[], double npwc[]){
    //B->Vll
    BdtoKstrll_obserr eo1; Bstophill_obserr eo2;

    double unv0[72]; for(int i=0;i<72;i++){unv0[i] = 0.0;}
    int iter=100, sdcl(int(iter*15.9/100)); int n(50);
    auto c1 = new TCanvas(); c1->SetGrid();
    auto gr = new TMultiGraph();
    auto leg = new TLegend(0.6,0.7,0.9,0.9); leg->SetHeader("Plot Legends", "C");
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->K*ll #######//
    //B0->K*mumu
    if(s=="AFB(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.AFB(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.AFB(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("AFB(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="FL(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.FL(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.FL(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("FL(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P1(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P1(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P1(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P1(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P2(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P2(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P2(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P2(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P4p(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P4p(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P4p(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P4p(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P5p(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P5p(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P5p(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P5p(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P6p(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P6p(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P6p(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P6p(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P8p(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P8p(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P8p(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P8p(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="Rmue(B0->K*ll)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.diffWidth(q2[i],eo1.mmu(unv0),smwc,npwc,unv0)/eo1.diffWidth(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.diffWidth(q2[i],eo1.mmu(unv),smwc,npwc,unv)/eo1.diffWidth(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("Rmue(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="Rtaumu(B0->K*ll)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.diffWidth(q2[i],eo1.mtau(unv0),smwc,npwc,unv0)/eo1.diffWidth(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.diffWidth(q2[i],eo1.mtau(unv),smwc,npwc,unv)/eo1.diffWidth(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("Rtaumu(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S1(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S1(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S1(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S1(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S2(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S2(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S2(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S2(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S3(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S3(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S3(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S3(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S4(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S4(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S4(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S4(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S5(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S5(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S5(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S5(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S6(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S6(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S6(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S6(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S7(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S7(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S7(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S7(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S8(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S8(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S8(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S8(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S9(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S9(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S9(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S9(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="dBR/dq2(B0->K*mumu)")
    {
        double qmin(4*pow(eo1.mmu(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.diffWidth(q2[i],eo1.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.diffWidth(q2[i],eo1.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("dBR/dq2(B0->K*mumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    //B0->K*tautau
    if(s=="AFB(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.AFB(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.AFB(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("AFB(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="FL(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.FL(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.FL(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("FL(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P1(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P1(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P1(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P1(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P2(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P2(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P2(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P2(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P4p(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P4p(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P4p(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P4p(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P5p(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P5p(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P5p(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P5p(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P6p(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P6p(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P6p(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P6p(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P8p(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P8p(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P8p(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P8p(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S1(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S1(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S1(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S1(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S2(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S2(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S2(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S2(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S3(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S3(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S3(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S3(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S4(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S4(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S4(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S4(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S5(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S5(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S5(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S5(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S6(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S6(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S6(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S6(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S7(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S7(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S7(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S7(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S8(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S8(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S8(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S8(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S9(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S9(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S9(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S9(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="dBR/dq2(B0->K*tautau)")
    {
        double qmin(4*pow(eo1.mtau(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.diffWidth(q2[i],eo1.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.diffWidth(q2[i],eo1.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("dBR/dq2(B0->K*tautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    //B0->K*ee
    else if(s=="AFB(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.AFB(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.AFB(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("AFB(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="FL(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.FL(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.FL(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("FL(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P1(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P1(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P1(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P1(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P2(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P2(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P2(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P2(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P4p(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P4p(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P4p(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P4p(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P5p(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P5p(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P5p(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P5p(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P6p(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P6p(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P6p(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P6p(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="P8p(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.P8p(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.P8p(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("P8p(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S1(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S1(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S1(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S1(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S2(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S2(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S2(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S2(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S3(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S3(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S3(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S3(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S4(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S4(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S4(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S4(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S5(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S5(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S5(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S5(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S6(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S6(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S6(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S6(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S7(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S7(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S7(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S7(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S8(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S8(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S8(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S8(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="S9(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.S9(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.S9(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("S9(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="dBR/dq2(B0->K*ee)")
    {
        double qmin(4*pow(eo1.me(unv0),2)), qmax(pow(eo1.mBd(unv0)-eo1.mKst(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo1.diffWidth(q2[i],eo1.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo1.mnd_default();}
                        data[j]=eo1.diffWidth(q2[i],eo1.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("dBR/dq2(B0->K*ee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    //Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}

void get_obsplot::BtoPll_obsplot(string s, double smwc[], double npwc[]){
    //B->Pll
    BdtoKll_obserr eo3;

    double unv0[72]; for(int i=0;i<72;i++){unv0[i] = 0.0;}
    int iter=100, sdcl(int(iter*15.9/100)); int n(50);
    auto c1 = new TCanvas(); c1->SetGrid();
    auto gr = new TMultiGraph();
    auto leg = new TLegend(0.6,0.7,0.9,0.9); leg->SetHeader("Plot Legends", "C");
    label1:
    cout << "Enter the observable:"; cin >> s;

//####### B0->Kll #######//
    //B0->Kmumu
    if(s=="AFB(B0->Kmumu)")
    {
        double qmin(4*pow(eo3.mmu(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffAFB(q2[i],eo3.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffAFB(q2[i],eo3.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("AFB(B0->Kmumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="FH(B0->Kmumu)")
    {
        double qmin(4*pow(eo3.mmu(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffFH(q2[i],eo3.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffFH(q2[i],eo3.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("FH(B0->Kmumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="Rmue(B0->Kll)")
    {
        double qmin(4*pow(eo3.mmu(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffBrnch(q2[i],eo3.mmu(unv0),smwc,npwc,unv0)/eo3.diffBrnch(q2[i],eo3.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffBrnch(q2[i],eo3.mmu(unv),smwc,npwc,unv)/eo3.diffBrnch(q2[i],eo3.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("Rmue(B0->Kll)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="Rtaumu(B0->Kll)")
    {
        double qmin(4*pow(eo3.mtau(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffBrnch(q2[i],eo3.mtau(unv0),smwc,npwc,unv0)/eo3.diffBrnch(q2[i],eo3.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffBrnch(q2[i],eo3.mtau(unv),smwc,npwc,unv)/eo3.diffBrnch(q2[i],eo3.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("Rtaumu(B0->Kll)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="dBR/dq2(B0->Kmumu)")
    {
        double qmin(4*pow(eo3.mmu(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffBrnch(q2[i],eo3.mmu(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffBrnch(q2[i],eo3.mmu(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("dBR/dq2(B0->Kmumu)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    //B0->Ktautau
    else if(s=="AFB(B0->Ktautau)")
    {
        double qmin(4*pow(eo3.mtau(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffAFB(q2[i],eo3.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffAFB(q2[i],eo3.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("AFB(B0->Ktautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="FH(B0->Ktautau)")
    {
        double qmin(4*pow(eo3.mtau(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffFH(q2[i],eo3.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffFH(q2[i],eo3.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("FH(B0->Ktautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="dBR/dq2(B0->Ktautau)")
    {
        double qmin(4*pow(eo3.mtau(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffBrnch(q2[i],eo3.mtau(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffBrnch(q2[i],eo3.mtau(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("dBR/dq2(B0->Ktautau)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    //B0->Kee
    else if(s=="AFB(B0->Kee)")
    {
        double qmin(4*pow(eo3.me(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffAFB(q2[i],eo3.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffAFB(q2[i],eo3.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("AFB(B0->Kee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="FH(B0->Kee)")
    {
        double qmin(4*pow(eo3.me(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffFH(q2[i],eo3.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffFH(q2[i],eo3.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("FH(B0->Kee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }

    else if(s=="dBR/dq2(B0->Kee)")
    {
        double qmin(4*pow(eo3.me(unv0),2)), qmax(pow(eo3.mBd(unv0)-eo3.mK(unv0),2.0));
        double q2[n], c_obs[n], pe_obs[n], me_obs[n];
        double data[iter]; double unv[72];
        for(int i=0; i<n; i++)
            {
                q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
                c_obs[i] = eo3.diffBrnch(q2[i],eo3.me(unv0),smwc,npwc,unv0);
                for(int j=0;j<iter;j++)
                    {
                        for(int k=0;k<72;k++) {unv[k] = eo3.mnd_default();}
                        data[j]=eo3.diffBrnch(q2[i],eo3.me(unv),smwc,npwc,unv);
                    }
                sort(data,data+iter);
                me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
            }
        auto cgr = new TGraph(n,q2,c_obs); auto mgr = new TGraph(n,q2,me_obs); auto pgr = new TGraph(n,q2,pe_obs);
        cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
        cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
        gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
        gr->GetXaxis()->SetTitle("q^2"); gr->GetYaxis()->SetTitle("dBR/dq2(B0->Kee)"); gr->Draw("AC");
        leg->AddEntry(cgr, "Central Prediction", "l"); leg->AddEntry(pgr, "68.2% Confidence Level", "l"); leg->Draw();
        c1->SaveAs("obsplot.pdf"); }


//Warning!!!
    else
        {cout << "Error!!" << endl;
        cout << "Press 'c' to continue or q to exit." << endl;
        cin >> s;
        if(s=="c") goto label1;
        else if(s=="q") cout << "The program has been terminated." << endl;}
}

