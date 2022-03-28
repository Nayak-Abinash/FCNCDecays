//ProjectFiles
#include "observables.h"
#include "observables_mc.h"

double zero_unv[] = {/*A0_a0*/0.0, /*A0_a1*/0.0, /*A0_a2*/0.0, /*A1_a0*/0.0, /*A1_a1*/0.0, /*A1_a2*/0.0, /*A12_a1*/0.0, /*A12_a2*/0.0, /*V_a0*/0.0, /*V_a1*/0.0,
                /*V_a2*/0.0, /*T1_a0*/0.0, /*T1_a1*/0.0, /*T1_a2*/0.0, /*T2_a1*/0.0, /*T2_a2*/0.0, /*T23_a0*/0.0, /*T23_a1*/0.0, /*T23_a2*/0.0, /*CKM_A*/0.0,
                /*CKM_lambda*/0.0, /*CKM_rhobr*/0.0, /*CKM_etabr*/0.0, /*C1*/0.0, /*C2*/0.0, /*C3*/0.0, /*C4*/0.0, /*C5*/0.0, /*C6*/0.0, /*C7effRe*/0.0,
                /*C7effIm*/0.0, /*C8eff*/0.0, /*C9Re*/0.0, /*C9Im*/0.0, /*C10Re*/0.0, /*C10Im*/0.0, /*C7RHRe*/0.0, /*C7RHIm*/0.0, /*C9RHRe*/0.0, /*C9RHIm*/0.0,
                /*C10RHRe*/0.0, /*C10RHIm*/0.0, /*CSRe*/0.0, /*CSIm*/0.0, /*CSRHRe*/0.0, /*CSRHIm*/0.0, /*CPRe*/0.0, /*CPIm*/0.0, /*CPRHRe*/0.0, /*CPRHIm*/0.0,
                /*CTRe*/0.0, /*CTIm*/0.0, /*CT5Re*/0.0, /*CT5Im*/0.0, /*md*/0.0, /*mc*/0.0, /*ms*/0.0, /*mb*/0.0, /*me*/0.0, /*mmu*/0.0, /*mtau*/0.0, /*mBd*/0.0,
                /*mBs*/0.0, /*mK*/0.0, /*mKst*/0.0, /*mphi*/0.0, /*tauBd*/0.0, /*tauBs*/0.0, /*DGamma_dbar*/0.0, /*DGamma_sbar*/0.0, /*fBd*/0.0, /*fBs*/0.0};

double smwc[] = {/*C_1*/-0.257, /*C_2*/1.009, /*C_3*/-0.005, /*C_4*/-0.078, /*C_5*/0.000, /*C_6*/0.001, /*C_7effRe*/-0.304, /*C_7effIm*/0.0,
                                    /*C_8eff*/-0.167, /*C_9Re*/4.211, /*C_9Im*/0.0, /*C_10Re*/-4.103, /*C_10Im*/0.0};

double npwc[] = {/*C_7RHRe*/-0.006, /*C_7RHIm*/0.0, /*C_9RHRe*/0.0, /*C_9RHIm*/0.0, /*C_10RHRe*/0.0, /*C_10RHIm*/0.0, /*C_SRe*/0.0,
                                    /*C_SIm*/0.0, /*C_SRHRe*/0.0, /*C_SRHIm*/0.0, /*C_PRe*/0.0, /*C_PIm*/0.0, /*C_PRHRe*/0.0, /*C_PRHIm*/0.0,
                                    /*C_TRe*/0.0, /*C_TIm*/0.0, /*C_T5Re*/0.0, /*C_T5Im*/0.0};

int main()
{
    BdtoKstrll_obs o1; BdtoKstrll_obserr eo1;
    double qmin = 4*pow(o1.mmu(),2);
    double qmax = pow(o1.mBd()-o1.mKst(),2.0);
///////////Error///////////////
////68.2%ConfidenceLevel (point counting from boundary values)////////////////
    auto c1 = new TCanvas();
//Boundary counting
    int n(100); double q2[n], c_obs[n], pe_obs[n], me_obs[n];
    int iter(1000), sdcl(int(iter*15.9/100)); double data[iter];
    for(int i=0; i<n; i++)
    {
        q2[i] = qmin + (i+1)*(qmax-qmin)/(n+1);
        c_obs[i] = o1.AFB(q2[i],o1.mmu());
        for(int j=0; j<iter; j++)
        {
            double unv[] = {/*A0_a0*/eo1.mnd_default(), /*A0_a1*/eo1.mnd_default(), /*A0_a2*/eo1.mnd_default(), /*A1_a0*/eo1.mnd_default(),
                            /*A1_a1*/eo1.mnd_default(), /*A1_a2*/eo1.mnd_default(), /*A12_a1*/eo1.mnd_default(), /*A12_a2*/eo1.mnd_default(),
                            /*V_a0*/eo1.mnd_default(), /*V_a1*/eo1.mnd_default(), /*V_a2*/eo1.mnd_default(), /*T1_a0*/eo1.mnd_default(),
                            /*T1_a1*/eo1.mnd_default(), /*T1_a2*/eo1.mnd_default(), /*T2_a1*/eo1.mnd_default(), /*T2_a2*/eo1.mnd_default(),
                            /*T23_a0*/eo1.mnd_default(), /*T23_a1*/eo1.mnd_default(), /*T23_a2*/eo1.mnd_default(), /*CKM_A*/eo1.mnd_default(),
                            /*CKM_lambda*/eo1.mnd_default(), /*CKM_rhobr*/eo1.mnd_default(), /*CKM_etabr*/eo1.mnd_default(), /*C1*/eo1.mnd_default(),
                            /*C2*/eo1.mnd_default(), /*C3*/eo1.mnd_default(), /*C4*/eo1.mnd_default(), /*C5*/eo1.mnd_default(),
                            /*C6*/eo1.mnd_default(), /*C7effRe*/eo1.mnd_default(), /*C7effIm*/eo1.mnd_default(), /*C8eff*/eo1.mnd_default(),
                            /*C9Re*/eo1.mnd_default(), /*C9Im*/eo1.mnd_default(), /*C10Re*/eo1.mnd_default(), /*C10Im*/eo1.mnd_default(),
                            /*C7RHRe*/eo1.mnd_default(), /*C7RHIm*/eo1.mnd_default(), /*C9RHRe*/eo1.mnd_default(), /*C9RHIm*/eo1.mnd_default(),
                            /*C10RHRe*/eo1.mnd_default(), /*C10RHIm*/eo1.mnd_default(), /*CSRe*/eo1.mnd_default(), /*CSIm*/eo1.mnd_default(),
                            /*CSRHRe*/eo1.mnd_default(), /*CSRHIm*/eo1.mnd_default(), /*CPRe*/eo1.mnd_default(), /*CPIm*/eo1.mnd_default(),
                            /*CPRHRe*/eo1.mnd_default(), /*CPRHIm*/eo1.mnd_default(), /*CTRe*/eo1.mnd_default(), /*CTIm*/eo1.mnd_default(),
                            /*CT5Re*/eo1.mnd_default(), /*CT5Im*/eo1.mnd_default(), /*md*/eo1.mnd_default(), /*mc*/eo1.mnd_default(),
                            /*ms*/eo1.mnd_default(), /*mb*/eo1.mnd_default(), /*me*/eo1.mnd_default(), /*mmu*/eo1.mnd_default(),
                            /*mtau*/eo1.mnd_default(), /*mBd*/eo1.mnd_default(), /*mBs*/eo1.mnd_default(), /*mK*/eo1.mnd_default(),
                            /*mKst*/eo1.mnd_default(), /*mphi*/eo1.mnd_default(), /*tauBd*/eo1.mnd_default(), /*tauBs*/eo1.mnd_default(),
                            /*DGamma_dbar*/eo1.mnd_default(), /*DGamma_sbar*/eo1.mnd_default(), /*fBd*/eo1.mnd_default(), /*fBs*/eo1.mnd_default()};
            data[j]= eo1.AFB(q2[i],eo1.mmu(unv),smwc,npwc,unv);
        }
        sort(data,data+iter);
        me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
        /*if(c_obs[i] >= 0)
        {
            me_obs[i] = data[sdcl]; pe_obs[i] = data[iter-sdcl-1];
        }
        else
        {
            me_obs[i] = data[iter-sdcl-1]; pe_obs[i] = data[sdcl];
        }*/
        //cout << me_obs[i] << ", " << c_obs[i] << ", " << pe_obs[i] << endl;
    }
    c1->SetGrid();
    auto gr = new TMultiGraph();
    auto cgr = new TGraph(n,q2,c_obs);
    auto pgr = new TGraph(n,q2,pe_obs);
    auto mgr = new TGraph(n,q2,me_obs);
    cgr->SetLineColor(2); pgr->SetLineColor(4); mgr->SetLineColor(4);
    cgr->SetLineWidth(2); pgr->SetLineWidth(2); mgr->SetLineWidth(2);
    //cgr->SetMarkerColor(2); pgr->SetMarkerColor(4); mgr->SetMarkerColor(4);
    //cgr->SetMarkerStyle(21); pgr->SetMarkerStyle(21); mgr->SetMarkerStyle(21);
    gr->Add(cgr); gr->Add(pgr); gr->Add(mgr);
    gr->GetXaxis()->SetTitle("q^2");
    gr->GetYaxis()->SetTitle("AFB(q^2)");
    gr->Draw("AC");
    //TCanvas::Update() draws the frame, after which one can change it
    c1->Update();
    c1->GetFrame()->SetBorderSize(12);
    c1->Modified();
    auto leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetHeader("Plot Legends", "C");
    leg->AddEntry(cgr, "Central Prediction", "l");
    leg->AddEntry(pgr, "68.2% Confidence Level", "l");
    leg->Draw();
    c1->SaveAs("obsplot.pdf");

    return 0;
}

