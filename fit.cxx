#include <iostream>
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TLegend.h"
#include "TF1.h"

void fit(){
    TFile *file=new TFile("result_o8.root");
    float momBins[4] = {1,2,3,5};

    Float_t fitRMin = 1.7;
    Float_t fitRMax = 2.;
    const float rotwthmin = 1.84; // peak mean fitting range
    const float rotwthmax = 1.89; //peak mean fitting range

    const int rebin = 5;

    TH1D *mass[3], *massBKG[3];
    TProfile *diFlowMass[3];

    for(int pT = 0; pT < 3; pT++){
        mass[pT] = (TH1D *)file->Get(Form("Mass_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]));
        massBKG[pT] = (TH1D *)file->Get(Form("Mass_BKG_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]));
        diFlowMass[pT] = (TProfile *)file->Get(Form("diFlowMass_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]));
    }

    //colors, rebinning
    for(int pT = 0; pT < 3; pT++){
        mass[pT]->SetMarkerColor(46);
        mass[pT]->SetLineColor(46);

        mass[pT]->Rebin(rebin);
        mass[pT]->GetXaxis()->SetRangeUser(fitRMin,fitRMax);

        massBKG[pT]->Rebin(rebin);
        massBKG[pT]->GetXaxis()->SetRangeUser(fitRMin,fitRMax);

    }

    TF1 *funUS = new TF1("funUS","pol1(0)+gaus(2)", fitRMin,fitRMax);
    funUS->SetParameters(1.,1.,1.,1.84,0.01);
    funUS->SetLineColor(2);
    funUS->SetLineStyle(1);
    funUS->SetParName(2,"height");
    funUS->SetParName(3,"mean");
    funUS->SetParName(4,"sigma");
    funUS->SetParLimits(3,rotwthmin,rotwthmax);
    funUS->SetParLimits(4,0.000001,0.1);
    funUS->SetLineColor(9);

    TF1 *funLS = new TF1("funLS","pol1(0)", fitRMin,fitRMax);
    funLS->SetParameters(1.,1.);
    funLS->SetLineColor(2);
    funLS->SetLineStyle(1);
    funLS->SetLineColor(9);



    TCanvas *c, *c2;


    TF1 *s, *s2, *b, *sb, *Nsig, *Nbkg, *spolu;
    double mean, sigma, norm, lin, cons;

    for(int pT = 0; pT < 3; pT++)
    {
        cout << "pt bin " << pT << endl;
        mass[pT]->Fit(funUS, "R");
        c = new TCanvas("c", "");
        mass[pT]->Draw();
        c->SaveAs(Form("Figures/Mass_pT_%.0f_%0.f.png", momBins[pT], momBins[pT+1]));
        massBKG[pT]->Fit(funLS, "L");
        c2 = new TCanvas("c2", "");
        massBKG[pT]->Draw();
        c2->SaveAs(Form("Figures/Mass_BKG_pT_%.0f_%0.f.png", momBins[pT], momBins[pT+1]));

        cons = funUS->GetParameter(0);
        lin = funUS->GetParameter(1);
        norm = funUS->GetParameter(2);
        mean = funUS->GetParameter(3);
        sigma = funUS->GetParameter(4);
        }

    //s = new TF1("s", "31*exp(-0.5*((x-1.86)/0.023)^2)", fitRMin,fitRMax);
    s = new TF1("s", "[0]*exp(-0.5*((x-[1])/[2])^2)", fitRMin,fitRMax);
    //s->SetParameters(norm, mean, sigma);
    s->FixParameter(0, norm);
    s->FixParameter(1, mean);
    s->FixParameter(2, sigma);


    cout << " MEAN    " << mean << endl;

    b = new TF1("b", "33 - 14*x", fitRMin,fitRMax);

    sb = new TF1("sb", "s+b");

    Nsig = new TF1("Nsig", "s/(s+b)");

    Nbkg = new TF1("Nbkg", "b/(s+b)");

    spolu = new TF1("spolu", "[0]*Nsig + ([1]+[2]*x)*Nbkg", fitRMin,fitRMax);
    spolu->SetParameters(0.02, 1, 1);

    new TCanvas;
    diFlowMass[2]->Fit(spolu, "R");
    diFlowMass[2]->Draw();


    /*
    TF1 *dflow = new TF1("dflow", "funUS*[0] + ([1]+[2]*x)*(1-funUS)", 1.75, 1.95);
    new TCanvas;
    diFlowMass[2]->Fit(dflow, "R");
    diFlowMass[2]->Draw();
    */
}