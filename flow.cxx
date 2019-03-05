#include <iostream>
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TLegend.h"

using namespace std;

void flow()
{
	/********************
	vytvorenie a nacitanie suborov z dat

	ako Q vektor oznacujem sumu cosinov/sinov pre forward/backward zlozku vsetkych nabitych hadronov (okrem kpi paru) v danom evente
	*********************/
 	TFile *file=new TFile("result.root");
 	TProfile *qVec[4]; //Q vektory (fcia multiplicity)
 	TProfile *refFlow; //referencny tok (fcia multiplicity)
    TProfile *dirFlow[5]; //diferencialny tok (fcia multiplicity a pT)
    TProfile *corrD[2][5];  //cos/sin D pre rozne mult a  pT
    TProfile *qVec2[4]; //Q vektory
    TProfile *refFlow2; //referencny tok 
    TProfile *dirFlow2; //diferencialny tok (fcia PT)
    TProfile *corrD2[2]; //cos/sin D pre rozne pT
    TProfile *ReQ;
    TProfile *ImQ;

    //bkg only
	TProfile *corrDBKG[2][5];
	TProfile *dirFlowBKG[5];
	TProfile *corrD2BKG[2];
	TProfile *dirFlow2BKG;
    
    //TProfile *c2;
    TH1D *c2; 
    TH1D *d2[5];
    TH1D *v2[5];
    
    //cumulanty, directed flow a vysledne v2
    //noC znamena, ze ide o vypocet bez zahrnutia korekcie na uniformnu akceptanciu detektoru

    TH1D *cum;
    TH1D *cum_noC;
    TH1D *d2_all_mult;
    TH1D *v2_all_mult;
    TH1D *d2_all_mult_noC;
    TH1D *v2_all_mult_noC;

	TH1D *d2_all_mult_noC_BKG;
	TH1D *v2_all_mult_noC_BKG;

	TH1D *v2_D0;

    //loading TProfiles
 	TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
  	float multBin[6] = {0,7,12,16,22,100};
  	for(int m = 0; m < 4; m++)
  	{
	  	TString aa = names[m];
	  	qVec[m] = (TProfile *)file->Get(aa.Data());
	  	aa += "_no_mult";
	  	qVec2[m] = (TProfile *)file->Get(aa.Data());
  	}
 	float momBins[4] = {1,2,3,5};
  	TString multBinNames[6] = {"0","7","12","16","22","100"};
  	for(int m = 0; m < 5; m++)
  	{
	  	corrD[0][m] = (TProfile *)file->Get(Form("cosD_%.0f_%.0f", multBin[m], multBin[m+1]));
		corrDBKG[0][m] = (TProfile *)file->Get(Form("cosD_%.0f_%.0f_BKG", multBin[m], multBin[m+1]));
		corrD[1][m] = (TProfile *)file->Get(Form("sinD_%.0f_%.0f", multBin[m], multBin[m+1]));
		corrDBKG[1][m] = (TProfile *)file->Get(Form("sinD_%.0f_%.0f_BKG", multBin[m], multBin[m+1]));
		dirFlow[m] = (TProfile *)file->Get(Form("dirFlow_%.0f_%.0f", multBin[m], multBin[m+1]));
		dirFlowBKG[m] = (TProfile *)file->Get(Form("dirFlow_%.0f_%.0f_BKG", multBin[m], multBin[m+1]));
	}
  	refFlow = (TProfile *)file->Get("refFlow");

  	corrD2[0] = (TProfile *)file->Get("cosD_no_mult");
  	corrD2[1] = (TProfile *)file->Get("sinD_no_mult");
	dirFlow2 = (TProfile *)file->Get("dirFlow_no_mult");
  	refFlow2 = (TProfile *)file->Get("refFlow_no_mult");

	corrD2BKG[0] = (TProfile *)file->Get("cosD_no_mult_BKG");
	corrD2BKG[1] = (TProfile *)file->Get("sinD_no_mult_BKG");
	dirFlow2BKG = (TProfile *)file->Get("dirFlow_no_mult_BKG");


    ReQ = (TProfile *)file->Get("ReQ");
    ImQ = (TProfile *)file->Get("ImQ");


  	//zacina samotny vypocet.... vzorceky z clanku arxiv:1010.0233, resp. proceedings wejcf, v kombinacii s Katkou (niektore veci ako eta gap) a s ana note v2 D0 v Au+Au

  	
	float d,r,c,s,cA,cB,sA,sB,rEr,cEr,sEr,final,errorC,error,cumulant,cumEr; //len pomocne premenne
	//c2 = new TProfile("c2","c_2{2}", 5, multBin);
	c2 = new TH1D("c2", "v2 ref = sqrt(c_2{2})", 5, multBin);
	cum = new TH1D("cum", "v2 ref = sqrt(c_2{2})", 1, 0, 100);
	cum_noC = new TH1D("cum_noC", "v2 ref = sqrt(c_2{2})", 1, 0, 100);
	for(int i = 1; i < 6; i++)
		{
			//computing c_2{2} via reference flow
			r = refFlow->GetBinContent(i);
			//cout << "multiplicity from " << multBinNames[i-1] << " to " << multBinNames[i] << "  value of c_2{2} " << r << endl;
			/************
			 *  dlhe zakomentovane casti
			 *  su korekcie na non-uniform veci
			 *  po konzultacii s KKG to idem prerobit
			 *  lebo som asi nemala dobre tie subory
			 *  (davala som tam Q vector len z jednej casti (forward, backward),
			 *  ale pravdepodobne tam mal byt sucet vsetkych)
			 ************/
			/*
			cA = qVec[0]->GetBinContent(i);
			sA = qVec[2]->GetBinContent(i);
			cB = qVec[1]->GetBinContent(i);
			sB = qVec[3]->GetBinContent(i);
			final = r - cA*cB - sA*sB;
			cout << "final " << final << " bin" << i << endl;
			*/
			final = r;
			c2->SetEntries(refFlow->GetEntries());
			c2->SetBinContent(i, (TMath::Sqrt(final)));
			c2->SetBinError(i, refFlow->GetBinError(i));
		}

	r = refFlow2->GetBinContent(1);
    //c = ReQ->GetBinContent(1);
    //s = ImQ->GetBinContent(1);
    /*
	c = qVec2[0]->GetBinContent(1);
	s = qVec2[2]->GetBinContent(1);
	final = r - c*c - s*s;
	rEr = refFlow2->GetBinError(1);
	cEr = qVec2[0]->GetBinError(1);
	sEr = qVec2[2]->GetBinError(1);
	cout << "r  " << r << "  c  " << c << "  s  " << s << endl;
	cout << "errors..... r " << rEr << "c " << cEr << "s " << sEr << endl;
	errorC = TMath::Sqrt( r*r*rEr*rEr + 4*c*c*cEr*cEr + 4*s*s*sEr*sEr  );
	cout << "error C   " << errorC << endl;
	cout << "value of fnal " << final << endl;
	cum->SetEntries(refFlow2->GetEntries());
	cum->SetBinContent(1, (TMath::Sqrt(final)));
	//cum->SetBinError(1, refFlow2->GetBinError(1));
	error = 0.5*errorC/(TMath::Sqrt(final)); 
	cout << "error   " << error << endl;
	cum->SetBinError(1, error);
	cum->Print("all");
	*/
    error = (0.5*refFlow2->GetBinContent(1))/TMath::Sqrt(r);
	cum_noC->SetEntries(refFlow2->GetEntries());
	cum_noC->SetBinContent(1, (TMath::Sqrt(r)));
	cum_noC->SetBinError(1, error);
	cout << "fancy error " <<  error << endl;
	cout << "error just from RF" << refFlow2->GetBinContent(1);

	printf("error just from ref flow %f \n", refFlow2->GetBinError(1));
	printf(" fancy error %f \n", error);

	cum_noC->Draw();

	for (int i = 0; i < 5; i++)
		{
			//computing v2 via directed flow
			d2[i] = new TH1D(TString::Format("d2_%d", i), "d_2{2}", 3, momBins);
			v2[i] = new TH1D(TString::Format("v2_%d", i), "v_{2};p_{T};v_{2}", 3, momBins);
			d2[i]->SetEntries(dirFlow[i]->GetEntries());
			v2[i]->SetEntries(dirFlow[i]->GetEntries());
			for(int j = 1; j < 4; j++)
			{
				d = dirFlow[i]->GetBinContent(j);
				c = (corrD[0][i]->GetBinContent(j))*(qVec[0]->GetBinContent(i+1));
				s = (corrD[1][i]->GetBinContent(j))*(qVec[2]->GetBinContent(i+1));
		    	//final = d - c - s;
				final = d;
				d2[i]->SetBinContent(j, final);
				d2[i]->SetBinError(j, dirFlow[i]->GetBinError(j));
				v2[i]->SetBinContent(j,(d2[i]->GetBinContent(j))/(c2->GetBinContent(i+1)));
				v2[i]->SetBinError(j,(d2[i]->GetBinError(j))/(c2->GetBinContent(i+1)));
	//			allm2->Fill(momBins[j-1]+0.5,v2[i]->GetBinContent(j),v2[i]->GetEntries());
			}
		}	

	d2_all_mult = new TH1D("d2_all_mult", "d_2{2}", 3, momBins);
	v2_all_mult = new TH1D("v2_all_mult", "v_{2};p_{T};v_{2}", 3, momBins);

	d2_all_mult_noC = new TH1D("d2_all_mult_noC", "d_2{2}", 3, momBins);
	v2_all_mult_noC = new TH1D("v2_all_mult_noC", "v_{2};p_{T};v_{2}", 3, momBins);

	d2_all_mult_noC_BKG = new TH1D("d2_all_mult_noC_BKG", "d_2{2}", 3, momBins);
	v2_all_mult_noC_BKG = new TH1D("v2_all_mult_noC_BKG", "v_{2};p_{T};v_{2}", 3, momBins);

	v2_D0 = new TH1D("v2_D0", "v_{2} of D0 meson;p_{T};v_{2}", 3, momBins);


	cumulant = cum_noC->GetBinContent(1);
	cumEr = cum_noC->GetBinError(1);

	double dirEr,dirF;

	/*
	float Nsig = 1;
	float Nbkg = 2;

	float fSig = Nsig/(Nsig + Nbkg);
	float fBkg = Nbkg/(Nsig + Nbkg);
	*/


	for(int j = 1; j < 4; j++)
			{
				d = dirFlow2->GetBinContent(j);
				c = (corrD2[0]->GetBinContent(j))*(qVec2[0]->GetBinContent(1));
				s = (corrD2[1]->GetBinContent(j))*(qVec2[2]->GetBinContent(1));
				//final = d - c - s;
				final = d;
				d2_all_mult->SetBinContent(j, final);
				d2_all_mult->SetBinError(j, dirFlow2->GetBinError(j));
				d2_all_mult_noC->SetBinContent(j, d);
				d2_all_mult_noC->SetBinError(j, dirFlow2->GetBinError(j));

				d2_all_mult_noC_BKG->SetBinContent(j, dirFlow2BKG->GetBinContent(j));
				d2_all_mult_noC_BKG->SetBinError(j, dirFlow2BKG->GetBinError(j));
				/*
				d2_all_mult->SetBinError(j,  TMath::Sqrt( dirFlow2->GetBinError(j)*dirFlow2->GetBinError(j) + corrD2[0]->GetBinContent(j)*corrD2[0]->GetBinContent(j)*qVec2[0]->GetBinError(1)*qVec2[0]->GetBinError(1) + qVec2[0]->GetBinContent(1)*qVec2[0]->GetBinContent(1)*corrD2[0]->GetBinError(j)*corrD2[0]->GetBinError(j) + corrD2[1]->GetBinContent(j)*corrD2[1]->GetBinContent(j)*qVec2[2]->GetBinError(1)*qVec2[2]->GetBinError(1) + qVec2[2]->GetBinContent(1)*qVec2[2]->GetBinContent(1)*corrD2[1]->GetBinError(j)*corrD2[1]->GetBinError(j) ) );
				printf("fancy error %f \n", d2_all_mult->GetBinError(j));
				printf("normal error from ref flow %f \n", dirFlow2->GetBinError(j));

				dirEr = d2_all_mult->GetBinError(j);
				dirF = d2_all_mult->GetBinContent(j);
				//error = TMath::Sqrt((dirEr*dirEr)/(cumulant*cumulant) +  (dirF*dirF)*(cumEr*cumEr)/(TMath::Power(cumulant,4)));

				v2_all_mult->SetBinContent(j,(d2_all_mult->GetBinContent(j))/(cum->GetBinContent(1)));
				v2_all_mult->SetEntries(dirFlow2->GetEntries());
				error = v2_all_mult->GetBinContent(j)*TMath::Sqrt( (cumEr/cumulant)*(cumEr/cumulant) + (dirEr/dirF)*(dirEr/dirF) );
				//v2_all_mult->SetBinError(j,(d2_all_mult->GetBinError(j))/(cum->GetBinContent(1)));
				v2_all_mult->SetBinError(j,error);
				*/
				v2_all_mult_noC->SetEntries(dirFlow2->GetEntries());
				v2_all_mult_noC->SetBinContent(j,(d2_all_mult_noC->GetBinContent(j))/(cum_noC->GetBinContent(1)));

				v2_all_mult_noC_BKG->SetEntries(dirFlow2BKG->GetEntries());
				v2_all_mult_noC_BKG->SetBinContent(j,(d2_all_mult_noC_BKG->GetBinContent(j))/(cum_noC->GetBinContent(1)));

                dirEr = d2_all_mult_noC->GetBinError(j);
                dirF = d2_all_mult_noC->GetBinContent(j);

				error = TMath::Sqrt( (dirEr*dirEr)/(cumulant*cumulant) + (dirF*dirF*cumEr*cumEr)/(cumulant*cumulant*cumulant*cumulant) );
				v2_all_mult_noC->SetBinError(j,error);

				dirEr = d2_all_mult_noC_BKG->GetBinError(j);
				dirF = d2_all_mult_noC_BKG->GetBinContent(j);

				error = TMath::Sqrt( (dirEr*dirEr)/(cumulant*cumulant) + (dirF*dirF*cumEr*cumEr)/(cumulant*cumulant*cumulant*cumulant) );
				v2_all_mult_noC_BKG->SetBinError(j,error);

				//cout << j << "    value      " << (v2_all_mult_noC->GetBinContent(j) - fBkg*(v2_all_mult_noC_BKG->GetBinContent(j)))/fSig << endl;
				//v2_D0->SetBinContent(j, (v2_all_mult_noC->GetBinContent(j) - fBkg*v2_all_mult_noC_BKG->GetBinContent(j))/fSig);
			}

	//lame pic is coming
	//TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
	//legend->SetHeader("v2 D0 analysis","C");

	//v2_all_mult_noC_BKG->Draw();
	//v2_all_mult_noC->SetLineColor(kRed);
	v2_all_mult_noC->Draw();

	//legend->AddEntry(v2_all_mult_noC_BKG,"Pure bkg (wrong sign)","f");
	//legend->AddEntry(v2_all_mult_noC,"Total (right sign)","f");
	//legend->Draw("same");

	//v2_D0->Draw();

	TFile *output = new TFile("v2_output.root", "recreate");
	c2->Write();
	for (int i = 0; i < 5; i++)
	{
		d2[i]->Write();
		v2[i]->Write();
	}
	//cum->Write();
	cum_noC->Write();
	//v2_all_mult->Write();
	//d2_all_mult->Write();
	v2_all_mult_noC->Write();
	d2_all_mult_noC->Write();

	v2_all_mult_noC_BKG->Write();
	d2_all_mult_noC_BKG->Write();

	//v2_all_mult->Draw();
	//v2_all_mult_noC->SetLineColor(kRed);
	//v2_all_mult_noC->Draw("same");
	
 }

 	