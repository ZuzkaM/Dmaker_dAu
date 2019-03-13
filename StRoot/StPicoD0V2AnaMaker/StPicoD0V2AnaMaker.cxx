#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoD0V2AnaMaker.h"
#include "TComplex.h"
ClassImp(StPicoD0V2AnaMaker)

float multBin[6] = {0,7,12,16,22,100};

// _________________________________________________________
StPicoD0V2AnaMaker::StPicoD0V2AnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName){
    // constructor
}

// _________________________________________________________
StPicoD0V2AnaMaker::~StPicoD0V2AnaMaker() {
    // destructor
}

// _________________________________________________________
int StPicoD0V2AnaMaker::InitHF() {
    DeclareHistograms();
    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");

    TFile *fileW=new TFile("/star/u/zuzana/zuzana/D0v2/Dmaker_dAu/StRoot/StPicoD0V2AnaMaker/weight.root");
    if((!fileW) || (!fileW->IsOpen())) printf("file does not exist");
    weights = (TH1D *)fileW->Get("hadron_phi");


    return kStOK;
}

// _________________________________________________________
void StPicoD0V2AnaMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoD0V2AnaMaker::FinishHF() {
    WriteHistograms();
    return kStOK;
}
// _________________________________________________________
int StPicoD0V2AnaMaker::MakeHF() {
    getHadronCorV2(1);
    return kStOK;
}

// _________________________________________________________
std::vector<int> StPicoD0V2AnaMaker::createCandidates() {

	std::vector<int> tracksofCand;

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++)  {
        StPicoTrack const* pion1 = mPicoDst->track(i);
        if (!mHFCuts -> isGoodPion(pion1)) continue;

        for(unsigned  int j=0;j<mPicoDst->numberOfTracks();j++)  {
            StPicoTrack const* kaon = mPicoDst->track(j);
            if (pion1->id() == kaon->id()) continue;
            if (!mHFCuts -> isGoodKaon(kaon)) continue;
            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), i, j, mPrimVtx, mBField, kTRUE); //the order (pion1, kaon) needs to stay same!
            if (!mHFCuts->isGoodSecondaryVertexPair(pair)) continue;
            if(pair->m() < 1.804 || pair->m() > 1.924 || pair->pt() < 1 || pair->pt() > 5) continue;
			if (!mHFCuts->isGoodSecondaryVertexPairPtBin(pair)) continue;
			if (pair->eta() < 0.075 && pair->eta() > -0.075) continue; //eta gap

			int charge = 0;
            if((kaon->charge() + pion1->charge() != 0) ) charge = 1;

			getCorV2(pair, 1, charge);


			tracksofCand.push_back(pion1->id());
			tracksofCand.push_back(kaon->id());


        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    return tracksofCand;
}

// _________________________________________________________
void StPicoD0V2AnaMaker::DeclareHistograms() {
    TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
    float multBin[6] = {0, 7, 12, 16, 22, 100};
    int nMultBins = sizeof(multBin)/sizeof(multBin[0])-1;

    float momBins[4] = {1,2,3,5};
    int nMomBins = sizeof(momBins)/sizeof(momBins[0])-1;

    for(int m = 0; m < 4; m++) {
        qVec[m] = new TProfile(names[m].Data(),"Q vector", nMultBins, multBin);
        qVec[m]->Sumw2();

        qVec2[m] = new TProfile((names[m]+"_no_mult").Data(),"Q vector", 1, 0, 100);
        qVec2[m]->Sumw2();
    }

    for(int m = 0; m < 5; m++) {
        corrD[0][m] = new TProfile(Form("cosD_%.0f_%.0f", multBin[m], multBin[m+1]), "", nMomBins, momBins);
        corrD[1][m] = new TProfile(Form("sinD_%.0f_%.0f", multBin[m], multBin[m+1]), "", nMomBins, momBins);
        dirFlow[m] = new TProfile(Form("dirFlow_%.0f_%.0f", multBin[m], multBin[m+1]), "", nMomBins, momBins);

        corrDBKG[0][m] = new TProfile(Form("cosD_%.0f_%.0f_BKG", multBin[m], multBin[m+1]), "", nMomBins, momBins);
        corrDBKG[1][m] = new TProfile(Form("sinD_%.0f_%.0f_BKG", multBin[m], multBin[m+1]), "", nMomBins, momBins);
        dirFlowBKG[m] = new TProfile(Form("dirFlow_%.0f_%.0f_BKG", multBin[m], multBin[m+1]), "", nMomBins, momBins);
    }

    corrD2[0] = new TProfile("cosD_no_mult","cos D", nMomBins, momBins);
    corrD2[1] = new TProfile("sinD_no_mult","sin D", nMomBins, momBins);
    dirFlow2 = new TProfile("dirFlow_no_mult","dir flow", nMomBins, momBins);

    corrD2BKG[0] = new TProfile("cosD_no_mult_BKG","cos D _BKG", nMomBins, momBins);
    corrD2BKG[1] = new TProfile("sinD_no_mult_BKG","sin D _BKG", nMomBins, momBins);
    dirFlow2BKG= new TProfile("dirFlow_no_mult_BKG","dir flow _BKG", nMomBins, momBins);

    refFlow = new TProfile("refFlow", "", nMultBins, multBin);
    refFlow_noC = new TProfile("refFlow_noC", "", nMultBins, multBin); //just to compare c22 with & withouth weights
    refFlow2 = new TProfile("refFlow_no_mult", "", 1, 0, 100);

    cosH = new TProfile("ReQ", "Re Q", 1, 0, 100);
    sinH = new TProfile("ImQ", "Im Q", 1, 0, 100);

    hadron_phi = new TH1D("hadron_phi", "Hadron phi", 2000, -5, 5);
    hadron_phi->Sumw2();

    hadron_check = new TH1D("hadron_check", "Hadron phi", 2000, -5, 5);
    hadron_check->Sumw2();

    D_phi = new TH1D("D_phi", "D phi", 2000, -5, 5);
    D_phi->Sumw2();

    hadron_phi_etaP = new TH1D("hadron_phi_etaP", "Hadron phi, eta > 0", 2000, -5, 5);
    D_phi_etaP = new TH1D("D_phi_etaP", "D phi, eta > 0", 2000, -5, 5);

    hadron_phi_etaN = new TH1D("hadron_phi_etaN", "Hadron phi, eta < 0", 2000, -5, 5);
    D_phi_etaN = new TH1D("D_phi_etaN", "D phi, eta < 0", 2000, -5, 5);

    diFlowMass = new TProfile("diFlowMass", "d2 D0 vs. m_inv; m_inv; d_2", 200, 1.75, 1.95);
    diFlowMassBKG = new TProfile("diFlowMassBKG", "d2 D0 vs. m_inv - wrong sign; m_inv; d_2", 200, 1.75, 1.95);

    for(int pT = 0; pT < 3; pT++){
        mass[pT] = new TH1D(Form("Mass_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]), "Mass of K pi pair - US", 2000, 0.4, 2.4);
        massBKG[pT] = new TH1D(Form("Mass_BKG_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]), "Mass of K pi pair - LS", 2000, 0.4, 2.4);
    }

    NtracksFvsCum = new TProfile("NtracksFvsCum", "NtracksFvsCum", 100, 0, 100);
    NtracksBvsCum = new TProfile("NtracksBvsCum", "NtracksBvsCum", 100, 0, 100);

    phiVsEta = new TH2D("phiVsEta", "phi vs. eta of charged hadrons", 2000, -5, 5,80, -2, 2);
}

// _________________________________________________________
void StPicoD0V2AnaMaker::WriteHistograms() {
    for(int m = 0; m < 4; m++) {
        qVec[m]->Write();
        qVec2[m]->Write();
    }
    refFlow->Write();
    refFlow_noC->Write();

    for(int m = 0; m < 5; m++) {
        corrD[0][m]->Write();
        corrD[1][m]->Write();
        dirFlow[m]->Write();

        corrDBKG[0][m]->Write();
        corrDBKG[1][m]->Write();
        dirFlowBKG[m]->Write();
    }

    corrD2[0]->Write();
    corrD2[1]->Write();
    dirFlow2->Write();
    refFlow2->Write();

    corrD2BKG[0]->Write();
    corrD2BKG[1]->Write();
    dirFlow2BKG->Write();

    hadron_phi->Write();
    hadron_check->Write();
    D_phi->Write();

    hadron_phi_etaP->Write();
    D_phi_etaP->Write();

    hadron_phi_etaN->Write();
    D_phi_etaN->Write();

    cosH->Write();
    sinH->Write();

    diFlowMass->Write();
    diFlowMassBKG->Write();

    for(int pT = 0; pT < 3; pT++){
        mass[pT]->Write();
        massBKG[pT]->Write();
    }

    NtracksBvsCum->Write();
    NtracksFvsCum->Write();

    phiVsEta->Write();
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::getHadronCorV2(int idxGap) {

    std::vector<int> tracksToRemove = createCandidates();

    double maxNentries = weights->GetMaximum();
    double weightHadron = 1;

    double QcosF[3] = {0};
    double QcosB[3] = {0};
    double QsinF[3] = {0};
    double QsinB[3] = {0};

    //no weights applied!
    double QcosF_noC[3] = {0};
    double QcosB_noC[3] = {0};
    double QsinF_noC[3] = {0};
    double QsinB_noC[3] = {0};

    double NtracksB = 0;
    double NtracksF = 0;

    TComplex QvectorF[3];
    TComplex QvectorB[3];

    //no weights applied!
    TComplex QvectorF_noC;
    TComplex QvectorB_noC;

	double etaGap[3] = {0,0.15,0.05};
    double mEtaGap = etaGap[idxGap];
    float hadronFill[7] = {0};
    float hadronFill_noC[7] = {0};
    float Qvec[3] = {0};
    const double reweight = 1;//mGRefMultCorrUtil->getWeight();
    // int centrality  = mGRefMultCorrUtil->getCentralityBin9();
    int mult = mPicoEvent->grefMult();


    //loop over all tracks
    for(unsigned int i=0;i<mPicoDst->numberOfTracks();++i) {
        StPicoTrack const* hadron = mPicoDst->track(i);
    
        if(!mHFCuts->isGoodTrack(hadron)) continue;
        if(!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
        float etaHadron = hadron->gMom().PseudoRapidity();
        float phiHadron = hadron->gMom().Phi();

        if(containsId(hadron->id(), tracksToRemove)) continue;

        weightHadron = maxNentries/(weights->GetBinContent( weights->FindBin(phiHadron) ));

        Qvec[0]++;
        Qvec[1] += cos(2 * phiHadron);
        Qvec[2] += sin(2 * phiHadron);

        if(etaHadron<-0.5*mEtaGap) {//backward sample
            hadronFill[0] += weightHadron;
            hadronFill[1] += weightHadron*sin(2 * phiHadron);
            hadronFill[2] += weightHadron*cos(2 * phiHadron);
            NtracksB++;
            hadronFill_noC[0] += 1;
            hadronFill_noC[1] += sin(2 * phiHadron);
            hadronFill_noC[2] += cos(2 * phiHadron);
            for(int ipow = 0; ipow < 3; ipow++)
            {
                QcosB[ipow] += TMath::Power(weightHadron, ipow)*TMath::Cos(2*phiHadron);
                QsinB[ipow] += TMath::Power(weightHadron, ipow)*TMath::Sin(2*phiHadron);
            }
        }

        if(etaHadron>0.5*mEtaGap) {//forward sample
            hadronFill[3] += weightHadron;
            hadronFill[4] += weightHadron*sin(2 * phiHadron);
            hadronFill[5] += weightHadron*cos(2 * phiHadron);
            hadronFill_noC[3] += 1;
            hadronFill_noC[4] += sin(2 * phiHadron);
            hadronFill_noC[5] += cos(2 * phiHadron);
            NtracksF++;
            for(int ipow = 0; ipow < 3; ipow++)
            {
                QcosF[ipow] += TMath::Power(weightHadron, ipow)*TMath::Cos(2*phiHadron);
                QsinF[ipow] += TMath::Power(weightHadron, ipow)*TMath::Sin(2*phiHadron);
            }
        }
        hadron_phi->Fill(phiHadron);
        hadron_check->Fill(phiHadron, weightHadron);
        if(etaHadron>0) hadron_phi_etaP->Fill(phiHadron);
        if(etaHadron<0) hadron_phi_etaN->Fill(phiHadron);
        phiVsEta->Fill(phiHadron, etaHadron);
    }


    //filling real and imaginary part of Q vector
    for(int ipow = 0; ipow < 3; ipow++)
    {
        QvectorB[ipow] = TComplex(QcosB[ipow], QsinB[ipow]);
        QvectorF[ipow] = TComplex(QcosF[ipow], QsinF[ipow]);
    }

    QvectorB_noC = TComplex(hadronFill_noC[2], hadronFill_noC[1]);
    QvectorF_noC = TComplex(hadronFill_noC[5], hadronFill_noC[4]);

    hadronFill[6] = mult;
    hadronFill[7] = reweight;
    //mHadronTuple->Fill(hadronFill);
    if(hadronFill[0]==0 || hadronFill[3]==0)
        return false;

    //Z code: reference flow creation: average sin/cos phi of a hadron in an event.... (no error!)

    if(idxGap==1)  {
        qVec[0]->Fill(mult, hadronFill[2]/hadronFill[0], reweight);
        qVec[1]->Fill(mult, hadronFill[5]/hadronFill[3], reweight);
        qVec[2]->Fill(mult, hadronFill[1]/hadronFill[0], reweight);
        qVec[3]->Fill(mult, hadronFill[4]/hadronFill[3], reweight);

        double c22 = (QvectorB[1]*(TComplex::Conjugate(QvectorF[1]))).Re();
        refFlow->Fill(mult, (c22/(hadronFill[0]*hadronFill[3])), reweight);
        //no weights
        c22 = (QvectorB_noC*(TComplex::Conjugate(QvectorF_noC))).Re();
        refFlow_noC->Fill(mult, (c22/(hadronFill_noC[0]*hadronFill_noC[3])), reweight);
        //no mult
        qVec2[0]->Fill(mult, hadronFill[2]/hadronFill[0], reweight);
        qVec2[1]->Fill(mult, hadronFill[5]/hadronFill[3], reweight);
        qVec2[2]->Fill(mult, hadronFill[1]/hadronFill[0], reweight);
        qVec2[3]->Fill(mult, hadronFill[4]/hadronFill[3], reweight);
        refFlow2->Fill(mult, (c22/(hadronFill[0]*hadronFill[3])), reweight);

        //for non-uniform acceptance
        cosH->Fill(mult, Qvec[1]/Qvec[0]);
        sinH->Fill(mult, Qvec[2]/Qvec[0]);

        //Ntracks vs cumulant
        NtracksBvsCum->Fill(NtracksB, c22);
        NtracksFvsCum->Fill(NtracksF, c22);
    }

    return true;
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::getCorV2(StHFPair *kp,double weight, int charge) {
    int mult = mPicoEvent->grefMult();

//    float multBin[6] = {0,7,12,16,22,100};
    double etaGap[3] = {0,0.15,0.05};
    float momBins[4] = {1,2,3,5};

    double maxNentries = weights->GetMaximum();

    double weightDcan = maxNentries/(weights->GetBinContent( weights->FindBin( kp->phi() ) ));

    TComplex QvecD;
    TComplex QvecHadrons;

    int k=0;
    double corFill[7] = {0};
    corFill[0] = 1 ;
    if(charge == 0) {
        corFill[1] = weightDcan*sin(2 * kp->phi());
        corFill[2] = weightDcan*cos(2 * kp->phi());

        QvecD = TComplex(corFill[2], corFill[1]);

        D_phi->Fill(kp->phi());
        if (kp->eta() > 0) D_phi_etaP->Fill(kp->phi());
        if (kp->eta() < 0) D_phi_etaN->Fill(kp->phi());

        for (unsigned int i = 0; i < mPicoDst->numberOfTracks(); i++) {
            StPicoTrack const *hadron = mPicoDst->track(i);
            if (!mHFCuts->isGoodTrack(hadron)) continue;
            if (!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
            if (i == kp->particle1Idx() || i == kp->particle2Idx()) continue;
            float etaHadron = hadron->gMom().PseudoRapidity();
            float phiHadron = hadron->gMom().Phi();
            double weightHadron = maxNentries / (weights->GetBinContent(weights->FindBin(phiHadron)));

            if (!isEtaGap(kp->eta(), etaGap[k], etaHadron)) continue;
            corFill[3] += weightHadron;
            corFill[4] += weightHadron * sin(2 * phiHadron);
            corFill[5] += weightHadron * cos(2 * phiHadron);
        }
        QvecHadrons = TComplex(corFill[5], corFill[4]);
        double dif22 = (QvecD*(TComplex::Conjugate(QvecHadrons))).Re();

            for (int m = 0; m < 5; m++) {
                if (mult >= multBin[m] && mult < multBin[m + 1]) {
                    corrD[0][m]->Fill(kp->pt(), corFill[2], weight);
                    corrD[1][m]->Fill(kp->pt(), corFill[1], weight);
                    dirFlow[m]->Fill(kp->pt(), dif22/(weightDcan*corFill[3]), weight);
                }
            }
            for (int pT = 0; pT < 3; pT++) {
                if(kp->pt() >= momBins[pT] && kp->pt() < momBins[pT+1]) mass[pT]->Fill( kp->m() );
            }
            corrD2[0]->Fill(kp->pt(), corFill[2], weight);
            corrD2[1]->Fill(kp->pt(), corFill[1], weight);
            dirFlow2->Fill(kp->pt(), dif22/(weightDcan*corFill[3]), weight);
            diFlowMass->Fill(kp->m(), dif22/(weightDcan*corFill[3]), weight);

        }
    else{
        corFill[1] = weightDcan*sin(2 * kp->phi());
        corFill[2] = weightDcan*cos(2 * kp->phi());

        QvecD = TComplex(corFill[2], corFill[1]);
        for (unsigned int i = 0; i < mPicoDst->numberOfTracks(); i++) {
            StPicoTrack const *hadron = mPicoDst->track(i);
            if (!mHFCuts->isGoodTrack(hadron)) continue;
            if (!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
            if (i == kp->particle1Idx() || i == kp->particle2Idx()) continue;
            float etaHadron = hadron->gMom().PseudoRapidity();
            float phiHadron = hadron->gMom().Phi();
            double weightHadron = maxNentries / (weights->GetBinContent(weights->FindBin(phiHadron)));
            if (!isEtaGap(kp->eta(), etaGap[k], etaHadron)) continue;
            corFill[3] += weightHadron;
            corFill[4] += weightHadron * sin(2 * phiHadron);
            corFill[5] += weightHadron * cos(2 * phiHadron);
            }

            QvecHadrons = TComplex(corFill[5], corFill[4]);
            double dif22 = (QvecD*(TComplex::Conjugate(QvecHadrons))).Re();

            for (int m = 0; m < 5; m++) {
                if (mult >= multBin[m] && mult < multBin[m + 1]) {
                    corrDBKG[0][m]->Fill(kp->pt(), corFill[2], weight);
                    corrDBKG[1][m]->Fill(kp->pt(), corFill[1], weight);
                    dirFlowBKG[m]->Fill(kp->pt(), dif22/(weightDcan*corFill[3]), weight);

                }
            }
            for (int pT = 0; pT < 3; pT++) {
                if(kp->pt() >= momBins[pT] && kp->pt() < momBins[pT+1]) massBKG[pT]->Fill( kp->m() );
            }
            corrD2BKG[0]->Fill(kp->pt(), corFill[2], weight);
            corrD2BKG[1]->Fill(kp->pt(), corFill[1], weight);
            dirFlow2BKG->Fill(kp->pt(), dif22/(weightDcan*corFill[3]), weight);
            diFlowMassBKG->Fill(kp->m(), dif22/(weightDcan*corFill[3]), weight);

    }

    return true;
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::isEtaGap(double dEta,double mGap,double hEta) {
    if(mGap == 0) return true;
    //double range =  2. - mGap*2;
    // if(dEta> (1.-2*mGap))
    //   return hEta<(dEta-mGap) && hEta>(dEta-mGap-range);
    // else if(dEta<(-1.+2*mGap))
    //   return hEta>(dEta+mGap) && hEta<(dEta+mGap+range);
    // else
    //   return (hEta>(dEta+mGap) || hEta<(dEta-mGap));
    if(dEta>0)
        return hEta<-0.5*mGap;
    else
        return hEta>0.5*mGap;
}

bool StPicoD0V2AnaMaker::containsId(int id, std::vector<int>& tracksToRemove)
{
	for(unsigned int i = 0; i < tracksToRemove.size(); i++)
	{
		if(id == tracksToRemove[i])
			return true;
	}
}
