#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoD0V2AnaMaker.h"
#include "TComplex.h"



ClassImp(StPicoD0V2AnaMaker)

const double etaGap[3] = {0,0.15,1.0};
const float multBin[6] = {0,7,12,16,22,100};
//set number of pT bins and binning
const int nptBins = 3;
const float momBins[nptBins+1] = {1,2,3,5};
//set means and sigmas from fit of the invariant mass -- !! must be in accordance with previous lines!! #of bins
float const meanFit[nptBins] = {1.866, 1.863, 1.864};
float const sigmaFit[nptBins] = {0.0137, 0.0131, 0.0234};
// SET harmonics (developed for elliptic flow v_2, therefore harmonics = 2)
const int harmonics = 2;


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

    //weight file - reweighting phi
    TFile *fileW=new TFile("/star/u/zuzana/zuzana/D0v2/Dmaker_dAu/StRoot/StPicoD0V2AnaMaker/weight.root");
    if((!fileW) || (!fileW->IsOpen())) printf("file does not exist");
    weights = (TH1D *)fileW->Get("hadron_phi");

    //TMVA init
    TMVA::Tools::Instance();

    TString dir    = "/star/u/zuzana/zuzana/D0v2/Dmaker_dAu/StRoot/weights/";
    TString prefix = "TMVAClassification";
    TString ptbin[nptBins] = {"12", "23", "35"};


    for (int pT = 0; pT < nptBins; pT++) {
        reader[pT] = new TMVA::Reader( "!Color:!Silent" );
        reader[pT]->AddVariable("k_dca", &k_dca[pT] );
        reader[pT]->AddVariable("pi1_dca", &pi1_dca[pT] );
        reader[pT]->AddVariable("dcaDaughters", &dcaDaughters[pT] );
        reader[pT]->AddVariable("cosTheta", &cosTheta[pT]  );
        reader[pT]->AddVariable("D_decayL", &D_decayL[pT] );
        reader[pT]->AddVariable("dcaD0ToPv", &dcaD0ToPv[pT] );

        TString methodName = "BDT method";
        TString weightfile = dir + prefix + TString("_BDT.weights.pt") + ptbin[pT] + TString(".xml");
        reader[pT]->BookMVA( methodName, weightfile );
    }


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
    getHadronCorV2(2);
    return kStOK;
}

// _________________________________________________________
std::vector<int> StPicoD0V2AnaMaker::createCandidates() {

	std::vector<int> tracksofCand; //to exclude daughter particles

    //tmva input cuts
    float const dcaV0ToPvCons = 0.05;
    float const decayLengthCons = 0.0005; //0.0005
    float const cosThetaCons = 0.5;
    float const dcaDaughtersCons = 0.02;
    float const kDca = 0.002;
    float const pDca = 0.002;
    float const minPt = 0.15;
    //from Lukas's ana
    //float const bdtCuts[nptBins] = {0.365, 0.299, 0.288};
    float const bdtCuts[nptBins] = {0.21, 0.2, 0.22}; //looser cuts

    //loop - all particles
    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++)  {
        StPicoTrack const* pion1 = mPicoDst->track(i);
        if (!mHFCuts -> isGoodPion(pion1)) continue;

        for(unsigned  int j=0;j<mPicoDst->numberOfTracks();j++)  {
            StPicoTrack const* kaon = mPicoDst->track(j);
            if (pion1->id() == kaon->id()) continue;
            if (!mHFCuts -> isGoodKaon(kaon)) continue;

            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), i, j, mPrimVtx, mBField, kTRUE); //the order (pion1, kaon) needs to stay same!

	    //analysis just in this interval
            if(pair->pt() < 1 || pair->pt() > 5) continue;
            if(pair->eta() < -1 || pair->eta() > 1) continue;

	    //unlike sign pairs < 2
            float flag = -99.;
            if(kaon->charge()<0 && pion1->charge()>0 ) flag=0.; // -+
            if(kaon->charge()>0 && pion1->charge()<0 ) flag=1.; // +-

            if(kaon->charge()<0 && pion1->charge()<0) flag=4.; // --
            if(kaon->charge()>0 && pion1->charge()>0) flag=5.; // ++

	    //find the correct pT bin
            int pTbin = 0;
            for (int pT = 0; pT < nptBins; pT++) {
                if(pair->pt() >= momBins[pT] && pair->pt() < momBins[pT+1]) pTbin = pT;
            }

            //assigning TMVA variables to those used in picos
            k_pt[pTbin] = kaon->gPt();
            pi1_pt[pTbin] = pion1->gPt();
            k_dca[pTbin] = pair->particle2Dca();
            pi1_dca[pTbin] = pair->particle1Dca();
            D_decayL[pTbin] = pair->decayLength();
            cosTheta[pTbin] = cos(pair->pointingAngle());
            dcaD0ToPv[pTbin] = pair->DcaToPrimaryVertex();
            dcaDaughters[pTbin] = pair->dcaDaughters();

	    //tmva input cuts - similar as during the training phase
            if  (k_pt[pTbin]>minPt && pi1_pt[pTbin]>minPt && D_decayL[pTbin]>decayLengthCons && D_decayL[pTbin]<0.2 &&
                 dcaDaughters[pTbin]<dcaDaughtersCons && k_dca[pTbin]>kDca && k_dca[pTbin]<0.2 &&
                 pi1_dca[pTbin]>pDca && pi1_dca[pTbin]<0.2 && dcaD0ToPv[pTbin] < dcaV0ToPvCons && cosTheta[pTbin] > cosThetaCons) {

                 //evaluate BDT, continue just pairs that have passed BDT cut
                 float valueMVA = reader[pTbin]->EvaluateMVA("BDT method");
                 if(valueMVA < bdtCuts[pTbin]) continue;

                 //filling plots of invariant mass for unlike and like sign pairs
                 if(flag < 2) mass[pTbin]->Fill( pair->m() );
                 else massBKG[pTbin]->Fill( pair->m() );

                 //v2 analysis
                 if(pair->m() < 1.2 || pair->m() > 2.4) continue; //cut for d22 vs. mInv
                 if(pair->eta() < 0.075 && pair->eta() > -0.075) continue; //eta gap

                 getCorV2(pair, 1, flag);

                 //excluding daughter tracks
                 if(pair->m() < meanFit[pTbin] - 3*sigmaFit[pTbin] || pair->m() < meanFit[pTbin] + 3*sigmaFit[pTbin]) continue;
                 tracksofCand.push_back(pion1->id());
                 tracksofCand.push_back(kaon->id());
            }
        }
    }

    return tracksofCand;
}

// _________________________________________________________
void StPicoD0V2AnaMaker::DeclareHistograms() {
    TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
    int nMultBins = sizeof(multBin)/sizeof(multBin[0])-1;
    int nMomBins = sizeof(momBins)/sizeof(momBins[0])-1;

    dirFlow2 = new TProfile("dirFlow_no_mult","dir flow", nMomBins, momBins);
    dirFlow2BKG= new TProfile("dirFlow_no_mult_BKG","dir flow _BKG", nMomBins, momBins);

    refFlow = new TProfile("refFlow", "", nMultBins, multBin);
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



    for(int pT = 0; pT < nptBins; pT++){
        mass[pT] = new TH1D(Form("Mass_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]), "Mass of K pi pair - US", 2000, 0.4, 2.4);
        mass[pT]->Sumw2();
        massBKG[pT] = new TH1D(Form("Mass_BKG_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]), "Mass of K pi pair - LS", 2000, 0.4, 2.4);
        massBKG[pT]->Sumw2();

        diFlowMass[pT] = new TProfile(Form("diFlowMass_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]), "d2 D0 vs. m_inv; m_inv; d_2", 2000, 0.4, 2.4);
        diFlowMassBKG[pT] = new TProfile(Form("diFlowMassBKG_pT_%.0f_%0.f", momBins[pT], momBins[pT+1]), "d2 D0 vs. m_inv - wrong sign; m_inv; d_2", 2000, 0.4, 2.4);
    }

    NtracksFvsCum = new TProfile("NtracksFvsCum", "NtracksFvsCum", 100, 0, 100);
    NtracksBvsCum = new TProfile("NtracksBvsCum", "NtracksBvsCum", 100, 0, 100);

    NtracksBvsF = new TH2D("NtracksBvsF", "", 100, 0, 100, 100, 0, 100);
    NtracksVsMult = new TH2D("NtracksVsMult", "", 100, 0, 100, 100, 0, 100);

    phiVsEta = new TH2D("phiVsEta", "phi vs. eta of charged hadrons", 1000, -5, 5,40, -2, 2);
    phiVsEtaDcand = new TH2D("phiVsEtaDcand", "phi vs. eta of D candidates", 1000, -5, 5,40, -2, 2);

    //PID capability plots
    TOF = new TH2D("TOF", "", 1000, 0, 3.5, 1000, 0, 2.5);
    TPC = new TH2D("TPC", "", 1000, 0, 3.5, 1000, 0, 6);

    //pT of daughter particles plots
    kPT = new TH1D("kPT", "", 1000, 0, 10);
    piPT = new TH1D("piPT", "", 100, 0, 10);

}

// _________________________________________________________
void StPicoD0V2AnaMaker::WriteHistograms() {

    const int nptBins = 3;

    refFlow->Write();
    dirFlow2->Write();
    refFlow2->Write();
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



    for(int pT = 0; pT < nptBins; pT++){
        mass[pT]->Write();
        massBKG[pT]->Write();
        diFlowMass[pT]->Write();
        diFlowMassBKG[pT]->Write();
    }

    NtracksBvsCum->Write();
    NtracksFvsCum->Write();

    NtracksBvsF->Write();
    NtracksVsMult->Write();

    phiVsEta->Write();
    phiVsEtaDcand->Write();

    TOF->Write();
    TPC->Write();
    kPT->Write();
    piPT->Write();

}

// _________________________________________________________
bool StPicoD0V2AnaMaker::getHadronCorV2(int idxGap) {

    std::vector<int> tracksToRemove = createCandidates();

    double maxNentries = weights->GetMaximum();
    double weightHadron = 1;

    double QcosF[harmonics+1][3] = {0};
    double QcosB[harmonics+1][3] = {0};
    double QsinF[harmonics+1][3] = {0};
    double QsinB[harmonics+1][3] = {0};

    double NtracksB = 0;
    double NtracksF = 0;

    TComplex QvectorF[harmonics+1][3];
    TComplex QvectorB[harmonics+1][3];

	double mEtaGap = etaGap[idxGap];

    int mult = mPicoEvent->grefMult();
    int Ntracks = 0;


    //loop over all tracks
    for(unsigned int i=0;i<mPicoDst->numberOfTracks();++i) {
        StPicoTrack const* hadron = mPicoDst->track(i);

    
        if(!mHFCuts->isGoodTrack(hadron)) continue;

        //PID capability plots
        TPC->Fill(hadron->gMom().Perp(), hadron->dEdx());

        float beta = mHFCuts->getTofBetaBase(hadron);
        if(beta > 0) TOF->Fill(hadron->gMom().Perp(), 1/beta);


        if(!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
        Ntracks++;

        float etaHadron = hadron->gMom().PseudoRapidity();
        float phiHadron = hadron->gMom().Phi();

        if(etaHadron < -1 || etaHadron > 1) continue;
        if(containsId(hadron->id(), tracksToRemove)) continue; //excluding daughter particles

        //control plots
        if(mHFCuts->isGoodKaon(hadron)) kPT->Fill(hadron->gMom().Perp());
        if(mHFCuts->isGoodPion(hadron)) piPT->Fill(hadron->gMom().Perp());

        //REWEIGHTING -- applying phi plots
        weightHadron = maxNentries/(weights->GetBinContent( weights->FindBin(phiHadron) ));

        //GENERIC FRAMEWORK
        //preparing Q vectors
        //
        //backward sample
        if(etaHadron<-0.5*mEtaGap) {
            NtracksB++;
            for(int harm = 0; harm < harmonics+1; harm++) {
                for (int ipow = 0; ipow < 3; ipow++) {
                    QcosB[harm][ipow] += TMath::Power(weightHadron, ipow) * TMath::Cos(harm * phiHadron);
                    QsinB[harm][ipow] += TMath::Power(weightHadron, ipow) * TMath::Sin(harm * phiHadron);
                }
            }
        }
        //forward sample
        if(etaHadron>0.5*mEtaGap) {
            NtracksF++;
            for(int harm = 0; harm < harmonics+1; harm++) {
                for (int ipow = 0; ipow < 3; ipow++) {
                    QcosF[harm][ipow] += TMath::Power(weightHadron, ipow) * TMath::Cos(harm * phiHadron);
                    QsinF[harm][ipow] += TMath::Power(weightHadron, ipow) * TMath::Sin(harm * phiHadron);
                }
            }
        }

        //control plots
        hadron_phi->Fill(phiHadron);
        hadron_check->Fill(phiHadron, weightHadron);
        if(etaHadron>0) hadron_phi_etaP->Fill(phiHadron);
        if(etaHadron<0) hadron_phi_etaN->Fill(phiHadron);
        phiVsEta->Fill(phiHadron, etaHadron);
    }

    //filling real and imaginary part of Q vector
    for(int harm = 0; harm < harmonics+1; harm++) {
        for (int ipow = 0; ipow < 3; ipow++) {
            QvectorB[harm][ipow] = TComplex(QcosB[harm][ipow], QsinB[harm][ipow]);
            QvectorF[harm][ipow] = TComplex(QcosF[harm][ipow], QsinF[harm][ipow]);
        }
    }

    if(NtracksB==0 || NtracksF==0)
        return false;

    double c22 = ((QvectorB[2][1]*(TComplex::Conjugate(QvectorF[2][1])))/(QvectorB[0][1]*QvectorF[0][1])).Re();
    refFlow->Fill(mult, c22);
    refFlow2->Fill(mult, c22);

    //Ntracks vs cumulant
    NtracksBvsCum->Fill(NtracksB, c22);
    NtracksFvsCum->Fill(NtracksF, c22);

    NtracksBvsF->Fill(NtracksB, NtracksF);
    NtracksVsMult->Fill(Ntracks, mult);

    return true;
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::getCorV2(StHFPair *kp,double weight, int flag) {
    int mult = mPicoEvent->grefMult();

    double maxNentries = weights->GetMaximum();

    double weightDcan = maxNentries/(weights->GetBinContent( weights->FindBin( kp->phi() ) ));

    TComplex QvecD[harmonics+1][3] = {0};
    TComplex QvecHadrons[harmonics+1][3] = {0};
    double QcosH[harmonics+1][3] = {0};
    double QsinH[harmonics+1][3] = {0};

    int pTbin = 0;
    for (int pT = 0; pT < nptBins; pT++) {
        if(pair->pt() >= momBins[pT] && pair->pt() < momBins[pT+1]) pTbin = pT;
    }

    double ntracks = 0;

    //charged hadrons loop, excluding daughter particles
    for (unsigned int i = 0; i < mPicoDst->numberOfTracks(); i++) {
        StPicoTrack const *hadron = mPicoDst->track(i);
        if (!mHFCuts->isGoodTrack(hadron)) continue;
        if (!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
        if (i == kp->particle1Idx() || i == kp->particle2Idx()) continue;
        float etaHadron = hadron->gMom().PseudoRapidity();
        float phiHadron = hadron->gMom().Phi();
        double weightHadron = maxNentries / (weights->GetBinContent(weights->FindBin(phiHadron)));

        if(etaHadron < -1 || etaHadron > 1) continue;

        if (!isEtaGap(kp->eta(), etaGap[1], etaHadron)) continue;
        for(int harm = 0; harm < harmonics+1; harm++) {
            for (int ipow = 0; ipow < 3; ipow++) {
                QcosH[harm][ipow] += TMath::Power(weightHadron, ipow) * TMath::Cos(harm * phiHadron);
                QsinH[harm][ipow] += TMath::Power(weightHadron, ipow) * TMath::Sin(harm * phiHadron);
            }
        }

        ntracks++;
    }

    if(ntracks < 10) return false;

    //creating Q vector of D candidate (p vector in GF paper) and Q vector of charged hadrons
    for(int harm = 0; harm < harmonics+1; harm++) {
        for (int ipow = 0; ipow < 3; ipow++) {
            QvecD[harm][ipow] = TComplex( TMath::Power(weightDcan, ipow) * TMath::Cos(harm * kp->phi() )  , TMath::Power(weightDcan, ipow) * TMath::Sin(harm * kp->phi() ) );
            QvecHadrons[harm][ipow] = TComplex(QcosH[harm][ipow], QsinH[harm][ipow]);
        }
    }

    double dif22 = ((QvecD[2][1]*(TComplex::Conjugate(QvecHadrons[2][1])))/(QvecD[0][1]*QvecHadrons[0][1]) ).Re();

    //US pairs
    if(flag < 2) {

        //control plots
        D_phi->Fill(kp->phi());
        if (kp->eta() > 0) D_phi_etaP->Fill(kp->phi());
        if (kp->eta() < 0) D_phi_etaN->Fill(kp->phi());
        phiVsEtaDcand->Fill(kp->phi(), kp->eta());

        diFlowMass[pTbin]->Fill(kp->m(), dif22);
        if(kp->m() > meanFit[pTbin] - 3*sigmaFit[pTbin] && kp->m() < meanFit[pTbin] + 3*sigmaFit[pTbin]) dirFlow2->Fill(kp->pt(), dif22);
    }

    //LS pairs
    else{
        diFlowMassBKG[pTbin]->Fill(kp->m(), dif22);
        if(kp->m() > meanFit[pTbin] - 3*sigmaFit[pTbin] && kp->m() < meanFit[pTbin] + 3*sigmaFit[pTbin]) dirFlow2BKG->Fill(kp->pt(), dif22);
    }

    return true;
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::isEtaGap(double dEta,double mGap,double hEta) {
    if(mGap == 0) return true;

    if(dEta > 0 && hEta < -0.5*mGap)
        return true;
    else if(dEta < 0 && hEta > -0.5*mGap)
        return true;
    else
        return false;
}

bool StPicoD0V2AnaMaker::containsId(int id, std::vector<int>& tracksToRemove)
{
	for(unsigned int i = 0; i < tracksToRemove.size(); i++)
	{
		if(id == tracksToRemove[i])
			return true;
	}
}
