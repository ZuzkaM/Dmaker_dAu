#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoD0V2AnaMaker.h"
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

			getCorV2(pair, 1);

			tracksofCand.push_back(pion1->id());
			tracksofCand.push_back(kaon->id());


        }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)

    return tracksofCand;
}

// _________________________________________________________
void StPicoD0V2AnaMaker::DeclareHistograms() {
    TString names[4] = {"cos_B", "cos_F", "sin_B", "sin_F"}; //backward and forward samples
//    float multBin[6] = {0, 7, 12, 16, 22, 100};
    int nMultBins = sizeof(multBin)/sizeof(multBin[0])-1;

    float momBins[7] = {0,1,2,3,4,5,10};
    int nMomBins = sizeof(momBins)/sizeof(momBins[0])-1;

    for(int m = 0; m < 4; m++) {
        qVec[m] = new TProfile(names[m].Data(),"Q vector", nMultBins, multBin);
        qVec[m]->Sumw2();

        qVecPow2[m] = new TProfile((names[m]+"Pow2").Data(),"Q vector", nMultBins, multBin);
        qVecPow2[m]->Sumw2();

        qVec2[m] = new TProfile((names[m]+"_no_mult").Data(),"Q vector", 1, 0, 100);
        qVec2[m]->Sumw2();
    }

    for(int m = 0; m < 5; m++) {
        corrD[0][m] = new TProfile(Form("cosD_%.0f_%.0f", multBin[m], multBin[m+1]), "", nMomBins, momBins);
        corrD[1][m] = new TProfile(Form("sinD_%.0f_%.0f", multBin[m], multBin[m+1]), "", nMomBins, momBins);
        dirFlow[m] = new TProfile(Form("dirFlow_%.0f_%.0f", multBin[m], multBin[m+1]), "", nMomBins, momBins);
    }

    corrD2[0] = new TProfile("cosD_no_mult","cos D", nMomBins, momBins);
    corrD2[1] = new TProfile("sinD_no_mult","sin D", nMomBins, momBins);
    dirFlow2 = new TProfile("dirFlow_no_mult","dir flow", nMomBins, momBins);
    refFlow = new TProfile("refFlow", "", nMultBins, multBin);
    refFlow2 = new TProfile("refFlow_no_mult", "", 1, 0, 100);

    hadron_phi = new TH1D("hadron_phi", "Hadron phi", 2000, -5, 5);
    D_phi = new TH1D("D_phi", "D phi", 2000, -5, 5);
}

// _________________________________________________________
void StPicoD0V2AnaMaker::WriteHistograms() {
    for(int m = 0; m < 4; m++) {
        qVec[m]->Write();
        qVecPow2[m]->Write();
        qVec2[m]->Write();
    }
    refFlow->Write();

    for(int m = 0; m < 5; m++) {
        corrD[0][m]->Write();
        corrD[1][m]->Write();
        dirFlow[m]->Write();
    }

    corrD2[0]->Write();
    corrD2[1]->Write();
    dirFlow2->Write();
    refFlow2->Write();

    hadron_phi->Write();
    D_phi->Write();
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::getHadronCorV2(int idxGap) {
	std::vector<int> tracksToRemove = createCandidates();
	
	double etaGap[3] = {0,0.15,0.05};
    double mEtaGap = etaGap[idxGap];
    float hadronFill[7] = {0};
    const double reweight = 1;//mGRefMultCorrUtil->getWeight();
    // int centrality  = mGRefMultCorrUtil->getCentralityBin9();
    int mult = mPicoEvent->grefMult();

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();++i) {
        StPicoTrack const* hadron = mPicoDst->track(i);
        if(containsId(hadron->id(), tracksToRemove)) continue;

        if(!mHFCuts->isGoodTrack(hadron)) continue;
        if(!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
        float etaHadron = hadron->gMom().PseudoRapidity();
        float phiHadron = hadron->gMom().Phi();

        if(etaHadron<-0.5*mEtaGap) {//backward sample
            hadronFill[0]++;
            hadronFill[1] += sin(2 * phiHadron);
            hadronFill[2] += cos(2 * phiHadron);
        }

        if(etaHadron>0.5*mEtaGap) {//forward sample
            hadronFill[3]++;
            hadronFill[4] += sin(2 * phiHadron);
            hadronFill[5] += cos(2 * phiHadron);
        }
        hadron_phi->Fill(phiHadron);
    }

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
        qVecPow2[0]->Fill(mult, (hadronFill[2]*hadronFill[2])/(hadronFill[0]*hadronFill[0]), reweight);
        qVecPow2[1]->Fill(mult, (hadronFill[5]*hadronFill[5])/(hadronFill[3]*hadronFill[3]), reweight);
        qVecPow2[2]->Fill(mult, (hadronFill[1]*hadronFill[1])/(hadronFill[0]*hadronFill[0]), reweight);
        qVecPow2[3]->Fill(mult, (hadronFill[4]*hadronFill[4])/(hadronFill[3]*hadronFill[3]), reweight);
        refFlow->Fill(mult, ((hadronFill[2]*hadronFill[5])/(hadronFill[0]*hadronFill[3])), reweight);
        //no mult
        qVec2[0]->Fill(mult, hadronFill[2]/hadronFill[0], reweight);
        qVec2[1]->Fill(mult, hadronFill[5]/hadronFill[3], reweight);
        qVec2[2]->Fill(mult, hadronFill[1]/hadronFill[0], reweight);
        qVec2[3]->Fill(mult, hadronFill[4]/hadronFill[3], reweight);
        refFlow2->Fill(mult, ((hadronFill[2]*hadronFill[5])/(hadronFill[0]*hadronFill[3])), reweight);
    }
    return true;
}

// _________________________________________________________
bool StPicoD0V2AnaMaker::getCorV2(StHFPair *kp,double weight) {
    int mult = mPicoEvent->grefMult();

    double hadronv2=1;
//    float multBin[6] = {0,7,12,16,22,100};
    double etaGap[3] = {0,0.15,0.05};

    int k=0;
    double corFill[7] = {0};
    corFill[0] = 1 ;
    corFill[1] = sin(2* kp->phi())/sqrt(hadronv2);
    corFill[2] = cos(2* kp->phi())/sqrt(hadronv2);

    D_phi->Fill(kp->phi());

    for(unsigned int i=0; i<mPicoDst->numberOfTracks();i++) {
        StPicoTrack const* hadron = mPicoDst->track(i);
        if(!mHFCuts->isGoodTrack(hadron)) continue;
        if(!mHFCuts->isGoodProton(hadron) && !mHFCuts->isGoodKaon(hadron) && !mHFCuts->isGoodPion(hadron)) continue;
        if(i==kp->particle1Idx() || i==kp->particle2Idx()) continue;
        float etaHadron = hadron->gMom().PseudoRapidity();
        float phiHadron = hadron->gMom().Phi();
        if(!isEtaGap(kp->eta(),etaGap[k],etaHadron))  continue;
        corFill[3]++;
        corFill[4] += sin(2*phiHadron)/sqrt(hadronv2);
        corFill[5] += cos(2*phiHadron)/sqrt(hadronv2);

        for(int m = 0; m < 5; m++) {
            if(mult >= multBin[m] && mult < multBin[m+1]) {
                corrD[0][m]->Fill(kp->pt(), corFill[2], weight);
                corrD[1][m]->Fill(kp->pt(), corFill[1], weight);
                dirFlow[m]->Fill(kp->pt(), corFill[2]*corFill[5]/corFill[3], weight);
            }
        }

        corrD2[0]->Fill(kp->pt(), corFill[2], weight);
        corrD2[1]->Fill(kp->pt(), corFill[1], weight);
        dirFlow2->Fill(kp->pt(), corFill[2]*corFill[5]/corFill[3], weight);
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
        return hEta<-mGap;
    else
        return hEta>mGap;
}

bool StPicoD0V2AnaMaker::containsId(int id, std::vector<int>& tracksToRemove)
{
	for(unsigned int i = 0; i < tracksToRemove.size(); i++)
	{
		if(id == tracksToRemove[i])
			return true;
	}
}