#ifndef StPicoD0V2AnaMaker_h
#define StPicoD0V2AnaMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
#include "TVector3.h"
#include "TComplex.h"

//#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"
#include "TProfile.h"
//#include "StPicoD0AnaHists.h"
#include <vector>
#include "TClonesArray.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
//#include "StPicoHFMaker/StHFTriplet.h"

#include "phys_constants.h"

#include "TH1F.h"
#include "TH3F.h"
#include <ctime>


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


using namespace TMVA;

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;

class StHFPair;
class StHFTriplet;
class StHFCuts;

class StPicoD0V2AnaMaker : public StPicoHFMaker
{
public:
    StPicoD0V2AnaMaker(char const*, StPicoDstMaker*, char const*);
    virtual ~StPicoD0V2AnaMaker();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();
    void DeclareHistograms();
    void WriteHistograms();
    bool isGoodHadron(StPicoTrack const*) const;

protected:

private:
    vector<int> createCandidates();

    TProfile *profV2[8][5][3];//i.S or B; j.flatten; k. differetn etaGap
    TH1D *hadronV2[5][3];
    TH1D *hadronV2_sum[5][3];
    TH1D *hadron_phi;
    TH1D *D_phi;
    TH1D *hadron_phi_etaP;
    TH1D *D_phi_etaP;
    TH1D *hadron_phi_etaN;
    TH1D *D_phi_etaN;
    TH2D *fitPhi[6];
    TH2D *massPt;
    TH2D *massPtLike;
    TH2D *massLike;
    TH2D *massLike2;
    TH2D *massUnlike;
    TH2D *v2Weight[8][3];
    TH2F *hPhiD[8][3];
    TH2F *hPhiHadron[8][3];
    TH2D *likeV2Mass[6][5];
    TH2D *likeV2Mass2[6][5];
    TH2D *unlikeV2Mass[6][5];
    TProfile *V2Mass[2][6][5];
    TProfile *candPt;

    TProfile *qVec[4];
    TProfile *refFlow;
    TProfile *dirFlow[5];
    TProfile *corrD[2][5];
    TProfile *dirFlowBKG[5];
    TProfile *corrDBKG[2][5];
    TProfile *qVec2[4];
    TProfile *refFlow2;
    TProfile *corrD2[2];
    TProfile *dirFlow2;
    TProfile *corrD2BKG[2];
    TProfile *dirFlow2BKG;
    TProfile *cosH;
    TProfile *sinH;
    TProfile *diFlowMass[3];
    TProfile *diFlowMassBKG[3];

    TH1D *weights;
    TH1D *hadron_check;

    TH1D *mass[3];
    TH1D *massBKG[3];

    TProfile *NtracksFvsCum;
    TProfile *NtracksBvsCum;
    TH2D *NtracksBvsF;
    TH2D *NtracksVsMult;

    TH2D *phiVsEta;
    TH2D *phiVsEtaDcand;
    TH2D *TOF;
    TH2D *TPC;
    TH1D *kPT;
    TH1D *piPT;

    TMVA::Reader *reader[3];
    Float_t k_pt[3], pi1_pt[3], k_dca[3], pi1_dca[3], dcaDaughters[3], cosTheta[3], D_decayL[3], dcaD0ToPv[3];




    bool getHadronCorV2(int );
    bool getCorV2(StHFPair *, double, int);
    bool isEtaGap(double, double ,double);
    bool containsId(int id, std::vector<int>& tracksToRemove);


    TString mOutFileBaseName;

    TFile* mOutFile;

    ClassDef(StPicoD0V2AnaMaker, 1) //set to 1
};

#endif
