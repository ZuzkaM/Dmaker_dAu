#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
<<<<<<< HEAD
#include "TProfile.h"
=======
#include "TVector3.h"

>>>>>>> f6b502a69ce0930c336534a90b74d11c6b2ffcaf
//#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"
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

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;

class StHFPair;
class StHFTriplet;
class StHFCuts;

class TProfile;

class StPicoD0AnaMaker : public StPicoHFMaker
{
public:
    StPicoD0AnaMaker(char const*, StPicoDstMaker*, char const*);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();

    void DeclareHistograms();
    void WriteHistograms();
    bool isGoodHadron(StPicoTrack const*) const;

protected:
    std::vector<unsigned short> mIdxPicoPions;
    std::vector<unsigned short> mIdxPicoKaons;

private:
    int createCandidates();
    int analyzeCandidates();

    TNtuple *ntp_DMeson_Signal;
    TNtuple *ntp_DMeson_Background;
<<<<<<< HEAD
//    TNtuple *ntp_kaon;
//    TNtuple *ntp_pion;

    TProfile *profV2[8][5][3];//i.S or B; j.flatten; k. differetn etaGap
    TH1D *hadronV2[5][3];
    TH1D *hadronV2_sum[5][3];
    TH1D *hadron_phi;
    TH1D *D_phi;
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
    TProfile *qVecPow2[4];
    TProfile *refFlow;
    TProfile *dirFlow[5];
    TProfile *corrD[2][5];
    TProfile *qVec2[4];
    TProfile *refFlow2;
    TProfile *corrD2[2];
    TProfile *dirFlow2; 

=======
>>>>>>> f6b502a69ce0930c336534a90b74d11c6b2ffcaf

    bool getHadronCorV2(int );
    bool getCorV2(StHFPair *, double);
    bool isEtaGap(double, double ,double);
    int mRunNumber;
    TString mOutFileBaseName;

    TFile* mOutFile;

    ClassDef(StPicoD0AnaMaker, 1) //set to 1
};

#endif
