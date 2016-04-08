#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCut.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"

#include <sstream>
#include <math.h>

#include <iostream>

#include <algorithm>


using namespace std;



//------------------------------------------------------------------------------
// Member functions
//------------------------------------------------------------------------------


int Wait() {
  cout << " Continue [<RET>|q]?  "; 
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0;
}

  
TH2F*   LoadHistogram(TString filename,
		      TString hname,
		      TString cname);
  
bool  IsIsolatedLepton(int k);
float ElectronIsolation(int k);
float MuonIsolation(int k);
bool  IsTightLepton(int k);
bool  IsLooseLepton(int k);
bool  IsLooseIsolatedLepton(int k);

Double_t GetPPFWeight(int a, int b, int s = 5);
Double_t GetPFFWeight(int a, int b, int s = 5);
Double_t GetFFFWeight(int a, int b, int s = 5);

Double_t GetPFWeight(int a, int b, int s = 5);
Double_t GetFFWeight(int a, int b, int s = 5);

Double_t GetFactor(TH2F*   h2,
		  Float_t leptonPt,
		  Float_t leptonEta,
		  Float_t leptonPtMax = -999.);

Double_t GetFactorError(TH2F*   h2,
		       Float_t leptonPt,
		       Float_t leptonEta,
		       Float_t leptonPtMax = -999.);

//------------------------------------------------------------------------------
// Define structures & GLB variables 
//------------------------------------------------------------------------------

const int nMuJetET = 7;
const int nEleJetET = 3;

TH2F*                         MuonFR[nMuJetET];
TH2F*                         ElecFR[nEleJetET];

TString muonJetPt[nMuJetET] = {"10","15", "20", "25", "30", "35", "45"};
enum { mu10, mu15, mu20, mu25, mu30, mu35, mu45};

TString elecJetPt[nEleJetET] = {"25", "35", "45"};

enum {ele25, ele35, ele45}; 

enum {ElUp, ElDown, MuUp, MuDown};

enum {Muon, Electron, Unknown};

enum {noTight, Tight};

int ch;
 
enum {mm=0,ee=1,em=2,me=3,  mmm=4, eee=5, mme=6, eem =7};

struct Lepton { 

  UInt_t         index;
  UInt_t         flavour;  // Muon, Electron
  UInt_t         type;    // Tight, Loose
  Float_t        charge;
  Float_t        pt;
  Float_t        eta;
  Float_t        fr;
  Float_t        pr;

};

std::vector<Lepton> AnalysisLeptons;
Float_t          run;
Float_t          evt;
Float_t         njet;
Float_t          jetRho;
Float_t          channel;
Float_t          baseW;
Float_t          trigger;
Float_t          mll;
vector<float>   *std_vector_lepton_isLooseLepton;
vector<float>   *std_vector_lepton_isTightLepton;
vector<float>   *std_vector_lepton_trackIso;
vector<float>   *std_vector_electron_ecalPFClusterIso;
vector<float>   *std_vector_electron_hcalPFClusterIso;
vector<float>   *std_vector_lepton_flavour;
vector<float>   *std_vector_lepton_pt;
vector<float>   *std_vector_lepton_eta;
vector<float>   *std_vector_lepton_phi;
vector<float>   *std_vector_lepton_isTightMuon; 
vector<float>   *std_vector_lepton_isMediumMuon; 
vector<float>   *std_vector_lepton_eleIdMedium; 
vector<float>   *std_vector_lepton_eleIdVeto; 
vector<float>   *std_vector_lepton_eleIdTight;
vector<float>   *std_vector_lepton_chargedHadronIso; 
vector<float>   *std_vector_lepton_photonIso;
vector<float>   *std_vector_lepton_neutralHadronIso;
vector<float>   *std_vector_lepton_sumPUPt;
vector<float>   *std_vector_electron_effectiveArea;
vector<float>   *std_vector_lepton_dz;
vector<float>   *std_vector_lepton_d0;

vector<float>   *std_vector_electron_dEtaIn;
vector<float>   *std_vector_electron_dPhiIn;
vector<float>   *std_vector_electron_full5x5_sigmaIetaIeta;
vector<float>   *std_vector_electron_hOverE;
vector<float>   *std_vector_electron_ooEmooP;
vector<float>   *std_vector_electron_expectedMissingInnerHits;
vector<float>   *std_vector_electron_passConversionVeto;
vector<float>   *std_vector_lepton_closejet_PartonFlavour;

//------------------------------------------------------------------------------
// MAIN function
//------------------------------------------------------------------------------

int addWJetsWeights(string input="latino_ll_LP_test.root",
		    float luminosity = 1.0,
		    string indir = "input/",
		    string outdir = "output/",
		    bool addWeight=false) {

  //addWJetsWeights("latino_Run2015D_05Oct2015_DoubleMuon_0.root", 1.0, "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/Data13TeV/21Oct/25ns/","output/", true, "WW"); 

  //  addWJetsWeights("latino_WJetsToLNu.root", 1.0, "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/21Oct_25ns_MC/", "output/", true, "WW"); 
 //  addWJetsWeights("latino_WJetsToLNu.root", 1.0, "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/21Oct_25ns_MC__l2sel__hadd/", "output/", true, "WW"); 

  printf("Start W+jets weights with %s, lumi=%.3f, indir=%s, outdir=%s, addWeight option %d\n",
	 input.c_str(), luminosity, indir.c_str(), outdir.c_str(), int(addWeight));
  
  
  // ------ root settings ---------
  //----------------------------------------------------------------------------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(0.80,"X");
  gStyle->SetTitleOffset(0.80,"Y");
  // ------------------------------ 


  string openMode = "read";


  // SF, FR, PR and trigger efficiencies histograms
  //----------------------------------------------------------------------------


  TH2F*                         MuonPR;
  TH2F*                         ElecPR;

 
  MuonPR = LoadHistogram("ZJets_MuPR_RunII_25ns_jet20_76X", "h_Muon_signal_pt_eta_bin", "MuonPR");

  ElecPR = LoadHistogram("ZJets_ElePR_RunII_25ns_jet35_76X", "h_Ele_signal_pt_eta_bin", "ElecPR");
  
  for ( int j=0; j<nMuJetET; j++) { 
    MuonFR[j] = LoadHistogram(Form("MuFR_RunII_25ns_jet%s_226pb_76X", muonJetPt[j].Data()), "FR_pT_eta_EWKcorr", Form("MuonFR_Jet%s", muonJetPt[j].Data()));
  }

  for ( int j=0; j<nEleJetET; j++) { 
    ElecFR[j] = LoadHistogram(Form("EGFR_RunII_25ns_jet%s_226pb_76X_New", elecJetPt[j].Data()), "FR_pT_eta_EWKcorr", Form("ElecFR_Jet%s", elecJetPt[j].Data()));
   

}
 

 //--- OPEN LOOSE-LOOSE FILE 
 //----------------------------------------------------------------------------


  TFile* file = new TFile(Form("%s%s", indir.c_str(),input.c_str()),openMode.c_str());

  TTree* tree = (TTree*)file->Get("latino");

  printf("input tree from %s%s: %d events\n", indir.c_str(),input.c_str(), (int)tree->GetEntries());

  tree->SetBranchAddress("evt", &evt);
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("njet", &njet);
 tree->SetBranchAddress("baseW  ", &baseW);
 tree->SetBranchAddress("jetRho", &jetRho);
 tree->SetBranchAddress("channel", &channel);
 tree->SetBranchAddress("trigger", &trigger);
 tree->SetBranchAddress("mll", &mll); 
tree->SetBranchAddress("std_vector_lepton_isLooseLepton", &std_vector_lepton_isLooseLepton);
 tree->SetBranchAddress("std_vector_lepton_isTightLepton", &std_vector_lepton_isTightLepton);
 tree->SetBranchAddress("std_vector_lepton_trackIso", &std_vector_lepton_trackIso);
 tree->SetBranchAddress("std_vector_electron_ecalPFClusterIso", &std_vector_electron_ecalPFClusterIso);
 tree->SetBranchAddress("std_vector_electron_hcalPFClusterIso", &std_vector_electron_hcalPFClusterIso);
 tree->SetBranchAddress("std_vector_lepton_pt", &std_vector_lepton_pt);
 tree->SetBranchAddress("std_vector_lepton_eta", &std_vector_lepton_eta);
 tree->SetBranchAddress("std_vector_lepton_phi", &std_vector_lepton_phi);
 tree->SetBranchAddress("std_vector_lepton_flavour", &std_vector_lepton_flavour);
 tree->SetBranchAddress("std_vector_lepton_isTightMuon", &std_vector_lepton_isTightMuon);
 tree->SetBranchAddress("std_vector_lepton_isMediumMuon", &std_vector_lepton_isMediumMuon);
 tree->SetBranchAddress("std_vector_lepton_eleIdMedium", &std_vector_lepton_eleIdMedium);
 tree->SetBranchAddress("std_vector_lepton_eleIdVeto", &std_vector_lepton_eleIdVeto);
 tree->SetBranchAddress("std_vector_lepton_eleIdTight", &std_vector_lepton_eleIdTight); 
 tree->SetBranchAddress("std_vector_lepton_chargedHadronIso", &std_vector_lepton_chargedHadronIso);
 tree->SetBranchAddress("std_vector_lepton_photonIso", &std_vector_lepton_photonIso);
 tree->SetBranchAddress("std_vector_lepton_neutralHadronIso", &std_vector_lepton_neutralHadronIso);
 tree->SetBranchAddress("std_vector_lepton_sumPUPt", &std_vector_lepton_sumPUPt);
 tree->SetBranchAddress("std_vector_electron_effectiveArea", &std_vector_electron_effectiveArea);
 tree->SetBranchAddress("std_vector_lepton_dz", &std_vector_lepton_dz);
 tree->SetBranchAddress("std_vector_lepton_d0", &std_vector_lepton_d0);
 tree->SetBranchAddress("std_vector_lepton_closejet_PartonFlavour", &std_vector_lepton_closejet_PartonFlavour);

  tree->SetBranchAddress("std_vector_electron_dEtaIn", &std_vector_electron_dEtaIn);
  tree->SetBranchAddress("std_vector_electron_dPhiIn", &std_vector_electron_dPhiIn);
  tree->SetBranchAddress("std_vector_electron_full5x5_sigmaIetaIeta", &std_vector_electron_full5x5_sigmaIetaIeta);
  tree->SetBranchAddress("std_vector_electron_hOverE",&std_vector_electron_hOverE);
  tree->SetBranchAddress("std_vector_electron_ooEmooP",&std_vector_electron_ooEmooP);
  tree->SetBranchAddress("std_vector_electron_expectedMissingInnerHits", &std_vector_electron_expectedMissingInnerHits);
  tree->SetBranchAddress("std_vector_electron_passConversionVeto", &std_vector_electron_passConversionVeto);


//--- New branches to be added
//----------------------------------------------------------------------------
 
 Float_t weightFF, weightPF, weightUp, weightDown; 
 Float_t weightFFF, weightPFF, weightPPF; 

 Float_t weight3lForWW ,weight3l, weight3lElUp, weight3lElDown, weight3lstatElUp, weight3lstatElDown, weight3lMuUp, weight3lMuDown, weight3lstatMuUp, weight3lstatMuDown; 
 Float_t weight2l0j, weight2l1j, weight2l2j, weight2l0jElUp, weight2l1jElUp, weight2l2jElUp, weight2l0jElDown, weight2l1jElDown, weight2l2jElDown, weight2l0jstatElUp, weight2l1jstatElUp, weight2l2jstatElUp, weight2l0jstatElDown, weight2l1jstatElDown, weight2l2jstatElDown,  weight2l0jMuUp, weight2l1jMuUp, weight2l2jMuUp, weight2l0jMuDown, weight2l1jMuDown, weight2l2jMuDown, weight2l0jstatMuUp, weight2l1jstatMuUp, weight2l2jstatMuUp, weight2l0jstatMuDown, weight2l1jstatMuDown, weight2l2jstatMuDown;

 //Float_t weight2lMu, weight2lEle;

 tree->SetBranchStatus("*", 1);
 // in case already exist
 tree->SetBranchStatus("fakeW3l",               0);
 tree->SetBranchStatus("fakeW3lElUp",             0);
 tree->SetBranchStatus("fakeW3lElDown",           0);
 tree->SetBranchStatus("fakeW3lstatElUp",         0);
 tree->SetBranchStatus("fakeW3lstatElDown",       0);
 tree->SetBranchStatus("fakeW3lMuUp",             0);
 tree->SetBranchStatus("fakeW3lMuDown",           0);
 tree->SetBranchStatus("fakeW3lstatMuUp",         0);
 tree->SetBranchStatus("fakeW3lstatMuDown",       0);
 tree->SetBranchStatus("fakeW2l0j",             0);
 //tree->SetBranchStatus("fakeW2lMu",             0);
 //tree->SetBranchStatus("fakeW2lEle",            0);
 tree->SetBranchStatus("fakeW2l1j",             0);
 tree->SetBranchStatus("fakeW2l2j",             0);
 tree->SetBranchStatus("fakeW2l0jstatUp",       0);
 tree->SetBranchStatus("fakeW2l1jstatUp",       0);
 tree->SetBranchStatus("fakeW2l2jstatUp",       0);
 tree->SetBranchStatus("fakeW2l0jstatDown",     0);
 tree->SetBranchStatus("fakeW2l1jstatDown",     0);
 tree->SetBranchStatus("fakeW2l2jstatDown",     0);
 tree->SetBranchStatus("fakeW2l0jUp",           0);
 tree->SetBranchStatus("fakeW2l1jUp",           0);
 tree->SetBranchStatus("fakeW2l2jUp",           0);
 tree->SetBranchStatus("fakeW2l0jDown",         0);
 tree->SetBranchStatus("fakeW2l1jDown",         0);
 tree->SetBranchStatus("fakeW2l2jDown",         0);

 //tree->SetBranchStatus("fakeWJet",        0);
 //tree->SetBranchStatus("fakeQCD",         0);

 /* to be added later 
 tree->SetBranchStatus("fakeWUp",         0);
 tree->SetBranchStatus("fakeWDown",       0); 
 tree->SetBranchStatus("fakeWMuUp",       0);
 tree->SetBranchStatus("fakeWMuDown",     0);
 tree->SetBranchStatus("fakeWElUp",       0);
 tree->SetBranchStatus("fakeWElDown",     0);
 tree->SetBranchStatus("fakeWMuStatUp",   0);
 tree->SetBranchStatus("fakeWMuStatDown", 0);
 tree->SetBranchStatus("fakeWElStatUp",   0);
 tree->SetBranchStatus("fakeWElStatDown", 0);
 */


 TFile *newfile;
 TTree *newtree;

 if(addWeight) { 

   gSystem->Exec(Form("mkdir %s", outdir.c_str()));
   newfile = new TFile(Form("%s/%s",outdir.c_str(),input.c_str()),"recreate");
   newfile->cd();
   // make a copy of the input tree
   newtree = tree->CloneTree(0);
   newtree->SetAutoSave(10000000);  // autosave when 10 Mbyte written
   newtree->Branch("fakeW3l",              &weight3l,                "fakeW3l/F");
   //newtree->Branch("fakeW2l",            &weight2l,                "fakeW2l/F");
   newtree->Branch("fakeW3lElUp",          &weight3lElUp,              "fakeW3lElUp/F");
   newtree->Branch("fakeW3lElDown",        &weight3lElDown,            "fakeW3lElDown/F");
   newtree->Branch("fakeW3lstatElUp",      &weight3lstatElUp,          "fakeW3lstatElUp/F");
   newtree->Branch("fakeW3lstatElDown",    &weight3lstatElDown,        "fakeW3lstatElDown/F");
   newtree->Branch("fakeW3lMuUp",          &weight3lMuUp,              "fakeW3lMuUp/F");
   newtree->Branch("fakeW3lMuDown",        &weight3lMuDown,            "fakeW3lMuDown/F");
   newtree->Branch("fakeW3lstatMuUp",      &weight3lstatMuUp,          "fakeW3lstatMuUp/F");
   newtree->Branch("fakeW3lstatMuDown",    &weight3lstatMuDown,        "fakeW3lstatMuDown/F");
   newtree->Branch("fakeW2l0j",            &weight2l0j,              "fakeW2l0j/F");
   newtree->Branch("fakeW2l1j",            &weight2l1j,              "fakeW2l1j/F");
   newtree->Branch("fakeW2l2j",            &weight2l2j,              "fakeW2l2j/F");
   //newtree->Branch("fakeW2lEle",         &weight2lEle,              "fakeW2lEle/F");
   //newtree->Branch("fakeW2lMu",          &weight2lMu,               "fakeW2lMu/F");
   newtree->Branch("fakeW2l0jElUp",        &weight2l0jElUp,            "fakeW2l0jElUp/F");
   newtree->Branch("fakeW2l1jElUp",        &weight2l1jElUp,            "fakeW2l1jElUp/F");
   newtree->Branch("fakeW2l2jElUp",        &weight2l2jElUp,            "fakeW2l2jElUp/F");
   newtree->Branch("fakeW2l0jElDown",      &weight2l0jElDown,          "fakeW2l0jElDown/F");
   newtree->Branch("fakeW2l1jElDown",      &weight2l1jElDown,          "fakeW2l1jElDown/F");
   newtree->Branch("fakeW2l2jElDown",      &weight2l2jElDown,          "fakeW2l2jElDown/F");
   newtree->Branch("fakeW2l0jstatElUp",    &weight2l0jstatElUp,        "fakeW2l0jstatElUp/F");
   newtree->Branch("fakeW2l1jstatElUp",    &weight2l1jstatElUp,        "fakeW2l1jstatElUp/F");
   newtree->Branch("fakeW2l2jstatElUp",    &weight2l2jstatElUp,        "fakeW2l2jstatElUp/F");
   newtree->Branch("fakeW2l0jstatElDown",  &weight2l0jstatElDown,      "fakeW2l0jstatElDown/F");
   newtree->Branch("fakeW2l1jstatElDown",  &weight2l1jstatElDown,      "fakeW2l1jstatElDown/F");
   newtree->Branch("fakeW2l2jstatElDown",  &weight2l2jstatElDown,      "fakeW2l2jstatElDown/F");



   newtree->Branch("fakeW2l0jMuUp",        &weight2l0jMuUp,            "fakeW2l0jMuUp/F");
   newtree->Branch("fakeW2l1jMuUp",        &weight2l1jMuUp,            "fakeW2l1jMuUp/F");
   newtree->Branch("fakeW2l2jMuUp",        &weight2l2jMuUp,            "fakeW2l2jMuUp/F");
   newtree->Branch("fakeW2l0jMuDown",      &weight2l0jMuDown,          "fakeW2l0jMuDown/F");
   newtree->Branch("fakeW2l1jMuDown",      &weight2l1jMuDown,          "fakeW2l1jMuDown/F");
   newtree->Branch("fakeW2l2jMuDown",      &weight2l2jMuDown,          "fakeW2l2jMuDown/F");
   newtree->Branch("fakeW2l0jstatMuUp",    &weight2l0jstatMuUp,        "fakeW2l0jstatMuUp/F");
   newtree->Branch("fakeW2l1jstatMuUp",    &weight2l1jstatMuUp,        "fakeW2l1jstatMuUp/F");
   newtree->Branch("fakeW2l2jstatMuUp",    &weight2l2jstatMuUp,        "fakeW2l2jstatMuUp/F");
   newtree->Branch("fakeW2l0jstatMuDown",  &weight2l0jstatMuDown,      "fakeW2l0jstatMuDown/F");
   newtree->Branch("fakeW2l1jstatMuDown",  &weight2l1jstatMuDown,      "fakeW2l1jstatMuDown/F");
   newtree->Branch("fakeW2l2jstatMuDown",  &weight2l2jstatMuDown,      "fakeW2l2jstatMuDown/F");

   //newtree->Branch("fakeWJet",         &weightPF,              "fakeWJet/F");
   //newtree->Branch("fakeQCD",          &weightFF,              "fakeQCD/F");

   /* to be added later 
    newtree->Branch("fakeWUp",          &weightUp,                 "fakeWUp/F");
    newtree->Branch("fakeWDown",        &weightDown,             "fakeWDown/F");
    newtree->Branch("fakeWMuUp",        &weightMuUp,             "fakeWMuUp/F");
    newtree->Branch("fakeWMuDown",      &weightMuDown,         "fakeWMuDown/F");
    newtree->Branch("fakeWElUp",        &weightElUp,             "fakeWElUp/F");
    newtree->Branch("fakeWElDown",      &weightElDown,         "fakeWElDown/F");
    newtree->Branch("fakeWMuStatUp",    &weightMuStatUp,     "fakeWMuStatUp/F");
    newtree->Branch("fakeWMuStatDown",  &weightMuStatDown, "fakeWMuStatDown/F");
    newtree->Branch("fakeWElStatUp",    &weightElStatUp,     "fakeWElStatUp/F");
    newtree->Branch("fakeWElStatDown",  &weightElStatDown, "fakeWElStatDown/F");
   */
 }


 // Loop over events
 //----------------------------------------------------------------------------
  
 double sumOfWeightsMuMu(0.);
 double sumOfWeightsMuEl(0.);
 double sumOfWeightsElMu(0.);
 double sumOfWeightsElEl(0.);

 double sumOfQCDWeightsMuMu(0.);
 double sumOfQCDWeightsMuEl(0.);
 double sumOfQCDWeightsElMu(0.);
 double sumOfQCDWeightsElEl(0.);

 int mmPP(0),mmPF(0),mmFF(0);
 int mePP(0),mePF(0),meFF(0);
 int emPP(0),emPF(0),emFF(0);
 int eePP(0),eePF(0),eeFF(0);

 Long64_t nentries = tree->GetEntries();

 unsigned int countEntriesByHand(0);
 unsigned int countEntriesByHand2(0);


  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    
    tree->GetEntry(jentry);
  
    weight2l0j = 0.0; weight2l1j = 0.0; weight2l2j = 0.0;
    weight2l0jElUp = 0.0; weight2l1jElUp = 0.0; weight2l2jElUp = 0.0;
    weight2l0jElDown = 0.0; weight2l1jElDown = 0.0; weight2l2jElDown = 0.0;
    weight2l0jstatElUp = 0.0; weight2l1jstatElUp = 0.0; weight2l2jstatElUp = 0.0;
    weight2l0jstatElDown = 0.0; weight2l1jstatElDown = 0.0; weight2l2jstatElDown = 0.0;
    weight2l0jMuUp = 0.0; weight2l1jMuUp = 0.0; weight2l2jMuUp = 0.0;
    weight2l0jMuDown = 0.0; weight2l1jMuDown = 0.0; weight2l2jMuDown = 0.0;
    weight2l0jstatMuUp = 0.0; weight2l1jstatMuUp = 0.0; weight2l2jstatMuUp = 0.0;
    weight2l0jstatMuDown = 0.0; weight2l1jstatMuDown = 0.0; weight2l2jstatMuDown = 0.0;
    weightPF = 0.0; 
    weightFF = 0.0; 

    
    weight3l = 0.0; weight3lElUp = 0.0; weight3lElDown = 0.0; weight3lMuUp = 0.0; weight3lMuDown = 0.0;
    weight3lstatElUp = 0.0; weight3lstatElDown = 0.0;  weight3lstatMuUp = 0.0; weight3lstatMuDown = 0.0;
    weightPPF = 0.0; 
    weightPFF = 0.0; 
    weightFFF = 0.0; 

    //weight2lMu = 0.0;  weight2lEle = 0.0;

    float sign = 1.; 
    //if ( input.Contains("Zg") ) sign = -1.; // flip sign to subtract non-fake contamination from Vg/Vg*
    
    if (jentry%10000==0) printf("Processing event %d \n", int(jentry));

    if ( trigger == 0) continue; 


    // Loop over leptons
    //--------------------------------------------------------------------------
    AnalysisLeptons.clear();
     
    int vsize = std_vector_lepton_pt->size();
 

    for (int i=0; i<vsize; i++) {
      
      // if (i< 2)  cout << i << "  "  << evt << "  " << run<< endl;

		   
      //if ( !IsLooseLepton(i)) continue;  
      //if ( !IsLooseIsolatedLepton(i) ) continue;   


      if (std_vector_lepton_isLooseLepton->at(i) != 1 ) continue;
     
      
      float pt  = std_vector_lepton_pt ->at(i);
      float eta = std_vector_lepton_eta->at(i);
      float phi = std_vector_lepton_phi->at(i);
      float id  = std_vector_lepton_flavour ->at(i);


      Lepton lep;
      
      lep.index  = i;
      lep.charge = id;
      lep.pt = pt; 
      lep.eta = eta;
      
      if (fabs(id) == 11) { 
	lep.flavour = Electron; 
      } else if (fabs(id) == 13) { 
	lep.flavour = Muon; 
      }

      if (fabs(id) == 11 && pt <= 13) continue;
      if (fabs(id) == 13 && pt <= 10) continue;       
      if (fabs(id) == 11 && fabs(eta) >= 2.5 ) continue;  
      if (fabs(id) == 13 && fabs(eta) >= 2.4 ) continue;  

      //if ( IsTightLepton(i) &&  IsIsolatedLepton(i) ) { 
      if (std_vector_lepton_isTightLepton->at(i) == 1) {
	lep.type = Tight;
      } else {
	lep.type =  noTight;
      }

           
     //------------------------------------------------------------------------
     //
     // Electron
     //
     //------------------------------------------------------------------------
      if (fabs(id) == 11)
	{
	  lep.pr = GetFactor(ElecPR,       pt, eta);
	}
      //------------------------------------------------------------------------
      //
      // Muon 
      //
      //------------------------------------------------------------------------
      else if (fabs(id) == 13)
	{
	  lep.pr = GetFactor(MuonPR,        pt, eta);  
	}
  
      AnalysisLeptons.push_back(lep);
 
    } // end loop on leptons


    // Data-driven estimates 2 lepton case
    //----------------------------------------------------------------------------
    
  
    if ( AnalysisLeptons.size() == 2)  {

      weightPF    = float(GetPFWeight(int(mu20), int(ele35))); 
      weightFF    = float(GetFFWeight(int(mu20), int(ele35))); 
      
      weight2l0j = weightPF + weightFF;

      weight2l1j = float(GetPFWeight(int(mu25), int(ele35))) + float(GetFFWeight(int(mu25), int(ele35)));
      weight2l2j = float(GetPFWeight(int(mu35), int(ele35))) + float(GetFFWeight(int(mu35), int(ele35)));
     
      weight2l0jElUp = float(GetPFWeight(int(mu20), int(ele45))) + float(GetFFWeight(int(mu20), int(ele45)));
      weight2l1jElUp = float(GetPFWeight(int(mu25), int(ele45))) + float(GetFFWeight(int(mu25), int(ele45)));
      weight2l2jElUp = float(GetPFWeight(int(mu35), int(ele45))) + float(GetFFWeight(int(mu35), int(ele45)));
    
      weight2l0jElDown = float(GetPFWeight(int(mu20), int(ele25))) + float(GetFFWeight(int(mu20), int(ele25)));
      weight2l1jElDown = float(GetPFWeight(int(mu25), int(ele25))) + float(GetFFWeight(int(mu25), int(ele25)));
      weight2l2jElDown = float(GetPFWeight(int(mu35), int(ele25))) + float(GetFFWeight(int(mu35), int(ele25)));
    
      weight2l0jstatElUp = float(GetPFWeight(int(mu20), int(ele35), int(ElUp))) + float(GetFFWeight(int(mu20), int(ele35), int(ElUp)));
      weight2l1jstatElUp = float(GetPFWeight(int(mu25), int(ele35), int(ElUp))) + float(GetFFWeight(int(mu25), int(ele35), int(ElUp)));
      weight2l2jstatElUp = float(GetPFWeight(int(mu35), int(ele35), int(ElUp))) + float(GetFFWeight(int(mu35), int(ele35), int(ElUp)));
    
      weight2l0jstatElDown = float(GetPFWeight(int(mu20), int(ele35), int(ElDown))) + float(GetFFWeight(int(mu20), int(ele35), int(ElDown)));
      weight2l1jstatElDown = float(GetPFWeight(int(mu25), int(ele35), int(ElDown))) + float(GetFFWeight(int(mu25), int(ele35), int(ElDown)));
      weight2l2jstatElDown = float(GetPFWeight(int(mu35), int(ele35), int(ElDown))) + float(GetFFWeight(int(mu35), int(ele35), int(ElDown)));

      weight2l0jMuUp = float(GetPFWeight(int(mu30), int(ele35))) + float(GetFFWeight(int(mu30), int(ele35)));
      weight2l1jMuUp = float(GetPFWeight(int(mu35), int(ele35))) + float(GetFFWeight(int(mu35), int(ele35)));
      weight2l2jMuUp = float(GetPFWeight(int(mu45), int(ele35))) + float(GetFFWeight(int(mu45), int(ele35)));
      
      weight2l0jMuDown = float(GetPFWeight(int(mu10), int(ele35))) + float(GetFFWeight(int(mu10), int(ele35)));
      weight2l1jMuDown = float(GetPFWeight(int(mu15), int(ele35))) + float(GetFFWeight(int(mu15), int(ele35)));
      weight2l2jMuDown = float(GetPFWeight(int(mu25), int(ele35))) + float(GetFFWeight(int(mu25), int(ele35)));
   
   weight2l0jstatMuUp = float(GetPFWeight(int(mu20), int(ele35), int(MuUp))) + float(GetFFWeight(int(mu20), int(ele35), int(MuUp)));
      weight2l1jstatMuUp = float(GetPFWeight(int(mu25), int(ele35), int(MuUp))) + float(GetFFWeight(int(mu25), int(ele35), int(MuUp)));
      weight2l2jstatMuUp = float(GetPFWeight(int(mu35), int(ele35), int(MuUp))) + float(GetFFWeight(int(mu35), int(ele35), int(MuUp)));

      weight2l0jstatMuDown = float(GetPFWeight(int(mu20), int(ele35), int(MuDown))) + float(GetFFWeight(int(mu20), int(ele35), int(MuDown)));
      weight2l1jstatMuDown = float(GetPFWeight(int(mu25), int(ele35), int(MuDown))) + float(GetFFWeight(int(mu25), int(ele35), int(MuDown)));
      weight2l2jstatMuDown = float(GetPFWeight(int(mu35), int(ele35), int(MuDown))) + float(GetFFWeight(int(mu35), int(ele35), int(MuDown)));
      
    }
    

    int nevet = 0;
    
    if ( AnalysisLeptons.size()  == 3 )  {
      nevet++;

      weightPPF   =  float(GetPPFWeight(int(mu20), int(ele35))); 
      weightPFF   =  float(GetPFFWeight(int(mu20), int(ele35))); 
      weightFFF   =  float(GetFFFWeight(int(mu20), int(ele35))); 
      weight3l     = weightPPF + weightPFF + weightFFF;
      
      weight3lElUp =  float(GetPPFWeight(int(mu20), int(ele45))) + float(GetPFFWeight(int(mu20), int(ele45))) + float(GetFFFWeight(int(mu20), int(ele45))); 
    
      weight3lElDown =  float(GetPPFWeight(int(mu20), int(ele25))) + float(GetPFFWeight(int(mu20), int(ele25))) + float(GetFFFWeight(int(mu20), int(ele25))); 

      weight3lMuUp =  float(GetPPFWeight(int(mu30), int(ele35))) + float(GetPFFWeight(int(mu30), int(ele35))) + float(GetFFFWeight(int(mu30), int(ele35))); 

      weight3lMuDown =  float(GetPPFWeight(int(mu10), int(ele35))) + float(GetPFFWeight(int(mu10), int(ele35))) + float(GetFFFWeight(int(mu10), int(ele35))); 

      weight3lstatElUp =  float(GetPPFWeight(int(mu20), int(ele35), int(ElUp))) + float(GetPFFWeight(int(mu20), int(ele35), int(ElUp))) + float(GetFFFWeight(int(mu20), int(ele35), int(ElUp))); 
   
      weight3lstatElDown =  float(GetPPFWeight(int(mu20), int(ele35), int(ElDown))) + float(GetPFFWeight(int(mu20), int(ele35), int(ElDown))) + float(GetFFFWeight(int(mu20), int(ele35), int(ElDown))); 

     weight3lstatMuUp =  float(GetPPFWeight(int(mu20), int(ele35), int(MuUp))) + float(GetPFFWeight(int(mu20), int(ele35), int(MuUp))) + float(GetFFFWeight(int(mu20), int(ele35), int(MuUp))); 

      weight3lstatMuDown =  float(GetPPFWeight(int(mu20), int(ele35), int(MuDown))) + float(GetPFFWeight(int(mu20), int(ele35), int(MuDown))) + float(GetFFFWeight(int(mu20), int(ele35), int(MuDown)));

    }    

 
    
    bool isTight1 = 0, isTight2 = 0,   isTight3 = 0; 
    if ( AnalysisLeptons[0].type == Tight)  {
      isTight1 = 1;
    } else if (AnalysisLeptons[0].type == noTight)  {
      isTight1 = 0;
    }

    if ( AnalysisLeptons[1].type == Tight) {
      isTight2 = 1;
    }  else if (AnalysisLeptons[1].type == noTight)  {
      isTight2 = 0;
    }
    if ( AnalysisLeptons[2].type == Tight) {
      isTight3 = 1;
    }  else if (AnalysisLeptons[2].type == noTight)  {
      isTight3 = 0;
    }



    //------ summing the weights only for those events passing the full selection
    float pt1  = AnalysisLeptons[0].pt; float pt2  = AnalysisLeptons[1].pt; float pt3  = AnalysisLeptons[2].pt;
    float ch1  = AnalysisLeptons[0].charge;  float ch2  = AnalysisLeptons[1].charge;
    int flav1 =  AnalysisLeptons[0].flavour;  int flav2 =  AnalysisLeptons[1].flavour;int flav3 =  AnalysisLeptons[2].flavour;

  
    
    bool commonSel =  ( AnalysisLeptons.size() == 2 && trigger == 1 && pt1 > 20 && pt2 > 20 && ch1*ch2 < 0 && mll>12);
 
 
    if(commonSel) {

      //cout << weight << endl;
      countEntriesByHand2++;

      if(flav1 == Muon     && flav2 == Muon) {
	sumOfWeightsMuMu += weightPF; 
	sumOfQCDWeightsMuMu += weightFF; 
      }

      if(flav1 == Muon     && flav2 == Electron) {
	sumOfWeightsMuEl += weightPF;
	sumOfQCDWeightsMuEl += weightFF; 
      }

      if(flav1 == Electron && flav2 == Muon) {
	sumOfWeightsElMu += weightPF; 
	sumOfQCDWeightsElMu += weightFF; 
      }      

      if(flav1 == Electron && flav2 == Electron) {
	sumOfWeightsElEl += weightPF; 
	sumOfQCDWeightsElEl += weightFF; 
      }
      
      char* label;
      if( isTight1 && isTight2 ) label="PP";
      if( !isTight1 && !isTight2) label="FF";
      if(( isTight1 && !isTight2) || ( !isTight1 &&  isTight2)) label="PF";
	

      if(label == "PP"){
	if(flav1 == Muon     && flav2 == Muon) mmPP++;
	if(flav1 == Muon     && flav2 == Electron) mePP++;
	if(flav1 == Electron && flav2 == Muon) emPP++;
	if(flav1 == Electron && flav2 == Electron) eePP++;
      }
      if(label == "PF"){
	if(flav1 == Muon     && flav2 == Muon) mmPF++;
	if(flav1 == Muon     && flav2 == Electron) mePF++;
	if(flav1 == Electron && flav2 == Muon) emPF++;
	if(flav1 == Electron && flav2 == Electron) eePF++;
      }
      if(label == "FF"){
	if(flav1 == Muon     && flav2 == Muon) mmFF++;
	if(flav1 == Muon     && flav2 == Electron) meFF++;
	if(flav1 == Electron && flav2 == Muon) emFF++;
	if(flav1 == Electron && flav2 == Electron) eeFF++;
      }

    } // end selection

 

    // normalise to a unit lumi 
    baseW *= 1./luminosity;
    countEntriesByHand++;


    //if(addWeight) newtree->Fill(); 
    if ( AnalysisLeptons.size() > 1  )  {  
      if(addWeight && trigger>0) newtree->Fill(); //in the momento we will have trigger.
    }
  } // end loop on entries



  //------------------------------------------------------------------------------
  // PRINTOUT 
  //------------------------------------------------------------------------------

  cout << "# total: " << countEntriesByHand << endl;
  cout << "# Prompt-Fake: passing WW 0-jet: " << countEntriesByHand2 << endl;
  
  cout << "sumOfWeights mm: " << sumOfWeightsMuMu << endl;
  cout << "sumOfWeights me: " << sumOfWeightsMuEl << endl;
  cout << "sumOfWeights em: " << sumOfWeightsElMu << endl;
  cout << "sumOfWeights ee: " << sumOfWeightsElEl << endl;
  
  cout << "  " << endl;
  
  cout << "# Fake-Fake: passing WW 0-jet: " << countEntriesByHand2 << endl;
  
  cout << "sumOfQCDWeights mm: " << sumOfQCDWeightsMuMu << endl;
  cout << "sumOfQCDWeights me: " << sumOfQCDWeightsMuEl << endl;
  cout << "sumOfQCDWeights em: " << sumOfQCDWeightsElMu << endl;
  cout << "sumOfQCDWeights ee: " << sumOfQCDWeightsElEl << endl;

  
  cout << "mm PP,PF,FF: " << mmPP << " , " << mmPF << " , " << mmFF << endl;
  cout << "me PP,PF,FF: " << mePP << " , " << mePF << " , " << meFF << endl;
  cout << "em PP,PF,FF: " << emPP << " , " << emPF << " , " << emFF << endl;
  cout << "ee PP,PF,FF: " << eePP << " , " << eePF << " , " << eeFF << endl;

  //====================================

  if (addWeight) { 
    newfile->Write();
    delete newfile;
  }
  delete file;
  return 0;

} // end main



//------------------------------------------------------------------------------
// Helped method for weight computation
//------------------------------------------------------------------------------

		     
//------------------------------------------------------------------------------
// GetPPFWeight
//------------------------------------------------------------------------------
Double_t GetPPFWeight(int mu, int ele, int stat) {
  Double_t promptProbability[3];
  Double_t fakeProbability[3];

  int nTight= 0;

  for (UInt_t i=0; i<3; i++) {
    
    Lepton lep = AnalysisLeptons[i];

    Double_t f = 0.0;
    Double_t fV = 0.0, fE = 0.0;
    
    if ( lep.flavour == Muon) { 
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 35.0);
      if (stat == MuUp) f = fV+fE; 
      else if (stat == MuDown) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 35.0);
      if (stat == ElUp) f = fV+fE; 
      else if (stat == ElDown) f = fV-fE;
      else f = fV; 
    }

    Double_t p = lep.pr;

    if (lep.type == Tight)
      {
	nTight++;
	promptProbability[i] = p * (1 - f);
	fakeProbability[i]   = f * (1 - p);
      }
    else if (lep.type == noTight ) 
      {
	promptProbability[i] = p * f;
	fakeProbability[i]   = p * f;
      }

    promptProbability[i] /= (p - f);
    fakeProbability[i]   /= (p - f);
  }

  Double_t PPF = promptProbability[0] * promptProbability[1] * fakeProbability[2];
  Double_t PFP = promptProbability[0] * fakeProbability[1]   * promptProbability[2];
  Double_t FPP = fakeProbability[0]   * promptProbability[1] * promptProbability[2];

  Double_t result = PPF + PFP + FPP;

    if (nTight == 1 || nTight == 3) result *= -1.;

  return result;
}



//------------------------------------------------------------------------------
// GetPFFWeight
//------------------------------------------------------------------------------
		     Double_t GetPFFWeight(int mu, int ele, int stat)
{
  Double_t promptProbability[3];
  Double_t fakeProbability[3];
 
  int nTight= 0;
  
  for (UInt_t i=0; i<3; i++) {

    Lepton lep = AnalysisLeptons[i];

    Double_t f = 0.0;
    Double_t fV = 0.0, fE = 0.0;
    
    if ( lep.flavour == Muon) { 
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 35.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 35.);
      if (stat == MuUp) f = fV+fE; 
      else if (stat == MuDown) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 35.);
      if (stat == ElUp) f = fV+fE; 
      else if (stat == ElDown) f = fV-fE;
      else f = fV; 
    }

    Double_t p = lep.pr;

    if (lep.type == Tight)
      {
	nTight++;
	promptProbability[i] = p * (1 - f);
	fakeProbability[i]   = f * (1 - p);
      }
    else if (lep.type == noTight )
      {
	promptProbability[i] = p * f;
	fakeProbability[i]   = p * f;
      }

    promptProbability[i] /= (p - f);
    fakeProbability[i]   /= (p - f);
  }

  Double_t PFF = promptProbability[0] * fakeProbability[1]   * fakeProbability[2];
  Double_t FPF = fakeProbability[0]   * promptProbability[1] * fakeProbability[2];
  Double_t FFP = fakeProbability[0]   * fakeProbability[1]   * promptProbability[2];

  Double_t result = PFF + FPF + FFP;

  if (nTight == 0 || nTight == 2) result *= -1.;

  return result;
}


//------------------------------------------------------------------------------
// GetFFFWeight
//------------------------------------------------------------------------------
		     Double_t GetFFFWeight(int mu, int ele, int stat)
{
  Double_t fakeProbability[3];
  int nTight= 0;
  
  for (UInt_t i=0; i<3; i++) {
    
    Lepton lep = AnalysisLeptons[i];
 
    Double_t f = 0.0;
    Double_t fV = 0.0, fE = 0.0;
    
    if ( lep.flavour == Muon) { 
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 35.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 35.);
      if (stat == MuUp) f = fV+fE; 
      else if (stat == MuDown) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 35.);
      if (stat == ElUp) f = fV+fE; 
      else if (stat== ElDown) f = fV-fE;
      else f = fV; 
    }

    Double_t p = lep.pr;

    if (lep.type == Tight)
      {
	nTight++;
	fakeProbability[i] = f * (1 - p);
      }
    else if (lep.type == noTight)
      {
	fakeProbability[i] = p * f;
      }

    fakeProbability[i] /= (p - f);
  }

  Double_t FFF = fakeProbability[0] * fakeProbability[1] * fakeProbability[2];

  if (nTight == 1 || nTight == 3) FFF *= -1.;

  return FFF;
}


//------------------------------------------------------------------------------
// GetPFWeight
//------------------------------------------------------------------------------
Double_t GetPFWeight(int mu, int ele, int stat)
{
  Double_t promptProbability[2];
  Double_t fakeProbability[2];

  int nTight = 0; 
 
  for (UInt_t i=0; i<2; i++) {
    
    Lepton lep = AnalysisLeptons[i];

    Double_t f = 0.0;
    Double_t fV = 0.0, fE = 0.0;
 
    if ( lep.flavour == Muon) { 
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 35.0);     

      if (stat == MuUp) f = fV+fE; 
      else if (stat == MuDown) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 35.0);


      if (stat == ElUp) f = fV+fE; 
      else if (stat == ElDown) f = fV-fE;
      else f = fV; 
    }

    Double_t p = lep.pr;

    if (lep.type == Tight)
      {
	nTight++;

	promptProbability[i] = p * (1 - f);
	fakeProbability[i]   = f * (1 - p);

	//cout << lep.flavour << "  "  << lep.pt << "  " << lep.eta << "  " << p << "  " << f << "  " << promptProbability[i] << "  " << fakeProbability[i]  << endl;

      }
    else if (lep.type == noTight ) 
      {
	promptProbability[i] = p * f;
	fakeProbability[i]   = p * f;
      }

    promptProbability[i] /= (p - f);
    fakeProbability[i]   /= (p - f);
  }

  Double_t PF = promptProbability[0] * fakeProbability[1];
  Double_t FP = promptProbability[1] * fakeProbability[0];

  Double_t result = PF + FP;

  
  if (nTight == 0 || nTight == 2 ) result *= -1.;

  return result;
}


//------------------------------------------------------------------------------
// GetFFWeight
//------------------------------------------------------------------------------
Double_t GetFFWeight(int mu, int ele, int stat)
{
  Double_t fakeProbability[2];

  int nTight= 0;

  for (UInt_t i=0; i<2; i++) {
    
    Lepton lep = AnalysisLeptons[i];

    Double_t f = 0.0;
    Double_t fV = 0.0, fE = 0.0;
    
    if ( lep.flavour == Muon) { 
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 35.0);
      if (stat == MuUp) f = fV+fE; 
      else if (stat == MuDown) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 35.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 35.);
      if (stat == ElUp) f = fV+fE; 
      else if (stat == ElDown) f = fV-fE;
      else f = fV; 
    }

    Double_t p = lep.pr;

    if (lep.type == Tight)
      {
	nTight++;
	fakeProbability[i] = f * (1 - p);
      }
    else if (lep.type == noTight)
      {
	fakeProbability[i] = p * f;
      }

    fakeProbability[i] /= (p - f);
  }

  Double_t FF = fakeProbability[0] * fakeProbability[1];

  if (nTight == 1) FF *= -1.;

  return FF;

}


//------------------------------------------------------------------------------
// LoadHistogram
//------------------------------------------------------------------------------


TH2F* LoadHistogram(TString filename,
		    TString hname,
		    TString cname) {

  TString path = "results/";

  TFile* inputfile = TFile::Open(path + filename + ".root");

  TH2F* hist = (TH2F*)inputfile->Get(hname)->Clone(cname);
  
  hist->SetDirectory(0);
  
  inputfile->Close();

  return hist;
}



//------------------------------------------------------------------------------
// IsTightLepton
//
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
// egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium
//------------------------------------------------------------------------------
bool IsTightLepton(int k)
{
 
bool is_tight_lepton = false;

  // Muon tight ID
  if (fabs(std_vector_lepton_flavour->at(k)) == 13)
    {
      float dxyCut = 0;
	
      if ( std_vector_lepton_pt->at(k) < 20 ) { 
	dxyCut = 0.01;
      }	else {
	dxyCut = 0.02;
      }

      is_tight_lepton = ( std_vector_lepton_isMediumMuon->at(k)              && 
			  fabs(std_vector_lepton_dz->at(k))  < 0.1       && 
			  fabs(std_vector_lepton_d0->at(k)) < dxyCut );
	}

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_flavour->at(k)) == 11)
    {
      is_tight_lepton = std_vector_lepton_eleIdTight->at(k);
    }
  
  return is_tight_lepton;


}

//------------------------------------------------------------------------------
// MuonIsolation
//------------------------------------------------------------------------------
float MuonIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_flavour->at(k);

  float relative_isolation = -999;

  if (fabs(id) != 13) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso->at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso->at(k) +
	      std_vector_lepton_neutralHadronIso->at(k) -
	      0.5*std_vector_lepton_sumPUPt->at(k)));

  relative_isolation /= pt;

 
  return relative_isolation;
}


//------------------------------------------------------------------------------
// ElectronIsolation
//------------------------------------------------------------------------------
float ElectronIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_flavour->at(k);

  float relative_isolation = -999;

  if (fabs(id) != 11) return relative_isolation;

  relative_isolation =
    std_vector_lepton_chargedHadronIso->at(k) +
    max(float(0.0),
	float(std_vector_lepton_photonIso->at(k) +
	      std_vector_lepton_neutralHadronIso->at(k) -
	      jetRho*std_vector_electron_effectiveArea->at(k)));
  
  relative_isolation /= pt;
  
  return relative_isolation;
}


//------------------------------------------------------------------------------
// IsIsolatedLepton
//------------------------------------------------------------------------------
bool IsIsolatedLepton(int k)
{
  float id = std_vector_lepton_flavour->at(k);
  float pt = std_vector_lepton_pt->at(k);
 
  bool is_isolated_lepton = false;

  if (fabs(id) == 11) {
    if ( std_vector_lepton_trackIso->at(k)/pt < 0.2            &&
         std_vector_electron_ecalPFClusterIso->at(k)/pt < 0.45 &&
         std_vector_electron_hcalPFClusterIso->at(k)/pt < 0.25 )  
         is_isolated_lepton = true; 
  }
  else if (fabs(id) == 13) {
    if ( MuonIsolation(k)   < 0.15  &&   
         std_vector_lepton_trackIso->at(k)/pt < 0.4 ) 
      is_isolated_lepton = true;
  }
  
  return is_isolated_lepton;

}



//------------------------------------------------------------------------------
// IsLooseLepton
//
//------------------------------------------------------------------------------
bool IsLooseLepton(int k)
{
 bool is_loose_lepton = false;
  
  // Muon tight ID
  if (fabs(std_vector_lepton_flavour->at(k)) == 13)
    {
      float dxyCut = 0;
      
      if ( std_vector_lepton_pt->at(k) < 20 ) { 
	dxyCut = 0.01;
      }	else {
	dxyCut = 0.02;
      }
     
      is_loose_lepton = ( std_vector_lepton_isMediumMuon->at(k)                   && 
			  fabs(std_vector_lepton_dz->at(k))  < 0.1      && 
			  fabs(std_vector_lepton_d0->at(k)) < dxyCut); 
    }
  

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_flavour->at(k)) == 11)
    {
      // Faltan los aislamientos!!!!!!

      float aeta = fabs(std_vector_lepton_eta->at(k));
      
      if (aeta <= 1.479) {
	
	if (fabs(std_vector_electron_dEtaIn->at(k)) < 0.01               &&
	    fabs(std_vector_electron_dPhiIn->at(k)) < 0.04               &&
	    std_vector_electron_full5x5_sigmaIetaIeta->at(k)    < 0.011              &&
	    std_vector_electron_hOverE->at(k)              < 0.08               &&
	    std_vector_electron_ooEmooP->at(k) < 0.01               && //????
	    std_vector_electron_expectedMissingInnerHits->at(k)<=2  &&
	    std_vector_lepton_d0->at(k)      < 0.1                &&
	    fabs(std_vector_lepton_dz->at(k))< 0.373              &&
	    std_vector_electron_passConversionVeto->at(k) )
	  {  
	    is_loose_lepton = true; 
	  }
      }
      else if (aeta > 1.479 && aeta < 2.5) { 

	if (fabs(std_vector_electron_dEtaIn->at(k)) < 0.01               &&
	    fabs(std_vector_electron_dPhiIn->at(k)) < 0.08               &&
	    std_vector_electron_full5x5_sigmaIetaIeta->at(k)    < 0.031              &&
	    std_vector_electron_hOverE->at(k)              < 0.08               &&
	    std_vector_electron_ooEmooP->at(k) < 0.01               &&
	    std_vector_electron_expectedMissingInnerHits->at(k)<=1  &&
	    std_vector_lepton_d0->at(k)      < 0.2                &&
	    fabs(std_vector_lepton_dz->at(k))< 0.602              &&
	    std_vector_electron_passConversionVeto->at(k) )
	  {  
	    is_loose_lepton = true; 
	  }
      }

    }
  
  return is_loose_lepton;
}



//------------------------------------------------------------------------------
// IsLooseIsolatedLepton
//------------------------------------------------------------------------------
bool IsLooseIsolatedLepton(int k)
{
  float id = std_vector_lepton_flavour->at(k);
  float pt = std_vector_lepton_pt->at(k);
 
  bool is_isolated_lepton = false;

  if (fabs(id) == 11) {
    if ( std_vector_lepton_trackIso->at(k)/pt < 0.2            &&
         std_vector_electron_ecalPFClusterIso->at(k)/pt < 0.45 &&
         std_vector_electron_hcalPFClusterIso->at(k)/pt < 0.25 )  
         is_isolated_lepton = true; 
  }
  else if (fabs(id) == 13) {
    if ( MuonIsolation(k)   < 0.4  &&   
         std_vector_lepton_trackIso->at(k)/pt < 0.4 ) 
      is_isolated_lepton = true;
  }
  
  return is_isolated_lepton;

}


//------------------------------------------------------------------------------
// GetFactor
//------------------------------------------------------------------------------
Double_t GetFactor(TH2F*   h2,
		  Float_t leptonPt,
		  Float_t leptonEta,
		  Float_t leptonPtMax)
{
  Float_t aeta = fabs(leptonEta);

  Int_t nbins = h2->GetNbinsX();

  Float_t ptMax = (leptonPtMax > 0) ? leptonPtMax : h2->GetXaxis()->GetBinCenter(nbins);

  Float_t factor = h2->GetBinContent(h2->FindBin(TMath::Min(leptonPt, ptMax), aeta));

  return factor;
}

//------------------------------------------------------------------------------
// GetFactor
//------------------------------------------------------------------------------
Double_t GetFactorError(TH2F*   h2,
		       Float_t leptonPt,
		       Float_t leptonEta,
		       Float_t leptonPtMax)
{
  Float_t aeta = fabs(leptonEta);

  Int_t nbins = h2->GetNbinsX();

  Float_t ptMax = (leptonPtMax > 0) ? leptonPtMax : h2->GetXaxis()->GetBinCenter(nbins);

  Float_t factor = h2->GetBinError(h2->FindBin(TMath::Min(leptonPt, ptMax), aeta));

  return factor;
}
