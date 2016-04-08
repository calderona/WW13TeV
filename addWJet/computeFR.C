#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
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
#include "TChain.h"
#include "TLorentzVector.h"

#include <sstream>
#include <math.h>
#include <iomanip>
#include <iostream>

#include <algorithm>
using namespace std;


//------------------------------------------------------------------------------
// Member functions
//------------------------------------------------------------------------------

bool  IsIsolatedLepton(int k);
float ElectronIsolation(int k);
float MuonIsolation(int k);
bool  IsTightLepton(int k);
bool  IsLooseLepton(int k);
bool  IsLooseIsolatedLepton(int k);
Double_t giveMT(int k, double MET, double MET_Phi);


//------------------------------------------------------------------------------
// Define structures & GLB variables 
//------------------------------------------------------------------------------

enum {Muon, Electron, Unknown};

enum {Tight, noTight};

enum ch {mm=0,ee=1,em=2,me=3};

struct Lepton { 

  UInt_t         index;
  UInt_t         flavor;  // Muon, Electron
  UInt_t         type;    // Tight, Loose
  Float_t        charge;
  Bool_t         fromW;
  Bool_t         iGen;
  TLorentzVector v;

};


struct Jet { 

  UInt_t         index;
  TLorentzVector v;

};

std::vector<Lepton> AnalysisLeptons;
std::vector<Lepton> AnalysisGenLeptons;                                                                   

std::vector<Jet> AnalysisJets;
std::vector<Jet> AnalysisLeptonJets;


Float_t          run;
Float_t          evt;
Float_t          jetRho;
Float_t          puW;
Float_t          channel;
Float_t          baseW;
Float_t          metPfType1;  
Float_t          metPfType1Phi;
Float_t          nvtx;
Float_t          triggerFakeRate;
Float_t          trigger;
Float_t          njet;

vector<float>   *std_vector_lepton_trackIso;
vector<float>   *std_vector_electron_ecalPFClusterIso;
vector<float>   *std_vector_electron_hcalPFClusterIso;

//vector<float>   *std_vector_lepton_isTightLepton;
//vector<float>   *std_vector_lepton_isLooseLepton;
vector<float>   *std_vector_lepton_flavour;
vector<float>   *std_vector_lepton_pt;
vector<float>   *std_vector_lepton_eta;
vector<float>   *std_vector_lepton_phi;
vector<float>   *std_vector_lepton_isTightMuon; 
vector<float>   *std_vector_lepton_isMediumMuon; 
vector<float>   *std_vector_lepton_eleIdMedium; 
vector<float>   *std_vector_lepton_eleIdTight; 
vector<float>   *std_vector_lepton_eleIdVeto;
vector<float>   *std_vector_lepton_chargedHadronIso; 
vector<float>   *std_vector_lepton_photonIso;
vector<float>   *std_vector_lepton_neutralHadronIso;
vector<float>   *std_vector_lepton_sumPUPt;
vector<float>   *std_vector_electron_effectiveArea;
vector<float>   *std_vector_lepton_dz;
vector<float>   *std_vector_lepton_d0;
vector<float>   *std_vector_lepton_SIP3D;


vector<float>   *std_vector_electron_dEtaIn;
vector<float>   *std_vector_electron_dPhiIn;
vector<float>   *std_vector_electron_full5x5_sigmaIetaIeta;
vector<float>   *std_vector_electron_hOverE;
vector<float>   *std_vector_electron_ooEmooP;
vector<float>   *std_vector_electron_expectedMissingInnerHits;
vector<float>   *std_vector_electron_passConversionVeto;


vector<float>   *std_vector_jet_cmvav2;

vector<float>   *std_vector_trigger;
vector<float>   *std_vector_trigger_prescale;

vector<float>   *std_vector_jet_pt;
vector<float>   *std_vector_jet_eta;
vector<float>   *std_vector_jet_phi;
vector<float>   *std_vector_lepton_closejet_pt;
vector<float>   *std_vector_lepton_closejet_eta;
vector<float>   *std_vector_lepton_closejet_phi;
vector<float>   *std_vector_lepton_closejet_PartonFlavour;
vector<float>   *std_vector_lepton_closejet_drlj;
vector<float>   *std_vector_jet_csvv2ivf;

Float_t  GEN_weight_SM = 1.0;
Float_t  negativeWeight = 1.0;

Float_t myBaseW = 1.0;

//------------------------------------------------------------------------------
// MAIN function
//------------------------------------------------------------------------------


void computeFR (TString theSample="WJets", bool doGen = false, int flavour = 0, float inputJetEt = 0) {

  // ------ root settings ---------
 
  TH1::SetDefaultSumw2();
  
  TString path = Form("/gpfs/csic_users/calderon/IFCA_2015/WW13TeV_V1/addWJet/rootfiles/");

  TFile* root_output;  
  TString jetFlavour; 

  if (flavour == 1 )   jetFlavour = "_D_"; 
  if (flavour == 2 )   jetFlavour = "_U_"; 
  if (flavour == 3 )   jetFlavour = "_S_"; 
  if (flavour == 4 )   jetFlavour = "_C_"; 
  if (flavour == 5 )   jetFlavour = "_B_"; 
  if (flavour == 21)   jetFlavour = "_G_"; 
  if (flavour == 0 )   jetFlavour = "_"; 

  TString jet = Form("%i", (int)inputJetEt);
  
  if (doGen) {
    root_output = new TFile(path + theSample + "_FR_FakeRate_226pb_76X_Normalize.root", "recreate");
  }  else {
    root_output =  new TFile(path + theSample + jetFlavour + "jet"+ jet +"_FakeRate_226pb_76X_NoBjets.root", "recreate");
  }

  //-- Define histograms 
  //----------------------------------------------------------------------------

  // Define pt bins for FR matrix
  int pTbin = 8;
  Double_t pTbins[9] = { 10., 15., 20., 25., 30., 35.,40., 45., 50.};

  // Define eta bins for FR matrix
  int etabin = 5;
  Double_t etabins[6] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

  //Define pt bins for FR matrix
  int dxypTbin = 7;
  Double_t dxypTbins[8] = {0,0.002,0.004,0.006,0.008,0.012,0.016,0.020};


  // MUON histograms

  TH1F  *h_triggerBits = new TH1F("h_triggerBits", "h_triggerBits", 80, 0, 80 ); 
  TH1F * h_LooseMuons = new TH1F("h_LooseMuons", "h_LooseMuons", 10, 0.1, 1.1);
  TH1F  *h_LooseMuons_Iso   = new TH1F("h_LooseMuons_Iso","h_LooseMuons_Iso", 10,0,1.0);
  TH1F  *h_LooseMuons_All  = new TH1F("h_LooseMuons_All","h_LooseMuons_All", 3,0,2);  

  TH1F  *h_Muon_fake_pT_bin   = new TH1F("h_Muon_fake_pT_bin","h_Muon_fake_pT_bin", pTbin,pTbins);
  TH1F  *h_Muon_signal_pT_bin = new TH1F("h_Muon_signal_pT_bin","h_Muon_signal_pT_bin",pTbin,pTbins);
  TH1F  *h_Muon_fake_eta_bin = new TH1F("h_Muon_fake_eta_bin","h_Muon_fake_eta_bin",etabin,etabins);
  TH1F  *h_Muon_signal_eta_bin = new TH1F("h_Muon_signal_eta_bin","h_Muon_signal_eta_bin",etabin,etabins);

  TH2F  *h_Muon_fake_pt_eta_bin = new TH2F("h_Muon_fake_pt_eta_bin","h_Muon_fake_pt_eta_bin",pTbin,pTbins,etabin,etabins);
  TH2F  *h_Muon_signal_pt_eta_bin = new TH2F("h_Muon_signal_pt_eta_bin","h_Muon_signal_pt_eta_bin",pTbin,pTbins,etabin,etabins);

  TH2F  *h_Muon_fake_ptIso_eta_bin = new TH2F("h_Muon_fake_ptIso_eta_bin","h_Muon_fake_ptIso_eta_bin",pTbin,pTbins,etabin,etabins);
  TH2F  *h_Muon_signal_ptIso_eta_bin = new TH2F("h_Muon_signal_ptIso_eta_bin","h_Muon_signal_ptIso_eta_bin",pTbin,pTbins,etabin,etabins);

  TH1F  *h_Muon_fake_nPV = new TH1F("h_Muon_fake_nPV","h_Muon_fake_nPV",10,0,30);
  TH1F  *h_Muon_signal_nPV = new TH1F("h_Muon_signal_nPV","h_Muon_signal_nPV",10,0,30);

  TH1F  *h_Muon_fake_iso   = new TH1F("h_Muon_fake_iso","h_Muon_fake_iso", 100,0,0.5);
  TH1F  *h_Muon_signal_iso   = new TH1F("h_Muon_signal_iso","h_Muon_signal_iso", 100,0,0.5);


  TH1F  *h_Muon_fake_ptRel   = new TH1F("h_Muon_fake_ptRel","h_Muon_fake_ptRel", 100,0,100);
  TH1F  *h_Muon_signal_ptRel   = new TH1F("h_Muon_signal_ptRel","h_Muon_signal_ptRel", 100,0,100);
 
 TH1F  *h_Muon_fake_ChIso   = new TH1F("h_Muon_fake_ChIso","h_Muon_fake_ChIso", 100,0,1);
  TH1F  *h_Muon_fake_NeuIso   = new TH1F("h_Muon_fake_NeuIso","h_Muon_fake_NeuIso", 100,0,0.05);
  TH1F  *h_Muon_fake_PhIso   = new TH1F("h_Muon_fake_PhIso","h_Muon_fake_PhIso", 100,0,0.05);

  TH1F  *h_Muon_signal_ChIso   = new TH1F("h_Muon_signal_ChIso","h_Muon_signal_ChIso", 100,0,1);
  TH1F  *h_Muon_signal_NeuIso   = new TH1F("h_Muon_signal_NeuIso","h_Muon_signal_NeuIso", 100,0,0.05);
  TH1F  *h_Muon_signal_PhIso   = new TH1F("h_Muon_signal_PhIso","h_Muon_signal_PhIso", 100,0,0.05);

  TH1F  *h_Muon_fake_btag   = new TH1F("h_Muon_fake_btag","h_Muon_fake_btag", 10,-1,1);
  TH1F  *h_Muon_signal_btag   = new TH1F("h_Muon_signal_btag","h_Muon_signal_btag", 10,-1,1);

  TH1F  *h_Muon_Mt1   = new TH1F("h_Muon_Mt1","h_Muon_Mt1", 100,0,200);
  TH1F  *h_Muon_Mt2  = new TH1F("h_Muon_Mt2","h_Muon_Mt2", 100,0,200);
  TH1F  *h_Muon_pt   = new TH1F("h_Muon_pt","h_Muon_pt", 100,15,200);

  TH1F  *h_Muon_fake_Mt   = new TH1F("h_Muon_fake_Mt","h_Muon_fake_Mt", 100,0,200);
  TH1F  *h_Muon_fake_Mt_Met20   = new TH1F("h_Muon_fake_Mt_Met20","h_Muon_fake_Mt_Met20", 100,0,200);
  TH1F  *h_Muon_fake_Met   = new TH1F("h_Muon_fake_Met","h_Muon_fake_Met", 100,0,200);
  TH1F  *h_Muon_fake_CRjetEt   = new TH1F("h_Muon_fake_CRjetEt","h_Muon_fake_CRjetEt", 100,0,200);
  TH1F  *h_Muon_fake_pt   = new TH1F("h_Muon_fake_pt","h_Muon_fake_pt", 100,15,200);
  TH1F  *h_Muon_jetEt   = new TH1F("h_Muon_jetEt","h_Muon_jetEt", 100,0,200);
  TH1F  *h_Muon_ptIso   = new TH1F("h_Muon_ptIso","h_Muon_ptIso", 100,0,200);
  TH2F  *h_Muon_jetEt_ptIso   = new TH2F("h_Muon_jetEt_ptIso","h_Muon_jetEt_ptIso",100,15,100,100,15,100); 
  TH2F  *h_Muon_jetEt_pt   = new TH2F("h_Muon_jetEt_pt","h_Muon_jetEt_pt",100,15,100,100,15,100); 
  TH2F  *h_Muon_jetEt_iso   = new TH2F("h_Muon_jetEt_iso","h_Muon_jetEt_iso",100,15,100,100,0,0.5); 

  TH1F  *h_Muon_signal_Mt   = new TH1F("h_Muon_signal_Mt","h_Muon_signal_Mt", 100,0,200);
  TH1F  *h_Muon_signal_Mt_Met20   = new TH1F("h_Muon_signal_Mt_Met20","h_Muon_signal_Mt_Met20", 100,0,200);
  TH1F  *h_Muon_signal_Met   = new TH1F("h_Muon_signal_Met","h_Muon_signal_Met", 100,0,200);
  TH1F  *h_Muon_signal_CRjetEt   = new TH1F("h_Muon_signal_CRjetEt","h_Muon_signal_CRjetEt", 100,0,200);
  TH1F  *h_Muon_signal_pt   = new TH1F("h_Muon_signal_pt","h_Muon_signal_pt", 100,15,200);
  
  TH1F  *h_Muon_fake_dxy_pt = new TH1F("h_Muon_fake_dxy_pt","h_Muon_fake_dxy_pt",50, 0, 0.3);
  TH1F  *h_Muon_signal_dxy_pt = new TH1F("h_Muon_signal_dxy_pt","h_Muon_fake_dxy_pt",50, 0, 0.3);

  TH1F  *h_Muon_fake_dxy = new TH1F("h_Muon_fake_dxy","h_Muon_fake_dxy",50, 0, 0.3);
  TH1F  *h_Muon_signal_dxy = new TH1F("h_Muon_signal_dxy","h_Muon_fake_dxy",50, 0, 0.3);

  TH1F  *h_Muon_fake_dxy_pt_bin = new TH1F("h_Muon_fake_dxy_pt_bin","h_Muon_fake_dxy_pt_bin",dxypTbin,dxypTbins);
  TH1F  *h_Muon_signal_dxy_pt_bin = new TH1F("h_Muon_signal_dxy_pt_bin","h_Muon_fake_dxy_pt_bin",dxypTbin,dxypTbins);

  TH1F  *h_Muon_fake_jetEt = new TH1F("h_Muon_fake_jetEt","h_Muon_fake_jetEt",20,10,150);
  TH1F  *h_Muon_signal_jetEt = new TH1F("h_Muon_signal_jetEt","h_Muon_signal_jetEt",20,10,150);

  TH1F  *h_Muon_fake_sip3d = new TH1F("h_Muon_fake_sip3d","h_Muon_fake_sip3d",10, 0, 40);
  TH1F  *h_Muon_signal_sip3d = new TH1F("h_Muon_signal_sip3d","h_Muon_fake_sip3d",10, 0, 40);

  TH1F  *h_Muon_fake_ptIso = new TH1F("h_Muon_fake_ptIso","h_Muon_fake_ptIso",pTbin,pTbins);
  TH1F  *h_Muon_signal_ptIso = new TH1F("h_Muon_signal_ptIso","h_Muon_fake_ptIso",pTbin,pTbins);

  TH1F  *h_Muon_fake_jetEta = new TH1F("h_Muon_fake_jetEta","h_Muon_fake_jetEta",50,-5,5);
  TH1F  *h_Muon_signal_jetEta = new TH1F("h_Muon_signal_jetEta","h_Muon_signal_jetEta",50,-5,5);

  TH1F  *h_Muon_fake_mll = new TH1F("h_Muon_fake_mll","h_Muon_fake_mll",50,60,120);
  TH1F  *h_Muon_signal_mll = new TH1F("h_Muon_signal_mll","h_Muon_signal_mll",50,60,120);


  // ELECTRON histograms

  TH1F  *h_Ele_fake_pT_bin   = new TH1F("h_Ele_fake_pT_bin","h_Ele_fake_pT_bin", pTbin,pTbins);
  TH1F  *h_Ele_signal_pT_bin = new TH1F("h_Ele_signal_pT_bin","h_Ele_signal_pT_bin",pTbin,pTbins);
  TH1F  *h_Ele_fake_eta_bin = new TH1F("h_Ele_fake_eta_bin","h_Ele_fake_eta_bin",etabin,etabins);
  TH1F  *h_Ele_signal_eta_bin = new TH1F("h_Ele_signal_eta_bin","h_Ele_signal_eta_bin",etabin,etabins);

  TH2F  *h_Ele_fake_pt_eta_bin = new TH2F("h_Ele_fake_pt_eta_bin","h_Ele_fake_pt_eta_bin",pTbin,pTbins,etabin,etabins);
  TH2F  *h_Ele_signal_pt_eta_bin = new TH2F("h_Ele_signal_pt_eta_bin","h_Ele_signal_pt_eta_bin",pTbin,pTbins,etabin,etabins);

  TH2F  *h_Ele_fake_ptIso_eta_bin = new TH2F("h_Ele_fake_ptIso_eta_bin","h_Ele_fake_ptIso_eta_bin",pTbin,pTbins,etabin,etabins);
  TH2F  *h_Ele_signal_ptIso_eta_bin = new TH2F("h_Ele_signal_ptIso_eta_bin","h_Ele_signal_ptIso_eta_bin",pTbin,pTbins,etabin,etabins);

  TH1F  *h_Ele_fake_nPV = new TH1F("h_Ele_fake_nPV","h_Ele_fake_nPV",10,0,30);
  TH1F  *h_Ele_signal_nPV = new TH1F("h_Ele_signal_nPV","h_Ele_signal_nPV",10,0,30);

  TH1F  *h_Ele_fake_iso   = new TH1F("h_Ele_fake_iso","h_Ele_fake_iso", 100,0,0.5);
  TH1F  *h_Ele_signal_iso   = new TH1F("h_Ele_signal_iso","h_Ele_signal_iso", 100,0,0.5);
  

  TH1F  *h_Ele_fake_ChIso   = new TH1F("h_Ele_fake_ChIso","h_Ele_fake_ChIso", 100,0,0.5);
  TH1F  *h_Ele_fake_NeuIso   = new TH1F("h_Ele_fake_NeuIso","h_Ele_fake_NeuIso", 100,0,0.5);
  TH1F  *h_Ele_fake_PhIso   = new TH1F("h_Ele_fake_PhIso","h_Ele_fake_PhIso", 100,0,0.5);

  TH1F  *h_Ele_signal_ChIso   = new TH1F("h_Ele_signal_ChIso","h_Ele_signal_ChIso", 100,0,0.5);
  TH1F  *h_Ele_signal_NeuIso   = new TH1F("h_Ele_signal_NeuIso","h_Ele_signal_NeuIso", 100,0,0.5);
  TH1F  *h_Ele_signal_PhIso   = new TH1F("h_Ele_signal_PhIso","h_Ele_signal_PhIso", 100,0,0.5);

  TH1F  *h_Ele_Mt1   = new TH1F("h_Ele_Mt1","h_Ele_Mt1", 100,0,200);
  TH1F  *h_Ele_Mt2  = new TH1F("h_Ele_Mt2","h_Ele_Mt2", 100,0,200);

  TH1F  *h_Ele_fake_Mt   = new TH1F("h_Ele_fake_Mt","h_Ele_fake_Mt", 100,0,200);
  TH1F  *h_Ele_fake_Met   = new TH1F("h_Ele_fake_Met","h_Ele_fake_Met", 100,0,200);
  TH1F  *h_Ele_fake_CRjetEt   = new TH1F("h_Ele_fake_CRjetEt","h_Ele_fake_CRjetEt", 100,0,200);
  TH1F  *h_Ele_fake_pt   = new TH1F("h_Ele_fake_pt","h_Ele_fake_pt", 100,15,200);
  TH1F  *h_Ele_jetEt   = new TH1F("h_Ele_jetEt","h_Ele_jetEt", 100,0,200);
  TH1F  *h_Ele_fake_Mt_Met20   = new TH1F("h_Ele_fake_Mt_Met20","h_Ele_fake_Mt_Met20", 100,0,200);
  TH1F  *h_Ele_pt   = new TH1F("h_Ele_pt","h_Ele_pt", 100,15,200);
  TH1F  *h_Ele_ptIso   = new TH1F("h_Ele_ptIso","h_Ele_ptIso", 100,0,200);
  TH2F  *h_Ele_jetEt_iso   = new TH2F("h_Ele_jetEt_iso","h_Ele_jetEt_iso",100,15,100,100,15,100); 


  TH1F  *h_Ele_fake_ptRel   = new TH1F("h_Ele_fake_ptRel","h_Ele_fake_ptRel", 100,0,0.5);

  TH1F  *h_Ele_signal_Mt   = new TH1F("h_Ele_signal_Mt","h_Ele_signal_Mt", 100,0,200);
  TH1F  *h_Ele_signal_Met   = new TH1F("h_Ele_signal_Met","h_Ele_signal_Met", 100,0,200);
  TH1F  *h_Ele_signal_CRjetEt   = new TH1F("h_Ele_signal_CRjetEt","h_Ele_signal_CRjetEt", 100,0,200);
  TH1F  *h_Ele_signal_pt   = new TH1F("h_Ele_signal_pt","h_Ele_signal_pt", 100,15,200);
  TH1F  *h_Ele_signal_Mt_Met20   = new TH1F("h_Ele_signal_Mt_Met20","h_Ele_signal_Mt_Met20", 100,0,200);
  
  TH1F  *h_Ele_fake_jetEt = new TH1F("h_Ele_fake_jetEt","h_Ele_fake_jetEt",20,10,150);
  TH1F  *h_Ele_signal_jetEt = new TH1F("h_Ele_signal_jetEt","h_Ele_signal_jetEt",20,10,150);

  TH1F  *h_Ele_fake_jetEta = new TH1F("h_Ele_fake_jetEta","h_Ele_fake_jetEta",50,-5,5);
  TH1F  *h_Ele_signal_jetEta = new TH1F("h_Ele_signal_jetEta","h_Ele_signal_jetEta",50,-5,5);

  TH1F  *h_Ele_fake_sip3d = new TH1F("h_Ele_fake_sip3d","h_Ele_fake_sip3d",10, 0, 40);
  TH1F  *h_Ele_signal_sip3d = new TH1F("h_Ele_signal_sip3d","h_Ele_fake_sip3d",10, 0, 40);

  TH1F  *h_Ele_fake_ptIso = new TH1F("h_Ele_fake_ptIso","h_Ele_fake_ptIso",pTbin,pTbins);
  TH1F  *h_Ele_signal_ptIso = new TH1F("h_Ele_signal_ptIso","h_Ele_fake_ptIso",pTbin,pTbins);

  TH1F  *h_Ele_fake_mll = new TH1F("h_Ele_fake_mll","h_Ele_fake_mll",50,60,120);
  TH1F  *h_Ele_signal_mll = new TH1F("h_Ele_signal_mll","h_Ele_signal_mll",50,60,120);

  //--- OPEN LOOSE-LOOSE FILE 
  //----------------------------------------------------------------------------
  
  TString filesPath1, filesPath2, filesPath3, filesPath4;
  filesPath1 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/Data13TeV/22Jan/";
  filesPath2 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/22Jan_25ns_mAODv2_MC/MCl1loose__EpTCorr/";//MCl1loose_EpTCorr/";
  //filesPath1 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/";
  
  TChain* tree = new TChain("latino", "latino");

  if (theSample == "WJets") {
    tree->Add(filesPath2 + "latino_WJetsToLNu_0000__part0.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0000__part1.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0000__part2.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0000__part3.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0000__part4.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0001__part0.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0001__part1.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0001__part2.root");
    tree->Add(filesPath2 + "latino_WJetsToLNu_0001__part3.root");

  }
  else if (theSample == "QCD2") {
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part1.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part2.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part3.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part4.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part5.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part6.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part7.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part8.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part9.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part10.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part11.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part12.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part13.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part14.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part15.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part16.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part17.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part18.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part19.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part20.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part21.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part22.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part23.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part24.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part25.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part26.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part27.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part28.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part29.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part30.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part31.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part32.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part33.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part34.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part35.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part36.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part37.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part38.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part39.root");
    tree->Add( filesPath2 + "latino_QCD_Pt-20toInf_MuEnrichedPt15__part40.root");  
}
  else if (theSample == "QCD1") {
    tree->Add(filesPath1 + "latino_QCD_Pt-15to20_MuEnrichedPt5.root");
    myBaseW = 1273190000/2362016;
  }
  else if (theSample == "QCDEM_0") {
    tree->Add(filesPath2 + "latino_QCD_Pt-15to20_EMEnriched.root");
  }
  else if (theSample == "QCDEM_1") {
    tree->Add(filesPath2 + "latino_QCD_Pt-20to30_EMEnriched.root");
  }
  else if (theSample == "QCDEM_2") {
    //tree->Add(filesPath1 + "latino_QCD_Pt-30to50_EMEnriched.root");  
    tree->Add(filesPath2 + "latino_QCD_Pt-30toInf_DoubleEMEnriched__part1.root");
    tree->Add(filesPath2 + "latino_QCD_Pt-30toInf_DoubleEMEnriched__part2.root");
    tree->Add(filesPath2 + "latino_QCD_Pt-30toInf_DoubleEMEnriched__part3.root");
    tree->Add(filesPath2 + "latino_QCD_Pt-30toInf_DoubleEMEnriched__part4.root");
    tree->Add(filesPath2 + "latino_QCD_Pt-30toInf_DoubleEMEnriched__part5.root");
  }
  else if (theSample == "QCDEM_3") {
    tree->Add(filesPath1 + "latino_QCD_Pt-50to80_EMEnriched.root");
  }
  else if (theSample == "ZJets_high") {
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0000__part0.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0000__part1.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0000__part2.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0000__part3.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0000__part4.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0001__part0.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0001__part1.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0001__part2.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0001__part3.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0001__part4.root");
    tree->Add(filesPath2 + "latino_DYJetsToLL_M-50_0002__part0.root");
  }
  else if (theSample == "ZJets_low") {
   tree->Add(filesPath2 + "latino_DYJetsToLL_M-10to50_0000__part0.root");
   tree->Add(filesPath2 + "latino_DYJetsToLL_M-10to50_0000__part1.root");
   tree->Add(filesPath2 + "latino_DYJetsToLL_M-10to50_0000__part2.root");
   tree->Add(filesPath2 + "latino_DYJetsToLL_M-10to50_0000__part3.root");
   tree->Add(filesPath2 + "latino_DYJetsToLL_M-10to50_0000__part4.root");
   tree->Add(filesPath2 + "latino_DYJetsToLL_M-10to50_0001__part0.root");
  }
  else if (theSample == "Top") {
    tree->Add(filesPath1 + "latino_TT.root");
    negativeWeight =  0.331658;
  }
  else if (theSample == "SingleTop") {
    tree->Add(filesPath1 + "latino_ST_t-channel.root");
    tree->Add(filesPath1 + "latino_ST_tW_antitop.root");
    tree->Add(filesPath1 + "latino_ST_tW_top.root");
  }
  else if (theSample == "dataDEG") {
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0000__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0000__part1.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0000__part2.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0000__part3.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0000__part4.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0001__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0001__part1.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0001__part2.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0001__part3.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0001__part4.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0002__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0002__part1.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0002__part2.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0002__part3.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0002__part4.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0003__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleEG_0003__part1.root");
    tree->Add(filesPath1 + "latino_Run2015C_16Dec2015_DoubleEG_0000__part0.root");
    tree->Add(filesPath1 + "latino_Run2015C_16Dec2015_DoubleEG_0000__part1.root");
  }
  else if (theSample == "dataDMu") {

    /*   tree->Add(filesPath1 + "03Mar_Run2015C_16Dec2015/l1loose__EpTCorr/latino_Run2015C_16Dec2015_DoubleMuon_0000__part0.root");
    tree->Add(filesPath1 + "03Mar_Run2015C_16Dec2015/l1loose__EpTCorr/latino_Run2015C_16Dec2015_DoubleMuon_0000__part1.root");
    tree->Add(filesPath1 + "03Mar_Run2015D_16Dec2015/l1loose__EpTCorr/latino_Run2015D_16Dec2015_DoubleMuon_0000__part0.root");
    tree->Add(filesPath1 + "03Mar_Run2015D_16Dec2015/l1loose__EpTCorr/latino_Run2015D_16Dec2015_DoubleMuon_0000__part1.root");
    */

    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0000__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0000__part1.root");  
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0000__part2.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0000__part3.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0000__part4.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0001__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0001__part1.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0001__part2.root");  
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0001__part3.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0001__part4.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0002__part0.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0002__part1.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0002__part2.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0002__part3.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0002__part4.root");
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0003__part0.root");  
    tree->Add(filesPath1 + "latino_Run2015D_16Dec2015_DoubleMuon_0003__part1.root");
    tree->Add(filesPath1 + "latino_Run2015C_16Dec2015_DoubleMuon_0000__part0.root");
    tree->Add(filesPath1 + "latino_Run2015C_16Dec2015_DoubleMuon_0000__part1.root");
}
  else {
    return;
  }

  tree->SetBranchAddress("evt", &evt);
  tree->SetBranchAddress("run", &run);
  tree->SetBranchAddress("baseW", &baseW);
  tree->SetBranchAddress("puW", &puW);
  tree->SetBranchAddress("jetRho", &jetRho);
  tree->SetBranchAddress("channel", &channel);
  tree->SetBranchAddress("metPfType1", &metPfType1);
  tree->SetBranchAddress("metPfType1Phi", &metPfType1Phi);
  tree->SetBranchAddress("nvtx",&nvtx);
  tree->SetBranchAddress("triggerFakeRate", &triggerFakeRate);
  tree->SetBranchAddress("njet",&njet);

  tree->SetBranchAddress("std_vector_lepton_trackIso", &std_vector_lepton_trackIso);
  tree->SetBranchAddress("std_vector_electron_ecalPFClusterIso", &std_vector_electron_ecalPFClusterIso);
  tree->SetBranchAddress("std_vector_electron_hcalPFClusterIso", &std_vector_electron_hcalPFClusterIso);
  //tree->SetBranchAddress("std_vector_lepton_isTightLepton", &std_vector_lepton_isTightLepton);
  //tree->SetBranchAddress("std_vector_lepton_isLooseLepton", &std_vector_lepton_isLooseLepton);

  tree->SetBranchAddress("std_vector_jet_cmvav2", &std_vector_jet_cmvav2);
  tree->SetBranchAddress("std_vector_lepton_pt", &std_vector_lepton_pt);
  tree->SetBranchAddress("std_vector_lepton_eta", &std_vector_lepton_eta);
  tree->SetBranchAddress("std_vector_lepton_phi", &std_vector_lepton_phi);
  tree->SetBranchAddress("std_vector_lepton_flavour", &std_vector_lepton_flavour);
  tree->SetBranchAddress("std_vector_lepton_isTightMuon", &std_vector_lepton_isTightMuon);
  tree->SetBranchAddress("std_vector_lepton_isMediumMuon", &std_vector_lepton_isMediumMuon); 
  tree->SetBranchAddress("std_vector_lepton_eleIdMedium", &std_vector_lepton_eleIdMedium);
  tree->SetBranchAddress("std_vector_lepton_eleIdTight", &std_vector_lepton_eleIdTight);
  tree->SetBranchAddress("std_vector_lepton_eleIdVeto", &std_vector_lepton_eleIdVeto);
  tree->SetBranchAddress("std_vector_lepton_chargedHadronIso", &std_vector_lepton_chargedHadronIso);
  tree->SetBranchAddress("std_vector_lepton_photonIso", &std_vector_lepton_photonIso);
  tree->SetBranchAddress("std_vector_lepton_neutralHadronIso", &std_vector_lepton_neutralHadronIso);
  tree->SetBranchAddress("std_vector_lepton_sumPUPt", &std_vector_lepton_sumPUPt);
  tree->SetBranchAddress("std_vector_electron_effectiveArea", &std_vector_electron_effectiveArea);
  tree->SetBranchAddress("std_vector_lepton_dz", &std_vector_lepton_dz);
  tree->SetBranchAddress("std_vector_lepton_d0", &std_vector_lepton_d0);
  tree->SetBranchAddress("std_vector_lepton_SIP3D", &std_vector_lepton_SIP3D);
 
  tree->SetBranchAddress("std_vector_electron_dEtaIn", &std_vector_electron_dEtaIn);
  tree->SetBranchAddress("std_vector_electron_dPhiIn", &std_vector_electron_dPhiIn);
  tree->SetBranchAddress("std_vector_electron_full5x5_sigmaIetaIeta", &std_vector_electron_full5x5_sigmaIetaIeta);
  tree->SetBranchAddress("std_vector_electron_hOverE",&std_vector_electron_hOverE);
  tree->SetBranchAddress("std_vector_electron_ooEmooP",&std_vector_electron_ooEmooP);
  tree->SetBranchAddress("std_vector_electron_expectedMissingInnerHits", &std_vector_electron_expectedMissingInnerHits);
 
  tree->SetBranchAddress("std_vector_electron_passConversionVeto", &std_vector_electron_passConversionVeto);


  tree->SetBranchAddress("std_vector_trigger", &std_vector_trigger);
  tree->SetBranchAddress("std_vector_trigger_prescale", &std_vector_trigger_prescale);
  tree->SetBranchAddress("std_vector_jet_csvv2ivf", &std_vector_jet_csvv2ivf);
  tree->SetBranchAddress("std_vector_jet_pt", &std_vector_jet_pt);
  tree->SetBranchAddress("std_vector_jet_eta", &std_vector_jet_eta);
  tree->SetBranchAddress("std_vector_jet_phi", &std_vector_jet_phi);
  tree->SetBranchAddress("std_vector_lepton_closejet_pt", &std_vector_lepton_closejet_pt);
  tree->SetBranchAddress("std_vector_lepton_closejet_eta", &std_vector_lepton_closejet_eta);
  tree->SetBranchAddress("std_vector_lepton_closejet_phi", &std_vector_lepton_closejet_phi);
  tree->SetBranchAddress("std_vector_lepton_closejet_PartonFlavour", &std_vector_lepton_closejet_PartonFlavour);
  tree->SetBranchAddress("std_vector_lepton_closejet_drlj", &std_vector_lepton_closejet_drlj);

  if ( tree->FindBranch("GEN_weight_SM") ) tree->SetBranchAddress("GEN_weight_SM", &GEN_weight_SM);


  //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------

  Long64_t nentries = tree->GetEntries();
  float nE =0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
    tree->GetEntry(jentry);
  
    if (jentry%10000==0) printf("Processing event %d \n", int(jentry));

    //Float_t luminosity = 1269; // pb-1
    //Float_t luminosity = 1567; // pb-1
    //Float_t luminosity = 2120.016; //pb-1
    Float_t luminosity = 2.263; //pb-1 --> New lumi norm tag gives 2.318


    Double_t totalW = 1; 

    if ( !theSample.Contains("data") ) {
      Double_t   genWeight   = 1.0;
      genWeight = GEN_weight_SM/abs(GEN_weight_SM);
      totalW  = genWeight * puW * baseW * luminosity;
    }


    nE++;

    //cout << "1:  " << nE << endl;

    // Loop over leptons
    //--------------------------------------------------------------------------
 
    AnalysisLeptons.clear();
     
    int vsize = std_vector_lepton_pt->size();
  
    for (int i=0; i<vsize; i++) {

      if ( !IsLooseLepton(i)) continue;  
      if ( !IsLooseIsolatedLepton(i) ) continue;   
  
      //      if (std_vector_lepton_isLooseLepton->at(i) == 0) continue;

     float pt  = std_vector_lepton_pt ->at(i);
     float eta = std_vector_lepton_eta->at(i);
     float phi = std_vector_lepton_phi->at(i);
     float id  = std_vector_lepton_flavour ->at(i); 


     if (fabs(id) == 11 && pt <= 13) continue;
     if (fabs(id) == 13 && pt <= 10) continue;       
     if (fabs(id) == 11 && fabs(eta) >= 2.5 ) continue;  
     if (fabs(id) == 13 && fabs(eta) >= 2.4 ) continue;  

     Lepton lep;
    
     lep.index  = i;
     lep.charge = id;
     
     if (fabs(id) == 11)
       { lep.flavor = Electron; }
     else if (fabs(id) == 13)
       { lep.flavor = Muon; }

 
     if (IsTightLepton(i) && IsIsolatedLepton(i) )  { 
     //if ( std_vector_lepton_isTightLepton->at(i) == 1) {
       lep.type = Tight;
     } else {
       lep.type = noTight;
     }
 
     lep.fromW = false;
     lep.iGen = false; 

     TLorentzVector tlv;                                                                                                                      
     tlv.SetPtEtaPhiM(pt, eta, phi, 0);                                                                         
     lep.v = tlv; 
     
    
     AnalysisLeptons.push_back(lep);

    } // end loop on leptons 

 
    // Compute Prompt Rates from MC
    //--------------------------------------------------------------------------  
    if ( theSample == "ZJets_high" && doGen) {


      for (int iMu1 = 0; iMu1 <  AnalysisLeptons.size(); iMu1++) {

	if ( AnalysisLeptons[iMu1].v.Perp()  < 10 ) continue; 
	if ( AnalysisLeptons[iMu1].type != Tight ) continue;
	if ( AnalysisLeptons[iMu1].flavor != Electron && AnalysisLeptons[iMu1].flavor != Muon) continue;

	for (int iMu2 = 0; iMu2 <  AnalysisLeptons.size(); iMu2++) {

	  if ( AnalysisLeptons[iMu1].flavor != AnalysisLeptons[iMu2].flavor) continue; 
	  if ( AnalysisLeptons[iMu1].index ==  AnalysisLeptons[iMu2].index) continue;

	  if ( AnalysisLeptons[iMu2].v.Perp()  < 10 ) continue; 

	  float Minv = (AnalysisLeptons[iMu1].v+AnalysisLeptons[iMu2].v).M();

	  if ( (AnalysisLeptons[iMu1].charge)*(AnalysisLeptons[iMu2].charge) > 0 ) continue; 
	  if ( Minv < 70 ||  Minv > 110 )     continue;
 
	  float pt = std_vector_lepton_pt->at(iMu2);
	  float eta = fabs(std_vector_lepton_eta->at(iMu2));
     
	  if (AnalysisLeptons[iMu1].flavor == Muon) {
   	   
	    h_Muon_fake_pT_bin->Fill(pt, totalW);
	    h_Muon_fake_eta_bin->Fill(eta, totalW);
	    h_Muon_fake_pt_eta_bin->Fill(pt, eta , totalW);

	    if ( AnalysisLeptons[iMu2].type == Tight ) {

	      h_Muon_signal_pT_bin->Fill(pt, totalW);
	      h_Muon_signal_eta_bin->Fill(eta, totalW);
	      h_Muon_signal_pt_eta_bin->Fill(pt, eta , totalW); 
	    }
	  }

	  else if ( AnalysisLeptons[iMu1].flavor == Electron) {

	    h_Ele_fake_pT_bin->Fill(pt, totalW);
	    h_Ele_fake_eta_bin->Fill(eta, totalW);
	    h_Ele_fake_pt_eta_bin->Fill(pt, eta , totalW); 

	    if ( AnalysisLeptons[iMu2].type == Tight ) {

	      h_Ele_signal_pT_bin->Fill(pt, totalW);
	      h_Ele_signal_eta_bin->Fill(eta, totalW);
	      h_Ele_signal_pt_eta_bin->Fill(pt, eta , totalW); 
	    }
	  }
	  

	}
      }
    }

   
  
    /*
    // For normalized the data. 
    //--------------------------------------------------------------------------  
    if ( doGen) {
      
      if (!passTrigger) continue;

      for (int iMu1 = 0; iMu1 <  AnalysisLeptons.size(); iMu1++) {

	if ( AnalysisLeptons[iMu1].v.Perp()  < 10 ) continue; 
	if ( AnalysisLeptons[iMu1].flavor != Electron && AnalysisLeptons[iMu1].flavor != Muon) continue;
	for (int iMu2 = 0; iMu2 <  AnalysisLeptons.size(); iMu2++) {

	  if ( AnalysisLeptons[iMu1].flavor != AnalysisLeptons[iMu2].flavor) continue; 
	  if ( AnalysisLeptons[iMu1].index ==  AnalysisLeptons[iMu2].index) continue;

	  if ( AnalysisLeptons[iMu2].v.Perp()  < 10 ) continue; 

	  float Minv = (AnalysisLeptons[iMu1].v+AnalysisLeptons[iMu2].v).M();

	  if ( (AnalysisLeptons[iMu1].charge)*(AnalysisLeptons[iMu2].charge) > 0 ) continue; 
	  if ( Minv < 70 ||  Minv > 110 )     continue;
 
	  float pt = std_vector_lepton_pt->at(iMu2);
	  float eta = fabs(std_vector_lepton_eta->at(iMu2));
     
	  if (AnalysisLeptons[iMu1].flavor == Muon) {
   	   
	    h_Muon_fake_mll->Fill(Minv);
	    if ( AnalysisLeptons[iMu1].type == Tight && AnalysisLeptons[iMu2].type == Tight ) {
	      h_Muon_signal_mll->Fill(Minv);
	    }
	  }

	  else if ( AnalysisLeptons[iMu1].flavor == Electron) {
	    h_Ele_fake_mll->Fill(Minv);
	    
	    if ( AnalysisLeptons[iMu1].type == Tight &&  AnalysisLeptons[iMu2].type == Tight  ) {

	     h_Ele_signal_mll->Fill(Minv);
	    }
	  }
	  
	}
      }
    }
    */

    // Perform analysis
    //--------------------------------------------------------------------------  

   if ( AnalysisLeptons.size() != 1) continue;

   // Selected triggers in data 
   //--------------------------------------------------------------------------  
    bool passTrigger = true;

    if ( theSample.Contains("dataDMu") && AnalysisLeptons[0].flavor == Muon) {
      
      passTrigger = false; 

      int iBit = std_vector_trigger->size();
      int nBits = 0; 
      
      for (int i=0; i < iBit; i++) {
	if ( std_vector_trigger->at(i) == 1) {
	  h_triggerBits->Fill(i); 	 
	  nBits ++ ;
	} 
      }

      //Selecting [17 - 20] HLT_Mu8, HLT_Mu17, HLT_Mu24, HLT_Mu34
      //Selecting [21 - 24] HLT_Mu8_TrkIsoVVL, HLT_Mu17_TrkIsoVVL, HLT_Mu24_TrkIsoVVL, HLT_Mu34_TrkIsoVVL
      if ( AnalysisLeptons[0].v.Perp() > 10   && 
	   AnalysisLeptons[0].v.Perp() <= 20 &&
	   (std_vector_trigger->at(22) == 1) ) {

	     passTrigger = true;
	     //totalW = totalW *(1567/1.702);  //Lumi HLT_Mu8: 0.57 pb; Lumi HLT_Mu8_TrkIsoVVL: 1.132 pb
	     totalW = totalW *(2263/1.386);  //Lumi HLT_Mu8: 0.674 pb; Lumi HLT_Mu8_TrkIsoVVL: 1.386 pb
      }

      /*if ( std_vector_trigger->at(3) != 1 &&
	   std_vector_trigger->at(4) != 1 &&
	   std_vector_trigger->at(5) != 1 &&
	   std_vector_trigger->at(6) != 1 
	   ) continue;*/

	if ( AnalysisLeptons[0].v.Perp() > 20  && 
	     (std_vector_trigger->at(23) == 1 )) {

	       passTrigger = true;
	       //totalW = totalW *(1567/271.288);  //Lumi HLT_Mu17: 135.644 pb; Lumi HLT_Mu17_TrkIsoVVL: 135.644 pb
	       totalW = totalW *(2263/201.951);  //Lumi HLT_Mu17: 190.107 pb; Lumi HLT_Mu17_TrkIsoVVL: 201.951 pb
      }

	
	// Selecting [48-51] HLT_Mu8_TrkIsoVVL, HLT_Mu17_TrkIsoVVL, HLT_Mu24_TrkIsoVVL, HLT_Mu34_TrkIsoVVL
	/*  if ( //std_vector_trigger->at(48) == 1 ||
	   //std_vector_trigger->at(49) == 1 ||
	   //std_vector_trigger->at(50) == 1 ||
	   std_vector_trigger->at(49) == 1 
	   ) passTrigger = true ;
     
      totalW =   totalW * (1.269 / 0.135); // HLT_Mu17_TrkIsoVVL_v2-->0.135 fb

	*/
    }

    if ( theSample.Contains("dataDEG") && AnalysisLeptons[0].flavor == Electron) {

      passTrigger = false;

      int iBit = std_vector_trigger->size();
      int nBits = 0; 

      for (int i=0; i < iBit; i++) {
	if ( std_vector_trigger->at(i) == 1) {
	  h_triggerBits->Fill(i); 	 
	  nBits ++ ;
	} 
      }
     
      if ( AnalysisLeptons[0].v.Perp() > 13   && 
	   AnalysisLeptons[0].v.Perp() <= 25  &&
	   std_vector_trigger->at(31) == 1 )  {

	passTrigger = true;
	//totalW = totalW *(1567/0.459);  //Lumi HLT_Ele8: 0.459 pb;
	//totalW = totalW *(2263/0.559);  //Lumi HLT_Ele8: 0.559 pb;
	totalW = totalW *(2263/11.204);  //Lumi HLT_Ele8_Iso: 11.204 pb;
      }

      /*if ( AnalysisLeptons[0].v.Perp() > 15   && 
	   AnalysisLeptons[0].v.Perp() <= 30  &&
	   (std_vector_trigger->at(27) == 1))  {
	
	passTrigger = true;
	//totalW = totalW *(1567/12.473);  //Lumi HLT_Ele12: 6.379 pb; HLT_Ele12_Iso: 6.094
	totalW = totalW *(2120.016/9.948);  //Lumi HLT_Ele12: 9.948 pb; HLT_Ele12_Iso: 8.876
	}*/

      if ( AnalysisLeptons[0].v.Perp() > 25  &&
	   (std_vector_trigger->at(33) == 1))  {
	
	passTrigger = true;
	//totalW = totalW *(1567/6.699);  //Lumi HLT_Ele23: 4.578 pb; HLT_Ele23_Iso: 2.121
	totalW = totalW *(2263/3.201);  //Lumi HLT_Ele23: 6.147  pb; HLT_Ele23_Iso: 3.201
      }

      /*if ( std_vector_trigger->at(18) == 1 
	   //std_vector_trigger->at(52) == 1 ||
	   //std_vector_trigger->at(53) == 1 ||
	   //std_vector_trigger->at(54) == 1 
	   ) passTrigger = true ;
      */
      //totalW =   totalW * std_vector_trigger_prescale->at(4);
    }
	  
    if (!passTrigger) continue;

    int index = AnalysisLeptons[0].index; 

    h_LooseMuons_All->Fill(1);
    h_LooseMuons_Iso->Fill(MuonIsolation(index));    

    double mt = giveMT(AnalysisLeptons[0].index, metPfType1, metPfType1Phi); 

    AnalysisJets.clear();
    AnalysisLeptonJets.clear();
   
    Float_t closeJetBtag = -999;

    // --------> away jet 
    int jetsize = std_vector_jet_pt->size();
        
    for ( int j=0; j<jetsize; j++) {
      
      Jet jet_; 
      
      jet_.index = j; 
      
      float pt = std_vector_jet_pt->at(j);
      float eta = std_vector_jet_eta->at(j);
      float phi = std_vector_jet_phi->at(j);
      float btag = std_vector_jet_csvv2ivf->at(j);     

      if ( pt <= 10.0) continue;  // follow Gui analysis... has sense.

      TLorentzVector jetvt;
      jetvt.SetPtEtaPhiM(pt, eta, phi, 0);
      
      jet_.v = jetvt;

      float dR =  jetvt.DeltaR(AnalysisLeptons[0].v);
      //cout << j << "  " << dR << endl;

      if ( dR <= 0.3 ) {
	closeJetBtag = std_vector_jet_cmvav2->at(j);
	if ( std_vector_jet_cmvav2->at(j) >= -0.715) { 
	  AnalysisLeptonJets.push_back(jet_);	
	}
      } else if ( dR > 1 )   {       
	AnalysisJets.push_back(jet_);
      }
         
    }

    if (AnalysisLeptonJets.size() > 0 ) continue;

    if (AnalysisJets.size() < 1 ) continue;
    if (AnalysisJets[0].v.Perp() <  inputJetEt) continue; 


    //int jetIndex = AnalysisJets[0].index;  
    //if ( std_vector_jet_pt->at(jetIndex) < 0 ) continue;
    
    
    
    if (doGen) continue;
      
    //if ( fabs(std_vector_lepton_closejet_PartonFlavour->at(0)) != 5) continue;
    //if ( std_vector_lepton_closejet_drlj->at(0) > 0.3 )   continue;

    if ( std_vector_lepton_pt->at(index) <= 10 ) continue; 


    if ( AnalysisLeptons[0].flavor == Muon ) { 
	  
      float pt = std_vector_lepton_pt->at(index);
      float eta = fabs(std_vector_lepton_eta->at(index));
      float dxy = fabs(std_vector_lepton_d0->at(index));
      float jetEt = std_vector_lepton_closejet_pt->at(index);
      float muSIP3D = std_vector_lepton_SIP3D->at(index);
      float elSIP3D =  std_vector_lepton_SIP3D->at(index);
      float ptIso = pt*(1+MuonIsolation(index));//pt*(1+min(0.0,(MuonIsolation(index)-0.15)));
      float drlj =  std_vector_lepton_closejet_drlj->at(index);
      //cout << std_vector_lepton_closejet_drlj->at(index) << endl;

      // control plots
      h_Muon_fake_Mt->Fill(mt, totalW);
      h_Muon_fake_Met->Fill(metPfType1, totalW);
      h_Muon_fake_pt->Fill(pt, totalW);
      h_Muon_fake_CRjetEt->Fill(pt, totalW);
      if ( metPfType1 > 20 )  h_Muon_fake_Mt_Met20->Fill(mt, totalW);     
      
      if ( AnalysisLeptons[0].type == Tight)  {
	
	h_Muon_signal_Mt->Fill(mt, totalW);
	h_Muon_signal_Met->Fill(metPfType1, totalW);
	h_Muon_signal_pt->Fill(pt, totalW);
	h_Muon_signal_CRjetEt->Fill(pt, totalW);
	if ( metPfType1 > 30 )  h_Muon_signal_Mt_Met20->Fill(mt, totalW);     
      
      }


     if ( mt < 20 && metPfType1 < 20 ) { 

	  h_Muon_fake_pT_bin->Fill(pt, totalW);
	  h_Muon_fake_eta_bin->Fill(eta, totalW);
	  h_Muon_fake_pt_eta_bin->Fill(pt, eta , totalW); 
	  h_Muon_fake_nPV->Fill(nvtx, totalW); 
	  h_Muon_fake_dxy_pt->Fill((dxy/pt)*100 , totalW);
	  h_Muon_fake_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
	  h_Muon_fake_dxy->Fill(dxy , totalW);
	  h_Muon_pt->Fill(pt, totalW);
	  h_Muon_fake_sip3d->Fill(muSIP3D, totalW);
	  h_Muon_fake_ptIso->Fill(ptIso, totalW);
	  h_Muon_ptIso ->Fill(ptIso, totalW);
	  h_Muon_fake_ptIso_eta_bin->Fill(ptIso, eta , totalW); 
	  h_Muon_fake_pt->Fill(pt, totalW);
	  h_Muon_fake_iso->Fill(MuonIsolation(index), totalW); 
	  h_Muon_fake_ChIso->Fill(std_vector_lepton_chargedHadronIso->at(index) , totalW); 
	  h_Muon_fake_PhIso->Fill(std_vector_lepton_photonIso->at(index), totalW); 
	  h_Muon_fake_NeuIso->Fill(std_vector_lepton_neutralHadronIso->at(index),   totalW); 
	  h_Muon_fake_btag->Fill(closeJetBtag, totalW); 
	 
	  if ( jetEt > 10 ) { 
	    h_Muon_fake_jetEt->Fill(jetEt, totalW);
	    h_Muon_jetEt->Fill(jetEt, totalW);
	    h_Muon_jetEt_pt->Fill(jetEt, pt, totalW); 
	    h_Muon_jetEt_iso->Fill(jetEt,MuonIsolation(index) , totalW); 
	    h_Muon_fake_ptRel->Fill(pt/jetEt ,  totalW); 
	  }


	  //if ( pt > 30 && pt < 35 ) 
	    //cout << "Run:  " << run << "; Event:  " << evt << "; Pt:  " << pt << endl;

	  if ( AnalysisLeptons[0].type == Tight)  {
	    
	    h_Muon_signal_pT_bin->Fill(pt, totalW);
	    h_Muon_signal_eta_bin->Fill(eta, totalW);
	    h_Muon_signal_pt_eta_bin->Fill(pt, eta , totalW); 
	    h_Muon_signal_nPV->Fill(nvtx, totalW); 
	    h_Muon_signal_dxy_pt->Fill((dxy/pt)*100 , totalW);
	    h_Muon_signal_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
	    h_Muon_signal_sip3d->Fill(muSIP3D, totalW);
	    h_Muon_signal_ptIso->Fill(ptIso, totalW);
	    h_Muon_signal_pt->Fill(pt, totalW);
	    h_Muon_signal_iso->Fill(MuonIsolation(index), totalW); 
	    h_Muon_signal_ptIso_eta_bin->Fill(ptIso, eta , totalW); 
	    h_Muon_signal_ChIso->Fill(std_vector_lepton_chargedHadronIso->at(index) , totalW); 
	    h_Muon_signal_PhIso->Fill(std_vector_lepton_photonIso->at(index), totalW); 
	    h_Muon_signal_NeuIso->Fill(std_vector_lepton_neutralHadronIso->at(index),   totalW); 

	    h_Muon_signal_btag->Fill(closeJetBtag, totalW); 
 
	    if ( jetEt > 10 ) { 
	      h_Muon_fake_jetEt->Fill(jetEt, totalW);
	      h_Muon_fake_ptRel->Fill(pt/jetEt ,  totalW); 
	    }
	  
	  }    
     }
    }
	
    if ( AnalysisLeptons[0].flavor == Electron ) { 
		 
      float pt = std_vector_lepton_pt->at(index);
      float eta = fabs(std_vector_lepton_eta->at(index));
      float dxy = fabs(std_vector_lepton_d0->at(index));
      float jetEt = std_vector_lepton_closejet_pt->at(index);
      float elSIP3D =  std_vector_lepton_SIP3D->at(index);
      float ptIso = pt*(1+MuonIsolation(index));
      float elptIso = pt*(1+ElectronIsolation(index));

     
      // control plots
      h_Ele_fake_Mt->Fill(mt, totalW);
      h_Ele_fake_Met->Fill(metPfType1, totalW);
      h_Ele_fake_pt->Fill(pt, totalW);
      h_Ele_fake_CRjetEt->Fill(pt, totalW);
      if ( metPfType1 > 20 )  h_Ele_fake_Mt_Met20->Fill(mt, totalW);     
      
      if ( AnalysisLeptons[0].type == Tight)  {

	h_Ele_signal_Mt->Fill(mt, totalW);
	h_Ele_signal_Met->Fill(metPfType1, totalW);
	h_Ele_signal_pt->Fill(pt, totalW);
	h_Ele_signal_CRjetEt->Fill(pt, totalW);
	if ( metPfType1 > 20 )  h_Ele_signal_Mt_Met20->Fill(mt, totalW);     
      

      }

      if ( mt < 20 && metPfType1 < 20 ) { 

	  h_Ele_fake_pT_bin->Fill(pt, totalW);
	  h_Ele_fake_eta_bin->Fill(eta, totalW);
	  h_Ele_fake_pt_eta_bin->Fill(pt, eta , totalW); 
	  h_Ele_fake_nPV->Fill(nvtx, totalW); 
	  h_Ele_jetEt->Fill(jetEt, totalW);
	  if ( pt > 13) h_Ele_fake_jetEt->Fill(jetEt, totalW);
	  h_Ele_pt->Fill(pt, totalW);
	  h_Ele_fake_sip3d->Fill(elSIP3D, totalW);
	  h_Ele_fake_ptIso->Fill(elptIso, totalW);
	  h_Ele_ptIso ->Fill(elptIso, totalW);
	  h_Ele_fake_pt->Fill(pt, totalW);
	  h_Ele_fake_iso->Fill(ElectronIsolation(1), totalW); 
	  h_Ele_fake_ptIso_eta_bin->Fill(elptIso, eta , totalW); 

	  if ( AnalysisLeptons[0].type == Tight)  {
	    
	    h_Ele_signal_pT_bin->Fill(pt, totalW);
	    h_Ele_signal_eta_bin->Fill(eta, totalW);
	    h_Ele_signal_pt_eta_bin->Fill(pt, eta , totalW); 
	    h_Ele_signal_nPV->Fill(nvtx, totalW); 
	    if ( pt > 13) h_Ele_signal_jetEt->Fill(jetEt, totalW);
	    h_Ele_signal_sip3d->Fill(elSIP3D, totalW);
	    h_Ele_signal_ptIso->Fill(elptIso, totalW);
	    h_Ele_signal_pt->Fill(pt, totalW);
	    h_Ele_signal_iso->Fill(ElectronIsolation(1), totalW); 
	    h_Ele_signal_ptIso_eta_bin->Fill(elptIso, eta , totalW); 

	  } 
	  }
    }
  
    

  } //end Loop
  

  // Save output
  //--------------------------------------------------------------------------
  root_output->cd();
 
  root_output->Write("", TObject::kOverwrite);
 
  root_output->Close();


} // end main function



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
	    std_vector_electron_passConversionVeto->at(k))
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
// Compute MT variable 
//------------------------------------------------------------------------------
Double_t giveMT(int k, double MET, double MET_Phi)
{
  TLorentzVector lep;
  lep.SetPtEtaPhiM(std_vector_lepton_pt->at(k), std_vector_lepton_eta->at(k), std_vector_lepton_phi->at(k), 0);

  // --> Transversal Mass ( lepton, MET) 
  Double_t pt = std_vector_lepton_pt->at(k);
  Double_t px = lep.Px();
  Double_t py = lep.Py();
  Double_t pxnu = cos(MET_Phi)*MET;
  Double_t pynu = sin(MET_Phi)*MET;
  Double_t ptnu = sqrt(pxnu*pxnu +  pynu*pynu);
	  
  Double_t MT = sqrt(fabs((pt+ptnu)*(pt+ptnu) -(px+pxnu)*(px+pxnu)-(py+pynu)*(py+pynu)));

  return MT;
}

