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
  TLorentzVector v;

};


struct Jet { 

  UInt_t         index;
  TLorentzVector v;

};


std::vector<Lepton> AnalysisLeptons;
std::vector<Jet> AnalysisJets;

Float_t          jetRho;
Float_t          puW;
Float_t          channel;
Float_t          baseW;
Float_t          pfType1Met;  
Float_t          pfType1Metphi;
Float_t          mth;
Float_t          dphilmet1;
Float_t          dphilmet2;
Float_t          nvtx;
Float_t          triggerFakeRate;
Float_t          trigger;
Float_t          njet;

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
vector<float>   *std_vector_lepton_BestTrackdz;
vector<float>   *std_vector_lepton_BestTrackdxy;
vector<float>   *std_vector_lepton_muSIP3D;
vector<float>   *std_vector_lepton_elSIP3D;

vector<float>   *std_vector_electron_dEtaIn;
vector<float>   *std_vector_electron_dPhiIn;
vector<float>   *std_vector_electron_full5x5_sigmaIetaIeta;
vector<float>   *std_vector_electron_hOverE;
vector<float>   *std_vector_electron_ooEmooP;
vector<float>   *std_vector_electron_expectedMissingInnerHits;
vector<float>   *std_vector_electron_d0;
vector<float>   *std_vector_electron_dz;
vector<float>   *std_vector_electron_passConversionVeto;

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

Float_t  GEN_weight_SM = 1.0;
Float_t  negativeWeight = 1.0;

Float_t myBaseW = 1.0;


//------------------------------------------------------------------------------
// MAIN function
//------------------------------------------------------------------------------

void computeLatinoFR (TString theSample="WJets", float inputJetEt = 0) {

  // ------ root settings ---------
 
  TH1::SetDefaultSumw2();
  
  TString path = Form("/gpfs/csic_users/calderon/IFCA_2015/WW13TeV_V1/addWJet/rootfiles/");

  TFile* root_output;  
  
  TString jet = Form("%i", (int)inputJetEt);
  
  root_output =  new TFile(path + theSample + "_jet"+ jet +"_FakeRate.root", "recreate");
  

  //------------------------------------------------------------------------------
  // Define histograms 
  //------------------------------------------------------------------------------

  // Define pt bins for FR matrix
  int pTbin = 7;
  Double_t pTbins[8] = { 15., 20., 25., 30., 35.,40., 45., 50.};

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


  //----------------------------------------------------------------------------
  //--- Open loose leptons files
  //----------------------------------------------------------------------------
  
  TString filesPath1, filesPath2, filesPath3, filesPath4;
  filesPath1 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/MC_Spring15/21Oct/25ns/";
  filesPath2 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/Data13TeV/21Oct/25ns/";
  filesPath3 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/MC_Spring15/17Sep/25ns/";
  filesPath4 = "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/MC_Spring15/05Aug/25ns/";

  TChain* tree = new TChain("latino", "latino");

  if (theSample == "WJets") {
    tree->Add(filesPath1 + "latino_WJetsToLNu.root");
    myBaseW = 61526.7/24184766;
    negativeWeight = 0.683938;
  }
  else if (theSample == "QCD2") {
    tree->Add( filesPath4 + "latino_QCD_Pt-20toInf_MuEnrichedPt15.root");
  }
  else if (theSample == "QCD1") {
    tree->Add(filesPath3 + "latino_QCD_Pt-15to20_MuEnrichedPt5.root");
    myBaseW = 1273190000/2362016;
  }
  else if (theSample == "QCDEM_0") {
    tree->Add(filesPath3 + "latino_QCD_Pt-15to20_EMEnriched.root");
  }
  else if (theSample == "QCDEM_1") {
    tree->Add(filesPath3 + "latino_QCD_Pt-20to30_EMEnriched.root");
  }
  else if (theSample == "QCDEM_2") {
    tree->Add(filesPath3 + "latino_QCD_Pt-30to50_EMEnriched.root");  
  }
  else if (theSample == "QCDEM_3") {
    tree->Add(filesPath1 + "latino_QCD_Pt-50to80_EMEnriched.root");
  }
  else if (theSample == "ZJets_high") {
    tree->Add(filesPath3 + "latino_DYJetsToLL_M-50.root");
    myBaseW = 6025.2/28579957;
    negativeWeight = 0.66998;
  }
  else if (theSample == "ZJets_low") {
   tree->Add(filesPath3 + "latino_DYJetsToLL_M-10to50.root");
   myBaseW = 18610/30021559;
   negativeWeight = 0.727601;
  }
  else if (theSample == "dataDEG") {

    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleEG_0.root");
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleEG_1.root");
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleEG_2.root");
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleEG_3.root");
    tree->Add(filesPath2 + "latino_Run2015D_05Oct2015_DoubleEG_0.root");
    tree->Add(filesPath2 + "latino_Run2015D_05Oct2015_DoubleEG_0.root");
  }
  else if (theSample == "dataDMu") {
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleMuon_0.root");
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleMuon_1.root");
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleMuon_2.root");
    tree->Add(filesPath2 + "latino_Run2015D_PromptReco_DoubleMuon_3.root");
    tree->Add(filesPath2 + "latino_Run2015D_05Oct2015_DoubleMuon_0.root");
    tree->Add(filesPath2 + "latino_Run2015D_05Oct2015_DoubleMuon_1.root");
  }
  else {
    return;
  }


  tree->SetBranchAddress("baseW", &baseW);
  tree->SetBranchAddress("puW", &puW);
  tree->SetBranchAddress("jetRho", &jetRho);
  tree->SetBranchAddress("channel", &channel);
  tree->SetBranchAddress("pfType1Met", &pfType1Met);
  tree->SetBranchAddress("pfType1Metphi", &pfType1Metphi);
  tree->SetBranchAddress("mth", &mth);
  tree->SetBranchAddress("dphilmet1", &dphilmet1);
  tree->SetBranchAddress("dphilmet2", &dphilmet2);
  tree->SetBranchAddress("nvtx",&nvtx);
  tree->SetBranchAddress("triggerFakeRate", &triggerFakeRate);
  tree->SetBranchAddress("njet",&njet);

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
  tree->SetBranchAddress("std_vector_lepton_BestTrackdz", &std_vector_lepton_BestTrackdz);
  tree->SetBranchAddress("std_vector_lepton_BestTrackdxy", &std_vector_lepton_BestTrackdxy);
  tree->SetBranchAddress("std_vector_lepton_muSIP3D", &std_vector_lepton_muSIP3D);
  tree->SetBranchAddress("std_vector_lepton_elSIP3D", &std_vector_lepton_elSIP3D);

  
  tree->SetBranchAddress("std_vector_electron_dEtaIn", &std_vector_electron_dEtaIn);
  tree->SetBranchAddress("std_vector_electron_dPhiIn", &std_vector_electron_dPhiIn);
  tree->SetBranchAddress("std_vector_electron_full5x5_sigmaIetaIeta", &std_vector_electron_full5x5_sigmaIetaIeta);
  tree->SetBranchAddress("std_vector_electron_hOverE",&std_vector_electron_hOverE);
  tree->SetBranchAddress("std_vector_electron_ooEmooP",&std_vector_electron_ooEmooP);
  tree->SetBranchAddress("std_vector_electron_expectedMissingInnerHits", &std_vector_electron_expectedMissingInnerHits);
  tree->SetBranchAddress("std_vector_electron_d0",&std_vector_electron_d0);
  tree->SetBranchAddress("std_vector_electron_dz", &std_vector_electron_dz);
  tree->SetBranchAddress("std_vector_electron_passConversionVeto", &std_vector_electron_passConversionVeto);

  tree->SetBranchAddress("std_vector_trigger", &std_vector_trigger);
  tree->SetBranchAddress("std_vector_trigger_prescale", &std_vector_trigger_prescale);

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

    Float_t luminosity = 1269; // pb-1

    if ( theSample == "QCD1") {
      baseW = 559921.3;
    }
    if ( theSample == "QCD2" ) {
      baseW = 161679.9;
    }
    if ( theSample == "QCDEM_0") {
      baseW = 599960.88236;
    }
    if ( theSample == "QCDEM_1") {
      baseW = 60767.86121;
    }
    if ( theSample == "QCDEM_2") {
      baseW = 28961.53312;
    }
    if ( theSample == "QCDEM_3") {
      baseW = 3788.10779;
    }

    Double_t totalW = 1; 

    if ( !theSample.Contains("data") ) {
      Double_t   genWeight   = 1.0;
      genWeight = GEN_weight_SM/abs(GEN_weight_SM)/negativeWeight;
      totalW  = genWeight * puW * myBaseW * luminosity;
    }

    nE++;

       
    //--------------------------------------------------------------------------
    // Loop over leptons
    //--------------------------------------------------------------------------
 
    AnalysisLeptons.clear();
     
    int vsize = std_vector_lepton_pt->size();
  
    for (int i=0; i<vsize; i++) {
     
     float pt  = std_vector_lepton_pt ->at(i);
     float eta = std_vector_lepton_eta->at(i);
     float phi = std_vector_lepton_phi->at(i);
     float id  = std_vector_lepton_flavour ->at(i);

     Lepton lep;
    
     lep.index  = i;
     lep.charge = id;
     
     if (fabs(id) == 11)
       { lep.flavor = Electron; }
     else if (fabs(id) == 13)
       { lep.flavor = Muon; }

     if (IsTightLepton(i) && IsIsolatedLepton(i) )  { 
       lep.type = Tight;
     } else {
       lep.type = noTight;
     }
  
     if ( !IsLooseLepton(i)) continue;  
     if ( !IsLooseIsolatedLepton(i) ) continue;   
    
     TLorentzVector tlv;                                                                                                                      
     tlv.SetPtEtaPhiM(pt, eta, phi, 0);                                                                         
     lep.v = tlv;   
   
     AnalysisLeptons.push_back(lep);

    } // end loop on leptons 
    

    //-------------------------------------------------------------------------- 
    // Specify the trigger selection
    //--------------------------------------------------------------------------  

    bool passTrigger = true;

    if ( theSample.Contains("dataDMu") ) {

      passTrigger = false; 

      int iBit = std_vector_trigger->size();
      int nBits = 0; 

      for (int i=0; i < iBit; i++) {
	if ( std_vector_trigger->at(i) == 1) {
	  h_triggerBits->Fill(i); 	 
	  nBits ++ ;
	} 
      }

      // Selecting HLT_Mu8, HLT_Mu17, HLT_Mu24, HLT_Mu34
      /*if ( std_vector_trigger->at(3) != 1 &&
	   std_vector_trigger->at(4) != 1 &&
	   std_vector_trigger->at(5) != 1 &&
	   std_vector_trigger->at(6) != 1 
	   ) continue;*/

      // Selecting [48-51] HLT_Mu8_TrkIsoVVL, HLT_Mu17_TrkIsoVVL, HLT_Mu24_TrkIsoVVL, HLT_Mu34_TrkIsoVVL
      if ( //std_vector_trigger->at(48) == 1 ||
	   //std_vector_trigger->at(49) == 1 ||
	   //std_vector_trigger->at(50) == 1 ||
	   std_vector_trigger->at(49) == 1 
	 ) passTrigger = true ;
     
      totalW =   totalW * (1.269 / 0.135); // HLT_Mu17_TrkIsoVVL_v2-->0.135 fb
    }

    if ( theSample.Contains("dataDEG") ) {

      passTrigger = false;

      int iBit = std_vector_trigger->size();
      int nBits = 0; 

      for (int i=0; i < iBit; i++) {
	if ( std_vector_trigger->at(i) == 1) {
	  h_triggerBits->Fill(i); 	 
	  nBits ++ ;
	} 
      }

      if ( std_vector_trigger->at(34) == 1 
	   //std_vector_trigger->at(52) == 1 ||
	   //std_vector_trigger->at(53) == 1 ||
	   //std_vector_trigger->at(54) == 1 
	   ) passTrigger = true ;
      
      //totalW =   totalW * std_vector_trigger_prescale->at(4);
    }
  

    if (!passTrigger) continue;  // require to pass the trigger 

    if ( AnalysisLeptons.size() != 1) continue;  // require one and only one lepton 

    int index = AnalysisLeptons[0].index; 

    h_LooseMuons_All->Fill(1);
    h_LooseMuons_Iso->Fill(MuonIsolation(index));    

    double mt = giveMT(AnalysisLeptons[0].index, pfType1Met, pfType1Metphi); 


    //--------------------------------------------------------------------------
    // loop over jets 
    //--------------------------------------------------------------------------

    AnalysisJets.clear();
  
    // --------> away jet 
    int jetsize = std_vector_jet_pt->size();    
    
    for ( int j=0; j<jetsize; j++) {
      
      Jet jet_; 
      
      jet_.index = j; 
      
      float pt = std_vector_jet_pt->at(j);
      float eta = std_vector_jet_eta->at(j);
      float phi = std_vector_jet_phi->at(j);
     
      if ( pt < 0) continue; 

      TLorentzVector jetvt;
      jetvt.SetPtEtaPhiM(pt, eta, phi, 0);
      
      jet_.v = jetvt;

      float dR =  jetvt.DeltaR(AnalysisLeptons[0].v);
      //cout << j << "  " << dR << endl;
      if ( dR < 1.0 ) continue;
      
      AnalysisJets.push_back(jet_);
      
    }
   
  

    if ( AnalysisJets[0].v.Perp() <  inputJetEt) continue; 

    //int jetIndex = AnalysisJets[0].index;  
    //if ( std_vector_jet_pt->at(jetIndex) < 0 ) continue;
    
    if ( std_vector_lepton_pt->at(index)  < 10 ) continue; 
    
   
    //if ( std_vector_lepton_closejet_drlj->at(0) > 0.5 )   continue;


    // Fill histograms for muon FR and PR      
    //--------------------------------------------------------------------------

    if ( AnalysisLeptons[0].flavor == Muon ) { 
	  
      float pt = std_vector_lepton_pt->at(index);
      float eta = fabs(std_vector_lepton_eta->at(index));
      float dxy = fabs(std_vector_lepton_BestTrackdxy->at(index));
      float jetEt = std_vector_lepton_closejet_pt->at(index);
      float muSIP3D = std_vector_lepton_muSIP3D->at(index);
      float elSIP3D = std_vector_lepton_elSIP3D->at(index);
      float ptIso = pt*(1+MuonIsolation(index));

      //cout << std_vector_lepton_closejet_drlj->at(index) << endl;

      // control plots
      h_Muon_fake_Mt->Fill(mt, totalW);
      h_Muon_fake_Met->Fill(pfType1Met, totalW);
      h_Muon_fake_pt->Fill(pt, totalW);
      h_Muon_fake_CRjetEt->Fill(pt, totalW);
      if ( pfType1Met > 20 )  h_Muon_fake_Mt_Met20->Fill(mt, totalW);     
      
      if ( AnalysisLeptons[0].type == Tight)  {
	
	h_Muon_signal_Mt->Fill(mt, totalW);
	h_Muon_signal_Met->Fill(pfType1Met, totalW);
	h_Muon_signal_pt->Fill(pt, totalW);
	h_Muon_signal_CRjetEt->Fill(pt, totalW);
	if ( pfType1Met > 30 )  h_Muon_signal_Mt_Met20->Fill(mt, totalW);     
      
      }


     if ( mt < 20 && pfType1Met < 20 ) { 

	  h_Muon_fake_pT_bin->Fill(pt, totalW);
	  h_Muon_fake_eta_bin->Fill(eta, totalW);
	  h_Muon_fake_pt_eta_bin->Fill(pt, eta , totalW); 
	  h_Muon_fake_nPV->Fill(nvtx, totalW); 
	  h_Muon_fake_dxy_pt->Fill((dxy/pt)*100 , totalW);
	  h_Muon_fake_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
	  h_Muon_jetEt->Fill(jetEt, totalW);
	  h_Muon_fake_dxy->Fill(dxy , totalW);
	  h_Muon_fake_jetEt->Fill(jetEt, totalW);
	  h_Muon_pt->Fill(pt, totalW);
	  h_Muon_fake_sip3d->Fill(muSIP3D, totalW);
	  h_Muon_fake_ptIso->Fill(ptIso, totalW);
	  h_Muon_ptIso ->Fill(ptIso, totalW);
	  h_Muon_fake_ptIso_eta_bin->Fill(ptIso, eta , totalW); 
	  h_Muon_jetEt_pt->Fill(jetEt, pt, totalW); 
	  h_Muon_fake_pt->Fill(pt, totalW);
	  h_Muon_fake_iso->Fill(MuonIsolation(1), totalW); 


	  if ( AnalysisLeptons[0].type == Tight)  {
	    
	    h_Muon_signal_pT_bin->Fill(pt, totalW);
	    h_Muon_signal_eta_bin->Fill(eta, totalW);
	    h_Muon_signal_pt_eta_bin->Fill(pt, eta , totalW); 
	    h_Muon_signal_nPV->Fill(nvtx, totalW); 
	    h_Muon_signal_dxy_pt->Fill((dxy/pt)*100 , totalW);
	    h_Muon_signal_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
	    h_Muon_signal_jetEt->Fill(jetEt, totalW);
	    h_Muon_signal_sip3d->Fill(muSIP3D, totalW);
	    h_Muon_signal_ptIso->Fill(ptIso, totalW);
	    h_Muon_signal_pt->Fill(pt, totalW);
	    h_Muon_signal_iso->Fill(MuonIsolation(1), totalW); 
	    h_Muon_signal_ptIso_eta_bin->Fill(ptIso, eta , totalW); 
	  
	  }    
     }
    }



    // Fill histograms for electron FR and PR      
    //--------------------------------------------------------------------------
	
    if ( AnalysisLeptons[0].flavor == Electron ) { 
		 
      float pt = std_vector_lepton_pt->at(index);
      float eta = fabs(std_vector_lepton_eta->at(index));
      float dxy = fabs(std_vector_lepton_BestTrackdxy->at(index));
      float jetEt = std_vector_lepton_closejet_pt->at(index);
      float elSIP3D = std_vector_lepton_elSIP3D->at(index);
      float ptIso = pt*(1+MuonIsolation(index));
      float elptIso = pt*(1+ElectronIsolation(index));

     
      // control plots
      h_Ele_fake_Mt->Fill(mt, totalW);
      h_Ele_fake_Met->Fill(pfType1Met, totalW);
      h_Ele_fake_pt->Fill(pt, totalW);
      h_Ele_fake_CRjetEt->Fill(pt, totalW);
      if ( pfType1Met > 20 )  h_Ele_fake_Mt_Met20->Fill(mt, totalW);     
      
      if ( AnalysisLeptons[0].type == Tight)  {

	h_Ele_signal_Mt->Fill(mt, totalW);
	h_Ele_signal_Met->Fill(pfType1Met, totalW);
	h_Ele_signal_pt->Fill(pt, totalW);
	h_Ele_signal_CRjetEt->Fill(pt, totalW);
	if ( pfType1Met > 20 )  h_Ele_signal_Mt_Met20->Fill(mt, totalW);     
      

      }

      if ( mt < 20 && pfType1Met < 20 ) { 

	  h_Ele_fake_pT_bin->Fill(pt, totalW);
	  h_Ele_fake_eta_bin->Fill(eta, totalW);
	  h_Ele_fake_pt_eta_bin->Fill(pt, eta , totalW); 
	  h_Ele_fake_nPV->Fill(nvtx, totalW); 
	  h_Ele_jetEt->Fill(jetEt, totalW);
	  h_Ele_fake_jetEt->Fill(jetEt, totalW);
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
	    h_Ele_signal_jetEt->Fill(jetEt, totalW);
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
			  fabs(std_vector_lepton_BestTrackdz->at(k))  < 0.1       && 
			  fabs(std_vector_lepton_BestTrackdxy->at(k)) < dxyCut );
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

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)   < 0.15);
  
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
			  fabs(std_vector_lepton_BestTrackdz->at(k))  < 0.1      && 
			  fabs(std_vector_lepton_BestTrackdxy->at(k)) < dxyCut); 
    }
  

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_flavour->at(k)) == 11)
    {
      // Faltan los aislamientos!!!!!!

      float aeta = fabs(std_vector_lepton_eta->at(k));
      /*
      if (aeta <= 1.479) {
	
	if (fabs(std_vector_electron_dEtaIn->at(k)) < 0.01               &&
	    fabs(std_vector_electron_dPhiIn->at(k)) < 0.04               &&
	    std_vector_electron_full5x5_sigmaIetaIeta->at(k)    < 0.011              &&
	    std_vector_electron_hOverE->at(k)              < 0.08               &&
	    std_vector_electron_ooEmooP->at(k) < 0.01               && //????
	    std_vector_electron_expectedMissingInnerHits->at(k)<=2  &&
	    std_vector_electron_d0->at(k)      < 0.1                &&
	    fabs(std_vector_electron_dz->at(k))< 0.373              &&
	    std_vector_electron_passConversionVeto )
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
	    std_vector_electron_d0->at(k)      < 0.2                &&
	    fabs(std_vector_electron_dz->at(k))< 0.602              &&
	    std_vector_electron_passConversionVeto )
	  {  
	    is_loose_lepton = true; 
	  }
	  }*/

    }
  
  return is_loose_lepton;
}



//------------------------------------------------------------------------------
// IsLooseIsolatedLepton
//------------------------------------------------------------------------------
bool IsLooseIsolatedLepton(int k)
{
  float id = std_vector_lepton_flavour->at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k) < 0.4);
  
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

