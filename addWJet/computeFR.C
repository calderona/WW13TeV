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
bool  comesFromW (int k);
bool  isGen (int k);
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

Float_t          jetRho;
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

vector<float>   *std_vector_lepton_flavour;
vector<float>   *std_vector_lepton_pt;
vector<float>   *std_vector_lepton_eta;
vector<float>   *std_vector_lepton_phi;
vector<float>   *std_vector_lepton_isTightMuon; 
vector<float>   *std_vector_lepton_eleIdMedium; 
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

vector<float>   *std_vector_leptonGen_pt;
vector<float>   *std_vector_leptonGen_eta;
vector<float>   *std_vector_leptonGen_phi;
vector<float>   *std_vector_leptonGen_pid;
vector<float>   *std_vector_leptonGen_status;
vector<float>   *std_vector_leptonGen_mpid;
vector<float>   *std_vector_leptonGen_mstatus;
vector<float>   *std_vector_trigger;

vector<float>   *std_vector_jet_pt;
vector<float>   *std_vector_jet_eta;
vector<float>   *std_vector_jet_phi;
vector<float>   *std_vector_lepton_closejet_pt;
vector<float>   *std_vector_lepton_closejet_eta;
vector<float>   *std_vector_lepton_closejet_phi;
vector<float>   *std_vector_lepton_closejet_PartonFlavour;
vector<float>   *std_vector_lepton_closejet_drlj;



//------------------------------------------------------------------------------
// MAIN function
//------------------------------------------------------------------------------


void computeFR (TString theSample="WJets", bool doGen = false) {

  // ------ root settings ---------
 
  TH1::SetDefaultSumw2();
  
  TString path = Form("/gpfs/csic_users/calderon/IFCA_2015/WW13TeV/addWJet/rootfiles/");

  TFile* root_output;  


  if (doGen) {
    root_output = new TFile(path + theSample + "_FR_FakeRate.root", "recreate");
  }  else {
    root_output =  new TFile(path + theSample + "_FakeRate.root", "recreate");
  }

  //-- Define histograms 
  //----------------------------------------------------------------------------

  // Define pt bins for FR matrix
  int pTbin = 7;
  Double_t pTbins[8] = { 15., 20., 25., 30., 35., 40., 45.0, 50.0};

  // Define eta bins for FR matrix
  int etabin = 5;
  Double_t etabins[6] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

  //Define pt bins for FR matrix
  int dxypTbin = 7;
  Double_t dxypTbins[8] = {0,0.002,0.004,0.006,0.008,0.012,0.016,0.020};


  // MUON histograms

  TH1F  *h_triggerBits = new TH1F("h_triggerBits", "h_triggerBits", 13, 0, 13 ); 

  TH1F  *h_Muon_fake_pT_bin   = new TH1F("h_Muon_fake_pT_bin","h_Muon_fake_pT_bin", pTbin,pTbins);
  TH1F  *h_Muon_signal_pT_bin = new TH1F("h_Muon_signal_pT_bin","h_Muon_signal_pT_bin",pTbin,pTbins);
  TH1F  *h_Muon_fake_eta_bin = new TH1F("h_Muon_fake_eta_bin","h_Muon_fake_eta_bin",etabin,etabins);
  TH1F  *h_Muon_signal_eta_bin = new TH1F("h_Muon_signal_eta_bin","h_Muon_signal_eta_bin",etabin,etabins);

  TH2F  *h_Muon_fake_pt_eta_bin = new TH2F("h_Muon_fake_pt_eta_bin","h_Muon_fake_pt_eta_bin",pTbin,pTbins,etabin,etabins);
  TH2F  *h_Muon_signal_pt_eta_bin = new TH2F("h_Muon_signal_pt_eta_bin","h_Muon_signal_pt_eta_bin",pTbin,pTbins,etabin,etabins);

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
  TH1F  *h_Muon_fake_pt   = new TH1F("h_Muon_fake_pt","h_Muon_fake_pt", 15,0,200);
  TH1F  *h_Muon_jetEt   = new TH1F("h_Muon_jetEt","h_Muon_jetEt", 100,0,200);

  TH1F  *h_Muon_signal_Mt   = new TH1F("h_Muon_signal_Mt","h_Muon_signal_Mt", 100,0,200);
  TH1F  *h_Muon_signal_Mt_Met20   = new TH1F("h_Muon_signal_Mt_Met20","h_Muon_signal_Mt_Met20", 100,0,200);
  TH1F  *h_Muon_signal_Met   = new TH1F("h_Muon_signal_Met","h_Muon_signal_Met", 100,0,200);
  TH1F  *h_Muon_signal_CRjetEt   = new TH1F("h_Muon_signal_CRjetEt","h_Muon_signal_CRjetEt", 100,0,200);
  TH1F  *h_Muon_signal_pt   = new TH1F("h_Muon_signal_pt","h_Muon_signal_pt", 15,0,200);
  
  TH1F  *h_Muon_fake_dxy_pt = new TH1F("h_Muon_fake_dxy_pt","h_Muon_fake_dxy_pt",50, 0, 0.3);
  TH1F  *h_Muon_signal_dxy_pt = new TH1F("h_Muon_signal_dxy_pt","h_Muon_fake_dxy_pt",50, 0, 0.3);

  TH1F  *h_Muon_fake_dxy = new TH1F("h_Muon_fake_dxy","h_Muon_fake_dxy",50, 0, 0.3);
  TH1F  *h_Muon_signal_dxy = new TH1F("h_Muon_signal_dxy","h_Muon_fake_dxy",50, 0, 0.3);

  TH1F  *h_Muon_fake_dxy_pt_bin = new TH1F("h_Muon_fake_dxy_pt_bin","h_Muon_fake_dxy_pt_bin",dxypTbin,dxypTbins);
  TH1F  *h_Muon_signal_dxy_pt_bin = new TH1F("h_Muon_signal_dxy_pt_bin","h_Muon_fake_dxy_pt_bin",dxypTbin,dxypTbins);

  TH1F  *h_Muon_fake_jetEt = new TH1F("h_Muon_fake_jetEt","h_Muon_fake_jetEt",10,10,150);
  TH1F  *h_Muon_signal_jetEt = new TH1F("h_Muon_signal_jetEt","h_Muon_signal_jetEt",10,10,150);

  TH1F  *h_Muon_fake_sip3d = new TH1F("h_Muon_fake_sip3d","h_Muon_fake_sip3d",50, 0, 40);
  TH1F  *h_Muon_signal_sip3d = new TH1F("h_Muon_signal_sip3d","h_Muon_fake_sip3d",50, 0, 40);

  TH1F  *h_Muon_fake_ptIso = new TH1F("h_Muon_fake_ptIso","h_Muon_fake_ptIso",pTbin,pTbins);
  TH1F  *h_Muon_signal_ptIso = new TH1F("h_Muon_signal_ptIso","h_Muon_fake_ptIso",pTbin,pTbins);


  // ELECTRON histograms
  TH1F  *h_Ele_fake_pT_bin   = new TH1F("h_Ele_fake_pT_bin","h_Ele_fake_pT_bin", pTbin,pTbins);
  TH1F  *h_Ele_signal_pT_bin = new TH1F("h_Ele_signal_pT_bin","h_Ele_signal_pT_bin",pTbin,pTbins);
  TH1F  *h_Ele_fake_eta_bin = new TH1F("h_Ele_fake_eta_bin","h_Ele_fake_eta_bin",etabin,etabins);
  TH1F  *h_Ele_signal_eta_bin = new TH1F("h_Ele_signal_eta_bin","h_Ele_signal_eta_bin",etabin,etabins);

  TH2F  *h_Ele_fake_pt_eta_bin = new TH2F("h_Ele_fake_pt_eta_bin","h_Ele_fake_pt_eta_bin",pTbin,pTbins,etabin,etabins);
  TH2F  *h_Ele_signal_pt_eta_bin = new TH2F("h_Ele_signal_pt_eta_bin","h_Ele_signal_pt_eta_bin",pTbin,pTbins,etabin,etabins);

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

  TH1F  *h_Ele_signal_Mt   = new TH1F("h_Ele_signal_Mt","h_Ele_signal_Mt", 100,0,200);
  TH1F  *h_Ele_signal_Met   = new TH1F("h_Ele_signal_Met","h_Ele_signal_Met", 100,0,200);
  TH1F  *h_Ele_signal_CRjetEt   = new TH1F("h_Ele_signal_CRjetEt","h_Ele_signal_CRjetEt", 100,0,200);
  TH1F  *h_Ele_signal_pt   = new TH1F("h_Ele_signal_pt","h_Ele_signal_pt", 100,15,200);
  TH1F  *h_Ele_signal_Mt_Met20   = new TH1F("h_Ele_signal_Mt_Met20","h_Ele_signal_Mt_Met20", 100,0,200);
  
  TH1F  *h_Ele_fake_jetEt = new TH1F("h_Ele_fake_jetEt","h_Ele_fake_jetEt",10,10,150);
  TH1F  *h_Ele_signal_jetEt = new TH1F("h_Ele_signal_jetEt","h_Ele_signal_jetEt",10,10,150);

  TH1F  *h_Ele_fake_sip3d = new TH1F("h_Ele_fake_sip3d","h_Ele_fake_sip3d",50, 0, 40);
  TH1F  *h_Ele_signal_sip3d = new TH1F("h_Ele_signal_sip3d","h_Ele_fake_sip3d",50, 0, 40);

  TH1F  *h_Ele_fake_ptIso = new TH1F("h_Ele_fake_ptIso","h_Ele_fake_ptIso",pTbin,pTbins);
  TH1F  *h_Ele_signal_ptIso = new TH1F("h_Ele_signal_ptIso","h_Ele_fake_ptIso",pTbin,pTbins);



  //--- OPEN LOOSE-LOOSE FILE 
  //----------------------------------------------------------------------------
  
  TString filesPath1, filesPath2, filesPath3;
  filesPath1 = "/gpfs/csic_projects/cms/piedra/latino/RunII/MC_Spring15/50ns/";
  filesPath2 = "/gpfs/csic_projects/cms/piedra/latino/RunII/Data13TeVRun2015B/dataRun2015B/";
  filesPath3 = "/gpfs/csic_projects/cms/piedra/latino/RunII/MC_Spring15/25ns/";

  TChain* tree = new TChain("latino", "latino");

  if (theSample == "WJets") {
    tree->Add(filesPath1 + "latino_WJetsToLNu.root");
  }
  else if (theSample == "QCD") {
    tree->Add(filesPath1 + "latino_QCD15to20Mu.root");
  }
  else if (theSample == "QCDEM") {
    tree->Add(filesPath1 + "latino_QCD15to20EM.root");
    tree->Add(filesPath1 + "latino_QCD20to30EM.root");
    tree->Add(filesPath1 + "latino_QCD30to50EM.root");  
    tree->Add(filesPath1 + "latino_QCD50to80EM.root");
    //tree->Add(filesPath1 + "latino_QCD30toInfDoubleEM.root");
  }
  else if (theSample == "ZJets") {
    //tree->Add(filesPath1 + "latino_DYJetsToLL1050.root");
    //tree->Add(filesPath1 + "latino_DYJetsToLL.root");
    tree->Add(filesPath1 + "latino_DYJetsToLL.root");
  }
  else if (theSample == "Top") {
    tree->Add(filesPath1 + "latino_TTJets.root");
  }
  else if (theSample == "SingleTop") {
    tree->Add(filesPath1 + "latino_ST_t-channel_top.root");
    tree->Add(filesPath1 + "latino_topTchannelAntitop.root");
  }
  else if (theSample == "dataEG") {
    tree->Add(filesPath2 + "latino_DoubleEG.root");
  }
  else {
    return;
  }


  tree->SetBranchAddress("baseW", &baseW);
  tree->SetBranchAddress("jetRho", &jetRho);
  tree->SetBranchAddress("channel", &channel);
  tree->SetBranchAddress("pfType1Met", &pfType1Met);
  tree->SetBranchAddress("pfType1Metphi", &pfType1Metphi);
  tree->SetBranchAddress("mth", &mth);
  tree->SetBranchAddress("dphilmet1", &dphilmet1);
  tree->SetBranchAddress("dphilmet2", &dphilmet2);
  tree->SetBranchAddress("nvtx",&nvtx);
  tree->SetBranchAddress("triggerFakeRate", &triggerFakeRate);

  tree->SetBranchAddress("std_vector_lepton_pt", &std_vector_lepton_pt);
  tree->SetBranchAddress("std_vector_lepton_eta", &std_vector_lepton_eta);
  tree->SetBranchAddress("std_vector_lepton_phi", &std_vector_lepton_phi);
  tree->SetBranchAddress("std_vector_lepton_flavour", &std_vector_lepton_flavour);
  tree->SetBranchAddress("std_vector_lepton_isTightMuon", &std_vector_lepton_isTightMuon);
  tree->SetBranchAddress("std_vector_lepton_eleIdMedium", &std_vector_lepton_eleIdMedium);
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


  tree->SetBranchAddress("std_vector_leptonGen_pt", &std_vector_leptonGen_pt);
  tree->SetBranchAddress("std_vector_leptonGen_eta", &std_vector_leptonGen_eta);
  tree->SetBranchAddress("std_vector_leptonGen_phi", &std_vector_leptonGen_phi);
  tree->SetBranchAddress("std_vector_leptonGen_pid", &std_vector_leptonGen_pid);
  tree->SetBranchAddress("std_vector_leptonGen_status", &std_vector_leptonGen_status);
  tree->SetBranchAddress("std_vector_leptonGen_mpid", &std_vector_leptonGen_mpid);
  tree->SetBranchAddress("std_vector_leptonGen_mstatus", &std_vector_leptonGen_mstatus);
  tree->SetBranchAddress("std_vector_trigger", &std_vector_trigger);

  tree->SetBranchAddress("std_vector_jet_pt", &std_vector_jet_pt);
  tree->SetBranchAddress("std_vector_jet_eta", &std_vector_jet_eta);
  tree->SetBranchAddress("std_vector_jet_phi", &std_vector_jet_phi);
  tree->SetBranchAddress("std_vector_lepton_closejet_pt", &std_vector_lepton_closejet_pt);
  tree->SetBranchAddress("std_vector_lepton_closejet_eta", &std_vector_lepton_closejet_eta);
  tree->SetBranchAddress("std_vector_lepton_closejet_phi", &std_vector_lepton_closejet_phi);
  tree->SetBranchAddress("std_vector_lepton_closejet_PartonFlavour", &std_vector_lepton_closejet_PartonFlavour);
  tree->SetBranchAddress("std_vector_lepton_closejet_drlj", &std_vector_lepton_closejet_drlj);


  //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------

  Long64_t nentries = tree->GetEntries();
  float nE =0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
    tree->GetEntry(jentry);
  
    if (jentry%10000==0) printf("Processing event %d \n", int(jentry));

    Float_t luminosity = 10.0; // 2fb-1

    Double_t totalW  = 1;//baseW * luminosity;

    nE++;

    // Help variables (only for deining loose+loose sample) 
    //--------------------------------------------------------------------------

    bool passLooseIDISO1 = IsLooseLepton(0) && IsLooseIsolatedLepton(0); 
    bool passLooseIDISO2 = IsLooseLepton(1) && IsLooseIsolatedLepton(1); 
    
    //if (!passLooseIDISO1) continue;
    //if (!passLooseIDISO2) continue;
    
  
 
    bool passIDISO1 = IsTightLepton(0) && IsIsolatedLepton(0); 
    bool passIDISO2 = IsTightLepton(1) && IsIsolatedLepton(1);

    
 
    // Loop over leptons
    //--------------------------------------------------------------------------

    if ( theSample == "Top" && doGen ) { 
     
      AnalysisGenLeptons.clear();        

      int vGenSize = std_vector_leptonGen_pt->size();                                                                  
    
      for (int i=0; i<vGenSize; i++) {  

	float pt  = std_vector_leptonGen_pt ->at(i);
	float eta = std_vector_leptonGen_eta->at(i);
	float phi = std_vector_leptonGen_phi->at(i);
	float id  = std_vector_leptonGen_pid ->at(i);      
     
	Lepton lep;                                           
	lep.index  = i;         
	lep.charge = id;        

	if( std_vector_leptonGen_status->at(i) != 1 ) continue;                 
	
	if (fabs(std_vector_leptonGen_pid->at(i)) != 11 && fabs(std_vector_leptonGen_pid->at(i) != 13 )) continue; 	
	if ( fabs(std_vector_leptonGen_mpid->at(i)) != 24 && fabs(std_vector_leptonGen_mpid->at(i)) != 15 ) continue; 

	if (fabs(id) == 11)          
	  { lep.flavor = Electron; }                          
        else if (fabs(id) == 13)        
	  { lep.flavor = Muon; } 
   
	lep.fromW = true;                      
	lep.iGen = true;      	        

	TLorentzVector tlv;   
	tlv.SetPtEtaPhiM(pt, eta, phi, 0);
	lep.v = tlv;                          

	AnalysisGenLeptons.push_back(lep);    

      }
    
      if ( AnalysisGenLeptons.size() !=1 ) continue; 
    }
    
    
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
     lep.fromW = false;
     lep.iGen = false; 

     TLorentzVector tlv;                                                                                                                      
     tlv.SetPtEtaPhiM(pt, eta, phi, 0);                                                                         
     lep.v = tlv; 
     
     if ( doGen )  { 
       //lep.fromW = comesFromW (i); 
       //lep.iGen = isGen(i);
       //float dR = tlv.DeltaR(AnalysisGenLeptons[0].v);     
       //if ( dR < 0.3 ) lep.fromW = true; 
     }
  
     AnalysisLeptons.push_back(lep);

    } // end loop on leptons 
    
  
    // do Gen information
    //--------------------------------------------------------------------------
   
 
    if ( (theSample == "WJets" || theSample == "SingleTop") && doGen ) {
	
      if (AnalysisLeptons.size() < 2 ) continue;

	  float pt = std_vector_lepton_pt->at(1);
	  float eta = fabs(std_vector_lepton_eta->at(1));
	  float dxy = fabs(std_vector_lepton_BestTrackdxy->at(1));
	  float jetEt = std_vector_lepton_closejet_pt->at(1);
	  float muSIP3D = std_vector_lepton_muSIP3D->at(1);
	  float elSIP3D = std_vector_lepton_elSIP3D->at(1);
	  float ptIso = pt*(1+MuonIsolation(1));

	  if (  AnalysisLeptons[0].type == Tight ) {
	    
	    if ( AnalysisLeptons[1].flavor == Muon ) { 
	    
	      h_Muon_fake_pT_bin->Fill(pt, totalW);
	      h_Muon_fake_dxy_pt->Fill((dxy/pt)*100 , totalW);
	      h_Muon_fake_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
	      h_Muon_jetEt->Fill(jetEt, totalW);
	      h_Muon_fake_dxy->Fill(dxy , totalW);
	      h_Muon_fake_jetEt->Fill(jetEt, totalW);
	      h_Muon_fake_eta_bin->Fill(eta, totalW);
	      h_Muon_fake_nPV->Fill(nvtx, totalW); 
	      h_Muon_fake_sip3d->Fill(muSIP3D, totalW);
	      h_Muon_fake_ptIso->Fill(ptIso, totalW);

	      if ( AnalysisLeptons[1].type == Tight)  {
	      
		h_Muon_signal_eta_bin->Fill(eta, totalW);
		h_Muon_signal_pT_bin->Fill(pt, totalW);
		h_Muon_signal_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
		h_Muon_signal_jetEt->Fill(jetEt, totalW);
		h_Muon_signal_nPV->Fill(nvtx, totalW); 
		h_Muon_signal_sip3d->Fill(muSIP3D, totalW);
		h_Muon_signal_ptIso->Fill(ptIso, totalW);
	      }
	    }

	    if ( AnalysisLeptons[1].flavor == Electron ) { 
	    
	      h_Ele_fake_pT_bin->Fill(pt, totalW);	  
	      h_Ele_jetEt->Fill(jetEt, totalW);
	      h_Ele_fake_jetEt->Fill(jetEt, totalW);
	      h_Ele_fake_eta_bin->Fill(eta, totalW);
	      h_Ele_fake_nPV->Fill(nvtx, totalW); 
	      h_Ele_fake_sip3d->Fill(elSIP3D, totalW);
	      h_Ele_fake_ptIso->Fill(ptIso, totalW);

	      if ( AnalysisLeptons[1].type == Tight)  {
	      
		h_Ele_signal_eta_bin->Fill(eta, totalW);
		h_Ele_signal_pT_bin->Fill(pt, totalW);
		h_Ele_signal_jetEt->Fill(jetEt, totalW);
		h_Ele_signal_nPV->Fill(nvtx, totalW); 
		h_Ele_signal_sip3d->Fill(elSIP3D, totalW);
		h_Ele_signal_ptIso->Fill(ptIso, totalW);

	      }
	    }	
	  }
	  
    }    
    //  }
    //} // END WJETS MC SAMPLE 


    if ( theSample == "Top"  && doGen) {

      if (AnalysisLeptons.size() < 2 ) continue;
	  
	  float pt = std_vector_lepton_pt->at(1);
	  float eta = fabs(std_vector_lepton_eta->at(1));
	  float dxy = fabs(std_vector_lepton_BestTrackdxy->at(1));
	  float jetEt = std_vector_lepton_closejet_pt->at(1);

	  if (  AnalysisLeptons[0].type == Tight ) {
	    
	    if ( AnalysisLeptons[1].flavor == Muon ) { 
	      
	      h_Muon_fake_pT_bin->Fill(pt, totalW);
	      h_Muon_fake_dxy_pt->Fill((dxy/pt)*100 , totalW);
	      h_Muon_fake_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
	      h_Muon_jetEt->Fill(jetEt, totalW);
	      h_Muon_fake_dxy->Fill(dxy , totalW);
	      h_Muon_fake_jetEt->Fill(jetEt, totalW);
	      h_Muon_fake_eta_bin->Fill(eta, totalW);
	      h_Muon_fake_nPV->Fill(nvtx, totalW); 
	      
	      if ( AnalysisLeptons[1].type == Tight)  {
	      
		h_Muon_signal_eta_bin->Fill(eta, totalW);
		h_Muon_signal_pT_bin->Fill(pt, totalW);
		h_Muon_signal_dxy_pt_bin->Fill((dxy/pt)*100 , totalW);
		h_Muon_signal_jetEt->Fill(jetEt, totalW);
		h_Muon_signal_nPV->Fill(nvtx, totalW); 

	      }
	    }	

	    if ( AnalysisLeptons[1].flavor == Electron ) { 
	    
	      h_Ele_fake_pT_bin->Fill(pt, totalW);	  
	      h_Ele_jetEt->Fill(jetEt, totalW);
	      h_Ele_fake_jetEt->Fill(jetEt, totalW);
	      h_Ele_fake_eta_bin->Fill(eta, totalW);
	      h_Ele_fake_nPV->Fill(nvtx, totalW); 
	      
	      if ( AnalysisLeptons[1].type == Tight)  {
	      
		h_Ele_signal_eta_bin->Fill(eta, totalW);
		h_Ele_signal_pT_bin->Fill(pt, totalW);
		h_Ele_signal_jetEt->Fill(jetEt, totalW);
		h_Ele_signal_nPV->Fill(nvtx, totalW); 

	      }
	    }
	  }
    } // END TOP MC SAMPLE 



    // Perform analysis
    //--------------------------------------------------------------------------  
    if ( theSample == "dataEG") {

      int iBit = std_vector_trigger->size();

      int nBits = 0; 

      for (int i=0; i < iBit; i++) {
	if ( std_vector_trigger->at(i) == 1) {
	  h_triggerBits->Fill(i);
	  nBits ++ ;
	  } 
      }

      //if ( std_vector_trigger->at(2) != 1 ) continue;
    }

  
   
    if ( AnalysisLeptons.size() != 1) continue;
 
    //if ( theSample == "QCD" && fabs(std_vector_lepton_closejet_PartonFlavour->at(0)) == 5) continue;

    int index = AnalysisLeptons[0].index; 

    double mt = giveMT(AnalysisLeptons[0].index, pfType1Met, pfType1Metphi); 

    AnalysisJets.clear();
    
    // --------> away jet 
    int jetsize = std_vector_jet_pt->size();
    
    int jetIndex = 0; 
    
    for ( int j=0; j<jetsize; j++) {
      
      Jet jet_; 
      
      jet_.index = j; 
      
      float pt = std_vector_jet_pt->at(j);
      float eta = std_vector_jet_eta->at(j);
      float phi = std_vector_jet_phi->at(j);
      
      if (pt < 0 ) continue;
      
      TLorentzVector jetvt;
      jetvt.SetPtEtaPhiM(pt, eta, phi, 0);
      
      jet_.v = jetvt;
      
      if ( jetvt.DeltaR(AnalysisLeptons[0].v) < 1 ) continue;
      
      AnalysisJets.push_back(jet_);
      
    }
    
      
    if ( AnalysisJets.size() < 1 ) continue;
    
    if ( std_vector_lepton_pt->at(index)  < 15 ) continue; 
    
    if (doGen) continue;
    
    //if ( fabs(std_vector_lepton_closejet_PartonFlavour->at(0)) != 5) continue;
    //if ( std_vector_lepton_closejet_drlj->at(0) > 0.3 ) continue;

    if ( AnalysisLeptons[0].flavor == Muon ) { 
	  
      float pt = std_vector_lepton_pt->at(index);
      float eta = fabs(std_vector_lepton_eta->at(index));
      float dxy = fabs(std_vector_lepton_BestTrackdxy->at(index));
      float jetEt = std_vector_lepton_closejet_pt->at(index);
      float muSIP3D = std_vector_lepton_muSIP3D->at(index);
      float elSIP3D = std_vector_lepton_elSIP3D->at(index);
      float ptIso = pt*(1+MuonIsolation(index));

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
	if ( pfType1Met > 20 )  h_Muon_signal_Mt_Met20->Fill(mt, totalW);     
      
      }


      //if ( mt < 20 && pfType1Met < 20 ) { 

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
	  }    
	  //}
    }
	
    if ( AnalysisLeptons[0].flavor == Electron ) { 
		 
      float pt = std_vector_lepton_pt->at(index);
      float eta = fabs(std_vector_lepton_eta->at(index));
      float dxy = fabs(std_vector_lepton_BestTrackdxy->at(index));
      float jetEt = std_vector_lepton_closejet_pt->at(index);
      float elSIP3D = std_vector_lepton_elSIP3D->at(index);
      float ptIso = pt*(1+MuonIsolation(index));

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

      //if ( mt < 20 && pfType1Met < 20 ) { 

	  h_Ele_fake_pT_bin->Fill(pt, totalW);
	  h_Ele_fake_eta_bin->Fill(eta, totalW);
	  h_Ele_fake_pt_eta_bin->Fill(pt, eta , totalW); 
	  h_Ele_fake_nPV->Fill(nvtx, totalW); 
	  h_Ele_jetEt->Fill(jetEt, totalW);
	  h_Ele_fake_jetEt->Fill(jetEt, totalW);
	  h_Ele_pt->Fill(pt, totalW);
	  h_Ele_fake_sip3d->Fill(elSIP3D, totalW);
	  h_Ele_fake_ptIso->Fill(ptIso, totalW);


	  if ( AnalysisLeptons[0].type == Tight)  {
	    
	    h_Ele_signal_pT_bin->Fill(pt, totalW);
	    h_Ele_signal_eta_bin->Fill(eta, totalW);
	    h_Ele_signal_pt_eta_bin->Fill(pt, eta , totalW); 
	    h_Ele_signal_nPV->Fill(nvtx, totalW); 
	    h_Ele_signal_jetEt->Fill(jetEt, totalW);
	    h_Ele_signal_sip3d->Fill(elSIP3D, totalW);
	    h_Ele_signal_ptIso->Fill(ptIso, totalW);
	  } 
	  //}
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

      is_tight_lepton = ( std_vector_lepton_isTightMuon->at(k)              && 
			  fabs(std_vector_lepton_BestTrackdz->at(k))  < 0.1       && 
			  fabs(std_vector_lepton_BestTrackdxy->at(k)) < dxyCut );
	}

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_flavour->at(k)) == 11)
    {
      is_tight_lepton = std_vector_lepton_eleIdMedium->at(k);
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
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)   < 0.12);
  
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
     
      is_loose_lepton = ( std_vector_lepton_isTightMuon->at(k)                   && 
			  fabs(std_vector_lepton_BestTrackdz->at(k))  < 0.1      && 
			  fabs(std_vector_lepton_BestTrackdxy->at(k)) < 0.2 ); 
    }


  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_flavour->at(k)) == 11)
    {
      is_loose_lepton = std_vector_lepton_eleIdVeto->at(k);
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
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k) < 0.5);
  
  return is_isolated_lepton;
}


//------------------------------------------------------------------------------
// Match to a gen lepton and identify if comes from a W. 
//------------------------------------------------------------------------------
bool isGen (int k) { 
   
  bool isGen = false;

  float dRmin = 999.;
  int genIndex = -1; 
  
  TLorentzVector lep;
  lep.SetPtEtaPhiM(std_vector_lepton_pt->at(k), std_vector_lepton_eta->at(k), std_vector_lepton_phi->at(k), 0);

  int vGenSize = std_vector_leptonGen_pt->size();
  
  for (int i=0; i<vGenSize; i++) {
    
    if( std_vector_leptonGen_status->at(i) != 1 ) continue;

    if ( fabs(std_vector_lepton_flavour ->at(k)) != fabs(std_vector_leptonGen_pid->at(i)) ) continue;
    
    TLorentzVector genLep;
    genLep.SetPtEtaPhiM(std_vector_leptonGen_pt->at(i), std_vector_leptonGen_eta->at(i), std_vector_leptonGen_phi->at(i), 0);
    
    float dR = lep.DeltaR(genLep); 
       
    if ( dR < 0.3 )  { 
      genIndex = i; 
      isGen = true; 
    }
   
  }

  return isGen;

}

//------------------------------------------------------------------------------
// Match to a gen lepton and identify if comes from a W. 
//------------------------------------------------------------------------------
bool comesFromW (int k) { 
  
  bool isW = false;
 
  float dRmin = 999.;
  int genIndex = -1; 
  
  TLorentzVector lep;
  lep.SetPtEtaPhiM(std_vector_lepton_pt->at(k), std_vector_lepton_eta->at(k), std_vector_lepton_phi->at(k), 0);

  int vGenSize = std_vector_leptonGen_pt->size();
  
  for (int i=0; i<vGenSize; i++) {

    if( std_vector_leptonGen_status->at(i) != 1 ) continue;

    if ( fabs(std_vector_lepton_flavour ->at(k)) != fabs(std_vector_leptonGen_pid->at(i)) ) continue;
    
    TLorentzVector genLep;
    genLep.SetPtEtaPhiM(std_vector_leptonGen_pt->at(i), std_vector_leptonGen_eta->at(i), std_vector_leptonGen_phi->at(i), 0);
    
    float dR = lep.DeltaR(genLep); 
       
    if ( dR < 0.3 )  genIndex = i; 
   
  }
  
  if ( genIndex > 0 && fabs(std_vector_leptonGen_mpid->at(genIndex)) == 24 ) isW = true;

  
  int tau_pid = -999;
  int tau_status = -999;

  if ( genIndex > 0 && fabs(std_vector_leptonGen_mpid->at(genIndex)) == 15 )  {
    tau_pid = std_vector_leptonGen_mpid->at(genIndex); // can be +15 or -15
    tau_status = std_vector_leptonGen_mstatus->at(genIndex);
  }

  bool tau_from_W = false;

  for (int j=0; j<vGenSize; j++) {
    if (std_vector_leptonGen_pid->at(j) != tau_pid) continue;
    if (std_vector_leptonGen_status->at(j) != tau_status) continue;
    tau_from_W = (fabs(std_vector_leptonGen_mpid->at(j)) == 24);
  }

  if ( tau_from_W ) isW = true;
  
  return isW;
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

