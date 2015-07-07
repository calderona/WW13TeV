#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>


const Double_t ZMASS = 91.1876;


//------------------------------------------------------------------------------
// Member functions
//------------------------------------------------------------------------------

bool  IsIsolatedLepton(int k);

float ElectronIsolation(int k);

float MuonIsolation(int k);

bool  IsTightLepton(int k);

bool  IsLooseLepton(int k);

bool IsLooseIsolatedLepton(int k);


//------------------------------------------------------------------------------
// Define structures & GLB variables 
//------------------------------------------------------------------------------

enum {Muon, Electron};

enum {Loose, Tight};

//std::vector<Lepton> AnalysisLeptons;

Float_t fakeW; Float_t puW;

Float_t          jetRho;
Float_t          channel;
Float_t          baseW;
Float_t          dataset;     
Float_t          effW;         
Float_t          triggW;      
Float_t          dphilljetjet; 
Float_t          pfType1Met;  
Float_t          dphill;      
Float_t          dphilljet;   
Float_t          ptll;   
Float_t          mll;   
Float_t          mth;   
Float_t          dphilmet1;
Float_t          dphilmet2;
Float_t          trkMet;
Float_t          drll;
Float_t          njet;
vector<float>   *std_vector_lepton_id;
vector<float>   *std_vector_lepton_pt;
vector<float>   *std_vector_lepton_eta;
vector<float>   *std_vector_lepton_phi;
vector<float>   *std_vector_lepton_isTightMuon; 
vector<float>   *std_vector_lepton_eleIdMedium; 
vector<float>   *std_vector_lepton_chargedHadronIso; 
vector<float>   *std_vector_lepton_photonIso;
vector<float>   *std_vector_lepton_neutralHadronIso;
vector<float>   *std_vector_lepton_sumPUPt;
vector<float>   *std_vector_electron_effectiveArea;
vector<float>   *std_vector_lepton_BestTrackdz;
vector<float>   *std_vector_lepton_BestTrackdxy;


//------------------------------------------------------------------------------
// LatinosTreeScriptLL
//------------------------------------------------------------------------------
void LatinosTreeScriptLL(Float_t luminosity,
			 Int_t   jetChannel,
		         TString flavorChannel,
		         TString theSample,
		         Bool_t  verbose)
{

  TH1::SetDefaultSumw2();

  TString path = Form("rootfiles/%djet/%s/", jetChannel, flavorChannel.Data());

  gSystem->mkdir(path, kTRUE);

  TFile* output = new TFile(path + theSample + ".root", "recreate");


  // Counting histograms
  //----------------------------------------------------------------------------
  TH1F* hWTrigger     = new TH1F("hWTrigger",     "", 3, 0, 3);
  TH1F* hWMetCut      = new TH1F("hWMetCut",      "", 3, 0, 3);
  TH1F* hWLowMinv     = new TH1F("hWLowMinv",     "", 3, 0, 3);
  TH1F* hWZVeto       = new TH1F("hWZVeto",       "", 3, 0, 3);
  TH1F* hWpMetCut     = new TH1F("hWpMetCut",     "", 3, 0, 3);
  TH1F* hWJetVeto     = new TH1F("hWJetVeto",     "", 3, 0, 3);
  TH1F* hWDeltaPhiJet = new TH1F("hWDeltaPhiJet", "", 3, 0, 3);
  TH1F* hWSoftMuVeto  = new TH1F("hWSoftMuVeto",  "", 3, 0, 3); 
  TH1F* hWExtraLepton = new TH1F("hWExtraLepton", "", 3, 0, 3);
  TH1F* hWPtll        = new TH1F("hWPtll",        "", 3, 0, 3);
  TH1F* hWTopTagging  = new TH1F("hWTopTagging",  "", 3, 0, 3);

  TH1F* hWeffTrigger     = new TH1F("hWeffTrigger",     "", 3, 0, 3);
  TH1F* hWeffMetCut      = new TH1F("hWeffMetCut",      "", 3, 0, 3);
  TH1F* hWeffLowMinv     = new TH1F("hWeffLowMinv",     "", 3, 0, 3);
  TH1F* hWeffZVeto       = new TH1F("hWeffZVeto",       "", 3, 0, 3);
  TH1F* hWeffpMetCut     = new TH1F("hWeffpMetCut",     "", 3, 0, 3);
  TH1F* hWeffJetVeto     = new TH1F("hWeffJetVeto",     "", 3, 0, 3);
  TH1F* hWeffDeltaPhiJet = new TH1F("hWeffDeltaPhiJet", "", 3, 0, 3);
  TH1F* hWeffSoftMuVeto  = new TH1F("hWeffSoftMuVeto",  "", 3, 0, 3); 
  TH1F* hWeffExtraLepton = new TH1F("hWeffExtraLepton", "", 3, 0, 3);
  TH1F* hWeffPtll        = new TH1F("hWeffPtll",        "", 3, 0, 3);
  TH1F* hWeffTopTagging  = new TH1F("hWeffTopTagging",  "", 3, 0, 3);


  // WW level histograms: Region A, from tight -to- loose method
  //----------------------------------------------------------------------------
  TH1F* hPtLepton1WWLevelLL       = new TH1F("hPtLepton1WWLevelLL",       "", 200, 0, 200);
  TH1F* hPtLepton2WWLevelLL       = new TH1F("hPtLepton2WWLevelLL",       "", 200, 0, 200);
  TH1F* hPtDiLeptonWWLevelLL      = new TH1F("hPtDiLeptonWWLevelLL",      "", 200, 0, 200);
  TH1F* hMinvWWLevelLL            = new TH1F("hMinvWWLevelLL",            "", 200, 0, 200);
  TH1F* hMtWWLevelLL              = new TH1F("hMtWWLevelLL",              "", 250, 0, 250);
  TH1F* hNJets30WWLevelLL         = new TH1F("hNJetsPF30WWLevelLL",       "",  10, 0,  10);
  TH1F* hpfMetWWLevelLL           = new TH1F("hpfMetWWLevelLL",           "", 150, 0, 150);
  TH1F* hppfMetWWLevelLL          = new TH1F("hppfMetWWLevelLL",          "", 150, 0, 150);
  TH1F* hchMetWWLevelLL           = new TH1F("hchMetWWLevelLL",           "", 150, 0, 150);
  TH1F* hpchMetWWLevelLL          = new TH1F("hpchMetWWLevelLL",          "", 150, 0, 150);
  TH1F* hpminMetWWLevelLL         = new TH1F("hpminMetWWLevelLL",         "", 150, 0, 150);
  TH1F* hDeltaRLeptonsWWLevelLL   = new TH1F("hDeltaRLeptonsWWLevelLL",   "",  50, 0,   5);
  TH1F* hDeltaPhiLeptonsWWLevelLL = new TH1F("hDeltaPhiLeptonsWWLevelLL", "",  32, 0, 3.2);
  TH1F* hDPhiPtllJetWWLevelLL     = new TH1F("hDPhiPtllJetWWLevelLL",     "",  32, 0, 3.2);

  TH1F* hFakeLeptonIsoWWLevelLL = new TH1F("hFakeLeptonIsoWWLevelLL", "", 100, 0, 1);


  // WW level histograms: Region A, Tight+Tight
  //----------------------------------------------------------------------------
  TH1F* hPtLepton1WWLevelA       = new TH1F("hPtLepton1WWLevelA",       "", 200, 0, 200);
  TH1F* hPtLepton2WWLevelA       = new TH1F("hPtLepton2WWLevelA",       "", 200, 0, 200);
  TH1F* hPtDiLeptonWWLevelA      = new TH1F("hPtDiLeptonWWLevelA",      "", 200, 0, 200);
  TH1F* hMinvWWLevelA            = new TH1F("hMinvWWLevelA",            "", 200, 0, 200);
  TH1F* hMtWWLevelA              = new TH1F("hMtWWLevelA",              "", 250, 0, 250);
  TH1F* hNJets30WWLevelA         = new TH1F("hNJetsPF30WWLevelA",       "",  10, 0,  10);
  TH1F* hpfMetWWLevelA           = new TH1F("hpfMetWWLevelA",           "", 150, 0, 150);
  TH1F* hppfMetWWLevelA          = new TH1F("hppfMetWWLevelA",          "", 150, 0, 150);
  TH1F* hchMetWWLevelA           = new TH1F("hchMetWWLevelA",           "", 150, 0, 150);
  TH1F* hpchMetWWLevelA          = new TH1F("hpchMetWWLevelA",          "", 150, 0, 150);
  TH1F* hpminMetWWLevelA         = new TH1F("hpminMetWWLevelA",         "", 150, 0, 150);
  TH1F* hDeltaRLeptonsWWLevelA   = new TH1F("hDeltaRLeptonsWWLevelA",   "",  50, 0,   5);
  TH1F* hDeltaPhiLeptonsWWLevelA = new TH1F("hDeltaPhiLeptonsWWLevelA", "",  32, 0, 3.2);
  TH1F* hDPhiPtllJetWWLevelA     = new TH1F("hDPhiPtllJetWWLevelA",     "",  32, 0, 3.2);

  TH1F* hFakeLeptonIsoWWLevelA = new TH1F("hFakeLeptonIsoWWLevelA", "", 100, 0, 1);


  // WW level histograms: Region B, Tight+NoTight
  //----------------------------------------------------------------------------
  TH1F* hPtLepton1WWLevelB       = new TH1F("hPtLepton1WWLevelB",       "", 200, 0, 200);
  TH1F* hPtLepton2WWLevelB       = new TH1F("hPtLepton2WWLevelB",       "", 200, 0, 200);
  TH1F* hPtDiLeptonWWLevelB      = new TH1F("hPtDiLeptonWWLevelB",      "", 200, 0, 200);
  TH1F* hMinvWWLevelB            = new TH1F("hMinvWWLevelB",            "", 200, 0, 200);
  TH1F* hMtWWLevelB              = new TH1F("hMtWWLevelB",              "", 250, 0, 250);
  TH1F* hNJets30WWLevelB         = new TH1F("hNJetsPF30WWLevelB",       "",  10, 0,  10);
  TH1F* hpfMetWWLevelB           = new TH1F("hpfMetWWLevelB",           "", 150, 0, 150);
  TH1F* hppfMetWWLevelB          = new TH1F("hppfMetWWLevelB",          "", 150, 0, 150);
  TH1F* hchMetWWLevelB           = new TH1F("hchMetWWLevelB",           "", 150, 0, 150);
  TH1F* hpchMetWWLevelB          = new TH1F("hpchMetWWLevelB",          "", 150, 0, 150);
  TH1F* hpminMetWWLevelB         = new TH1F("hpminMetWWLevelB",         "", 150, 0, 150);
  TH1F* hDeltaRLeptonsWWLevelB   = new TH1F("hDeltaRLeptonsWWLevelB",   "",  50, 0,   5);
  TH1F* hDeltaPhiLeptonsWWLevelB = new TH1F("hDeltaPhiLeptonsWWLevelB", "",  32, 0, 3.2);
  TH1F* hDPhiPtllJetWWLevelB     = new TH1F("hDPhiPtllJetWWLevelB",     "",  32, 0, 3.2);

  TH1F* hFakeLeptonIsoWWLevelB = new TH1F("hFakeLeptonIsoWWLevelB", "", 100, 0, 1);
  TH1F* hpfMetWWLevelBNoMET = new TH1F("hppfMetWWLevelBNoMET",          "", 150, 0, 150);

  // WW level histograms: Region C, NoTight+NoTight
  //----------------------------------------------------------------------------
  TH1F* hPtLepton1WWLevelC       = new TH1F("hPtLepton1WWLevelC",       "", 200, 0, 200);
  TH1F* hPtLepton2WWLevelC       = new TH1F("hPtLepton2WWLevelC",       "", 200, 0, 200);
  TH1F* hPtDiLeptonWWLevelC      = new TH1F("hPtDiLeptonWWLevelC",      "", 200, 0, 200);
  TH1F* hMinvWWLevelC            = new TH1F("hMinvWWLevelC",            "", 200, 0, 200);
  TH1F* hMtWWLevelC              = new TH1F("hMtWWLevelC",              "", 250, 0, 250);
  TH1F* hNJets30WWLevelC         = new TH1F("hNJetsPF30WWLevelC",       "",  10, 0,  10);
  TH1F* hpfMetWWLevelC           = new TH1F("hpfMetWWLevelC",           "", 150, 0, 150);
  TH1F* hppfMetWWLevelC          = new TH1F("hppfMetWWLevelC",          "", 150, 0, 150);
  TH1F* hchMetWWLevelC           = new TH1F("hchMetWWLevelC",           "", 150, 0, 150);
  TH1F* hpchMetWWLevelC          = new TH1F("hpchMetWWLevelC",          "", 150, 0, 150);
  TH1F* hpminMetWWLevelC         = new TH1F("hpminMetWWLevelC",         "", 150, 0, 150);
  TH1F* hDeltaRLeptonsWWLevelC   = new TH1F("hDeltaRLeptonsWWLevelC",   "",  50, 0,   5);
  TH1F* hDeltaPhiLeptonsWWLevelC = new TH1F("hDeltaPhiLeptonsWWLevelC", "",  32, 0, 3.2);
  TH1F* hDPhiPtllJetWWLevelC     = new TH1F("hDPhiPtllJetWWLevelC",     "",  32, 0, 3.2);

  TH1F* hFakeLeptonIsoWWLevelC = new TH1F("hFakeLeptonIsoWWLevelC", "", 100, 0, 1);


   //----------------------------------------------------------------------------
  // Input files
  //----------------------------------------------------------------------------
  TString filesPath;

  filesPath = "/gpfs/csic_projects/cms/piedra/latino/25ns/";
		 

  TChain* tree = new TChain("latino", "latino");

  double mybaseW = 1.0;

  if (theSample == "WJetsFakes_Total") {
    tree->Add("addWJet/output/latino_WJetsToLNu.root");
  }
  else if (theSample == "WWTo2L2Nu_pow") {
    tree->Add(filesPath + "latino_WWTo2L2Nu.root");
  }
  else if (theSample == "ZZ") {
    tree->Add(filesPath + "latino_ZZ.root");
  }
  else {
    return;
  }


  // Declaration of leaf types
  //----------------------------------------------------------------------------
  tree->SetBranchAddress("baseW"       , &baseW);
  tree->SetBranchAddress("channel"     , &channel);
  tree->SetBranchAddress("dataset"     , &dataset);
  tree->SetBranchAddress("effW"        , &effW);
  tree->SetBranchAddress("triggW"      , &triggW);
  tree->SetBranchAddress("dphilljetjet", &dphilljetjet);
  tree->SetBranchAddress("pfType1Met"  , &pfType1Met);
  tree->SetBranchAddress("dphill"      , &dphill);
  tree->SetBranchAddress("dphilljet"   , &dphilljet);
  tree->SetBranchAddress("jetRho"      , &jetRho);
  tree->SetBranchAddress("channel"     , &channel);
  tree->SetBranchAddress("baseW"       , &baseW);
  tree->SetBranchAddress("ptll"        , &ptll);
  tree->SetBranchAddress("mll"         , &mll);
  tree->SetBranchAddress("mth"         , &mth);
  tree->SetBranchAddress("dphilmet1"   , &dphilmet1);
  tree->SetBranchAddress("dphilmet2"   , &dphilmet2);
  tree->SetBranchAddress("trkMet"      , &trkMet);
  tree->SetBranchAddress("njet"      , &njet);
 

  tree->SetBranchAddress("std_vector_lepton_pt", &std_vector_lepton_pt);
  tree->SetBranchAddress("std_vector_lepton_eta", &std_vector_lepton_eta);
  tree->SetBranchAddress("std_vector_lepton_phi", &std_vector_lepton_phi);
  tree->SetBranchAddress("std_vector_lepton_id", &std_vector_lepton_id);
  tree->SetBranchAddress("std_vector_lepton_isTightMuon", &std_vector_lepton_isTightMuon);
  tree->SetBranchAddress("std_vector_lepton_eleIdMedium", &std_vector_lepton_eleIdMedium);
  tree->SetBranchAddress("std_vector_lepton_chargedHadronIso", &std_vector_lepton_chargedHadronIso);
  tree->SetBranchAddress("std_vector_lepton_photonIso", &std_vector_lepton_photonIso);
  tree->SetBranchAddress("std_vector_lepton_neutralHadronIso", &std_vector_lepton_neutralHadronIso);
  tree->SetBranchAddress("std_vector_lepton_sumPUPt", &std_vector_lepton_sumPUPt);
  tree->SetBranchAddress("std_vector_electron_effectiveArea", &std_vector_electron_effectiveArea);
  tree->SetBranchAddress("std_vector_lepton_BestTrackdz", &std_vector_lepton_BestTrackdz);
  tree->SetBranchAddress("std_vector_lepton_BestTrackdxy", &std_vector_lepton_BestTrackdxy);

  tree->SetBranchAddress("fakeW", &fakeW);
  tree->SetBranchAddress("puW", &puW);
  

  // Set the channel
  //----------------------------------------------------------------------------
  Float_t SelectedChannel = -999;

  if      (flavorChannel == "MuMu") SelectedChannel =  0;
  else if (flavorChannel == "EE"  ) SelectedChannel =  1;
  else if (flavorChannel == "EMu" ) SelectedChannel =  2;
  else if (flavorChannel == "MuE" ) SelectedChannel =  3;
  else if (flavorChannel == "OF" )  SelectedChannel =  4;
  else if (flavorChannel == "SF" )  SelectedChannel =  5;
  
 
  
  //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------
  for (int ievent=0; ievent<tree->GetEntries(); ievent++) {
    
    tree->GetEntry(ievent);
    
    Double_t efficiencyW = effW * triggW;
    Double_t totalW      = -999;

    totalW = 1.0;


    bool isSemi = 0; 
  
    /*
    
    if (theSample.Contains("TopFakes")) {
      
      if ( leptonGenpid1 > 0 &&
	   leptonGenpid2 < 0  ) isSemi = 1; 
      
	   }*/	   
	    
   


   
    // Help variables
    //--------------------------------------------------------------------------
  
    float pt1  = std_vector_lepton_pt ->at(0); float pt2  = std_vector_lepton_pt ->at(1);
    float eta1 = std_vector_lepton_eta->at(0); float eta2 = std_vector_lepton_eta->at(1);
    float ch1  = std_vector_lepton_id->at(0);  float ch2  = std_vector_lepton_id->at(1);
    bool passLooseIDISO1 = IsLooseLepton(0) && IsLooseIsolatedLepton(0); 
    bool passLooseIDISO2 = IsLooseLepton(1) && IsLooseIsolatedLepton(1); 
    bool passIDISO1 = IsTightLepton(0) && IsIsolatedLepton(0); 
    bool passIDISO2 = IsTightLepton(1) && IsIsolatedLepton(1); 

    Float_t dphimin = (min(dphilmet1,dphilmet2));
    Float_t fullpmet = 0;
    Float_t trkpmet = 0;

    if (dphimin < TMath::Pi() / 2)
      fullpmet = pfType1Met * sin(dphimin);
    else 
      fullpmet = pfType1Met;

    if (dphimin < TMath::Pi() / 2)
      trkpmet = trkMet * sin(dphimin);
    else
      trkpmet = trkMet;
    
    Float_t mpmet = min(trkpmet,fullpmet);

 

    // The selection begins here
    //--------------------------------------------------------------------------
    
    if (pt2 <= 10) continue; // increase the pt of the leptons to further reduce Wjets 
    if (pt1 <= 10) continue; // increase the pt of the leptons to further reduce Wjets 
    //if (ch1*ch2 > 0) continue;
    if (!passLooseIDISO1) continue;
    if (!passLooseIDISO2) continue;

    if ( (SelectedChannel == -1)                                   || 
	 (channel == SelectedChannel)                              || 
	 (SelectedChannel == 4 && (channel == 2 || channel == 3) ) || 
	 (SelectedChannel == 5 && (channel == 0 || channel == 1) ) 
	 ) {
     

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      //
      // Main analisis
      //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
 
      //int preSel = (nextra == 0) && mll > 12 && (zveto==1 || !sameflav) && (mpmet > 20  && dyMVA ) && (dphiv || !sameflav) && bveto_mu &&  ( ptll>30 && (!sameflav || ptll>45)) && (jetbin == jetChannel) && (bveto_ip==1 &&  nbjettche==0) ;


     
      if ( true) { 

     	//if ( pfmet > 20 )  {	  

	  hPtLepton1WWLevelLL      ->Fill(pt1,          totalW*fakeW);
	  hPtLepton2WWLevelLL      ->Fill(pt2,          totalW*fakeW);
	  hPtDiLeptonWWLevelLL     ->Fill(ptll,         totalW*fakeW);
	  hMinvWWLevelLL           ->Fill(mll,          totalW*fakeW);
	  hMtWWLevelLL             ->Fill(mth,          totalW*fakeW);
	  hNJets30WWLevelLL        ->Fill(njet,         totalW*fakeW);
	  hpfMetWWLevelLL          ->Fill(pfType1Met,   totalW*fakeW);
	  //hppfMetWWLevelLL         ->Fill(ppfmet,       totalW*fakeW);
	  //hchMetWWLevelLL          ->Fill(chmet,        totalW*fakeW);
	  //hpchMetWWLevelLL         ->Fill(pchmet,       totalW*fakeW);
	  hpminMetWWLevelLL        ->Fill(mpmet,        totalW*fakeW);
	  hDeltaRLeptonsWWLevelLL  ->Fill(drll,         totalW*fakeW);
	  hDeltaPhiLeptonsWWLevelLL->Fill(dphill,       totalW*fakeW);
	  hDPhiPtllJetWWLevelLL    ->Fill(dphilljet,    totalW*fakeW);
	  //}

	  if ( passIDISO1 && passIDISO2 ) { //&&  isSemi == 1) { //&& pfmet > 20 )  {
     
	    hPtLepton1WWLevelA      ->Fill(pt1,       totalW);
	    hPtLepton2WWLevelA      ->Fill(pt2,       totalW);
	    hPtDiLeptonWWLevelA     ->Fill(ptll,      totalW);
	    hMinvWWLevelA           ->Fill(mll,       totalW);
	    hMtWWLevelA             ->Fill(mth,       totalW);
	    hNJets30WWLevelA        ->Fill(njet,      totalW);
	    hpfMetWWLevelA          ->Fill(pfType1Met,     totalW);
	    //hppfMetWWLevelA         ->Fill(ppfmet,    totalW);
	    //hchMetWWLevelA          ->Fill(chmet,     totalW);
	    //hpchMetWWLevelA         ->Fill(pchmet,    totalW);
	    hpminMetWWLevelA        ->Fill(mpmet,     totalW);
	    hDeltaRLeptonsWWLevelA  ->Fill(drll,      totalW);
	    hDeltaPhiLeptonsWWLevelA->Fill(dphill,    totalW);
	    hDPhiPtllJetWWLevelA    ->Fill(dphilljet, totalW);
	    
	  }
	
	  if ( (passIDISO1 && !passIDISO2) || (!passIDISO1 && passIDISO2) ) {
	    
	    hpfMetWWLevelBNoMET          ->Fill(pfType1Met,     totalW); 
	    
	    hPtLepton1WWLevelB      ->Fill(pt1,       totalW);
	    hPtLepton2WWLevelB      ->Fill(pt2,       totalW);
	    hPtDiLeptonWWLevelB     ->Fill(ptll,      totalW);
	    hMinvWWLevelB           ->Fill(mll,       totalW);
	    hMtWWLevelB             ->Fill(mth,       totalW);
	    hNJets30WWLevelB        ->Fill(njet,      totalW);
	    hpfMetWWLevelB          ->Fill(pfType1Met,     totalW);
	    //hppfMetWWLevelB         ->Fill(ppfmet,    totalW);
	    //hchMetWWLevelB          ->Fill(chmet,     totalW);
	    //hpchMetWWLevelB         ->Fill(pchmet,    totalW);
	    hpminMetWWLevelB        ->Fill(mpmet,     totalW);
	    hDeltaRLeptonsWWLevelB  ->Fill(drll,      totalW);
	    hDeltaPhiLeptonsWWLevelB->Fill(dphill,    totalW);
	    hDPhiPtllJetWWLevelB    ->Fill(dphilljet, totalW);
	  }	
	
	  if ( !passIDISO1 && !passIDISO2)  { // && pfmet > 20 ) { 
	    
	    hPtLepton1WWLevelC      ->Fill(pt1,       totalW);
	    hPtLepton2WWLevelC      ->Fill(pt2,       totalW);
	    hPtDiLeptonWWLevelC     ->Fill(ptll,      totalW);
	    hMinvWWLevelC           ->Fill(mll,       totalW);
	    hMtWWLevelC             ->Fill(mth,       totalW);
	    hNJets30WWLevelC        ->Fill(njet,      totalW);
	    hpfMetWWLevelC          ->Fill(pfType1Met,     totalW);
	    //hppfMetWWLevelC         ->Fill(ppfmet,    totalW);
	    //hchMetWWLevelC          ->Fill(chmet,     totalW);
	    //hpchMetWWLevelC         ->Fill(pchmet,    totalW);
	    hpminMetWWLevelC        ->Fill(mpmet,     totalW);
	    hDeltaRLeptonsWWLevelC  ->Fill(drll,      totalW);
	    hDeltaPhiLeptonsWWLevelC->Fill(dphill,    totalW);
	    hDPhiPtllJetWWLevelC    ->Fill(dphilljet, totalW);	    
	  
	  }
		
      }
      
    }
    
  }


	  
  verbose = true;

  // Print
  //----------------------------------------------------------------------------
  if (verbose) {

    float norm = 1933235;
  
  }
 
   

  // Save the histograms
  //----------------------------------------------------------------------------
  output->cd();
  output->Write("", TObject::kOverwrite);
  output->Close();
  
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
  if (fabs(std_vector_lepton_id->at(k)) == 13)
    {
      bool dxyCut = 0;
	
      if ( std_vector_lepton_pt->at(k) < 20 ) { 
	dxyCut = 0.01;
      }	else {
	dxyCut = 0.02;
      }

      is_tight_lepton = ( std_vector_lepton_isTightMuon->at(k)              && 
			  std_vector_lepton_BestTrackdz->at(k)  < 0.1       && 
			  std_vector_lepton_BestTrackdxy->at(k) < dxyCut );
	}

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_id->at(k)) == 11)
    {
      is_tight_lepton = true;//std_vector_lepton_eleIdMedium->at(k);
    }
  
  return is_tight_lepton;
}


//------------------------------------------------------------------------------
// MuonIsolation
//------------------------------------------------------------------------------
float MuonIsolation(int k)
{
  float pt = std_vector_lepton_pt->at(k);
  float id = std_vector_lepton_id->at(k);

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
  float id = std_vector_lepton_id->at(k);

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
  float id = std_vector_lepton_id->at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.12);
  
  return is_isolated_lepton;
}


//------------------------------------------------------------------------------
// IsLooseIsolatedLepton
//------------------------------------------------------------------------------
bool IsLooseIsolatedLepton(int k)
{
  float id = std_vector_lepton_id->at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.4);
  
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
  if (fabs(std_vector_lepton_id->at(k)) == 13)
    {
     
      is_loose_lepton = ( std_vector_lepton_isTightMuon->at(k)              && 
			  std_vector_lepton_BestTrackdz->at(k)  < 0.1       && 
			  std_vector_lepton_BestTrackdxy->at(k) < 0.2 ); 
    }

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_id->at(k)) == 11)
    {
      is_loose_lepton = true;//std_vector_lepton_eleIdMedium->at(k);
    }

  return is_loose_lepton;
}
