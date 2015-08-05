#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>


const Double_t ZMASS = 91.1876;

const UInt_t numberMetCuts = 5;
const UInt_t numberDYMVACuts = 5;

Double_t MetCut[] = {20, 25, 30, 45, 1000};

Double_t DYMVACut_0j[numberDYMVACuts] = {-0.9, -0.86, -0.6, 0.88, 1000};

Double_t DYMVACut_1j[numberDYMVACuts] = {-0.9, -0.86, -0.6, 0.84, 1000};

Bool_t runAtOviedo = false;
Bool_t runAtIfca   = !runAtOviedo;


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
Float_t          trigger;
vector<float>   *std_vector_lepton_id;
vector<float>   *std_vector_lepton_pt;
vector<float>   *std_vector_lepton_eta;
vector<float>   *std_vector_lepton_phi;
vector<float>   *std_vector_lepton_isTightMuon; 
vector<float>   *std_vector_lepton_isMediumMuon;
vector<float>   *std_vector_lepton_eleIdMedium; 
vector<float>   *std_vector_lepton_eleIdVeto;
vector<float>   *std_vector_lepton_chargedHadronIso; 
vector<float>   *std_vector_lepton_photonIso;
vector<float>   *std_vector_lepton_neutralHadronIso;
vector<float>   *std_vector_lepton_sumPUPt;
vector<float>   *std_vector_electron_effectiveArea;
vector<float>   *std_vector_lepton_BestTrackdz;
vector<float>   *std_vector_lepton_BestTrackdxy;




//gROOT->LoadMacro("../utils/TResultsTable.C+");

//------------------------------------------------------------------------------
// LatinosTreeScript
//------------------------------------------------------------------------------
void LatinosTreeScript(Float_t luminosity,
		       Int_t   jetChannel,
		       TString flavorChannel,
		       TString theSample,
		       Bool_t  verbose = 0)
{
  TH1::SetDefaultSumw2();

  TString path = Form("rootfilesTT/%djet/%s/", jetChannel, flavorChannel.Data());

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


  // WW level histograms
  //----------------------------------------------------------------------------
  TH1F* hPtLepton1WWLevel       = new TH1F("hPtLepton1WWLevel",       "", 200, 0, 200);
  TH1F* hPtLepton2WWLevel       = new TH1F("hPtLepton2WWLevel",       "", 200, 0, 200);
  TH1F* hPtDiLeptonWWLevel      = new TH1F("hPtDiLeptonWWLevel",      "", 200, 0, 200);
  TH1F* hMinvWWLevel            = new TH1F("hMinvWWLevel",            "", 200, 0, 200);
  TH1F* hMinvWWLevelHM          = new TH1F("hMinvWWLevelHM",            "", 1000, 0, 3000);
  TH1F* hMtWWLevel              = new TH1F("hMtWWLevel",              "", 250, 0, 250);
  TH1F* hNJets30WWLevel         = new TH1F("hNJetsPF30WWLevel",       "",  10, 0,  10);
  TH1F* hpfMetWWLevel           = new TH1F("hpfMetWWLevel",           "", 150, 0, 150);
  TH1F* hppfMetWWLevel          = new TH1F("hppfMetWWLevel",          "", 150, 0, 150);
  TH1F* hchMetWWLevel           = new TH1F("hchMetWWLevel",           "", 150, 0, 150);
  TH1F* hpchMetWWLevel          = new TH1F("hpchMetWWLevel",          "", 150, 0, 150);
  TH1F* hpminMetWWLevel         = new TH1F("hpminMetWWLevel",         "", 150, 0, 150);
  TH1F* hDeltaRLeptonsWWLevel   = new TH1F("hDeltaRLeptonsWWLevel",   "",  50, 0,   5);
  TH1F* hDeltaPhiLeptonsWWLevel = new TH1F("hDeltaPhiLeptonsWWLevel", "",  32, 0, 3.2);
  TH1F* hDPhiPtllJetWWLevel     = new TH1F("hDPhiPtllJetWWLevel",     "",  32, 0, 3.2);

  TH1F* h_WWLevel_TightFailEvents = new TH1F("h_WWLevel_TightFailEvents", "", 3, 0 , 3);
  TH1F* h_WWLevel_TightTightEvents = new TH1F("h_WWLevel_TightTightEvents", "", 3, 0 , 3);
  TH1F* h_WWLevel_TightLooseEvents = new TH1F("h_WWLevel_TightLooseEvents", "", 3, 0 , 3);


  // TwoLeptons level histograms
  //----------------------------------------------------------------------------
  TH1F* hPtLepton1TwoLeptonsLevel       = new TH1F("hPtLepton1TwoLeptonsLevel",       "", 200, 0, 200);
  TH1F* hPtLepton2TwoLeptonsLevel       = new TH1F("hPtLepton2TwoLeptonsLevel",       "", 200, 0, 200);
  TH1F* hPtDiLeptonTwoLeptonsLevel      = new TH1F("hPtDiLeptonTwoLeptonsLevel",      "", 200, 0, 200);
  TH1F* hMinvTwoLeptonsLevel            = new TH1F("hMinvTwoLeptonsLevel",            "", 200, 0, 200);
  TH1F* hMinvTwoLeptonsLevelHM          = new TH1F("hMinvTwoLeptonsLevelHM",            "", 1000, 0, 3000);
  TH1F* hMtTwoLeptonsLevel              = new TH1F("hMtTwoLeptonsLevel",              "", 250, 0, 250);
  TH1F* hNJets30TwoLeptonsLevel         = new TH1F("hNJetsPF30TwoLeptonsLevel",       "",  10, 0,  10);
  TH1F* hpfMetTwoLeptonsLevel           = new TH1F("hpfMetTwoLeptonsLevel",           "", 150, 0, 150);
  TH1F* hppfMetTwoLeptonsLevel          = new TH1F("hppfMetTwoLeptonsLevel",          "", 150, 0, 150);
  TH1F* hchMetTwoLeptonsLevel           = new TH1F("hchMetTwoLeptonsLevel",           "", 150, 0, 150);
  TH1F* hpchMetTwoLeptonsLevel          = new TH1F("hpchMetTwoLeptonsLevel",          "", 150, 0, 150);
  TH1F* hpminMetTwoLeptonsLevel         = new TH1F("hpminMetTwoLeptonsLevel",         "", 150, 0, 150);
  TH1F* hDeltaRLeptonsTwoLeptonsLevel   = new TH1F("hDeltaRLeptonsTwoLeptonsLevel",   "",  50, 0,   5);
  TH1F* hDeltaPhiLeptonsTwoLeptonsLevel = new TH1F("hDeltaPhiLeptonsTwoLeptonsLevel", "",  32, 0, 3.2);
  TH1F* hDPhiPtllJetTwoLeptonsLevel     = new TH1F("hDPhiPtllJetTwoLeptonsLevel",     "",  32, 0, 3.2);

  TH1F* h_TwoLeptons_TightFailEvents = new TH1F("h_TwoLeptons_TightFailEvents", "", 3, 0 , 3);
  TH1F* h_TwoLeptons_TightTightEvents = new TH1F("h_TwoLeptons_TightTightEvents", "", 3, 0 , 3);
  TH1F* h_TwoLeptons_TightLooseEvents = new TH1F("h_TwoLeptons_TightLooseEvents", "", 3, 0 , 3);





  // Data-driven methods: Z+jets
  //----------------------------------------------------------------------------
  TH1F* hNinZevents     [numberDYMVACuts];
  TH1F* hNoutZevents    [numberDYMVACuts];
  TH1F* hNinLooseZevents[numberDYMVACuts];
  TH1F* hMassInZevents  [numberDYMVACuts];
  TH1F* hMassOutZevents [numberDYMVACuts];

  for (Int_t nC=0; nC<numberDYMVACuts; nC++) {
    hNinZevents     [nC] = new TH1F(Form("hNinZevents%.1i",      nC+1 ), "",   3, 0,   3);
    hNoutZevents    [nC] = new TH1F(Form("hNoutZevents%.1i",     nC+1 ), "",   3, 0,   3);
    hNinLooseZevents[nC] = new TH1F(Form("hNinLooseZevents%.1i", nC+1 ), "",   3, 0,   3);
    hMassInZevents  [nC] = new TH1F(Form("hMassInZevents%.1i",   nC+1 ), "", 200, 0, 200);
    hMassOutZevents [nC] = new TH1F(Form("hMassOutZevents%.1i",  nC+1 ), "", 200, 0, 200);

   


  }


  // Data-driven methods: Top
  //----------------------------------------------------------------------------
  TH1F* hTopTaggedEvents            = new TH1F("hTopTaggedEvents",            "", 3, 0, 3);
  TH1F* hNTopControlRegion          = new TH1F("hNTopControlRegion",          "", 3, 0, 3);
  TH1F* hNTopTaggedTopControlRegion = new TH1F("hNTopTaggedTopControlRegion", "", 3, 0, 3);

  TH1F* hbTagDisTopTaggedEvents            = new TH1F("hbTagDisTopTaggedEvents",            "", 300, -10, 20);
  TH1F* hbTagDisNTopControlRegion          = new TH1F("hbTagDisNTopControlRegion",          "", 300, -10, 20);
  TH1F* hbTagDisNTopTaggedTopControlRegion = new TH1F("hbTagDisNTopTaggedTopControlRegion", "", 300, -10, 20);


  //----------------------------------------------------------------------------
  // Input files
  //----------------------------------------------------------------------------
  TString filesPath;

  filesPath = "/gpfs/csic_projects/cms/piedra/latino/RunII/MC_Spring15/50ns/";//"/gpfs/csic_projects/cms/piedra/latino/RunII/MC_Spring15/";

  TChain* tree = new TChain("latino", "latino");

  if (theSample == "data") {
    tree->Add("/gpfs/csic_projects/cms/piedra/latino/RunII/Data13TeVRun2015B/dataRun2015B/latino_MuonEG.root");
    tree->Add("/gpfs/csic_projects/cms/piedra/latino/RunII/Data13TeVRun2015B/dataRun2015B/latino_SingleElectron.root");
    tree->Add("/gpfs/csic_projects/cms/piedra/latino/RunII/Data13TeVRun2015B/dataRun2015B/latino_SingleMuon.root");
  }
  else if (theSample == "WJets") {
    tree->Add(filesPath +"latino_WJetsToLNu.root");
  }
  else if (theSample == "WWTo2L2Nu_pow") {
    tree->Add(filesPath + "latino_WWTo2L2Nu.root");
  }
  else if (theSample == "Top") {
    tree->Add(filesPath + "latino_TTJets.root");
  }
  else if (theSample == "WZ") {
    tree->Add(filesPath + "latino_WZ.root");
  }
  else if (theSample == "ZZ") {
    tree->Add(filesPath + "latino_ZZ.root");
  }
  else if (theSample == "DY") {
    tree->Add(filesPath + "latino_DYJetsToLL.root");
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
 
  tree->SetBranchAddress("trigger", &trigger);
  tree->SetBranchAddress("std_vector_lepton_pt", &std_vector_lepton_pt);
  tree->SetBranchAddress("std_vector_lepton_eta", &std_vector_lepton_eta);
  tree->SetBranchAddress("std_vector_lepton_phi", &std_vector_lepton_phi);
  tree->SetBranchAddress("std_vector_lepton_id", &std_vector_lepton_id);
  tree->SetBranchAddress("std_vector_lepton_isTightMuon", &std_vector_lepton_isTightMuon);
  tree->SetBranchAddress("std_vector_lepton_isMediumMuon", &std_vector_lepton_isMediumMuon);
  tree->SetBranchAddress("std_vector_lepton_eleIdMedium", &std_vector_lepton_eleIdMedium);
  tree->SetBranchAddress("std_vector_lepton_eleIdVeto", &std_vector_lepton_eleIdVeto);
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
  
  //  else if (flavorChannel == "All" ) SelectedChannel = -1;


  //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------
  for (int ievent=0; ievent<tree->GetEntries(); ievent++) {
    
    tree->GetEntry(ievent);
    
    Double_t efficiencyW = effW * triggW;
    Double_t totalW      = -999;

    if (theSample.Contains("data")) {
      totalW = 1.0;
    }   
    else { 
      totalW = baseW*luminosity;
    } 

    // Help variables
    //--------------------------------------------------------------------------
       
    float pt1  = std_vector_lepton_pt ->at(0); float pt2  = std_vector_lepton_pt ->at(1);
    float eta1 = std_vector_lepton_eta->at(0); float eta2 = std_vector_lepton_eta->at(1);
    float ch1  = std_vector_lepton_id->at(0);  float ch2  = std_vector_lepton_id->at(1);
    bool passLooseIDISO1 = IsLooseLepton(0) && IsLooseIsolatedLepton(0); 
    bool passLooseIDISO2 = IsLooseLepton(1) && IsLooseIsolatedLepton(1); 
    bool passIDISO1 = IsTightLepton(0) && IsIsolatedLepton(0); 
    bool passIDISO2 = IsTightLepton(1) && IsIsolatedLepton(1); 

    Float_t jetbin = njet;

    if (njet == 3) jetbin = 2;  



    // The selection begins here
    //--------------------------------------------------------------------------
    //if (theSample == "DY" && mctruth == 2 && (flavorChannel == "EMu" || flavorChannel == "MuE")) continue;
    if (trigger != 1) continue;
    if (pt2 <= 20) continue; // increase the pt of the leptons to further reduce Wjets 
    if (pt1 <= 20) continue; // increase the pt of the leptons to further reduce Wjets 
    if (ch1*ch2 > 0) continue;

    if (!passIDISO1) continue;
    if (!passIDISO2) continue;


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

        
  

      hWTrigger   ->Fill(1, totalW); 
      hWeffTrigger->Fill(1, efficiencyW);
      


      hPtLepton1TwoLeptonsLevel      ->Fill(pt1,       totalW);
      hPtLepton2TwoLeptonsLevel      ->Fill(pt2,       totalW);
      hPtDiLeptonTwoLeptonsLevel     ->Fill(ptll,      totalW);
      hMinvTwoLeptonsLevel           ->Fill(mll,       totalW);
      hMinvTwoLeptonsLevelHM         ->Fill(mll,       totalW);
      hMtTwoLeptonsLevel             ->Fill(mth,       totalW);
      hNJets30TwoLeptonsLevel        ->Fill(njet,      totalW);
      //hpfMetTwoLeptonsLevel          ->Fill(pfmet,     totalW);
      //hppfMetTwoLeptonsLevel         ->Fill(ppfmet,    totalW);
      //hchMetTwoLeptonsLevel          ->Fill(chmet,     totalW);
      //hpchMetTwoLeptonsLevel         ->Fill(pchmet,    totalW);
      //hpminMetTwoLeptonsLevel        ->Fill(mpmet,     totalW);
      hDeltaRLeptonsTwoLeptonsLevel  ->Fill(drll,      totalW);
      hDeltaPhiLeptonsTwoLeptonsLevel->Fill(dphill,    totalW);
      hDPhiPtllJetTwoLeptonsLevel    ->Fill(dphilljet, totalW);
      
      /*
      if (nextra == 0) {
	
	hWExtraLepton->Fill(1, totalW);
	hWeffExtraLepton->Fill(1, efficiencyW);
	
	if (pfmet > 20 ) { // removed for differential xsec
	  
	  hWMetCut->Fill(1, totalW);
	  hWeffMetCut->Fill(1, efficiencyW);
	  
	  if (mll > 12) {
	    
	    hWLowMinv->Fill(1, totalW);
	    hWeffLowMinv->Fill(1, efficiencyW);
	    
	    if (zveto==1 || !sameflav) {
	      
	      hWZVeto->Fill(1, totalW); 
	      hWeffZVeto->Fill(1, efficiencyW); 
	      
	      if (mpmet > 20  && dyMVA ) { 
		
		
		hWpMetCut->Fill(1, totalW);
		hWeffpMetCut->Fill(1, efficiencyW);		  
		
		
		if (dphiv || !sameflav) {
		  
		  
		  hWDeltaPhiJet->Fill(1, totalW);
		  hWeffDeltaPhiJet->Fill(1, efficiencyW);
		  
		  
		  if (bveto_mu) { //--> third soft lepton with pt> 3GeV 
		    
		    hWSoftMuVeto->Fill(1, totalW); 
		    hWeffSoftMuVeto->Fill(1, efficiencyW);
		    
		    
		    
		    if ( ptll>30 && (!sameflav || ptll>45) ) {
		      
		      hWPtll->Fill(1, totalW);
		      hWeffPtll->Fill(1, efficiencyW);
		      
		      
		      if (jetbin == jetChannel) {
			
			hWJetVeto->Fill(1, totalW);
			hWeffJetVeto->Fill(1, efficiencyW);
			
			
			if (bveto_ip==1 &&  nbjettche==0) {
			  
			  hWTopTagging->Fill(1, totalW);
			  hWeffTopTagging->Fill(1, efficiencyW);
			  
			  
			  
			  hPtLepton1WWLevel      ->Fill(pt1,       totalW);
			  hPtLepton2WWLevel      ->Fill(pt2,       totalW);
			  hPtDiLeptonWWLevel     ->Fill(ptll,      totalW);
			  hMinvWWLevel           ->Fill(mll,       totalW);
			  hMinvWWLevelHM         ->Fill(mll,       totalW);
			  hMtWWLevel             ->Fill(mth,       totalW);
			  hNJets30WWLevel        ->Fill(njet,      totalW);
			  hpfMetWWLevel          ->Fill(pfmet,     totalW);
			  hppfMetWWLevel         ->Fill(ppfmet,    totalW);
			  hchMetWWLevel          ->Fill(chmet,     totalW);
			  hpchMetWWLevel         ->Fill(pchmet,    totalW);
			  hpminMetWWLevel        ->Fill(mpmet,     totalW);
			  hDeltaRLeptonsWWLevel  ->Fill(drll,      totalW);
			  hDeltaPhiLeptonsWWLevel->Fill(dphill,    totalW);
			  hDPhiPtllJetWWLevel    ->Fill(dphilljet, totalW);
			  
			  hPtLepton1WWLevel_Diff  ->Fill(pt1,       totalW);
			  hEtaLepton1WWLevel_Diff ->Fill(eta1,      totalW);
			  hDileptonWWLevel_Diff   ->Fill(ptll,      totalW);
			  hMinvWWLevel_Diff       ->Fill(mll,       totalW);
			  hDeltaPhiWWLevel_Diff   ->Fill(dphill,    totalW);


			}
		      }
		    }		      
		  }
		}
	      }
	    }
	  }
	}
      }
      */
      
    }
  
  }
	  
  verbose = true;

  // Print
  //----------------------------------------------------------------------------
  if (verbose) {
    
    if (!theSample.Contains("Data")) {
      cout << endl;
      cout << " Normalized to " << luminosity << " 1/fb" << endl;
      cout << " ------------------+-----------" << endl;
      cout << " trigger           | " << hWTrigger    ->GetSumOfWeights() << " +- " << 
	hWTrigger    ->GetBinError(2)<< endl;
      cout << " extra lepton veto | " << hWExtraLepton->GetSumOfWeights() << endl;
      cout << " MET cut           | " << hWMetCut     ->GetSumOfWeights() << endl;
      cout << " low minv cut      | " << hWLowMinv    ->GetSumOfWeights() << endl;
      cout << " Z veto            | " << hWZVeto      ->GetSumOfWeights() << endl;
      cout << " projected MET cut | " << hWpMetCut    ->GetSumOfWeights() << endl;
      cout << " DeltaPhiJet veto  | " << hWDeltaPhiJet->GetSumOfWeights() << endl;     
      cout << " soft muon veto    | " << hWSoftMuVeto ->GetSumOfWeights() << endl;
      cout << " jet veto          | " << hWJetVeto    ->GetSumOfWeights() << endl;
      cout << " top tagging       | " << hWTopTagging ->GetSumOfWeights() << endl; 

      cout << endl;

    }
 
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
      float dxyCut = 0;
	
      if ( std_vector_lepton_pt->at(k) < 20 ) { 
	dxyCut = 0.01;
      }	else {
	dxyCut = 0.02;
      }

      is_tight_lepton = ( std_vector_lepton_isTightMuon->at(k)         && 
			  fabs(std_vector_lepton_BestTrackdz->at(k))  < 0.1       && 
      			  fabs(std_vector_lepton_BestTrackdxy->at(k)) < dxyCut );
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
  if (fabs(std_vector_lepton_id->at(k)) == 13)
    {
     
      is_loose_lepton = ( std_vector_lepton_isTightMuon->at(k)                   && 
			  fabs(std_vector_lepton_BestTrackdz->at(k))  < 0.1      && 
			  fabs(std_vector_lepton_BestTrackdxy->at(k)) < 0.2 ); 
    }

  // Electron cut based medium ID
  else if (fabs(std_vector_lepton_id->at(k)) == 11)
    {
      is_loose_lepton = true;//std_vector_lepton_eleIdVeto->at(k);
    }
 
  return is_loose_lepton;
}



//------------------------------------------------------------------------------
// IsLooseIsolatedLepton
//------------------------------------------------------------------------------
bool IsLooseIsolatedLepton(int k)
{
  float id = std_vector_lepton_id->at(k);

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k) < 0.5);
  
  return is_isolated_lepton;
}
