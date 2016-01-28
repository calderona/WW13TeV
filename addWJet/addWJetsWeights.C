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

Double_t GetPPFWeight(int a, int b, int s = 2);
Double_t GetPFFWeight(int a, int b, int s = 2);
Double_t GetFFFWeight(int a, int b, int s = 2);

Double_t GetPFWeight(int a, int b, int s = 2);
Double_t GetFFWeight(int a, int b, int s = 2);

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

enum {up, down};

enum {Muon, Electron, Unknown};

enum {noTight, Tight};

int ch;
 
enum {mm=0,ee=1,em=2,me=3};

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

Float_t         njet;
Float_t          jetRho;
Float_t          channel;
Float_t          baseW;
Float_t          trigger;
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
vector<float>   *std_vector_lepton_BestTrackdz;
vector<float>   *std_vector_lepton_BestTrackdxy;

vector<float>   *std_vector_electron_dEtaIn;
vector<float>   *std_vector_electron_dPhiIn;
vector<float>   *std_vector_electron_full5x5_sigmaIetaIeta;
vector<float>   *std_vector_electron_hOverE;
vector<float>   *std_vector_electron_ooEmooP;
vector<float>   *std_vector_electron_expectedMissingInnerHits;
vector<float>   *std_vector_electron_d0;
vector<float>   *std_vector_electron_dz;
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

 
  MuonPR = LoadHistogram("ZJets_MuPR_RunII_25ns_08Jan_2120pb", "h_Muon_signal_pt_eta_bin", "MuonPR");

  ElecPR = LoadHistogram("EGPR_RunII_25ns_08Jan_2120pb", "h_Ele_signal_pt_eta_bin", "ElecPR");
  
  for ( int j=0; j<nMuJetET; j++) { 
  
    MuonFR[j] = LoadHistogram(Form("MuFR_RunII_25ns_jet%s_08Jan", muonJetPt[j].Data()), "FR_pT_eta_EWKcorr", Form("MuonFR_Jet%s", muonJetPt[j].Data()));

  }

  for ( int j=0; j<nEleJetET; j++) { 
    ElecFR[j] = LoadHistogram(Form("EGFR_RunII_25ns_jet%s_08Jan", elecJetPt[j].Data()), "FR_pT_eta_EWKcorr", Form("ElecFR_Jet%s", elecJetPt[j].Data()));
 }
 

 //--- OPEN LOOSE-LOOSE FILE 
 //----------------------------------------------------------------------------


  TFile* file = new TFile(Form("%s%s", indir.c_str(),input.c_str()),openMode.c_str());

  TTree* tree = (TTree*)file->Get("latino");

  printf("input tree from %s%s: %d events\n", indir.c_str(),input.c_str(), (int)tree->GetEntries());
  
  tree->SetBranchAddress("njet", &njet);
 tree->SetBranchAddress("baseW  ", &baseW);
 tree->SetBranchAddress("jetRho", &jetRho);
 tree->SetBranchAddress("channel", &channel);
 tree->SetBranchAddress("trigger", &trigger);
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
 tree->SetBranchAddress("std_vector_lepton_BestTrackdz", &std_vector_lepton_BestTrackdz);
 tree->SetBranchAddress("std_vector_lepton_BestTrackdxy", &std_vector_lepton_BestTrackdxy);
 tree->SetBranchAddress("std_vector_lepton_closejet_PartonFlavour", &std_vector_lepton_closejet_PartonFlavour);

  tree->SetBranchAddress("std_vector_electron_dEtaIn", &std_vector_electron_dEtaIn);
  tree->SetBranchAddress("std_vector_electron_dPhiIn", &std_vector_electron_dPhiIn);
  tree->SetBranchAddress("std_vector_electron_full5x5_sigmaIetaIeta", &std_vector_electron_full5x5_sigmaIetaIeta);
  tree->SetBranchAddress("std_vector_electron_hOverE",&std_vector_electron_hOverE);
  tree->SetBranchAddress("std_vector_electron_ooEmooP",&std_vector_electron_ooEmooP);
  tree->SetBranchAddress("std_vector_electron_expectedMissingInnerHits", &std_vector_electron_expectedMissingInnerHits);
  tree->SetBranchAddress("std_vector_electron_d0",&std_vector_electron_d0);
  tree->SetBranchAddress("std_vector_electron_dz", &std_vector_electron_dz);
  tree->SetBranchAddress("std_vector_electron_passConversionVeto", &std_vector_electron_passConversionVeto);


//--- New branches to be added
//----------------------------------------------------------------------------
 
 Float_t weightFF, weightPF, weightUp, weightDown; 
 Float_t weightFFF, weightPFF, weightPPF; 

 Float_t weight3l, weight3lUp, weight3lDown, weight3lstatUp, weight3lstatDown; 
 Float_t weight2l0j, weight2l1j, weight2l2j, weight2l0jUp, weight2l1jUp, weight2l2jUp, weight2l0jDown, weight2l1jDown, weight2l2jDown, weight2l0jstatUp, weight2l1jstatUp, weight2l2jstatUp, weight2l0jstatDown, weight2l1jstatDown, weight2l2jstatDown;


 tree->SetBranchStatus("*", 1);
 // in case already exist
 tree->SetBranchStatus("fakeW3l",               0);
 tree->SetBranchStatus("fakeW3lUp",             0);
 tree->SetBranchStatus("fakeW3lDown",           0);
 tree->SetBranchStatus("fakeW3lstatUp",         0);
 tree->SetBranchStatus("fakeW3lstatDown",       0);
 tree->SetBranchStatus("fakeW2l0j",             0);
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

   newtree->Branch("fakeW3l",            &weight3l,                "fakeW3l/F");
   //newtree->Branch("fakeW2l",            &weight2l,                "fakeW2l/F");
   newtree->Branch("fakeW3lUp",          &weight3lUp,              "fakeW3lUp/F");
   newtree->Branch("fakeW3lDown",        &weight3lDown,            "fakeW3lDown/F");
   newtree->Branch("fakeW3lstatUp",      &weight3lstatUp,          "fakeW3lstatUp/F");
   newtree->Branch("fakeW3lstatDown",    &weight3lstatDown,        "fakeW3lstatDown/F");
   newtree->Branch("fakeW2l0j",          &weight2l0j,              "fakeW2l0j/F");
   newtree->Branch("fakeW2l1j",          &weight2l1j,              "fakeW2l1j/F");
   newtree->Branch("fakeW2l2j",          &weight2l2j,              "fakeW2l2j/F");
   newtree->Branch("fakeW2l0jUp",        &weight2l0jUp,            "fakeW2l0jUp/F");
   newtree->Branch("fakeW2l1jUp",        &weight2l1jUp,            "fakeW2l1jUp/F");
   newtree->Branch("fakeW2l2jUp",        &weight2l2jUp,            "fakeW2l2jUp/F");
   newtree->Branch("fakeW2l0jDown",      &weight2l0jDown,          "fakeW2l0jDown/F");
   newtree->Branch("fakeW2l1jDown",      &weight2l1jDown,          "fakeW2l1jDown/F");
   newtree->Branch("fakeW2l2jDown",      &weight2l2jDown,          "fakeW2l2jDown/F");
   newtree->Branch("fakeW2l0jstatUp",    &weight2l0jstatUp,        "fakeW2l0jstatUp/F");
   newtree->Branch("fakeW2l1jstatUp",    &weight2l1jstatUp,        "fakeW2l1jstatUp/F");
   newtree->Branch("fakeW2l2jstatUp",    &weight2l2jstatUp,        "fakeW2l2jstatUp/F");
   newtree->Branch("fakeW2l0jstatDown",  &weight2l0jstatDown,      "fakeW2l0jstatDown/F");
   newtree->Branch("fakeW2l1jstatDown",  &weight2l1jstatDown,      "fakeW2l1jstatDown/F");
   newtree->Branch("fakeW2l2jstatDown",  &weight2l2jstatDown,      "fakeW2l2jstatDown/F");


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
    weight2l0jUp = 0.0; weight2l1jUp = 0.0; weight2l2jUp = 0.0;
    weight2l0jDown = 0.0; weight2l1jDown = 0.0; weight2l2jDown = 0.0;
    weight2l0jstatUp = 0.0; weight2l1jstatUp = 0.0; weight2l2jstatUp = 0.0;
    weight2l0jstatDown = 0.0; weight2l1jstatDown = 0.0; weight2l2jstatDown = 0.0;
    weightPF = 0.0; 
    weightFF = 0.0; 

    weight3l = 0.0; weight3lUp = 0.0; weight3lDown = 0.0;
    weight3lstatUp = 0.0; weight3lstatDown = 0.0;
    weightPPF = 0.0; 
    weightPFF = 0.0; 
    weightFFF = 0.0; 


    float sign = 1.; 
    //if ( input.Contains("Zg") ) sign = -1.; // flip sign to subtract non-fake contamination from Vg/Vg*
    
    if (jentry%10000==0) printf("Processing event %d \n", int(jentry));

    if ( trigger == 0) continue; 


    // Loop over leptons
    //--------------------------------------------------------------------------
    AnalysisLeptons.clear();
     
    int vsize = std_vector_lepton_pt->size();
    
    for (int i=0; i<vsize; i++) {

      if ( !IsLooseLepton(i) ) continue; 
      if ( !IsLooseIsolatedLepton(i) ) continue; 
 
      //if (fabs(std_vector_lepton_closejet_PartonFlavour->at(i)) == 5) continue;
      
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
    
     if (pt <=10 ) continue;
     if (fabs(id) == 11 && fabs(eta) >= 2.5 ) continue;  
     if (fabs(id) == 13 && fabs(eta) >= 2.4 ) continue;  

     if ( IsTightLepton(i) &&  IsIsolatedLepton(i) ) { 
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
    

    
    if ( AnalysisLeptons.size() == 2 )  {
      weightPF    = float(GetPFWeight(int(mu20), int(ele35))); 
      weightFF    = float(GetFFWeight(int(mu20), int(ele35))); 
      weight2l0j = weightPF + weightFF;
      weight2l1j = float(GetPFWeight(int(mu25), int(ele35))) + float(GetFFWeight(int(mu25), int(ele35)));
      weight2l2j = float(GetPFWeight(int(mu35), int(ele35))) + float(GetFFWeight(int(mu35), int(ele35)));
      weight2l0jUp = float(GetPFWeight(int(mu30), int(ele45))) + float(GetFFWeight(int(mu30), int(ele45)));
      weight2l1jUp = float(GetPFWeight(int(mu35), int(ele45))) + float(GetFFWeight(int(mu35), int(ele45)));
      weight2l2jUp = float(GetPFWeight(int(mu45), int(ele45))) + float(GetFFWeight(int(mu45), int(ele45)));
      weight2l0jDown = float(GetPFWeight(int(mu10), int(ele25))) + float(GetFFWeight(int(mu10), int(ele25)));
      weight2l1jDown = float(GetPFWeight(int(mu15), int(ele25))) + float(GetFFWeight(int(mu15), int(ele25)));
      weight2l2jDown = float(GetPFWeight(int(mu25), int(ele25))) + float(GetFFWeight(int(mu25), int(ele25)));
      weight2l0jstatUp = float(GetPFWeight(int(mu20), int(ele35), int(up))) + float(GetFFWeight(int(mu30), int(ele45), int(up)));
      weight2l1jstatUp = float(GetPFWeight(int(mu25), int(ele35), int(up))) + float(GetFFWeight(int(mu25), int(ele35), int(up)));
      weight2l2jstatUp = float(GetPFWeight(int(mu35), int(ele35), int(up))) + float(GetFFWeight(int(mu35), int(ele35), int(up)));
      weight2l0jstatDown = float(GetPFWeight(int(mu20), int(ele35), int(down))) + float(GetFFWeight(int(mu30), int(ele45), int(down)));
      weight2l1jstatDown = float(GetPFWeight(int(mu25), int(ele35), int(down))) + float(GetFFWeight(int(mu25), int(ele35), int(down)));
      weight2l2jstatDown = float(GetPFWeight(int(mu35), int(ele35), int(down))) + float(GetFFWeight(int(mu35), int(ele35), int(down)));
   
    }
    
    if ( AnalysisLeptons.size() == 3 )  {
      weightPPF   =  float(GetPPFWeight(int(mu20), int(ele35))); 
      weightPFF   =  float(GetPFFWeight(int(mu20), int(ele35))); 
      weightFFF   =  float(GetFFFWeight(int(mu20), int(ele35))); 
      weight3l     = weightPPF + weightPFF + weightFFF;
      weight3lUp =  float(GetPPFWeight(int(mu30), int(ele45))) + float(GetPFFWeight(int(mu30), int(ele45))) + float(GetFFFWeight(int(mu30), int(ele45))); 
      weight3lDown =  float(GetPPFWeight(int(mu10), int(ele25))) + float(GetPFFWeight(int(mu10), int(ele25))) + float(GetFFFWeight(int(mu10), int(ele25))); 
      weight3lstatUp =  float(GetPPFWeight(int(mu20), int(ele35), int(up))) + float(GetPFFWeight(int(mu20), int(ele35), int(up))) + float(GetFFFWeight(int(mu20), int(ele35), int(up))); 
      weight3lstatDown =  float(GetPPFWeight(int(mu20), int(ele35), int(down))) + float(GetPFFWeight(int(mu20), int(ele35), int(down))) + float(GetFFFWeight(int(mu20), int(ele35), int(down))); 
    }
     


    bool isTight1 = 0, isTight2 = 0; 
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

    //------ summing the weights only for those events passing the full selection
    float pt1  = AnalysisLeptons[0].pt; float pt2  = AnalysisLeptons[1].pt;
    float ch1  = AnalysisLeptons[0].charge;  float ch2  = AnalysisLeptons[1].charge;
    int flav1 =  AnalysisLeptons[0].flavour;  int flav2 =  AnalysisLeptons[1].flavour;

    
    if (flav1 == Muon     && flav2 == Muon )    ch = mm; 
    if( flav1 == Electron && flav2 == Electron) ch = ee; 
    if (flav1 == Muon     && flav2 == Electron) ch = me; 
    if (flav1 == Electron && flav2 == Muon)     ch = em; 
    
    
    bool commonSel =  ( AnalysisLeptons.size() == 2 && trigger == 1 && pt1 > 10 && pt2 > 10 && ch1*ch2 < 0);

    if(commonSel) {
      //cout << weight << endl;
      countEntriesByHand2++;

      if(ch == mm) {
	sumOfWeightsMuMu += weightPF; 
	sumOfQCDWeightsMuMu += weightFF; 
      }

      if(ch == me) {
	sumOfWeightsMuEl += weightPF;
	sumOfQCDWeightsMuEl += weightFF; 
      }

      if(ch == em) {
	sumOfWeightsElMu += weightPF; 
	sumOfQCDWeightsElMu += weightFF; 
      }      

      if(ch == ee) {
	sumOfWeightsElEl += weightPF; 
	sumOfQCDWeightsElEl += weightFF; 
      }
      
      char* label;
      if( isTight1 && isTight2 ) label="PP";
      if( !isTight1 && !isTight2) label="FF";
      if(( isTight1 && !isTight2) || ( !isTight1 &&  isTight2)) label="PF";
	

      if(label == "PP"){
	if(ch == mm) mmPP++;
	if(ch == me) mePP++;
	if(ch == em) emPP++;
	if(ch == ee) eePP++;
      }
      if(label == "PF"){
	if(ch == mm) mmPF++;
	if(ch == me) mePF++;
	if(ch == em) emPF++;
	if(ch == ee) eePF++;
      }
      if(label == "FF"){
	if(ch == mm) mmFF++;
	if(ch == me) meFF++;
	if(ch == em) emFF++;
	if(ch == ee) eeFF++;
      }

    } // end selection

    // normalise to a unit lumi 
    baseW *= 1./luminosity;
    countEntriesByHand++;

    //if(addWeight) newtree->Fill(); 
    if ( AnalysisLeptons.size() >= 2 )  {
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
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 30.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 30.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
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
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 30.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 30.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
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
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 30.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 30.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
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
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 30.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 30.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
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
      fV = GetFactor(MuonFR[mu],  lep.pt, lep.eta, 30.);
      fE = GetFactorError(MuonFR[mu],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
      else f = fV;

    } else if ( lep.flavour == Electron) {
      fV = GetFactor(ElecFR[ele],  lep.pt, lep.eta, 30.0);
      fE = GetFactorError(ElecFR[ele],  lep.pt, lep.eta, 30.);
      if (stat == up) f = fV+fE; 
      else if (stat==down) f = fV-fE;
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
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k) < 0.15);
  
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
      
      if (aeta <= 1.479) {
	
	if (fabs(std_vector_electron_dEtaIn->at(k)) < 0.01               &&
	    fabs(std_vector_electron_dPhiIn->at(k)) < 0.04               &&
	    std_vector_electron_full5x5_sigmaIetaIeta->at(k)    < 0.011              &&
	    std_vector_electron_hOverE->at(k)              < 0.08               &&
	    std_vector_electron_ooEmooP->at(k) < 0.01               && //????
	    std_vector_electron_expectedMissingInnerHits->at(k)<=2  &&
	    std_vector_electron_d0->at(k)      < 0.1                &&
	    fabs(std_vector_electron_dz->at(k))< 0.373              &&
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
	    std_vector_electron_d0->at(k)      < 0.2                &&
	    fabs(std_vector_electron_dz->at(k))< 0.602              &&
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

  bool is_isolated_lepton = false;

  if      (fabs(id) == 11) is_isolated_lepton = true;  //(ElectronIsolation(k) < 0.15);
  else if (fabs(id) == 13) is_isolated_lepton = (MuonIsolation(k)     < 0.4);
  
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
