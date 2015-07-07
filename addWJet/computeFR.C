#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>

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
  Float_t        fr;
  Float_t        pr;

};

std::vector<Lepton> AnalysisLeptons;


Float_t          jetRho;
Float_t          channel;
Float_t          baseW;
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


vector<float>   *std_vector_leptonGen_pt;
vector<float>   *std_vector_leptonGen_eta;
vector<float>   *std_vector_leptonGen_phi;
vector<float>   *std_vector_leptonGen_id;
vector<float>   *std_vector_leptonGen_status;
vector<float>   *std_vector_leptonGen_mpid;
vector<float>   *std_vector_leptonGen_mstatus;




//------------------------------------------------------------------------------
// MAIN function
//------------------------------------------------------------------------------


int computeFR () {

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


  //--- OPEN LOOSE-LOOSE FILE 
  //----------------------------------------------------------------------------
  string openMode = "read";
  
  TString filesPath;
  filesPath = "/gpfs/csic_projects/cms/piedra/latino/25ns/";
		 
  TChain* tree = new TChain("latino", "latino");

  tree->Add(filesPath + "latino_WJetsToLNu.root");

 
  tree->SetBranchAddress("baseW  ", &baseW);
  tree->SetBranchAddress("jetRho", &jetRho);
  tree->SetBranchAddress("channel", &channel);
  
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

  tree->SetBranchAddress("std_vector_leptonGen_pt", &std_vector_leptonGen_pt);
  tree->SetBranchAddress("std_vector_leptonGen_eta", &std_vector_leptonGen_eta);
  tree->SetBranchAddress("std_vector_leptonGen_phi", &std_vector_leptonGen_phi);
  tree->SetBranchAddress("std_vector_leptonGen_id", &std_vector_leptonGen_id);
  tree->SetBranchAddress("std_vector_leptonGen_status", &std_vector_leptonGen_status);
  tree->SetBranchAddress("std_vector_leptonGen_mpid", &std_vector_leptonGen_mpid);
  tree->SetBranchAddress("std_vector_leptonGen_mstatus", &std_vector_leptonGen_mstatus);




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
      is_tight_lepton = std_vector_lepton_isTightMuon->at(k);
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
