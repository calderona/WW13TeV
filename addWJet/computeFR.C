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
  Float_t        fr;
  Float_t        pr;
  TLorentzVector v;

};

std::vector<Lepton> AnalysisLeptons;


Float_t          jetRho;
Float_t          channel;
Float_t          baseW;
Float_t          pfType1Met;  
Float_t          mth;
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
vector<float>   *std_vector_leptonGen_pid;
vector<float>   *std_vector_leptonGen_status;
vector<float>   *std_vector_leptonGen_mpid;
vector<float>   *std_vector_leptonGen_mstatus;




//------------------------------------------------------------------------------
// MAIN function
//------------------------------------------------------------------------------


void computeFR (TString theSample="WJets") {

  // ------ root settings ---------
 
  TH1::SetDefaultSumw2();
 
  TFile* root_output = new TFile(theSample + "_FakeRate.root", "recreate");


  //-- Define histograms 
  //----------------------------------------------------------------------------

  // Define pt bins for FR matrix
  int pTbin = 18;
  Double_t pTbins[19] = {10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
 
  // Define eta bins for FR matrix
  int etabin = 4;
  Double_t etabins[5] = {0,1.0,1.479,2.0,2.5};

  TH1F  *h_Muon_fake_pT_bin   = new TH1F("h_Muon_fake_pT_bin","h_Muon_fake_pT_bin", pTbin,pTbins);
  TH1F  *h_Muon_signal_pT_bin = new TH1F("h_Muon_signal_pT_bin","h_Muon_signal_pT_bin",pTbin,pTbins);
  TH1F  *h_Muon_fake_eta_bin = new TH1F("h_Muon_fake_eta_bin","h_Muon_fake_eta_bin",etabin,etabins);
  TH1F  *h_Muon_signal_eta_bin = new TH1F("h_Muon_signal_eta_bin","h_Muon_signal_eta_bin",etabin,etabins);

  TH1F  *h_Muon_fake_nPV = new TH1F("h_Muon_fake_nPV","h_Muon_fake_nPV",30,0,30);
  TH1F  *h_Muon_signal_nPV = new TH1F("h_Muon_signal_nPV","h_Muon_signal_nPV",30,0,30);

  TH1F  *h_Muon_fake_iso   = new TH1F("h_Muon_fake_iso","h_Muon_fake_iso", 100,0,0.5);
  TH1F  *h_Muon_signal_iso   = new TH1F("h_Muon_signal_iso","h_Muon_signal_iso", 100,0,0.5);



  //--- OPEN LOOSE-LOOSE FILE 
  //----------------------------------------------------------------------------
  
  TString filesPath;
  filesPath = "/gpfs/csic_projects/cms/piedra/latino/25ns/";
		 
  TChain* tree = new TChain("latino", "latino");

  if (theSample == "WJets") {
    tree->Add(filesPath + "latino_WJetsToLNu.root");
  }
  else {
    return;
  }


  tree->SetBranchAddress("baseW  ", &baseW);
  tree->SetBranchAddress("jetRho", &jetRho);
  tree->SetBranchAddress("channel", &channel);
  tree->SetBranchAddress("pfType1Met", &pfType1Met);
  tree->SetBranchAddress("mth", &mth);

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
  tree->SetBranchAddress("std_vector_leptonGen_pid", &std_vector_leptonGen_pid);
  tree->SetBranchAddress("std_vector_leptonGen_status", &std_vector_leptonGen_status);
  tree->SetBranchAddress("std_vector_leptonGen_mpid", &std_vector_leptonGen_mpid);
  tree->SetBranchAddress("std_vector_leptonGen_mstatus", &std_vector_leptonGen_mstatus);


 
  //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------

  Long64_t nentries = tree->GetEntries();
  float nE =0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
    tree->GetEntry(jentry);
  
    if (jentry%10000==0) printf("Processing event %d \n", int(jentry));

    Double_t totalW      = 1.0;

    nE++;

    // Help variables (only for deining loose+loose sample) 
    //--------------------------------------------------------------------------
    
    bool passLooseIDISO1 = IsLooseLepton(0) && IsLooseIsolatedLepton(0); 
    bool passLooseIDISO2 = IsLooseLepton(1) && IsLooseIsolatedLepton(1); 
    
    //if (!passLooseIDISO1) continue;
    //if (!passLooseIDISO2) continue;
    
    if ( pfType1Met >= 20 ) continue;
    //if ( mth >= 20) continue;    

    bool passIDISO1 = IsTightLepton(0) && IsIsolatedLepton(0); 
    bool passIDISO2 = IsTightLepton(1) && IsIsolatedLepton(1);


    // Loop over leptons
    //--------------------------------------------------------------------------
    AnalysisLeptons.clear();
     
    int vsize = std_vector_lepton_pt->size();
   
    for (int i=0; i<vsize; i++) {
     
     float pt  = std_vector_lepton_pt ->at(i);
     float eta = std_vector_lepton_eta->at(i);
     float phi = std_vector_lepton_phi->at(i);
     float id  = std_vector_lepton_id ->at(i);

     Lepton lep;
    
     lep.index  = i;
     lep.charge = id;
     
     if (fabs(id) == 11)
       { lep.flavor = Electron; }
     else if (fabs(id) == 13)
       { lep.flavor = Muon; }
     else lep.flavor = Unknown;


     if (IsTightLepton(i) && IsIsolatedLepton(i) )  { 
       lep.type = Tight;
     } else {
       lep.type = noTight;
     }

     //if ( !IsLooseLepton(i)) continue;  
     //if ( !IsLooseIsolatedLepton(i) ) continue; 
   

     TLorentzVector tlv;
     tlv.SetPtEtaPhiM(pt, eta, phi, 0);
     lep.v = tlv;

     lep.fromW = comesFromW (i);
   
     AnalysisLeptons.push_back(lep);

    } // end loop on leptons 


    //if ( AnalysisLeptons.size() > 2 ) continue;


    // Perform analysis
    //--------------------------------------------------------------------------
    if ( theSample == "WJets" ) {

      cout << nE << "  "  << AnalysisLeptons[0].fromW << "  " << AnalysisLeptons[1].fromW << 
	"  "  << AnalysisLeptons[2].fromW << "  "  << AnalysisLeptons[3].fromW <<  "  " << AnalysisLeptons[4].fromW << "  "  << AnalysisLeptons[5].fromW <<  AnalysisLeptons[6].fromW << "  "  << AnalysisLeptons[7].fromW <<  AnalysisLeptons[8].fromW << "  "  << AnalysisLeptons[9].fromW << endl;

      if (AnalysisLeptons[0].fromW)  { 
	h_Muon_signal_iso->Fill ( MuonIsolation(0), totalW);
      } else if (!AnalysisLeptons[0].fromW) {
	h_Muon_fake_iso->Fill ( MuonIsolation(0), totalW);
      }

      if (AnalysisLeptons[1].fromW)  { 
	h_Muon_signal_iso->Fill ( MuonIsolation(1), totalW);
      } else if (!AnalysisLeptons[1].fromW) {
	h_Muon_fake_iso->Fill ( MuonIsolation(1), totalW);
      }



      /*if ( !AnalysisLeptons[0].fromW ||  !AnalysisLeptons[1].fromW ) { 
	
	if ( !AnalysisLeptons[0].fromW && AnalysisLeptons[0].flavor == Muon) {
	  
	  h_Muon_fake_pT_bin->Fill(std_vector_lepton_pt->at(0), totalW);
	
	  if ( passIDISO1 ) h_Muon_signal_pT_bin->Fill(std_vector_lepton_pt->at(0), totalW);
	

	} else if ( !AnalysisLeptons[1].fromW && AnalysisLeptons[1].flavor == Muon) { 

	  h_Muon_fake_pT_bin->Fill(std_vector_lepton_pt->at(1), totalW);

	  if ( passIDISO2 ) h_Muon_signal_pT_bin->Fill(std_vector_lepton_pt->at(1), totalW);
	
	}
	}*/
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

    if( std_vector_leptonGen_pid->at(i)<0 ) continue;

    if( std_vector_leptonGen_status->at(i) != 1 ) continue;

    if ( fabs(std_vector_lepton_id ->at(k)) != fabs(std_vector_leptonGen_pid->at(i)) ) continue;
    
    TLorentzVector genLep;
    genLep.SetPtEtaPhiM(std_vector_leptonGen_pt->at(i), std_vector_leptonGen_eta->at(i), std_vector_leptonGen_phi->at(i), 0);
    
    float dR = lep.DeltaR(genLep); 
       
    if ( dR < 0.1 )   genIndex = i; 
   
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
     
