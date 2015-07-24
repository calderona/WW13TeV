#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"


const UInt_t nProcesses = 2;

enum {iWJets, iQCD};


TFile* input[nProcesses];

TString process[nProcesses];

process[iWJets]  = "WJets_FR";
process[iQCD]  = "QCDEM";


Color_t color[nProcesses];
color[iWJets]  = kRed+1;
color[iQCD]    = kAzure-5;



void doCompare(TString lep = "Muon" ) {



// LOADING MACROS. 
  gROOT->LoadMacro("../utils/TResultsTable.C+");
  gInterpreter->LoadMacro   ("../utils/draw.C+");
  gInterpreter->LoadMacro   ("../utils/utils.C+"); 
  gInterpreter->ExecuteMacro("../utils/ChargeRatioStyle.C");
  gStyle      ->SetOptStat  (0);
  gStyle      ->SetPalette  (1);


  bool doFRvsPt_bin = 1;
  bool doFRvsEta_bin = 1;
  bool doFRvsPV = 1;


//------------------------------------------------------------------------------
// Files and variables
//------------------------------------------------------------------------------


  TString path = "rootfiles/";
  
  for (UInt_t ip=0; ip<nProcesses; ip++) 
    input[ip] = new TFile(path + process[ip] + "_FakeRate.root", "read");
  
 
  TH1F* hist_fake_pt[nProcesses];
  TH1F* hist_signal_pt[nProcesses];
  TH1F* hist_FR_pt[nProcesses];

  TH1F* hist_fake_eta[nProcesses];
  TH1F* hist_signal_eta[nProcesses];
  TH1F* hist_FR_eta[nProcesses];

  TH1F* hist_fake_pv[nProcesses];
  TH1F* hist_signal_pv[nProcesses];
  TH1F* hist_FR_pv[nProcesses];


 
  for (UInt_t ip=0; ip<nProcesses; ip++) {

    if ( lep == "Muon" ) {
      
      hist_fake_pt[ip] = (TH1F*)input[ip]->Get("h_Muon_fake_pT_bin");
      hist_signal_pt[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_pT_bin");
      hist_FR_pt[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_pT_bin");
      
      hist_fake_eta[ip] = (TH1F*)input[ip]->Get("h_Muon_fake_eta_bin");
      hist_signal_eta[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_eta_bin");
      hist_FR_eta[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_eta_bin");
      
      hist_fake_pv[ip] = (TH1F*)input[ip]->Get("h_Muon_fake_nPV");
      hist_signal_pv[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_nPV");
      hist_FR_pv[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_nPV");
      
     
    }     else {

      hist_fake_pt[ip] = (TH1F*)input[ip]->Get("h_Ele_fake_pT_bin");
      hist_signal_pt[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_pT_bin");
      hist_FR_pt[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_pT_bin");
      
      hist_fake_eta[ip] = (TH1F*)input[ip]->Get("h_Ele_fake_eta_bin");
      hist_signal_eta[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_eta_bin");
      hist_FR_eta[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_eta_bin");
      
      hist_fake_pv[ip] = (TH1F*)input[ip]->Get("h_Ele_fake_nPV");
      hist_signal_pv[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_nPV");
      hist_FR_pv[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_nPV");
    }


    hist_FR_pt[ip]->Divide(hist_signal_pt[ip],hist_fake_pt[ip],1.,1.,"");
    hist_FR_eta[ip]->Divide(hist_signal_eta[ip],hist_fake_eta[ip],1.,1.,"");
    hist_FR_pv[ip]->Divide(hist_signal_pv[ip],hist_fake_pv[ip],1.,1.,"");
  }



//+++++++++++++++++++++++++++++ Fake rate plots vs Pt 


 if (doFRvsPt_bin) 
   {

     TCanvas * fake_pT_bin = new TCanvas("fake_pT_bin", "fake_pT_bin", 750, 750);
     fake_pT_bin->cd();


     DrawTH1F(fake_pT_bin, hist_FR_pt[iQCD],"","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_pt[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_FR_pt[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_pT_bin, hist_FR_pt[iWJets],"sameE1","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_pt[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_pt[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_pt[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_pt[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_pt[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_pt[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_pT_Mu9->SetAxisRange(0,1,"Y");
   }




 if (doFRvsEta_bin) 
   {

     TCanvas * fake_eta_bin = new TCanvas("fake_eta_bin", "fake_eta_bin", 750, 750);
     fake_eta_bin->cd();


     DrawTH1F(fake_eta_bin, hist_FR_eta[iQCD],"","","eta (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_eta[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_FR_eta[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_eta_bin, hist_FR_eta[iWJets],"sameE1","","eta (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_eta[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_eta[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_eta[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_eta[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_eta[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_eta[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_et_Mu9->SetAxisRange(0,1,"Y");
   }



 if (doFRvsPV) 
   {

     TCanvas * fake_pv_bin = new TCanvas("fake_pv_bin", "fake_pv_bin", 750, 750);
     fake_pv_bin->cd();


     DrawTH1F(fake_pv_bin, hist_FR_pv[iQCD],"","","pv (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_pv[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_FR_pv[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_pv_bin, hist_FR_pv[iWJets],"sameE1","","pv (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_pv[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_pv[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_pv[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_pv[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_pv[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_pv[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_et_Mu9->SetAxisRange(0,1,"Y");
   }

}
