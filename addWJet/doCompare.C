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
process[iQCD]  = "QCD";


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


  bool doFRvsPt_bin = 0;
  bool doFRvsEta_bin = 0;
  bool doFRvsPV = 0;
  bool doFRvsptIso = 0;
  bool doFRvssip3d = 1;
  bool doFRvsJetEt = 0;
  bool doJetEt = 0; 

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
  
  TH1F* hist_fake_jetet[nProcesses];
  TH1F* hist_signal_jetet[nProcesses];
  TH1F* hist_FR_jetet[nProcesses];

  TH1F* hist_fake_ptIso[nProcesses];
  TH1F* hist_signal_ptIso[nProcesses];
  TH1F* hist_FR_ptIso[nProcesses];

  TH1F* hist_fake_sip3d[nProcesses];
  TH1F* hist_signal_sip3d[nProcesses];
  TH1F* hist_FR_sip3d[nProcesses];

  TH1F* hist_jetEt[nProcesses];

 
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
      
      hist_fake_jetet[ip] = (TH1F*)input[ip]->Get("h_Muon_fake_jetEt");
      hist_signal_jetet[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_jetEt");
      hist_FR_jetet[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_jetEt");
      
      hist_jetEt[ip] =  (TH1F*)input[ip]->Get("h_Muon_jetEt");
      hist_jetEt[ip]->Scale(1./hist_jetEt[ip]->Integral());

      hist_fake_ptIso[ip] = (TH1F*)input[ip]->Get("h_Muon_fake_ptIso");
      hist_signal_ptIso[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_ptIso");
      hist_FR_ptIso[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_ptIso");

      hist_fake_sip3d[ip] = (TH1F*)input[ip]->Get("h_Muon_fake_sip3d");
      hist_signal_sip3d[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_sip3d");
      hist_FR_sip3d[ip] = (TH1F*)input[ip]->Get("h_Muon_signal_sip3d");

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

      hist_fake_jetet[ip] = (TH1F*)input[ip]->Get("h_Ele_fake_jetEt");
      hist_signal_jetet[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_jetEt");
      hist_FR_jetet[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_jetEt");

      hist_fake_ptIso[ip] = (TH1F*)input[ip]->Get("h_Ele_fake_ptIso");
      hist_signal_ptIso[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_ptIso");
      hist_FR_ptIso[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_ptIso");

      hist_fake_sip3d[ip] = (TH1F*)input[ip]->Get("h_Ele_fake_sip3d");
      hist_signal_sip3d[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_sip3d");
      hist_FR_sip3d[ip] = (TH1F*)input[ip]->Get("h_Ele_signal_sip3d");
    }


    hist_FR_pt[ip]->Divide(hist_signal_pt[ip],hist_fake_pt[ip],1.,1.,"");
    hist_FR_eta[ip]->Divide(hist_signal_eta[ip],hist_fake_eta[ip],1.,1.,"");
    hist_FR_pv[ip]->Divide(hist_signal_pv[ip],hist_fake_pv[ip],1.,1.,"");
    hist_FR_jetet[ip]->Divide(hist_signal_jetet[ip],hist_fake_jetet[ip],1.,1.,"");
    hist_FR_ptIso[ip]->Divide(hist_signal_ptIso[ip],hist_fake_ptIso[ip],1.,1.,"");
  
    hist_FR_sip3d[ip]->Divide(hist_signal_sip3d[ip],hist_fake_sip3d[ip],1.,1.,"");
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


 if (doFRvsptIso) 
   {

     TCanvas * fake_ptIso = new TCanvas("fake_ptIso", "fake_ptIso", 750, 750);
     fake_ptIso->cd();


     DrawTH1F(fake_ptIso, hist_FR_ptIso[iQCD],"","","pT*(1+Iso)","#Tight / #Loose");
     LineOpt(hist_FR_ptIso[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_FR_ptIso[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_ptIso, hist_FR_ptIso[iWJets],"sameE1","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_ptIso[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_ptIso[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_ptIso[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_ptIso[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_ptIso[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_ptIso[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_pT_Mu9->SetAxisRange(0,1,"Y");
   }


 if (doFRvssip3d) 
   {

     TCanvas * fake_sip3d = new TCanvas("fake_sip3d", "fake_sip3d", 750, 750);
     fake_sip3d->cd();


     DrawTH1F(fake_sip3d, hist_FR_sip3d[iQCD],"","","SIP3D","#Tight / #Loose");
     LineOpt(hist_FR_sip3d[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_FR_sip3d[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_sip3d, hist_FR_sip3d[iWJets],"sameE1","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_sip3d[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_sip3d[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_sip3d[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_sip3d[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_sip3d[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_sip3d[iQCD]->SetTitleSize(0.05,"X");
     
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
     
     
     DrawTH1F(fake_pv_bin, hist_FR_pv[iWJets],"sameE1","","pv ","#Tight / #Loose");
     LineOpt(hist_FR_pv[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_pv[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_pv[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_pv[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_pv[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_pv[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_et_Mu9->SetAxisRange(0,1,"Y");
   }

 if (doFRvsJetEt) 
   {

     TCanvas * fake_jetet_bin = new TCanvas("fake_jetet_bin", "fake_jetet_bin", 750, 750);
     fake_jetet_bin->cd();


     DrawTH1F(fake_jetet_bin, hist_FR_jetet[iQCD],"","","jetet (GeV/c)","#Tight / #Loose");
     LineOpt(hist_FR_jetet[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_FR_jetet[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_jetet_bin, hist_FR_jetet[iWJets],"sameE1","","jetet ","#Tight / #Loose");
     LineOpt(hist_FR_jetet[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_FR_jetet[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.25,0.8,hist_FR_jetet[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.25,0.7,hist_FR_jetet[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_FR_jetet[iQCD]->SetTitleSize(0.05,"Y");
     hist_FR_jetet[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_et_Mu9->SetAxisRange(0,1,"Y");
   }


 if (doJetEt) 
   {

     TCanvas * fake_jetEt = new TCanvas("fake_jetEt", "fake_jetEt", 750, 750);
     fake_jetEt->cd();


     DrawTH1F(fake_jetEt, hist_jetEt[iQCD],"hist","","E_{T} (GeV/c)","#Tight / #Loose");
     LineOpt(hist_jetEt[iQCD],kBlack,2,kSolid);
     MarkerOpt(hist_jetEt[iQCD],kBlack,1,kFullCircle);
     
     
     DrawTH1F(fake_jetEt, hist_jetEt[iWJets],"histsame","","E_{T} (GeV/c)","#Tight / #Loose");
     LineOpt(hist_jetEt[iWJets],kRed,2,kSolid);
     MarkerOpt(hist_jetEt[iWJets],kRed,1,kFullCircle);
     
     
     DrawLegend(0.75,0.8,hist_jetEt[iQCD] , " QCD", "LP",0.050,0.24, 0.10);
     DrawLegend(0.75,0.7,hist_jetEt[iWJets] , " W+Jets ", "LP",0.050,0.24, 0.10);
     
     hist_jetEt[iQCD]->SetTitleSize(0.05,"Y");
     hist_jetEt[iQCD]->SetTitleSize(0.05,"X");
     
     //Muon_FR_et_Mu9->SetAxisRange(0,1,"Y");
   }

}
