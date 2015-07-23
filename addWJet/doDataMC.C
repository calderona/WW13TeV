#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"

const UInt_t nProcesses = 3;

enum {iZJets, iWJets, iQCD};


TFile* input[nProcesses];

TString process[nProcesses];

process[iZJets]  = "ZJets";
process[iWJets]  = "WJets";
process[iQCD]  = "QCD";



Color_t color[nProcesses];

color[iZJets]  = kGreen+2;
color[iWJets]  = kGray+1;
color[iQCD]    = kAzure-5;

Bool_t   _setLogy    = false;


//----------------------------------------------------------------
// drawDistributions
//------------------------------------------------------------------------------
void doDataMC() {

 // LOADING MACROS. 
  gROOT->LoadMacro("../utils/TResultsTable.C+");
  gInterpreter->LoadMacro   ("../utils/draw.C+");
  gInterpreter->LoadMacro   ("../utils/utils.C+"); 
  gInterpreter->ExecuteMacro("../utils/ChargeRatioStyle.C");
  gStyle      ->SetOptStat  (0);
  gStyle      ->SetPalette  (1);

  
  gStyle->SetHatchesLineWidth(1.00);
  gStyle->SetHatchesSpacing  (0.55);



  // Read input files
  //----------------------------------------------------------------------------
  //TString path = Form("rootfiles/%djet/%s/", _njet, _channel.Data());
  
  TString path = "rootfiles/";
  
  for (UInt_t ip=0; ip<nProcesses; ip++) {
    input[ip] = new TFile(path + process[ip] + "_FakeRate.root", "read");
  }


  // Distributions
  //----------------------------------------------------------------------------
 
  if (1) { 
    DrawHistogram("h_Muon_Mt", "M_{T}", 10, 0, "GeV", 0, 200);
    DrawHistogram("h_Muon_Met", "MET", 10, 0, "GeV", 0, 200);
    DrawHistogram("h_Muon_pt",  "p_{T}", 10, 0, "GeV", 15, 200);
    DrawHistogram("h_Muon_jetEt",  "E_{T}", 10, 0, "GeV", 15, 200);
  }

}

//------------------------------------------------------------------------------
// DrawHistogram
//------------------------------------------------------------------------------
void DrawHistogram(TString  hname,
		   TString  xtitle,
		   Int_t    ngroup       = -1,
		   Int_t    precision    = 1,
		   TString  units        = "NULL",
		   Double_t xmin         = -999,
		   Double_t xmax         =  999,
		   Bool_t   moveOverflow = true)
{
  
  TCanvas* canvas;
  TPad* pad1;
  
  canvas = new TCanvas(hname, hname, 550, 600);
  pad1 = (TPad*)canvas->GetPad(0);
  

  //----------------------------------------------------------------------------
  // pad1
  //----------------------------------------------------------------------------
  pad1->cd();

  pad1->SetLogy(_setLogy);

  THStack* hstack = new THStack(hname, hname);

  TH1F* hist[nProcesses];

  for (UInt_t ip=0; ip<nProcesses; ip++) {

    hist[ip] = (TH1F*)input[ip]->Get(hname);
    hist[ip]->SetName(hname + process[ip]);

    if (ngroup > 0) hist[ip]->Rebin(ngroup);

    if (moveOverflow) MoveOverflowBins  (hist[ip], xmin, xmax);
    else              ZeroOutOfRangeBins(hist[ip], xmin, xmax);
    
    
    hist[ip]->SetFillColor(color[ip]);
    hist[ip]->SetFillStyle(1001);
    hist[ip]->SetLineColor(color[ip]);
     
      hstack->Add(hist[ip]);
  }





  // Draw
  //----------------------------------------------------------------------------
 
  hstack     ->Draw("hist");


  // Legend
  //----------------------------------------------------------------------------
  Double_t yoffset = 0.048;
  Double_t x0      = 0.400; 

  DrawLegend(0.73, 0.74 + 2.*(yoffset+0.001), hist[iQCD], " QCD ",    "lp", 0.035, 0.2, yoffset);
  DrawLegend(0.73, 0.74 + 1.*(yoffset+0.001), hist[iWJets],   " W+Jets",      "f",  0.035, 0.2, yoffset);
  DrawLegend(0.73, 0.74,                      hist[iZJets],   " Z+Jets",      "f",  0.035, 0.2, yoffset);
 

}
