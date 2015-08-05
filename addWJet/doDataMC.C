#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"

const UInt_t nProcesses = 4;

enum {iData, iZJets, iWJets, iQCD};


TFile* input[nProcesses];

TString process[nProcesses];

process[iData]  = "dataEG";
process[iZJets]  = "ZJets";
process[iWJets]  = "WJets";
process[iQCD]  = "QCD";



Color_t color[nProcesses];

color[iData] = kBlack;
color[iZJets]  = kGreen+2;
color[iWJets]  = kGray+1;
color[iQCD]    = kAzure-5;

Bool_t   _setLogy    = true;


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

  
  



  // Read input files
  //----------------------------------------------------------------------------
  //TString path = Form("rootfiles/%djet/%s/", _njet, _channel.Data());
  
  TString path = "rootfiles/";
  
  for (UInt_t ip=0; ip<nProcesses; ip++) 
    input[ip] = new TFile(path + process[ip] + "_FakeRate.root", "read");
  



  // Distributions
  //----------------------------------------------------------------------------
 
  if (0) { 
    DrawHistogram("h_Muon_fake_Mt", "M_{T}", 1, 0, "GeV", 0, 200);
    DrawHistogram("h_Muon_fake_Mt_Met20", "M_{T}", 1, 0, "GeV", 0, 200);
    DrawHistogram("h_Muon_fake_Met", "MET", 1, 0, "GeV", 0, 200);
    DrawHistogram("h_Muon_fake_pt",  "p_{T}", 1, 0, "GeV", 15, 200);
    DrawHistogram("h_Muon_fake_CRjetEt",  "E_{T}", 1, 0, "GeV", 15, 200);
  }

 if (1) { 
    DrawHistogram("h_Ele_fake_Mt", "M_{T}", 1, 0, "GeV", 0, 200);
    DrawHistogram("h_Ele_fake_Mt_Met20", "M_{T}", 1, 0, "GeV", 0, 200);
    DrawHistogram("h_Ele_fake_Met", "MET", 1, 0, "GeV", 0, 200);
    DrawHistogram("h_Ele_fake_pt",  "p_{T}", 1, 0, "GeV", 15, 200);
    DrawHistogram("h_Ele_fake_CRjetEt",  "E_{T}", 1, 0, "GeV", 15, 200);
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
    
    if (ip == iData) {
      hist[ip]->SetMarkerStyle(kFullCircle);
    } 
    else {

      hist[ip]->Scale(40.03/1000);

      hist[ip]->SetFillColor(color[ip]);
      hist[ip]->SetFillStyle(1001);
      hist[ip]->SetLineColor(color[ip]);
     
      hstack->Add(hist[ip]);
    }
  }
  


  // Axis labels
  //----------------------------------------------------------------------------
  TAxis* xaxis = hist[iData]->GetXaxis();
  TAxis* yaxis = hist[iData]->GetYaxis();

  TString ytitle = Form("entries / %s.%df", "%", precision);

  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(Form(ytitle.Data(), hist[iData]->GetBinWidth(0)));
  yaxis->SetTitleOffset(1.6);

  if (!units.Contains("NULL")) {
    
    xaxis->SetTitle(Form("%s [%s]", xaxis->GetTitle(), units.Data()));
    yaxis->SetTitle(Form("%s %s",   yaxis->GetTitle(), units.Data()));
  }


  // Draw
  //----------------------------------------------------------------------------
  xaxis->SetRangeUser(xmin, xmax);

  hist[iData]->Draw("ep");
  hstack     ->Draw("hist, same");
  hist[iData]->Draw("ep, same");


  // Legend
  //----------------------------------------------------------------------------
  Double_t yoffset = 0.048;
  Double_t x0      = 0.400; 

  //DrawLegend(0.73, 0.74 + 2.*(yoffset+0.001), hist[iQCD], " QCD ",    "f", 0.035, 0.2, yoffset);
  DrawLegend(0.73, 0.74 + 1.*(yoffset+0.001), hist[iWJets],   " W+Jets",      "f",  0.035, 0.2, yoffset);
  DrawLegend(0.73, 0.74,                      hist[iZJets],   " Z+Jets",      "f",  0.035, 0.2, yoffset);
 

}

//------------------------------------------------------------------------------
// MoveOverflowBins
//------------------------------------------------------------------------------
void MoveOverflowBins(TH1* h,
		      Double_t xmin,
		      Double_t xmax) const
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  Int_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax != -999) ? axis->FindBin(xmax) : nbins;

  Double_t firstVal = 0;
  Double_t firstErr = 0;

  Double_t lastVal = 0;
  Double_t lastErr = 0;

  for (UInt_t i=0; i<=nbins+1; i++) {

    if (i <= firstBin) {
      firstVal += h->GetBinContent(i);
      firstErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i >= lastBin) {
      lastVal += h->GetBinContent(i);
      lastErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i < firstBin || i > lastBin) {
      h->SetBinContent(i, 0);
      h->SetBinError  (i, 0);
    }
  }

  firstErr = sqrt(firstErr);
  lastErr  = sqrt(lastErr);

  h->SetBinContent(firstBin, firstVal);
  h->SetBinError  (firstBin, firstErr);

  h->SetBinContent(lastBin, lastVal);
  h->SetBinError  (lastBin, lastErr);
}


//------------------------------------------------------------------------------
// ZeroOutOfRangeBins
//------------------------------------------------------------------------------
void ZeroOutOfRangeBins(TH1* h, Double_t xmin, Double_t xmax) const
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  Int_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax != -999) ? axis->FindBin(xmax) : nbins;

  for (UInt_t i=0; i<=nbins+1; i++) {

    if (i < firstBin || i > lastBin) {
      h->SetBinContent(i, 0);
      h->SetBinError  (i, 0);
    }
  }
}
