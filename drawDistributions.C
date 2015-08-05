#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TTree.h"


const UInt_t nProcesses = 7;

//enum {iData, iWW, iWZ, iZZ, iWg, itt, itW, iWj, iDY, iDYtau, iH125};
//enum {iData, iWW, iqqWW, itt, itW, iWZ, iZZ, iWg, iWgS, iWj, iDY, iDYtau,  iZg, iVVV, iH125};


enum {iData, iWW, itt, iWZ, iZZ, iWj, iDY}; 

TFile* input[nProcesses];


TString process[nProcesses];

process[iData]  = "data";
process[itt]    = "Top";
//process[itW]    = "TW";
process[iWW]    = "WWTo2L2Nu_pow";
//process[iqqWW]  = "WWTo2L2Nu_pow_nnll_smear";
process[iWZ]    = "WZ";
process[iZZ]    = "ZZ";
//process[iWg]    = "WgammaNoStar";
//process[iWgS]   = "WgammaStar";
process[iWj]    = "WJets";
process[iDY]    = "DY";
//process[iDYtau] = "DYtautau";
//process[iZg]    = "Zgamma";
//process[iVVV]   = "VVV";
//process[iH125]  = "HWW125";



Color_t color[nProcesses];

color[iData]  = kBlack;
color[itt]    = kYellow;
//color[itW]    = kYellow;
color[iWW]    = kAzure-9;
//color[iqqWW]  = kAzure-9;
color[iWZ]    = kAzure-2;
color[iZZ]    = kAzure-2;
//color[iWg]    = kAzure-2;
//color[iWgS]   = kAzure-2;
color[iWj]    = kGray+1;
color[iDY]    = kGreen+2;
//color[iDYtau] = kGreen+2;
//color[iZg]    = kAzure-2;
//color[iVVV]   = kAzure-2;
//color[iH125]  = kRed;



// Relative errors directly from data cards (in %) 
//------------------------------------------------------------------------------

enum backgrounds {
   WW,
   ggWW,
   Vg,
   Wjets,
   VVV,
   DYtt,
   VV,
   ggH,
   VgS
};
 

enum myChannels {
  jet0_SF=0*9,
  jet0_OF=1*9,
  jet1_SF=2*9,
  jet1_OF=3*9
};

const UInt_t Nsyst = 36;

//[*] WW / ggWW/  Vg/ Wjets/ VVV/ DY/ VV/ ggH/ VgS

Double_t Systematics [Nsyst] = { 6.9,  30.4,  31.3,  27.4,  50.6,   2.6,   5.9,  32.5,  34.1,   // jet0_SF
				 6.5,  30.4,  30.0,  28.3,  50.4,  30.7,   6.2,  28.7,  30.3,   // jet0_OF
				 13.7,  32.1,  33.3,  23.9,  50.3,   2.0,   6.4,  41.0,  58.3,  // jet1_SF 
				 13.3,  31.8,  30.1,  26.4,  50.3,  30.4,   5.6,  35.5,  31.2};  // jet1_OF
				 //9.5,	43.0,  43.4,  39.4,  71.4,  30.8,   8.6,  43.4,	45.6,   // jet0
				 //19.1,	45.2,  44.9,  35.6,  71.1,  30.5,   8.5,  54.2,	66.1 }; // jet1


Double_t systError[nProcesses];

systError[iData] = 0.0;





// Settings
//------------------------------------------------------------------------------
TString  _channel;
Int_t    _njet;
Double_t _luminosity;
TString  _format;
Bool_t   _drawRatio;
Bool_t   _dataDriven;
Bool_t   _savePlots;
Bool_t   _setLogy;


// Scale factors
//------------------------------------------------------------------------------
Double_t ttScale[] = {1.07, 1.09, 1.09, 1.09};
Double_t tWScale[] = {1.07, 1.09, 1.09, 1.09};
Double_t WWScale[] = {1.00, 1.00, 1.00, 1.00};
Double_t ZjScale[] = {1.00, 1.00, 1.00, 1.00} ;//{1.00, 4.20, 1.80, 4.00};


//------------------------------------------------------------------------------
// drawDistributions
//------------------------------------------------------------------------------
void drawDistributions(Int_t    njet       = 0,
		       TString  channel    = "OF",
		       Double_t luminosity = 40.03,
		       TString  format     = "pdf",
		       Bool_t   drawRatio  = false,
		       Bool_t   dataDriven = true,
		       Bool_t   savePlots  = false,
		       Bool_t   setLogy    = false)
{
  _channel    = channel;
  _njet       = njet;
  _luminosity = luminosity;
  _format     = format;
  _drawRatio  = drawRatio;
  _dataDriven = dataDriven;
  _savePlots  = savePlots;
  _setLogy    = setLogy;

  gSystem->mkdir(_format, kTRUE);

  gStyle->SetHatchesLineWidth(1.00);
  gStyle->SetHatchesSpacing  (0.55);

  // LOADING MACROS. 
  
  gInterpreter->LoadMacro   ("utils/draw.C+");
  gInterpreter->LoadMacro   ("utils/utils.C+"); 
  gInterpreter->ExecuteMacro("utils/ChargeRatioStyle.C");
  gStyle      ->SetOptStat  (0);
  gStyle      ->SetPalette  (1);



  // Read input files
  //----------------------------------------------------------------------------
  //TString path = Form("rootfiles/%djet/%s/", _njet, _channel.Data());

  TString path = Form("rootfilesTT/%djet/%s/", _njet, _channel.Data());

  for (UInt_t ip=0; ip<nProcesses; ip++)
    input[ip] = new TFile(path + process[ip] + ".root", "read");



  // Differential distributions
  //----------------------------------------------------------------------------
  if (0) {

    DrawHistogram("hPtLepton1WWLevel_Diff",  "p_{T}^{max}",           1, 0, "GeV");
    DrawHistogram("hDileptonWWLevel_Diff",  "p_{T} (ll) ",           1, 0, "GeV");
    DrawHistogram("hDeltaPhiWWLevel_Diff",  "#Delta#phi_{ll}", 1, 0, "rad");
    DrawHistogram("hMinvWWLevel_Diff",  "m_{#font[12]{ll}}",           1, 0, "GeV");
  }

  // DY distributions
  //----------------------------------------------------------------------------
  if (0) {
    DrawHistogram("hMassInZevents45.0",  "m_{#font[12]{ll}}", 1, 0, "GeV", 76, 106);
    DrawHistogram("hMassOutZevents45.0", "m_{#font[12]{ll}}", 5, 0, "GeV");
    
    DrawHistogram("hMassInZevents30.0",  "m_{#font[12]{ll}}", 1, 0, "GeV", 76, 106);
    DrawHistogram("hMassOutZevents30.0", "m_{#font[12]{ll}}", 5, 0, "GeV");
  }


  // Top distributions
  //----------------------------------------------------------------------------
  if (0) {
    DrawHistogram("hbTagDisNTopTaggedTopControlRegion", "2^{nd} jet TCHE", 5, 1, "NULL", -999, -999, false);
    DrawHistogram("hbTagDisNTopControlRegion",          "2^{nd} jet TCHE", 5, 1, "NULL", -999, -999, false);
    DrawHistogram("hbTagDisTopTaggedEvents",            "2^{nd} jet TCHE", 5, 1, "NULL", -999, -999, false);
  }


  // PAS distributions
  //----------------------------------------------------------------------------
  if (0) {
    DrawHistogram("hPtLepton1WWLevel",  "p_{T}^{max}",           10, 0, "GeV", 20, 200);
    DrawHistogram("hPtLepton2WWLevel",  "p_{T}^{min}",           5, 0, "GeV", 20, 100);
    DrawHistogram("hPtDiLeptonWWLevel", "p_{T}^{#font[12]{ll}}", 5, 0, "GeV");
    DrawHistogram("hMinvWWLevel",       "m_{#font[12]{ll}}",     5, 0, "GeV");
    //DrawHistogram("hNJetsPF30WWLevel",  "nJets",                 1, 0, "");
    //DrawHistogram("hDeltaRLeptonsWWLevel",   "#DeltaR_{ll}",     5, 0, "rad");
    DrawHistogram("hDeltaPhiLeptonsWWLevel", "#Delta#phi_{ll}",  1, 0, "rad");
    //DrawHistogram("hpminMetWWLevel",         "min.proj.E_{T}",   6, 0, "GeV", 20, 200);
      
  }

  if(1){
    _dataDriven = false;
    systError[itt]    = 0.; // 10%
    //systError[itW]    = 0.; // 10%
    systError[iWW]    = 0.; // 7% 
    systError[iWZ]    = 0.; // 8%;
    systError[iZZ]    = 0.; // 8%;
    //systError[iWg]    = 0.; // 32%
    systError[iWj]    = 0.; // 36%
    systError[iDY]    = 0.; // 37%
    //systError[iDYtau] = 0.; // ?
    //systError[iH125]  = 0.;// 7.4%
    _setLogy = false;

    DrawHistogram("hPtLepton1TwoLeptonsLevel",  "p_{T}^{max}",           4, 0, "GeV", 20, 200);
    DrawHistogram("hPtLepton2TwoLeptonsLevel",  "p_{T}^{min}",           4, 0, "GeV", 20, 100);
    DrawHistogram("hPtDiLeptonTwoLeptonsLevel", "p_{T}^{#font[12]{ll}}", 4, 0, "GeV", 20, 200);
    DrawHistogram("hMinvTwoLeptonsLevel",       "m_{#font[12]{ll}}",     4, 0, "GeV");
    DrawHistogram("hNJetsPF30TwoLeptonsLevel",  "nJets",                 1, 0, "");
    DrawHistogram("hDeltaRLeptonsTwoLeptonsLevel",   "#DeltaR_{ll}",     4, 0, "rad");
    DrawHistogram("hDeltaPhiLeptonsTwoLeptonsLevel", "#Delta#phi_{ll}",  4, 0, "rad");
    DrawHistogram("hpminMetTwoLeptonsLevel",         "min.proj.E_{T}",   10, 0, "GeV", 0, 200);
    DrawHistogram("hppfMetTwoLeptonsLevel",          "proj.PFE_{T}",     10, 0, "GeV", 0, 200);

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
  TPad* pad2;

  if (_drawRatio) {
    canvas = new TCanvas(hname, hname, 550, 1.2*600);

    pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetTopMargin   (0.05);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3); 
    pad2->SetTopMargin   (0.08);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
  }
  else {
    canvas = new TCanvas(hname, hname, 550, 600);
    pad1 = (TPad*)canvas->GetPad(0);
  }


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
      hist[ip]->SetFillColor(color[ip]);
      hist[ip]->SetFillStyle(1001);
      hist[ip]->SetLineColor(color[ip]);

      //if(ip != iData)
      //hist[ip]->Scale(1.25);

      //      if (ip == iH125) hist[ip]->SetFillStyle(3354);

      if (_dataDriven && ip == itt)    hist[ip]->Scale(ttScale[_njet]);
      //if (_dataDriven && ip == itW)    hist[ip]->Scale(tWScale[_njet]);
      if (_dataDriven && ip == iWW)    hist[ip]->Scale(WWScale[_njet]);
      //if (_dataDriven && ip == iqqWW)    hist[ip]->Scale(WWScale[_njet]);
      if (_dataDriven && ip == iDY)    hist[ip]->Scale(ZjScale[_njet]);
      //if (_dataDriven && ip == iDYtau) hist[ip]->Scale(ZjScale[_njet]);

      hstack->Add(hist[ip]);
    }
  }


  // Include relative systematics 
  //----------------------------------------------------------------------------

  Int_t defineChannel = -1; 
  
  if ( _njet == 0 && _channel == "SF")  defineChannel  = jet0_SF;
  if ( _njet == 0 && _channel == "OF")  defineChannel  = jet0_OF;
  if ( _njet == 1 && _channel == "SF")  defineChannel  = jet1_SF;
  if ( _njet == 1 && _channel == "OF")  defineChannel  = jet1_OF;

  systError[itt]     = 0.10; // 10%
  //systError[itW]     = 0.10; // 10%
  
  //[*] WW / ggWW/  Vg/ Wjets/ VVV/ DY/ VV/ ggH/ VgS

  systError[iDY]     = Systematics [DYtt+defineChannel]/ 1e2;     
  systError[iWW]     = Systematics [ggWW+defineChannel]/ 1e2;
  //systError[iqqWW]   = Systematics [WW+defineChannel]  / 1e2;
  systError[iWZ]     = Systematics [VV+defineChannel]  / 1e2;
  systError[iZZ]     = Systematics [VV+defineChannel]  / 1e2;
  //systError[iWg]     = Systematics [Vg+defineChannel]  / 1e2;
  //systError[iWgS]    = Systematics [VgS+defineChannel] / 1e2;
  systError[iWj]     = Systematics [Wjets+defineChannel] / 1e2;
  //systError[iDYtau]  = Systematics [DYtt+defineChannel]/ 1e2;     
  //systError[iZg]     = Systematics [Vg+defineChannel] / 1e2;
  //systError[iVVV]    = Systematics [VVV+defineChannel]/ 1e2;
  //systError[iH125]   = Systematics [ggH+defineChannel]/ 1e2;


  // All MC
  //----------------------------------------------------------------------------
  TH1F* allmc = (TH1F*)hist[iData]->Clone("allmc");

  allmc->SetFillColor  (kGray+2);
  allmc->SetFillStyle  (   3345);
  allmc->SetLineColor  (kGray+2);
  allmc->SetMarkerColor(kGray+2);
  allmc->SetMarkerSize (      0);

  for (UInt_t ibin=1; ibin<=allmc->GetNbinsX(); ibin++) {

    Double_t binValue = 0;
    Double_t binError = 0;

    for (UInt_t ip=0; ip<nProcesses; ip++) {

      if (ip == iData) continue;

      Double_t binContent = hist[ip]->GetBinContent(ibin);
      
      binValue += binContent;
      binError += (hist[ip]->GetBinError(ibin) * hist[ip]->GetBinError(ibin));

    

      /*
      if (_dataDriven) { 
	//binError += (systError[ip]*binContent * systError[ip]*binContent);
	//}
    
	if ( nProcesses != itW &&  nProcesses != itt) { 
	  binError += (systError[ip]*binContent * systError[ip]*binContent);
	  // luminosity syst
	  binError += (0.026*binContent * 0.026*binContent);
	} else {
	  binError += (systError[ip]*binContent*ttScale[_njet] * systError[ip]*binContent*ttScale[_njet]);
	  // luminosity syst
	  binError += (0.026*binContent *ttScale[_njet]* 0.026*binContent*tScale[_njet]);
	}
	}*/

      binError = sqrt(binError);

      allmc->SetBinContent(ibin, binValue);
      allmc->SetBinError  (ibin, binError);
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
  hstack     ->Draw("hist,same");
  allmc      ->Draw("e2,same");
  hist[iData]->Draw("ep,same");


  // Adjust scale
  //----------------------------------------------------------------------------
  Float_t theMax   = GetMaximumIncludingErrors(hist[iData], xmin, xmax);
  Float_t theMaxMC = GetMaximumIncludingErrors(allmc,       xmin, xmax);

  if (theMaxMC > theMax) theMax = theMaxMC;

  if (pad1->GetLogy()) {

    theMax = TMath::Power(10, TMath::Log10(theMax) + 4.0);//2.7

    hist[iData]->SetMinimum(0.05);
  }
  else theMax *= 1.55;

  hist[iData]->SetMaximum(theMax);


  // Legend
  //----------------------------------------------------------------------------
  Double_t yoffset = (_drawRatio) ? 0.054 : 0.048;
  Double_t x0      = (_drawRatio) ? 0.370 : 0.400; 

  TString allmcTitle = (_dataDriven) ? " stat #oplus syst" : " #sigma_{stat}";

  DrawLegend(0.23, 0.74 + 2.*(yoffset+0.001), hist[iData], " data",    "lp", 0.035, 0.2, yoffset);
  DrawLegend(0.23, 0.74 + 1.*(yoffset+0.001), hist[iWW],   " WW",      "f",  0.035, 0.2, yoffset);
  DrawLegend(0.23, 0.74,                      hist[iWZ],   " VV",      "f",  0.035, 0.2, yoffset);
  //DrawLegend(0.23, 0.74 - 1.*(yoffset+0.001), hist[iH125], " H125",    "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74 + 2.*(yoffset+0.001), hist[iDY],   " Z+jets",  "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74 + 1.*(yoffset+0.001), hist[iWj],   " W+jets",  "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74,                      hist[itt],   " top",     "f",  0.035, 0.2, yoffset);
  DrawLegend(  x0, 0.74 - 1.*(yoffset+0.001), allmc,       allmcTitle, "f",  0.035, 0.2, yoffset);


  // Additional titles
  //----------------------------------------------------------------------------
  Double_t deltaY = (_drawRatio) ? 0.02 : 0.0;

  DrawTLatex(0.9, 0.860 + deltaY, 0.04, "CMS preliminary");
  DrawTLatex(0.9, 0.815 + deltaY, 0.03, Form("L = %.3f pb^{-1}", _luminosity/1e3));

  TString channelLabel; 

  if (_njet == 0) channelLabel = "0-jet";

  //CJ
  DrawTLatex(0.9, 0.725 + deltaY, 0.03, _channel);
  DrawTLatex(0.9, 0.695 + deltaY, 0.03, channelLabel);
  
  if (_dataDriven) {
    DrawTLatex(0.9, 0.770 + deltaY, 0.03, "data-driven normalization");
  }


  //----------------------------------------------------------------------------
  // pad2
  //----------------------------------------------------------------------------
  if (_drawRatio) {

    pad2->cd();
    
    TH1F* ratio       = hist[iData]->Clone("ratio");
    TH1F* uncertainty = allmc->Clone("uncertainty");

    for (UInt_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {

      Double_t mcValue = allmc->GetBinContent(ibin);
      Double_t mcError = allmc->GetBinError  (ibin);

      Double_t dtValue = ratio->GetBinContent(ibin);
      Double_t dtError = ratio->GetBinError  (ibin);

      Double_t ratioValue       = (mcValue > 0) ? dtValue/mcValue : 0.0;
      Double_t ratioError       = (mcValue > 0) ? dtError/mcValue : 0.0;
      Double_t uncertaintyError = (mcValue > 0) ? mcError/mcValue : 0.0;

      ratio->SetBinContent(ibin, ratioValue);
      ratio->SetBinError  (ibin, ratioError);

      uncertainty->SetBinContent(ibin, 1.0);
      uncertainty->SetBinError  (ibin, uncertaintyError);
    }


    TAxis* uaxis = (TAxis*)uncertainty->GetXaxis();
    
    uaxis->SetRangeUser(xmin, xmax);
    
    
    uncertainty->Draw("e2");
    ratio      ->Draw("ep,same");

    uncertainty->GetYaxis()->SetRangeUser(0.4, 1.6);

    Pad2TAxis(uncertainty, hist[iData]->GetXaxis()->GetTitle(), "data / prediction");
  }


  // Save
  //----------------------------------------------------------------------------
  pad1->cd();
  pad1->GetFrame()->DrawClone();
  pad1->RedrawAxis();

  if (_drawRatio) {
    pad2->cd();
    pad2->GetFrame()->DrawClone();
    pad2->RedrawAxis();
  }

  if (_savePlots) {
    canvas->cd();
    canvas->SaveAs(Form("%s/%s_%djet_%s.%s",
			_format.Data(),
			hname.Data(),
			_njet,
			_channel.Data(),
			_format.Data())
		   );
  }
}


//------------------------------------------------------------------------------
// GetMaximumIncludingErrors
//------------------------------------------------------------------------------
Float_t GetMaximumIncludingErrors(TH1F*    h,
				  Double_t xmin,
				  Double_t xmax)
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  Int_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
  Int_t lastBin  = (xmax != -999) ? axis->FindBin(xmax) : nbins;

  Float_t maxWithErrors = 0;

  for (Int_t i=firstBin; i<=lastBin; i++) {

    Float_t binHeight = h->GetBinContent(i) + h->GetBinError(i);

    if (binHeight > maxWithErrors) maxWithErrors = binHeight;
  }

  return maxWithErrors;
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


//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t x, Double_t y, Double_t tsize, const char* text)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC();
  tl->SetTextAlign(   32);
  tl->SetTextFont (   42);
  tl->SetTextSize (tsize);

  tl->Draw("same");
}


//------------------------------------------------------------------------------
// Pad2TAxis
//------------------------------------------------------------------------------
void Pad2TAxis(TH1* hist, TString xtitle, TString ytitle)
{
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.025);
  xaxis->SetLabelSize  (   0.1);
  xaxis->SetNdivisions (   505);
  xaxis->SetTitle      (xtitle);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.35);
  xaxis->SetTitleSize  (  0.11);
  
  yaxis->CenterTitle   (      );
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset(  0.02);
  yaxis->SetLabelSize  (   0.1);
  yaxis->SetNdivisions (   505);
  yaxis->SetTitle      (ytitle);
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  0.75);
  yaxis->SetTitleSize  (  0.11);
}
