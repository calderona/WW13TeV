#include "draw.h"


//------------------------------------------------------------------------------
// GetTH1F
//------------------------------------------------------------------------------
TH1F* GetTH1F(TString filename, TString name)
{
  TFile* file = new TFile(filename);
  
  TH1F* hist = (TH1F*)file->Get(name);

  if (!hist) {ExitPrint("GetTH1F()"); return NULL;}
  
  return hist;
}


//------------------------------------------------------------------------------
// GetTH2F
//------------------------------------------------------------------------------
TH2F* GetTH2F(TString filename, TString name)
{
  TFile* file = new TFile(filename);
  
  TH2F* hist = (TH2F*)file->Get(name);

  if (!hist) {ExitPrint("GetTH2F()"); return NULL;}
  
  return hist;
}


//------------------------------------------------------------------------------
// DrawTH1F
//------------------------------------------------------------------------------
void DrawTH1F(TPad*    canvas,
	      TH1*     hist,
	      TString  option,
	      TString  title,
	      TString  xtitle,
	      TString  ytitle,
	      Float_t  xoffset,
	      Float_t  yoffset)
{
  if (!hist) {ExitPrint("DrawTH1F()"); return;}

  TString name = hist->GetName();

  canvas ? canvas->cd() : canvas = new TCanvas(name, name);

   hist->SetDirectory(     0);
  hist->Draw        (option);

  if (!option.Contains("same")) {
    hist->SetTitle(title);
    AxisTH1F(hist, xtitle, ytitle, xoffset, yoffset);
  }
}


//------------------------------------------------------------------------------
// DrawTH1F
//------------------------------------------------------------------------------
void DrawTH1F(TVirtualPad*    canvas,
	      TH1*     hist,
	      TString  option,
	      TString  title,
	      TString  xtitle,
	      TString  ytitle,
	      Float_t  xoffset,
	      Float_t  yoffset)
{
  if (!hist) {ExitPrint("DrawTH1F()"); return;}

  TString name = hist->GetName();

  canvas ? canvas->cd() : canvas = new TCanvas(name, name);

   hist->SetDirectory(     0);
  hist->Draw        (option);

  if (!option.Contains("same")) {
    hist->SetTitle(title);
    AxisTH1F(hist, xtitle, ytitle, xoffset, yoffset);
  }
}



//------------------------------------------------------------------------------
// DrawTH2F
//------------------------------------------------------------------------------
void DrawTH2F(TPad*    canvas,
	      TH2F*    hist,
	      TString  option,
	      TString  title,
	      TString  xtitle,
	      TString  ytitle,
	      Float_t  xoffset,
	      Float_t  yoffset)
{
  if (!hist) {ExitPrint("DrawTH2F()"); return;}

  TString name = hist->GetName();

  canvas ? canvas->cd() : canvas = new TCanvas(name, name);

  hist->SetDirectory(     0);
  hist->Draw        (option);

  if (!option.Contains("same")) {
    hist->SetTitle(title);
    AxisTH2F(hist, xtitle, ytitle, xoffset, yoffset);
  }
}



//------------------------------------------------------------------------------
// DrawTHStack
//------------------------------------------------------------------------------
void DrawTHStack(TPad*    canvas,
	      THStack*     hist,
	      TString  option,
	      TString  title,
	      TString  xtitle,
	      TString  ytitle,
	      Float_t  xoffset,
	      Float_t  yoffset)
{
  if (!hist) {ExitPrint("DrawTHStack()"); return;}

  TString name = hist->GetName();

  canvas ? canvas->cd() : canvas = new TCanvas(name, name);

  //  hist->SetDirectory(     0);
  hist->Draw        (option);

  if (!option.Contains("same")) {
    hist->SetTitle(title);
    AxisTHStack(hist, xtitle, ytitle, xoffset, yoffset);
  }
}



//------------------------------------------------------------------------------
// Axis
//------------------------------------------------------------------------------
void Axis(TAxis*  axis,
	  TString title,
	  Float_t offset,
	  Bool_t  center)
{
  axis->CenterTitle   (center);
  axis->SetTitle      ( title);
  axis->SetTitleOffset(offset);
}


//------------------------------------------------------------------------------
// AxisTH1F
//------------------------------------------------------------------------------
void AxisTH1F(TH1*    hist,
	      TString xtitle,
	      TString ytitle,
	      Float_t xoffset,
	      Float_t yoffset)
{
  if (!hist) {ExitPrint("AxisTH1F()"); return;}

  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  Axis(xaxis, xtitle, xoffset, kFALSE);
  Axis(yaxis, ytitle, yoffset, kTRUE );
}


//------------------------------------------------------------------------------
// AxisTH2F
//------------------------------------------------------------------------------
void AxisTH2F(TH2F*   hist,
	      TString xtitle,
	      TString ytitle,
	      Float_t xoffset,
	      Float_t yoffset)
{
  if (!hist) {ExitPrint("AxisTH2F()"); return;}
  
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  Axis(xaxis, xtitle, xoffset, kFALSE);
  Axis(yaxis, ytitle, yoffset, kTRUE );
}


//------------------------------------------------------------------------------
// AxisTHStack
//------------------------------------------------------------------------------
void AxisTHStack(THStack* hist,
		 TString  xtitle,
		 TString  ytitle,
		 Float_t  xoffset,
		 Float_t  yoffset)
{
  if (!hist) {ExitPrint("AxisTHStack()"); return;}
  
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  Axis(xaxis, xtitle, xoffset, kFALSE);
  Axis(yaxis, ytitle, yoffset, kTRUE );
}


//------------------------------------------------------------------------------
// FillOpt
//------------------------------------------------------------------------------
void FillOpt(TH1F*   hist,
	     Color_t fcolor,
	     Style_t fstyle)
{
  if (!hist) {ExitPrint("FillOpt()"); return;}

  hist->SetFillColor(fcolor);
  hist->SetFillStyle(fstyle);
}


//------------------------------------------------------------------------------
// LineOpt
//------------------------------------------------------------------------------
void LineOpt(TH1F*   hist,
	     Color_t lcolor,
	     Width_t lwidth,
	     Style_t lstyle)
{
  if (!hist) {ExitPrint("LineOpt()"); return;}

  hist->SetLineColor(lcolor);
  hist->SetLineStyle(lstyle);
  hist->SetLineWidth(lwidth);
}

// LineOpt
//------------------------------------------------------------------------------
void LineOpt(TH1*   hist,
	     Color_t lcolor,
	     Width_t lwidth,
	     Style_t lstyle)
{
  if (!hist) {ExitPrint("LineOpt()"); return;}

  hist->SetLineColor(lcolor);
  hist->SetLineStyle(lstyle);
  hist->SetLineWidth(lwidth);
}




//------------------------------------------------------------------------------
// MarkerOpt
//------------------------------------------------------------------------------
void MarkerOpt(TH1F*   hist,
	       Color_t mcolor,
	       Size_t  msize,
	       Style_t mstyle)
{
  if (!hist) {ExitPrint("MarkerOpt()"); return;}

  hist->SetMarkerColor(mcolor);
  hist->SetMarkerSize (msize );
  hist->SetMarkerStyle(mstyle);
}

//------------------------------------------------------------------------------
// MarkerOpt
//------------------------------------------------------------------------------
void MarkerOpt(TH1*   hist,
	       Color_t mcolor,
	       Size_t  msize,
	       Style_t mstyle)
{
  if (!hist) {ExitPrint("MarkerOpt()"); return;}

  hist->SetMarkerColor(mcolor);
  hist->SetMarkerSize (msize );
  hist->SetMarkerStyle(mstyle);
}


//------------------------------------------------------------------------------
// ExitPrint
//------------------------------------------------------------------------------
void ExitPrint(TString method)
{
  printf(" draw::%-11s exit\n", method.Data());
}


//------------------------------------------------------------------------------
// MoveStats
//------------------------------------------------------------------------------
void MoveStats(TCanvas* canvas,
	       Float_t  xoffset)
{
  if (!canvas) {ExitPrint("RemoveStats()"); return;}

  canvas->Update();  // Needed by GetPrimitive

  TPad*       pad = (TPad*)canvas->GetPad(0);
  TPaveStats* tps = (TPaveStats*)pad->GetPrimitive("stats");
   
  tps->SetX1NDC(tps->GetX1NDC() + xoffset);
  tps->SetX2NDC(tps->GetX2NDC() + xoffset);

  canvas->Modified();  // Needed by SetX1NDC
}


//------------------------------------------------------------------------------
// SetTF1
//------------------------------------------------------------------------------
void SetTF1(TF1*    f,
	    TString option,
	    Color_t lcolor,
	    Style_t lstyle)
{
  f->SetLineColor(lcolor);
  f->SetLineStyle(lstyle);

  f->Draw(option.Data());
}


//------------------------------------------------------------------------------
// SetTGraph
//------------------------------------------------------------------------------
void SetTGraph(TGraph* g,
	       TString option,
	       Color_t mcolor,
	       Style_t mstyle,
	       TString title,
	       TString xtitle,
	       TString ytitle)
{
  g->Draw(option);

  g->SetMarkerColor(mcolor);
  g->SetMarkerStyle(mstyle);
  g->SetMarkerSize (   0.9);
  g->SetLineColor  (kBlack);
  g->SetLineStyle  (     1);
  g->SetLineWidth  (     1);

  g->SetTitle(title);

  TAxis* xaxis = g->GetHistogram()->GetXaxis();
  TAxis* yaxis = g->GetHistogram()->GetYaxis();

  Axis(xaxis, xtitle, default_xoffset, kFALSE);
  Axis(yaxis, ytitle, default_yoffset, kTRUE );
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
TLegend* DrawLegend(Float_t x1,
		    Float_t y1,
		    TH1F*   hist,
		    TString label,
		    TString option,
		    Float_t tsize,
		    Float_t xoffset,
		    Float_t yoffset)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + xoffset,
				y1 + yoffset);
  
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (tsize);

  legend->AddEntry(hist, label.Data(), option.Data());
  legend->Draw();

  return legend;
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
TLegend* DrawLegend(Float_t x1,
		    Float_t y1,
		    TH1*   hist,
		    TString label,
		    TString option,
		    Float_t tsize,
		    Float_t xoffset,
		    Float_t yoffset)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + xoffset,
				y1 + yoffset);
  
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (tsize);

  legend->AddEntry(hist, label.Data(), option.Data());
  legend->Draw();

  return legend;
}



//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
TLegend* DrawLegendG (Float_t x1,
		      Float_t y1,
		      TGraph*   hist,
		      TString label,
		      TString option,
		      Float_t tsize,
		      Float_t xoffset,
		      Float_t yoffset)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + xoffset,
				y1 + yoffset);
  
  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (tsize);

  legend->AddEntry(hist, label.Data(), option.Data());
  legend->Draw();

  return legend;
}


//------------------------------------------------------------------------------
// DrawLogo
//------------------------------------------------------------------------------
void DrawLogo(TCanvas* canvas,
	      TString  image,
	      Float_t  xlow,
	      Float_t  yup)
{
  canvas->cd();

  TPad *pad = new TPad("pad", "pad", xlow, yup-0.07, xlow+0.07, yup);

  pad->SetBorderMode  (0);
  pad->SetBorderSize  (0);
  pad->SetFillColor   (0);
  pad->SetFillStyle   (0);

  pad->SetBottomMargin(0);
  pad->SetLeftMargin  (0);
  pad->SetRightMargin (0);
  pad->SetTopMargin   (0);
  
  pad->Draw();
  pad->cd();

  TImage* mark = (TImage*)TImage::Open(image.Data());

  mark->SetConstRatio();
  mark->Draw();

  canvas->cd();
}


//------------------------------------------------------------------------------
// SetCanvas
//------------------------------------------------------------------------------
TCanvas* SetCanvas(TString title, Long_t iCanvas)
{
  if (iCanvas > -1) title += iCanvas;

  TCanvas* canvas = new TCanvas(title, title);

  return canvas;
}


//------------------------------------------------------------------------------
// SaveCanvas
//------------------------------------------------------------------------------
void SaveCanvas(TCanvas* canvas,
		TString  title,
		TString  extension,
		Long_t   iCanvas)
{
  canvas->cd();

  if (iCanvas > -1) title += iCanvas;

  title += extension;

  canvas->Update();
  canvas->GetFrame()->DrawClone();
  canvas->SaveAs(title.Data());
}


//------------------------------------------------------------------------------
// Ratio
//------------------------------------------------------------------------------
Float_t Ratio(Float_t a, Float_t b)
{
  if (b == 0) return 0.0;

  return a / b;
}


//------------------------------------------------------------------------------
// RatioError
//------------------------------------------------------------------------------
Float_t RatioError(Float_t a,
		   Float_t b,
		   Float_t aErr,
		   Float_t bErr)
{
  if (b == 0) return 0.0;

  if (aErr < 0) aErr = sqrt(a);
  if (bErr < 0) bErr = sqrt(b);

  return Ratio(a,b) * sqrt((aErr/a)*(aErr/a) + (bErr/b)*(bErr/b));
}


//------------------------------------------------------------------------------
// GetRatio
//------------------------------------------------------------------------------
TH1F* GetRatio(TH1F* hnum, TH1F* hden, Float_t mean)
{
  TH1F* hratio = (TH1F*)hnum->Clone("hratio");

  hratio->Reset();

  Int_t nbins = hnum->GetNbinsX();

  // bin = 0;       underflow bin
  // bin = 1;       first bin with low-edge xlow INCLUDED
  // bin = nbins;   last bin with upper-edge xup EXCLUDED
  // bin = nbins+1; overflow bin

  for (Int_t bin=0; bin<=nbins+1; bin++) {
    
    Float_t numBinContent = hnum->GetBinContent(bin);
    Float_t numBinError   = hnum->GetBinError  (bin);

    Float_t denBinContent = hden->GetBinContent(bin);
    Float_t denBinError   = hden->GetBinError  (bin);

    Float_t ratioBinContent = Ratio     (numBinContent, denBinContent);
    Float_t ratioBinError   = RatioError(numBinContent, denBinContent,
					 numBinError,   denBinError);

    if (mean > -999) ratioBinContent /= mean;
    if (mean > -999) ratioBinError   /= mean;

    hratio->SetBinContent(bin, ratioBinContent);
    hratio->SetBinError  (bin, ratioBinError  );
  }

  return hratio;
}


//------------------------------------------------------------------------------
// GetPtRatio
//------------------------------------------------------------------------------
void GetPtRatio(TH1F*& hratio, TH1F* hnum, TH1F* hden, Float_t mean)
{
  hratio->Reset();

  Int_t nbins = hratio->GetNbinsX();

  // bin = 0;       underflow bin
  // bin = 1;       first bin with low-edge xlow INCLUDED
  // bin = nbins;   last bin with upper-edge xup EXCLUDED
  // bin = nbins+1; overflow bin

  for (Int_t bin=0; bin<=nbins; bin++) {
    
    Float_t numBinContent = hnum->GetBinContent(bin);
    Float_t numBinError   = hnum->GetBinError  (bin);

    Float_t denBinContent = hden->GetBinContent(bin);
    Float_t denBinError   = hden->GetBinError  (bin);

    Float_t ratioBinContent = Ratio     (numBinContent, denBinContent);
    Float_t ratioBinError   = RatioError(numBinContent, denBinContent,
					 numBinError,   denBinError);

    if (mean > -999) ratioBinContent /= mean;
    if (mean > -999) ratioBinError   /= mean;

    hratio->SetBinContent(bin, ratioBinContent);
    hratio->SetBinError  (bin, ratioBinError  );
  }
}


//------------------------------------------------------------------------------
// DrawRatioBinContent
//------------------------------------------------------------------------------
void DrawRatioBinContent(TH1F* hnum, TH1F* hden)
{
  TH1F* hratio = GetRatio(hnum, hden);

  TLatex* num;
  TLatex* den;

  Double_t maxError = -999;

  Char_t title[100];

  Int_t nbins = hnum->GetNbinsX();

  for (Int_t bin=1; bin<=nbins; bin++)
    if (hratio->GetBinError(bin) > maxError) maxError = hratio->GetBinError(bin);

  for (Int_t bin=1; bin<=nbins; bin++) {

    Double_t x      = hnum  ->GetBinCenter (bin);
    Double_t y      = hratio->GetBinContent(bin);
    Double_t yDelta = 1.5 * maxError;
      
    sprintf(title, "%.0f #pm %.0f", hnum->GetBinContent(bin), hnum->GetBinError(bin));

    num = new TLatex(x, y+yDelta, title);
      
    sprintf(title, "%.0f #pm %.0f", hden->GetBinContent(bin), hden->GetBinError(bin));
      
    den = new TLatex(x, y-yDelta, title);
      
    SetLatex(num, 12, 90, kRed  -7, 0.025);
    SetLatex(den, 32, 90, kAzure-7, 0.025);
  }
}


//------------------------------------------------------------------------------
// SetLatex
//------------------------------------------------------------------------------
void SetLatex(TLatex* latex,
	      Short_t align,
	      Float_t tangle,
	      Color_t tcolor,
	      Float_t tsize)
{

  latex->SetTextAlign( align);
  latex->SetTextAngle(tangle);
  latex->SetTextColor(tcolor);
  latex->SetTextFont (    42);
  latex->SetTextSize ( tsize);
  
  latex->Draw("same");
}



//------------------------------------------------------------------------------
// SetLatex
//------------------------------------------------------------------------------
void SetLatex(TLatex* latex, char* tex)
{

  latex = new TLatex(3.570061,23.08044,tex);

  latex->SetTextAlign( 13);
  latex->SetNDC();
  latex->SetTextFont (    42);
  latex->SetTextSize ( 20);
  latex->SetX(0.597855);
  latex->SetY(0.971928);



  latex->Draw("same");

 
}


//------------------------------------------------------------------------------
// ProjectY
//------------------------------------------------------------------------------
void ProjectY(TH2F* hist)
{
  if (!hist) {ExitPrint("ProjectY()"); return;}

  Char_t name[100];

  TH1F* hpx = (TH1F*)hist->ProjectionX();

  for (Int_t bin=1; bin<=hpx->GetNbinsX(); bin++) {

    Float_t xl = hpx->GetBinLowEdge(bin);
    Float_t xh = hpx->GetBinLowEdge(bin) + hpx->GetBinWidth(bin);

    sprintf(name, "%s_%.0f_%.0f", hist->GetName(), xl, xh);

    TH1F* hpy = (TH1F*)hist->ProjectionY(name, bin, bin);

    TCanvas* c1 = SetCanvas(name);
    DrawTH1F(c1, hpy, "hist", name);
  }
}


//------------------------------------------------------------------------------
// SetTH1F
//-----------------------------------------------------------------------------
TH1F* SetTH1F(TString  name,
	      Int_t    nbinsx,
	      Double_t xmin,
	      Double_t xmax,
	      Int_t    icharge)
{
  TH1F* h = new TH1F((name + charge[icharge]).Data(),
		     (name + charge[icharge]).Data(),
		     nbinsx, xmin, xmax);

  return h;
}


//------------------------------------------------------------------------------
// SetTH1F
//-----------------------------------------------------------------------------
TH1F* SetTH1F(TString         name,
	      Int_t           nbinsx,
	      const Double_t* xbins,
	      Int_t           icharge)
{
  TH1F* h = new TH1F((name + charge[icharge]).Data(),
		     (name + charge[icharge]).Data(),
		     nbinsx, xbins);

  return h;
}


//------------------------------------------------------------------------------
// SetTH2F
//-----------------------------------------------------------------------------
TH2F* SetTH2F(TString  name,
	      Int_t    nbinsx,
	      Double_t xmin,
	      Double_t xmax,
	      Int_t    nbinsy,
	      Double_t ymin,
	      Double_t ymax,
	      Int_t    icharge)
{
  TH2F* h = new TH2F((name + charge[icharge]).Data(),
		     (name + charge[icharge]).Data(),
		     nbinsx, xmin, xmax,
		     nbinsy, ymin, ymax);

  return h;
}


//------------------------------------------------------------------------------
// SetTH2F
//-----------------------------------------------------------------------------
TH2F* SetTH2F(TString         name,
	      Int_t           nbinsx,
	      const Double_t* xbins,
	      Int_t           nbinsy,
	      Double_t        ymin,
	      Double_t        ymax,
	      Int_t           icharge)
{
  TH2F* h = new TH2F((name + charge[icharge]).Data(),
		     (name + charge[icharge]).Data(),
		     nbinsx, xbins, nbinsy, ymin, ymax);

  return h;
}


//------------------------------------------------------------------------------
// DeltaPhi
//------------------------------------------------------------------------------
Double_t DeltaPhi(Double_t phi1, Double_t phi2)
{
  if (phi1 < 0) phi1 += 2*acos(-1);
  if (phi2 < 0) phi2 += 2*acos(-1);

  return (phi1 - phi2);
}


//------------------------------------------------------------------------------
// DeltaR
//------------------------------------------------------------------------------
Double_t DeltaR(Double_t eta1, Double_t eta2,
		Double_t phi1, Double_t phi2)
{
  return sqrt((eta1-eta2) * (eta1-eta2)
	      + DeltaPhi(phi1, phi2) * DeltaPhi(phi1, phi2));
}


//------------------------------------------------------------------------------
// SetTPaveText
//------------------------------------------------------------------------------
TPaveText* SetTPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TPaveText* t = new TPaveText(x1, y1, x2, y2, "ndc");

  t->SetBorderSize(   0);
  t->SetFillStyle (   0);
  t->SetTextAlign (  12);
  t->SetTextFont  (  42);
  t->SetTextSize  (0.03);

  return t;
}


//------------------------------------------------------------------------------
// Chi-squared of two histograms
//------------------------------------------------------------------------------
TH1F* Chi2Rebin(TH1F*     hist1,
		TH1F*     hist2,
		Double_t& chi2,
		Int_t&    NDF,
		Double_t  low,
		Double_t  high,
		Int_t     verbose,
		Int_t     threshold,
		Int_t     fitParams)
{
  TH1F* hChi2 = (TH1F*)hist1->Clone("hChi2");

  Double_t chi2Sum = 0.0;

  Double_t nEvts1     = 0;
  Double_t nEvts2     = 0;
  Int_t    nBins      = 0;
  Int_t    binStarted = 1;

  Double_t evtErr1    = 0;
  Double_t evtErr2    = 0;
  Double_t binChi2    = 0;

  Int_t minBin;
  Int_t maxBin; 
  
  if (low == 0 && high == 0) {
    minBin = 1;
    maxBin = hist1->GetNbinsX(); 
  }
  else {
    minBin = hist1->FindBin(low );
    maxBin = hist1->FindBin(high);
  }
  
  for (Int_t i=minBin; i<=maxBin; i++) {
    nEvts1 += (Double_t)hist1->GetBinContent(i);
    nEvts2 += (Double_t)hist2->GetBinContent(i);
    
    evtErr1 += (Double_t)hist1->GetBinError(i)*hist1->GetBinError(i);
    evtErr2 += (Double_t)hist2->GetBinError(i)*hist2->GetBinError(i);
    
    if (nEvts1 < threshold && i != hist1->GetNbinsX()) continue;
    
    binChi2 = (nEvts1 - nEvts2) * (nEvts1 - nEvts2);
    
    if (evtErr1 != 0 || evtErr2 != 0) binChi2 /= (evtErr1 + evtErr2);

    chi2Sum += binChi2;

    hChi2->SetBinContent(i, binChi2);

    if (verbose >= 1)
      printf("%3d)[%3d,%3d] data: %7.2f th: %7.2f  chi2: %7.2f total: %8.2f\n",
             nBins, binStarted, i, nEvts1, nEvts2, binChi2, chi2Sum); 
 
    binStarted = i+1;
    nEvts1     = 0.0;
    nEvts2     = 0.0;
    evtErr1    = 0.0;
    evtErr2    = 0.0;
    nBins++;
  }
  
  chi2 = chi2Sum;
  NDF  = nBins - 1 - fitParams;
  
  return hChi2;
}


//------------------------------------------------------------------------------
// fitGaus
//------------------------------------------------------------------------------
void fitGaus(TH1*      hist,
	     TF1*      f,
	     Double_t& mean,
	     Double_t& wide,
	     Double_t& meanErr,
	     Double_t& wideErr)
{
  f->SetParameter(0, hist->GetMaximum());
  f->SetParameter(1, hist->GetMean());
  f->SetParameter(2, hist->GetRMS());

  hist->Fit(f, "qr");

  mean    = f->GetParameter(1);
  meanErr = f->GetParError (1);

  wide    = f->GetParameter(2);
  wideErr = f->GetParError (2);
  
  if (wide < 0.0) {

    std::cout << "\n [fitHisto] Negative width" << std::endl;

    f->SetParameter(2, TMath::Abs(wide));

    hist->Fit(f, "r");

    wide    = f->GetParameter(2);
    wideErr = f->GetParError (2);
  }
}


//------------------------------------------------------------------------------
// spillOver
//------------------------------------------------------------------------------
TGraph* spillOver(TH2D*     matrix,
		  Double_t& result,
		  Bool_t    verbose)
{
  Int_t ncbins = matrix->GetNbinsX();

  Double_t offDiagonalFraction[ncbins];

  Double_t offDiagonalAverage     = 0.0;
  Double_t offDiagonalSumFraction = 0.0;

  TH1D* proj = matrix->ProjectionX();

  TGraph* g = new TGraph(ncbins);

  for (Int_t i=1; i<=ncbins; i++) {

    Double_t diagonal    = 0.0;
    Double_t offDiagonal = 0.0;

    for (Int_t j=1; j<=ncbins; j++) {
      
      if (j == i) diagonal     = matrix->GetBinContent(i,j);
      else        offDiagonal += matrix->GetBinContent(i,j);
    }

    offDiagonalFraction[i-1] = offDiagonal / (diagonal + offDiagonal);

    offDiagonalAverage += offDiagonalFraction[i-1];

    g->SetPoint(i-1, proj->GetBinCenter(i), offDiagonalFraction[i-1]);

  }


  // Cosmetics
  //----------------------------------------------------------------------------
  g->GetYaxis()->SetTitle("spill-over fraction");

  g->GetYaxis()->CenterTitle();
  g->GetYaxis()->SetTitleOffset(1.9);

  g->SetMarkerStyle(kFullCircle);
  g->SetTitle("");


  // Result
  //----------------------------------------------------------------------------
  offDiagonalAverage /= ncbins;

  for (Int_t i=0; i<ncbins; i++) {

    //    if (verbose) printf(" [draw] The spill-over for bin %2d is %f\n", i, offDiagonalFraction[i]);

    if (offDiagonalFraction[i] > offDiagonalAverage) offDiagonalSumFraction += (offDiagonalFraction[i] - offDiagonalAverage);
  }

  result = offDiagonalAverage + offDiagonalSumFraction;

  if (verbose) printf(" [draw] The spill-over (average + delta_i > average) for %s is %f\n",
		      matrix->GetName(), result);

  return g;
}


//------------------------------------------------------------------------------
// isShaftMuon
//------------------------------------------------------------------------------
Bool_t isShaftMuon(Float_t vx, Float_t vz)
{
  Bool_t thruShaft1 = InCircumference(vz, vx, -30, 0, 30);

  //  Bool_t thruShaft2 = InCircumference(vz, vx, 20, 32, 5);

  Bool_t thruShaft2        = InCircumference(vz, vx, Z_PM54/1e3,  X_PM54/1e3, R_PM54_WALL/1e3);
  Bool_t thruShaft2_mirror = InCircumference(vz, vx, Z_PM54/1e3, -X_PM54/1e3, R_PM54_WALL/1e3);
  
  return (thruShaft1 || thruShaft2 || thruShaft2_mirror);
}


//------------------------------------------------------------------------------
// SetHistGraphColors
//------------------------------------------------------------------------------
void SetHistGraphColors(TH1*    hist,
			TGraph* g,
			Color_t color,
			Style_t mstyle)
{
  hist->SetLineColor  (color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(mstyle);

  g->SetLineColor  (color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(mstyle);
}


//------------------------------------------------------------------------------
// ProduceGraph
//------------------------------------------------------------------------------
TGraphErrors* ProduceGraph(TH1D*       hist,
			   TH2D*       matrix,
			   const char* name,
			   Bool_t      verbose)
{
  if (verbose) printf("\n[draw::ProduceGraph] Producing graph from %s\n", hist->GetName());

  Int_t ncbins = hist->GetNbinsX();
  Int_t npbins = ncbins / 2;

  Double_t *x, *ex, *y, *ey;

  x  = new Double_t[npbins];
  ex = new Double_t[npbins];
  y  = new Double_t[npbins];
  ey = new Double_t[npbins];

  Double_t numer;
  Double_t denom;
  Double_t sigNum;
  Double_t sigDen;
  Double_t ratio;
  Double_t correlPart;
  
  Int_t iBin, jBin;
  
  for (Int_t i=1; i<=npbins; i++) {
    
    x [i-1] = -2.0 / (hist->GetBinLowEdge(i+1) + hist->GetBinLowEdge(i));

    ex[i-1] = 0.0;

    denom = hist->GetBinContent(i);
    numer = hist->GetBinContent(ncbins-i+1);

    sigDen = hist->GetBinError(i);
    sigNum = hist->GetBinError(ncbins-i+1);

    ratio = (denom > 0) ? numer / denom : 0.0;

    correlPart = 0;

    if (matrix) {

      for (Int_t k=1; k<=ncbins; k++) {
	
	iBin = matrix->GetBin(k, i);
	jBin = matrix->GetBin(k, ncbins-i+1);
	
	correlPart += 
	  matrix->GetBinContent(iBin) *
	  matrix->GetBinContent(jBin) *
	  hist->GetBinError(k) * hist->GetBinError(k);
      }
      
      correlPart *= -(2.0 / numer / denom);
    }
    
    //    correlPart = 0;

    y[i-1] = ratio;

    if (denom > 0)
      ey[i-1] = TMath::Abs(ratio) * TMath::Sqrt(sigNum*sigNum/numer/numer + sigDen*sigDen/denom/denom + correlPart);
    else
      ey[i-1] = 0.0;
    
    if (verbose)
      printf("[%d] %.0f +- %.0f / %.0f +- %.0f = %.3f +- %.3f\n",
	     i, numer, sigNum, denom, sigDen, ratio, ey[i-1]);
  }

  if (verbose) printf("\n");

  TGraphErrors* tmpGraph = new TGraphErrors(npbins, x, y, ex, ey);

  tmpGraph->SetName       (name);
  tmpGraph->SetLineColor  (kBlack);
  tmpGraph->SetMarkerColor(kBlack);
  tmpGraph->SetMarkerStyle(kFullCircle);

  return tmpGraph;
}


//------------------------------------------------------------------------------
// GetWeight1D
//------------------------------------------------------------------------------
Float_t GetWeight1D(Float_t x, TH1D* th1)
{
  if (x < th1->GetXaxis()->GetXmin() || x >= th1->GetXaxis()->GetXmax()) return 0;

  Int_t ibinx = th1->GetXaxis()->FindBin(x);

  return th1->GetBinContent(ibinx);
}


//------------------------------------------------------------------------------
// GetWeight2D
//------------------------------------------------------------------------------
Float_t GetWeight2D(Float_t x, Float_t y, TH2D* th2)
{
  if (x < th2->GetXaxis()->GetXmin() || x >= th2->GetXaxis()->GetXmax()) return 0;
  if (y < th2->GetXaxis()->GetXmin() || y >= th2->GetXaxis()->GetXmax()) return 0;

  Int_t ibinx = th2->GetXaxis()->FindBin(x);
  Int_t ibiny = th2->GetYaxis()->FindBin(y);

  return th2->GetBinContent(ibinx, ibiny);
}


//------------------------------------------------------------------------------
// InCircumference
//------------------------------------------------------------------------------
Bool_t InCircumference(Float_t       vz,
		       Float_t       vx,
		       const Float_t cz,
		       const Float_t cx,
		       const Float_t radius)
{
  Bool_t result = ((vz-cz)*(vz-cz) + (vx-cx)*(vx-cx) < (radius*radius));
  
  return result;
}
