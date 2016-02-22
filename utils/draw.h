#ifndef _DRAW_H
#define _DRAW_H (1)

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TImage.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TString.h"


const Color_t kDefaultColor   = kRed-7;
const Float_t default_xoffset =   1.30;
const Float_t default_yoffset =   1.80;
const Float_t default_xlow    =  0.195;
const Float_t default_yup     =  0.905;


// From /data/2b/kypreos/cosmicCheck/generating/CMSSW_3_1_1/src/GeneratorInterface/CosmicMuonGenerator/src/Point5terialMap.cc
//------------------------------------------------------------------------------
const Float_t Z_PX56 = -14000;
const Float_t Z_PM54 = Z_PX56 + 30100;
const Float_t Z_PM56 = Z_PX56 - 24227;

const Float_t X_PX56 = 0;
const Float_t X_PM54 = 26600;
const Float_t X_PM56 = 18220;

const Float_t R_PX56_AIR = 10250;
const Float_t R_PM54_AIR =  6050;
const Float_t R_PM56_AIR =  3550;

const Float_t R_PX56_WALL = 12400;
const Float_t R_PM54_WALL = 6050+2150/10250*6050 + 1800;
const Float_t R_PM56_WALL = 3550+2150/10250*3550 + 1800;


// Old
//------------------------------------------------------------------------------
const TString charge[]  = {"", "_positive", "_negative"};

const TString prefix    = "/home/piedra/plots/";
const TString craftFile = "/data/2a/piedra/cosmics_v1_4_2/out/craft_histograms.root";
const TString mcFile    = "/data/2a/piedra/cosmics_v1_4_2/out/craft_histograms_mc.root";

const Color_t cmbColor  = kDefaultColor;
const Color_t staColor  = kBlack;
const Color_t trkColor  = kAzure-7;
const Color_t segColor  = kGreen+2;

const Style_t cmbStyle  = kFullCircle;
const Style_t staStyle  = kFullTriangleDown;
const Style_t trkStyle  = kOpenSquare;
const Style_t segStyle  = kOpenTriangleUp;


// Functions
//------------------------------------------------------------------------------
TH1F*         GetTH1F            (TString            filename,
				  TString            name);

TH2F*         GetTH2F            (TString            filename,
				  TString            name);

void          DrawTH1F           (TPad*              canvas,
				  TH1*               hist,
				  TString            option,
				  TString            title   = "",
				  TString            xtitle  = "",
				  TString            ytitle  = "",
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          DrawTH1F           (TVirtualPad*              canvas,
				  TH1*               hist,
				  TString            option,
				  TString            title   = "",
				  TString            xtitle  = "",
				  TString            ytitle  = "",
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          DrawTH2F           (TPad*              canvas,
				  TH2F*              hist,
				  TString            option,
				  TString            title   = "",
				  TString            xtitle  = "",
				  TString            ytitle  = "",
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          DrawTHStack        (TPad*              canvas,
	                          THStack*           hist,
	                          TString            option,
	                          TString            title,
	                          TString            xtitle,
	                          TString            ytitle,
	                          Float_t            xoffset = default_xoffset,
	                          Float_t            yoffset = default_yoffset);


void          Axis               (TAxis*             axis,
				  TString            title,
				  Float_t            offset,
				  Bool_t             center);

void          AxisTH1F           (TH1*               hist,
				  TString            xtitle,
				  TString            ytitle,
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          AxisTH2F           (TH2F*              hist,
				  TString            xtitle,
				  TString            ytitle,
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          AxisTHStack        (THStack*           hist,
				  TString            xtitle,
				  TString            ytitle,
				  Float_t            xoffset = default_xoffset,
				  Float_t            yoffset = default_yoffset);

void          FillOpt            (TH1F*              hist,
				  Color_t            fcolor = kDefaultColor,
				  Style_t            fstyle = 1001);

void          LineOpt            (TH1F*              hist,
				  Color_t            lcolor = kDefaultColor,
				  Width_t            lwidth = 0,
				  Style_t            lstyle = kSolid);

void          LineOpt            (TH1*              hist,
				  Color_t            lcolor = kDefaultColor,
				  Width_t            lwidth = 0,
				  Style_t            lstyle = kSolid);

void          MarkerOpt          (TH1F*              hist,
				  Color_t            mcolor = kDefaultColor,
				  Size_t             msize  = 0,
				  Style_t            mstyle = kFullCircle);

void          MarkerOpt          (TH1*              hist,
				  Color_t            mcolor = kDefaultColor,
				  Size_t             msize  = 0,
				  Style_t            mstyle = kFullCircle);

void          ExitPrint          (TString            method);

void          MoveStats          (TCanvas*           canvas,
				  Float_t            xoffset = 2.0);

void          SetTF1             (TF1*               f,
				  TString            option,
				  Color_t            lcolor,
				  Style_t            lstyle);

void          SetTGraph          (TGraph*            g,
				  TString            option,
				  Color_t            mcolor,
				  Style_t            mstyle,
				  TString            title  = "",
				  TString            xtitle = "",
				  TString            ytitle = "");

TLegend*      DrawLegend         (Float_t            x1,
				  Float_t            y1,
				  TH1F*              hist,
				  TString            label,
				  TString            option,
				  Float_t            tsize   = 0.04,
				  Float_t            xoffset = 0.25,
				  Float_t            yoffset = 0.07);



TLegend*      DrawLegend         (Float_t            x1,
				  Float_t            y1,
				  TH1*              hist,
				  TString            label,
				  TString            option,
				  Float_t            tsize   = 0.04,
				  Float_t            xoffset = 0.25,
				  Float_t            yoffset = 0.07);

TLegend*      DrawLegendG         (Float_t            x1,
				   Float_t            y1,
				   TGraph*            hist,
				   TString            label,
				   TString            option,
				   Float_t            tsize   = 0.04,
				   Float_t            xoffset = 0.25,
				   Float_t            yoffset = 0.07);


void          DrawLogo           (TCanvas*           canvas,
				  TString            image,
				  Float_t            xlow = default_xlow,
				  Float_t            yup  = default_yup);

TCanvas*      SetCanvas          (TString            title,
				  Long_t             iCanvas = -1);

void          SaveCanvas         (TCanvas*           canvas,
				  TString            title,
				  TString            extension = ".gif",
				  Long_t             iCanvas   = -1);

Float_t       Ratio              (Float_t            a,
				  Float_t            b);

Float_t       RatioError         (Float_t            a,
				  Float_t            b,
				  Float_t            aErr = -999,
				  Float_t            bErr = -999);

TH1F*         GetRatio           (TH1F*              hnum,
				  TH1F*              hden,
				  Float_t            mean = -999);

void          GetPtRatio         (TH1F*&             hratio,
				  TH1F*              hnum,
				  TH1F*              hden,
				  Float_t            mean = -999);

void          DrawRatioBinContent(TH1F*              hnum,
				  TH1F*              hden);

void          SetLatex           (TLatex*            latex,
				  Short_t            align,
				  Float_t            tangle,
				  Color_t            tcolor,
				  Float_t            tsize);

void          SetLatex           (TLatex*           latex, 
				  char*              tex);

void          ProjectY           (TH2F*              hist);

TH1F*         SetTH1F            (TString            name,
				  Int_t              nbinsx,
				  Double_t           xmin,
				  Double_t           xmax,
				  Int_t              icharge = 0);

TH1F*         SetTH1F            (TString            name,
				  Int_t              nbinsx,
				  const Double_t*    xbins,
				  Int_t              icharge = 0);

TH2F*         SetTH2F            (TString            name,
				  Int_t              nbinsx,
				  Double_t           xmin,
				  Double_t           xmax,
				  Int_t              nbinsy,
				  Double_t           ymin,
				  Double_t           ymax,
				  Int_t              icharge = 0);

TH2F*         SetTH2F            (TString            name,
				  Int_t              nbinsx,
				  const Double_t*    xbins,
				  Int_t              nbinsy,
				  Double_t           ymin,
				  Double_t           ymax,
				  Int_t              icharge = 0);

Double_t      DeltaPhi           (Double_t           phi1,
				  Double_t           phi2);

Double_t      DeltaR             (Double_t           eta1,
				  Double_t           eta2,
				  Double_t           phi1,
				  Double_t           phi2);

TPaveText*    SetTPaveText       (Double_t           x1,
				  Double_t           y1,
				  Double_t           x2,
				  Double_t           y2);

TH1F*         Chi2Rebin          (TH1F*              hist1,
				  TH1F*              hist2,
				  Double_t&          chi2,
				  Int_t&             NDF,
				  Double_t           low       =  0.0,
				  Double_t           high      =  0.0,
				  Int_t              verbose   =    0,
				  Int_t              threshold = -999,
				  Int_t              fitParams =    0);

void          fitGaus            (TH1*               hist,
				  TF1*               f,
				  Double_t&          mean,
				  Double_t&          wide,
				  Double_t&          meanErr,
				  Double_t&          wideErr);

TGraph*       spillOver          (TH2D*              matrix,
				  Double_t&          result,
				  Bool_t             verbose = false);

Bool_t        isShaftMuon        (Float_t            vx,
				  Float_t            vz);

void          SetHistGraphColors (TH1*               hist,
				  TGraph*            g,
				  Color_t            color,
				  Style_t            mstyle);

TGraphErrors* ProduceGraph       (TH1D*              hist,
				  TH2D*              matrix,
				  const char*        name,
				  Bool_t             verbose);

Float_t       GetWeight1D        (Float_t            x,
				  TH1D*              th1);

Float_t       GetWeight2D        (Float_t            x,
				  Float_t            y,
				  TH2D*              th2);

Bool_t        InCircumference    (Float_t            z,
				  Float_t            x,
				  const Float_t      cz,
				  const Float_t      cx,
				  const Float_t      radius);


#endif
