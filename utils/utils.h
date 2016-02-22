#ifndef _UTILS_H
#define _UTILS_H (1)

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


// Functions
//------------------------------------------------------------------------------

TH2F*  divide2DHistos( TH2F *h_den, 
		       TH2F *h_num );

void  print2D( TH2F *h); 

void print2D_PosEta(TH2F *h);

TH2F* minus2DHisto( TH2F *h_signal, 
		    TH2F *h_ewk);

TH1F *getPtProjection (TH2F *h_input,
                       TH1F* h_out);

TH1F *getEtaProjection(TH2F *h_input, 
		       TH1F* h_out);

TGraph* GetROC(TH1F* h_Sig, 
               TH1F* h_Back, 
               Style_t mstyle, 
               Color_t color);

#endif
