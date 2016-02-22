#include "utils.h"



//------------------------------------------------------------------------------
// do efficiency plot
//------------------------------------------------------------------------------


TGraph* GetROC(TH1F* h_Sig, TH1F* h_Back, Style_t mstyle, Color_t color)
{

 
  Int_t nbins =  h_Sig->GetNbinsX();

  TGraph* roc = new TGraph(nbins+1);

  for (Int_t i=1; i<=nbins+1; i++) {

    Float_t x = 1e2 * h_Back->Integral(1,i) / h_Back->Integral();
    Float_t y = 1e2 * h_Sig->Integral(1,i) / h_Sig->Integral();


    roc->SetPoint(i-1, x, y);
  }

  roc->SetMarkerColor(color);
  roc->SetLineColor(color);
  roc->SetMarkerStyle(mstyle);

  return roc;
}



//------------------------------------------------------------------------------
// get eta projection 
//------------------------------------------------------------------------------


TH1F *getEtaProjection(TH2F *h_input, TH1F* h_out)

{

  double p_bin = h_input->GetNbinsX();
  double e_bin = h_input->GetNbinsY();


  for (int  e = 1; e < (e_bin+1); e++) {
    
    double auxSum = 0;

    for (int p = 1; p < (p_bin+1); p++) {
      
      auxSum = auxSum +  h_input->GetBinContent(p,e); 
  
    }

    h_out->SetBinContent(e, auxSum);
    h_out->SetBinError(e, sqrt(auxSum));
 }

  return h_out;

}






//------------------------------------------------------------------------------
// get pt projection 
//------------------------------------------------------------------------------


TH1F *getPtProjection(TH2F *h_input, TH1F* h_out)

{

  double p_bin = h_input->GetNbinsX();
  double e_bin = h_input->GetNbinsY();


  for (int p = 1; p < (p_bin+1); p++) {
    
    double auxSum = 0;

    for (int e = 1; e < (e_bin+1); e++) {
      
      auxSum = auxSum +  h_input->GetBinContent(p,e); 
  
    }

    h_out->SetBinContent(p, auxSum);
    h_out->SetBinError(p, sqrt(auxSum));
 }

  return h_out;

}



//------------------------------------------------------------------------------
// Divide 2 TH2D histos and give the binomial error
//------------------------------------------------------------------------------

TH2F *divide2DHistos(TH2F *h_den, TH2F *h_num ) 
{

  TH2F *h_out = (TH2F*) h_den->Clone();

  double p_bin = h_den->GetNbinsX();
  double e_bin = h_den->GetNbinsY();

 
  for (int p = 1; p < (p_bin+1); p++) {
    
    for (int e = 1; e < (e_bin+1); e++) {


      double den = h_den->GetBinContent(p,e);   
      double num = h_num->GetBinContent(p,e);
 
      if (num == 0) cout << "Numerator without entries" << endl; 
      if (den == 0) cout << "Denominator without entries" << endl; 

      double fr = num / den;
      

      //++ Compute binomial error per bin
      double fr_err = sqrt(num*(1-fr)) / den;

      //++ Compute normal propagation errors: Poisson statistics 
      //double fr_err = sqrt(num*(1+fr)) / den;  



      //++Fill 2D final FR histo
      h_out->SetBinContent(p,e,fr);
      h_out->SetBinError(p,e,fr_err);
    
    }
  }

  return h_out;

}


//------------------------------------------------------------------------------
// Print out 2D histo content 
//------------------------------------------------------------------------------

void print2D(TH2F *h) 
{

  double p_bin = h->GetNbinsX();
  double e_bin = h->GetNbinsY();

  cout << "---------------------------------" <<endl;  
  for (int p = 1; p < (p_bin+1); p++) {
    
   double low_p =  h->ProjectionX()->GetBinLowEdge(p);
   double high_p = h->ProjectionX()->GetBinLowEdge(p+1);

   cout << "Pt bin: " << low_p << "-" << high_p << endl; 

    for (int e = 1; e < (e_bin+1); e++) {

      double low_e =  h->ProjectionY()->GetBinLowEdge(e);
      double high_e = h->ProjectionY()->GetBinLowEdge(e+1);

      cout << "  Eta bin:\t" << low_e << "-" << high_e << "\t" << h->GetBinContent(p,e) << "\t" << h->GetBinError(p,e)  << endl;

    }
  }
  cout << "---------------------------------" <<endl;

}




//------------------------------------------------------------------------------
// Print out 2D histo content (only positive eta) 
//------------------------------------------------------------------------------

void print2D_PosEta(TH2F *h) 
{

  double p_bin = h->GetNbinsX();
  double e_bin = h->GetNbinsY();

  cout << "---------------------------------" <<endl;  
  
  for (int p = 1; p < (p_bin+1); p++) {
    
   double low_p =  h->ProjectionX()->GetBinLowEdge(p);
   double high_p = h->ProjectionX()->GetBinLowEdge(p+1);
   
   cout << "Pt bin: " << low_p << "-" << high_p << endl; 
   
    for (int e = 1; e < (e_bin+1); e++) {

      double low_e =  h->ProjectionY()->GetBinLowEdge(e);
      double high_e = h->ProjectionY()->GetBinLowEdge(e+1);

      // Double_t meanEta = (h->GetBinContent(p,e) + h->GetBinContent(p,(e_bin+1-e)))/2.;

      cout << "  Eta bin:\t" << low_e << "-" << high_e << "\t" << h->GetBinContent(p,e) << "\t" << h->GetBinError(p,e)  << endl;

      }
  }
  cout << "---------------------------------" <<endl;

}






//------------------------------------------------------------------------------
// Divide 2 TH2D histos and give the binomial error
//------------------------------------------------------------------------------

TH2F *minus2DHisto(TH2F *h_signal, TH2F *h_ewk) 
{

  TH2F *h_out = (TH2F*) h_signal->Clone();

  double p_bin = h_signal->GetNbinsX();
  double e_bin = h_signal->GetNbinsY();

 
  for (int p = 1; p < (p_bin+1); p++) {
    
    for (int e = 1; e < (e_bin+1); e++) {

      double signal = h_signal->GetBinContent(p,e);   
      double ewk = h_ewk->GetBinContent(p,e);

      double signal_err = h_signal->GetBinError(p,e);   
      double ewk_err = h_ewk->GetBinError(p,e);

 
      if (signal == 0) cout << "Signal without entries" << endl; 
      if (ewk  == 0) cout << "EWK without entries" << endl; 

      double minus = signal - ewk;
      
      //++ Compute errors bin-by-bin
      //----- double minus_err = sqrt(signal*ewk*(signal+ewk)); // OLD!!! 
      double minus_err = sqrt(signal_err*signal_err + ewk_err*ewk_err);


      //++Fill 2D final FR histo
      h_out->SetBinContent(p,e,minus);
      h_out->SetBinError(p,e,minus_err);
    
    }
  }

  return h_out;

}
