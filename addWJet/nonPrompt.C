#include "TFile.h"                                         
#include "TH1F.h"                  



//------------------------------------------------------------------------------                                                               
// Top                                                                                                                                         
//------------------------------------------------------------------------------                                                               
void nonPrompt(Int_t     njet = 0,
	       TString   channel = "OF", 
	       TString   directory = "rootfiles") { 


 TString path = Form("%s/%djet/%s/", directory.Data(), njet, channel.Data());    

 //TFile* inputWj   = new TFile(path + "fakes_total.root");            
 TFile* inputWj   = new TFile(path + "TopFakes_Total.root");            

 TH1F* hWjetsTT_pt1               = (TH1F*) inputWj->Get("hPtLepton1WWLevelA");                     
 TH1F* hWjetsLL_pt1               = (TH1F*) inputWj->Get("hPtLepton1WWLevelLL");       


 hWjetsTT_pt1->RebinX(5);
 hWjetsLL_pt1->RebinX(5);

 hWjetsLL_pt1->SetLineColor(2);

 hWjetsTT_pt1->Draw("E1");
 hWjetsLL_pt1->Draw("histsame");




}
