void doPlots_LeptonPR (TString theSample="MC", TString lep ="Muon"){


  // LOADING MACROS. 
  gROOT->LoadMacro("../utils/TResultsTable.C+");
  gInterpreter->LoadMacro   ("../utils/draw.C+");
  gInterpreter->LoadMacro   ("../utils/utils.C+"); 
  gInterpreter->ExecuteMacro("../utils/ChargeRatioStyle.C");
  gStyle      ->SetOptStat  (0);
  gStyle      ->SetPalette  (1);



  /// PRINTING !!! /// 
  cout << ">> Loading macro Table.C..." << endl;
  gROOT->LoadMacro("../utils/TResultsTable.C+");
  
  cout << "Computing muon PR ......" << endl;


  //------------------------------------------------------------------------------
  // Define plots to print
  //------------------------------------------------------------------------------
 
  // If you want to print out the results
  bool printout = 1;
  
  // If you want to do PU plots
  bool doPUplots = 1; 

 
  double intLumi = 5;
  double effLumi = 1.2;




  //------------------------------------------------------------------------------
  // Files and variables
  //------------------------------------------------------------------------------


  TString pwd1 = "";
  
  TFile* p_data_jetX = 0;
  TFile* p_Wjets= 0;
  TFile* p_Zjets= 0;

  // DATA -- samples
 
  if (theSample == "data") { 
    if ( lep == "Muon" ) p_data_jetX = TFile::Open("rootfiles/dataDMu_jet20_FakeRate.root");
    if ( lep == "Electron" ) p_data_jetX = TFile::Open("rootfiles/dataDEG_jet20_FakeRate.root");

  } else if ( theSample == "MC" ) {

    p_Zjets = TFile::Open("rootfiles/ZJets_Muons_jet20_FakeRate.root");
  }


  
  //------------------------------------------------------------------------------
  // Create output file
  //------------------------------------------------------------------------------
   
    if ( lep == "Muon" )      TFile *MuFR_Data_file = new TFile("ZJets_MuPR_RunII_25ns_jet2_21Oct.root", "RECREATE");
    if ( lep == "Electron" )  TFile *EleFR_Data_file = new TFile("EGPR_RunII_25ns_jet2_21Oct.root", "RECREATE");
 

  
  //------------------------------------------------------------------------------
  // LEGENDS
  //------------------------------------------------------------------------------


  TLatex *tex1 = new TLatex(3.570061,23.08044,"CMS Preliminary");
  //TLatex *tex1 = new TLatex(3.570061,23.08044,"InclusiveMu15");
  tex1->SetNDC();
  tex1->SetTextAlign(13);
  tex1->SetX(0.172922);
  tex1->SetY(0.971928);
  //tex1->SetTextFont(40);
  tex1->SetTextSizePixels(18);
  
  TLatex *tex2 = new TLatex(3.570061,23.08044,"Wjets ~ fb^{-1} at #sqrt{s}= 13 TeV, 25ns");
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetX(0.597855);
  tex2->SetY(0.971928);
  //tex2->SetTextFont(40);
  tex2->SetTextSizePixels(12);


  TLatex *tex3 = new TLatex(3.570061,23.08044,"");
  //TLatex *tex1 = new TLatex(3.570061,23.08044,"InclusiveMu15");
  tex3->SetNDC();
  tex3->SetTextAlign(13);
  tex3->SetX(0.172922);
  tex3->SetY(0.971928);
  //tex1->SetTextFont(40);
  tex3->SetTextSizePixels(18);


  TLatex *tex4 = new TLatex(3.570061,23.08044,"");
  //TLatex *tex1 = new TLatex(3.570061,23.08044,"InclusiveMu15");
  tex4->SetNDC();
  tex4->SetTextAlign(13);
  tex4->SetX(0.172922);
  tex4->SetY(0.971928);
  //tex1->SetTextFont(40);
  tex4->SetTextSizePixels(18);
  



  //------------------------------------------------------------------------------
  // Clone Histograms
  //------------------------------------------------------------------------------
 

  if (  lep == "Muon" ) {

    if (theSample == "data") { 
      p_data_jetX->cd();
    } else if (theSample == "MC") {
      p_Zjets->cd();
    }


    TH1F *Muon_fake_pT_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();
    TH1F *Muon_signal_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();
    TH1F *Muon_FR_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();
 
    
    TH1F *Muon_fake_eta_jet15 = (TH1F*) h_Muon_fake_eta_bin->Clone();
    TH1F *Muon_signal_eta_jet15 = (TH1F*) h_Muon_signal_eta_bin->Clone();
    TH1F *Muon_FR_eta_jet15 = (TH1F*) h_Muon_signal_eta_bin->Clone();

    TH1F *Muon_fake_nPV_jet15 = (TH1F*) h_Muon_fake_nPV->Clone();
    TH1F *Muon_signal_nPV_jet15 = (TH1F*) h_Muon_signal_nPV->Clone();
    TH1F *Muon_FR_nPV_jet15 = (TH1F*) h_Muon_signal_nPV->Clone();
    
    TH1F *Muon_fake_jetEt_jet15 = (TH1F*) h_Muon_fake_jetEt->Clone();
    TH1F *Muon_signal_jetEt_jet15 = (TH1F*) h_Muon_signal_jetEt->Clone();
    TH1F *Muon_FR_jetEt_jet15 = (TH1F*) h_Muon_signal_jetEt->Clone();
    

    Muon_FR_pT_jet15->Divide(Muon_signal_pT_jet15,Muon_fake_pT_jet15,1.,1.,"B");
    Muon_FR_eta_jet15->Divide(Muon_signal_eta_jet15,Muon_fake_eta_jet15,1.,1.,"B");
    Muon_FR_nPV_jet15->Divide(Muon_signal_nPV_jet15,Muon_fake_nPV_jet15,1.,1.,"B");
    Muon_FR_jetEt_jet15->Divide(Muon_signal_jetEt_jet15,Muon_fake_jetEt_jet15,1.,1.,"B");

    TH2F *Muon_fake_pT_eta_jet15 = (TH2F*) h_Muon_fake_pt_eta_bin->Clone();
    TH2F *Muon_signal_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
    TH2F *Muon_FR_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();

      
  }

  if (lep == "Electron" ) {
   
    if (theSample == "data") { 
      p_data_jetX->cd();
    } else if (theSample == "MC") {
      p_Zjets->cd();
    }


    TH1F *Ele_fake_pT_jet15 = (TH1F*) h_Ele_fake_pT_bin->Clone();
    TH1F *Ele_signal_pT_jet15= (TH1F*) h_Ele_signal_pT_bin->Clone();
    TH1F *Ele_FR_pT_jet15= (TH1F*) h_Ele_signal_pT_bin->Clone();
 
    
    TH1F *Ele_fake_eta_jet15 = (TH1F*) h_Ele_fake_eta_bin->Clone();
    TH1F *Ele_signal_eta_jet15 = (TH1F*) h_Ele_signal_eta_bin->Clone();
    TH1F *Ele_FR_eta_jet15 = (TH1F*) h_Ele_signal_eta_bin->Clone();

    TH1F *Ele_fake_nPV_jet15 = (TH1F*) h_Ele_fake_nPV->Clone();
    TH1F *Ele_signal_nPV_jet15 = (TH1F*) h_Ele_signal_nPV->Clone();
    TH1F *Ele_FR_nPV_jet15 = (TH1F*) h_Ele_signal_nPV->Clone();

    TH1F *Ele_fake_jetEt_jet15 = (TH1F*) h_Ele_fake_jetEt->Clone();
    TH1F *Ele_signal_jetEt_jet15 = (TH1F*) h_Ele_signal_jetEt->Clone();
    TH1F *Ele_FR_jetEt_jet15 = (TH1F*) h_Ele_signal_jetEt->Clone();
   

    Ele_FR_pT_jet15->Divide(Ele_signal_pT_jet15,Ele_fake_pT_jet15,1.,1.,"B");
    Ele_FR_eta_jet15->Divide(Ele_signal_eta_jet15,Ele_fake_eta_jet15,1.,1.,"B");
    Ele_FR_nPV_jet15->Divide(Ele_signal_nPV_jet15,Ele_fake_nPV_jet15,1.,1.,"B");
    Ele_FR_jetEt_jet15->Divide(Ele_signal_jetEt_jet15,Ele_fake_jetEt_jet15,1.,1.,"B");

    
    TH2F *Ele_fake_pT_eta_jet15 = (TH2F*) h_Ele_fake_pt_eta_bin->Clone();
    TH2F *Ele_signal_pT_eta_jet15 = (TH2F*) h_Ele_signal_pt_eta_bin->Clone();
    TH2F *Ele_FR_pT_eta_jet15 = (TH2F*) h_Ele_signal_pt_eta_bin->Clone();
    
  }



  //------------------------------------------------------------------------------
  // Make Plots
  //------------------------------------------------------------------------------
  
    if ( lep == "Muon" ) {

      TCanvas * fake_pT_bin = new TCanvas("fake_pT_bin", "fake_pT_bin", 750, 750);
      fake_pT_bin->cd();
             
      DrawTH1F(fake_pT_bin, Muon_FR_pT_jet15,"E1","","pT (GeV/c)","#Tight / #Loose");
      LineOpt(Muon_FR_pT_jet15,kBlack,2,kSolid);
      MarkerOpt(Muon_FR_pT_jet15,kBlack,1,kFullCircle);
     
      DrawLegend(0.25,0.7,Muon_FR_pT_jet15 , "Z+Jets", "LP",0.030,0.12, 0.10);
      
      Muon_FR_pT_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();

 
      TCanvas * fake_eta_bin = new TCanvas("fake_eta_bin", "fake_eta_bin", 750, 750);
      fake_eta_bin->cd();
      
      DrawTH1F(fake_eta_bin, Muon_FR_eta_jet15,"E1","","eta (GeV/c)","#Tight / #Loose");
      LineOpt(Muon_FR_eta_jet15,kBlack,2,kSolid);
      MarkerOpt(Muon_FR_eta_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Muon_FR_eta_jet15 , "Data SingleMu", "LP",0.030,0.12, 0.10);
      
      Muon_FR_eta_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();
     

      TCanvas * fake_nPV = new TCanvas("fake_nPV", "fake_nPV", 750, 750);
      fake_nPV->cd();
      
      DrawTH1F(fake_nPV, Muon_FR_nPV_jet15,"E1","","nPV (GeV/c)","#Tight / #Loose");
      LineOpt(Muon_FR_nPV_jet15,kBlack,2,kSolid);
      MarkerOpt(Muon_FR_nPV_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Muon_FR_nPV_jet15 , "Z+Jets", "LP",0.030,0.12, 0.10);
      
      Muon_FR_nPV_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();


      TCanvas * fake_jetEt_bin = new TCanvas("fake_jetEt_bin", "fake_jetEt_bin", 750, 750);
      fake_jetEt_bin->cd();
      
      DrawTH1F(fake_jetEt_bin, Muon_FR_jetEt_jet15,"E1","","jetEt (GeV/c)","#Tight / #Loose");
      LineOpt(Muon_FR_jetEt_jet15,kBlack,2,kSolid);
      MarkerOpt(Muon_FR_jetEt_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Muon_FR_jetEt_jet15 , "Z+Jets", "LP",0.030,0.12, 0.10);
      
      Muon_FR_jetEt_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();
      


      if (printout) {

	MuFR_Data_file ->cd();

	Muon_FR_pT_eta_jet15->Divide(Muon_signal_pT_eta_jet15, Muon_fake_pT_eta_jet15, 1., 1., "B");

	cout << "-------- Results without MC corrections -------- " << endl;
	print2D_PosEta(Muon_FR_pT_eta_jet15);     

	Muon_FR_pT_eta_jet15->Write();

      
	// Output latex table: 
	
	int nrows =  6;
	int ncolumns = 4;
	
	cout << ">> Creating table..." << endl;
	TResultsTable t(nrows, ncolumns, true);

	cout << ">> Setting column titles..." << endl; 
	t.SetColumnTitle(0, "0 < /eta <= 1.0");
	t.SetColumnTitle(1, "1.0 < /eta <= 1.479");
	t.SetColumnTitle(2, "1.479 < /eta <= 2.0");
	t.SetColumnTitle(3, "2.0 < /eta <= 2.5");
	
	cout << ">> Setting row titles..." << endl;
	t.SetRowTitle(0, "10 < pt <= 15");
	t.SetRowTitle(1, "15 < pt <= 20");
	t.SetRowTitle(2, "20 < pt <= 25"); 
	t.SetRowTitle(3, "25 < pt <= 30"); 
	t.SetRowTitle(4, "30 < pt <= 35");
	t.SetRowTitle(5, "35 < pt <= 40");
	

	cout << ">> Filling table..." << endl;
	
	for (unsigned int i = 0; i < nrows ; i++) {
	
	  for (unsigned int j = 0; j < ncolumns; j++) {

	    t[i][j] = Muon_FR_pT_eta_jet15->GetBinContent(i+1,j+1);
	    t[i][j].SetError(Muon_FR_pT_eta_jet15->GetBinError(i+1,j+1));       

	  }
	}

	cout << ">> Saving to somefile.tex..." << endl;
	t.SaveAs("MuFR_Moriond13_jet05.tex");
	t.Print();
	
      } // END PRINTOUT
      
    }


    if ( lep == "Electron" ) {

      TCanvas * fake_pT_bin = new TCanvas("fake_pT_bin", "fake_pT_bin", 750, 750);
      fake_pT_bin->cd();
             
      DrawTH1F(fake_pT_bin, Ele_FR_pT_jet15,"E1","","pT (GeV/c)","#Tight / #Loose");
      LineOpt(Ele_FR_pT_jet15,kBlack,2,kSolid);
      MarkerOpt(Ele_FR_pT_jet15,kBlack,1,kFullCircle);
     
      DrawLegend(0.25,0.7,Ele_FR_pT_jet15 , "DataMu 2012 ", "LP",0.030,0.12, 0.10);
      
      Ele_FR_pT_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();

 
      TCanvas * fake_eta_bin = new TCanvas("fake_eta_bin", "fake_eta_bin", 750, 750);
      fake_eta_bin->cd();
      
      DrawTH1F(fake_eta_bin, Ele_FR_eta_jet15,"E1","","eta (GeV/c)","#Tight / #Loose");
      LineOpt(Ele_FR_eta_jet15,kBlack,2,kSolid);
      MarkerOpt(Ele_FR_eta_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Ele_FR_eta_jet15 , "Data SingleMu", "LP",0.030,0.12, 0.10);
      
      Ele_FR_eta_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();

      TCanvas * fake_nPV = new TCanvas("fake_nPV", "fake_nPV", 750, 750);
      fake_nPV->cd();
      
      DrawTH1F(fake_nPV, Ele_FR_nPV_jet15,"E1","","nPV (GeV/c)","#Tight / #Loose");
      LineOpt(Ele_FR_nPV_jet15,kBlack,2,kSolid);
      MarkerOpt(Ele_FR_nPV_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Ele_FR_nPV_jet15 , "Data SingleMu", "LP",0.030,0.12, 0.10);
      
      Ele_FR_nPV_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();
     

      TCanvas * fake_jetEt_bin = new TCanvas("fake_jetEt_bin", "fake_jetEt_bin", 750, 750);
      fake_jetEt_bin->cd();
      
      DrawTH1F(fake_jetEt_bin, Ele_FR_jetEt_jet15,"E1","","jetEt (GeV/c)","#Tight / #Loose");
      LineOpt(Ele_FR_jetEt_jet15,kBlack,2,kSolid);
      MarkerOpt(Ele_FR_jetEt_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Ele_FR_jetEt_jet15 , "Data SingleMu", "LP",0.030,0.12, 0.10);
      
      Ele_FR_jetEt_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();
      

      if (printout) {

	EleFR_Data_file ->cd();

	Ele_FR_pT_eta_jet15->Divide(Ele_signal_pT_eta_jet15, Ele_fake_pT_eta_jet15, 1., 1., "B");

	cout << "-------- Results without MC corrections -------- " << endl;
	print2D_PosEta(Ele_FR_pT_eta_jet15);     

	Ele_FR_pT_eta_jet15->Write();
      }
    }

  
       


} // END
