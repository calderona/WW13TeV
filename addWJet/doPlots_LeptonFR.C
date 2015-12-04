void doPlots_LeptonFR (bool correction, TString theSample="WJets", TString lep ="Muon"){


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
  
  cout << "Computing muon FR ......" << endl;
  
  

  //------------------------------------------------------------------------------
  // Define plots to print
  //------------------------------------------------------------------------------
 
  // If you want to print out the results
  bool printout = 1;
   if (theSample == "data") { 
      if ( lep == "Muon" ) p_data_jetX = TFile::Open("rootfiles/dataDMu_jet20_FakeRate.root");
      if ( lep == "Electron" ) p_data_jetX = TFile::Open("rootfiles/dataDEG_jet20_FakeRate.root");
      if ( correction ) {
	p_Wjets = TFile::Open("rootfiles/WJets_Muons_jet20_FakeRate.root");
	p_Zjets = TFile::Open("rootfiles/ZJets_Muons_jet20_FakeRate.root");
      }
    }

  


  
  //------------------------------------------------------------------------------
  // Create output file
  //------------------------------------------------------------------------------
 
 
    if ( lep == "Muon" )      TFile *MuFR_Data_file = new TFile("MuFR_RunII_25ns_jet25_21Oct.root", "RECREATE");
    if ( lep == "Electron" )  TFile *EleFR_Data_file = new TFile("EGFR_RunII_25ns_jet25_21Oct.root", "RECREATE");
  
 
  
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
  tex1->SetTextSizePixels(15);
  
  TLatex *tex2 = new TLatex(3.570061,23.08044,"Wjets ~ 1.269 fb^{-1} at #sqrt{s}= 13 TeV, 25ns");
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetX(0.597855);
  tex2->SetY(0.971928);
  //tex2->SetTextFont(40);
  tex2->SetTextSizePixels(15);


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
 

  if ( lep == "Muon" ) {

    p_data_jetX->cd();

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
    

    Muon_FR_pT_jet15->Divide(Muon_signal_pT_jet15,Muon_fake_pT_jet15,1.,1.,"");
    Muon_FR_eta_jet15->Divide(Muon_signal_eta_jet15,Muon_fake_eta_jet15,1.,1.,"");
    Muon_FR_nPV_jet15->Divide(Muon_signal_nPV_jet15,Muon_fake_nPV_jet15,1.,1.,"");
    Muon_FR_jetEt_jet15->Divide(Muon_signal_jetEt_jet15,Muon_fake_jetEt_jet15,1.,1.,"");

    TH2F *Muon_fake_pT_eta_jet15 = (TH2F*) h_Muon_fake_pt_eta_bin->Clone();
    TH2F *Muon_signal_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
    TH2F *Muon_FR_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();

    if ( correction) {

      p_Wjets->cd();
      
      TH1F *wjets_fake_pT_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *wjets_signal_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();
      TH1F *wjets_FR_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();
      
      TH1F *WZjets_fake_pT_jet15= (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *WZjets_signal_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();

      //WZjets_fake_pT_jet15->Scale(0.9504*effLumi/intLumi);
      //WZjets_signal_pT_jet15->Scale(0.9504*effLumi/intLumi);

      TH1F *wjets_fake_eta_jet15 = (TH1F*) h_Muon_fake_eta_bin->Clone();
      TH1F *wjets_signal_eta_jet15 = (TH1F*) h_Muon_signal_eta_bin->Clone();
      TH1F *wjets_FR_eta_jet15 = (TH1F*) h_Muon_signal_eta_bin->Clone();

      TH2F *wjets_fake_pT_eta_jet15 = (TH2F*) h_Muon_fake_pt_eta_bin->Clone();
      TH2F *wjets_signal_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
      TH2F *wjets_FR_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();


      //wjets_fake_pT_eta_jet15->Scale(effLumi/intLumi);
      //wjets_signal_pT_eta_jet15->Scale(effLumi/intLumi);

      p_Zjets->cd();
      
      TH1F *zjets_fake_pT_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *zjets_signal_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();
      TH1F *zjets_FR_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();
      
      TH1F *WZjets_fake_pT_jet15= (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *WZjets_signal_pT_jet15= (TH1F*) h_Muon_signal_pT_bin->Clone();

      //WZjets_fake_pT_jet15->Scale(0.9504*effLumi/intLumi);
      //WZjets_signal_pT_jet15->Scale(0.9504*effLumi/intLumi);

      TH1F *zjets_fake_eta_jet15 = (TH1F*) h_Muon_fake_eta_bin->Clone();
      TH1F *zjets_signal_eta_jet15 = (TH1F*) h_Muon_signal_eta_bin->Clone();
      TH1F *zjets_FR_eta_jet15 = (TH1F*) h_Muon_signal_eta_bin->Clone();

      TH2F *zjets_fake_pT_eta_jet15 = (TH2F*) h_Muon_fake_pt_eta_bin->Clone();
      TH2F *zjets_signal_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
      TH2F *zjets_FR_pT_eta_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();

      //zjets_fake_pT_eta_jet15->Scale(effLumi/intLumi);
      //zjets_signal_pT_eta_jet15->Scale(effLumi/intLumi);

      p_data_jetX->cd();

      TH2F *Muon_fake_pT_eta_corrW_jet15 = (TH2F*) h_Muon_fake_pt_eta_bin->Clone(); 
      TH2F *Muon_signal_pT_eta_corrW_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
      TH2F *Muon_FR_pT_eta_corrW_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
      
      TH2F *Muon_fake_pT_eta_corrWZ_jet15 = (TH2F*) h_Muon_fake_pt_eta_bin->Clone();
      TH2F *Muon_signal_pT_eta_corrWZ_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
      TH2F *Muon_FR_pT_eta_corrWZ_jet15 = (TH2F*) h_Muon_signal_pt_eta_bin->Clone();
      
      TH1F *Muon_FR_pT_corrW_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *Muon_FR_eta_corrW_jet15 = (TH1F*) h_Muon_fake_eta_bin->Clone();
      TH1F *Muon_FR_pT_corrWZ_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *Muon_FR_eta_corrWZ_jet15 = (TH1F*) h_Muon_fake_eta_bin->Clone();
     
      TH1F *Muon_RelCorrection_signal_pT_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();
      TH1F *Muon_RelCorrection_fake_pT_jet15 = (TH1F*) h_Muon_fake_pT_bin->Clone();

    }
  }

  if (  lep == "Electron" ) {

    p_data_jetX->cd();

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
   

    Ele_FR_pT_jet15->Divide(Ele_signal_pT_jet15,Ele_fake_pT_jet15,1.,1.,"");
    Ele_FR_eta_jet15->Divide(Ele_signal_eta_jet15,Ele_fake_eta_jet15,1.,1.,"");
    Ele_FR_nPV_jet15->Divide(Ele_signal_nPV_jet15,Ele_fake_nPV_jet15,1.,1.,"");
    Ele_FR_jetEt_jet15->Divide(Ele_signal_jetEt_jet15,Ele_fake_jetEt_jet15,1.,1.,"");

    
    TH2F *Ele_fake_pT_eta_jet15 = (TH2F*) h_Ele_fake_pt_eta_bin->Clone();
    TH2F *Ele_signal_pT_eta_jet15 = (TH2F*) h_Ele_signal_pt_eta_bin->Clone();
    TH2F *Ele_FR_pT_eta_jet15 = (TH2F*) h_Ele_signal_pt_eta_bin->Clone();
    
  }



  //------------------------------------------------------------------------------
  // Make Plots
  //------------------------------------------------------------------------------
  
  if(!correction) {
    
    if (  lep == "Muon" ) {

      TCanvas * fake_pT_bin = new TCanvas("fake_pT_bin", "fake_pT_bin", 750, 750);
      fake_pT_bin->cd();
             
      DrawTH1F(fake_pT_bin, Muon_FR_pT_jet15,"E1","","pT (GeV/c)","#Tight / #Loose");
      LineOpt(Muon_FR_pT_jet15,kBlack,2,kSolid);
      MarkerOpt(Muon_FR_pT_jet15,kBlack,1,kFullCircle);
     
      DrawLegend(0.25,0.7,Muon_FR_pT_jet15 , "DataMu 2012 ", "LP",0.030,0.12, 0.10);
      
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
      
      DrawLegend(0.25,0.7,Muon_FR_nPV_jet15 , "Data SingleMu", "LP",0.030,0.12, 0.10);
      
      Muon_FR_nPV_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();


      TCanvas * fake_jetEt_bin = new TCanvas("fake_jetEt_bin", "fake_jetEt_bin", 750, 750);
      fake_jetEt_bin->cd();
      
      DrawTH1F(fake_jetEt_bin, Muon_FR_jetEt_jet15,"E1","","jetEt (GeV/c)","#Tight / #Loose");
      LineOpt(Muon_FR_jetEt_jet15,kBlack,2,kSolid);
      MarkerOpt(Muon_FR_jetEt_jet15,kBlack,1,kFullCircle);
      
      DrawLegend(0.25,0.7,Muon_FR_jetEt_jet15 , "Data SingleMu", "LP",0.030,0.12, 0.10);
      
      Muon_FR_jetEt_jet15->SetAxisRange(0,1,"Y");
      
      tex1->Draw();
      tex2->Draw();
      


      if (printout) {

	MuFR_Data_file ->cd();

	Muon_FR_pT_eta_jet15->Divide(Muon_signal_pT_eta_jet15, Muon_fake_pT_eta_jet15, 1., 1., "");

	cout << "-------- Results without MC corrections -------- " << endl;
	print2D_PosEta(Muon_FR_pT_eta_jet15);     

	Muon_FR_pT_eta_jet15->Write();

      
	// Output latex table: 
	
	int nrows =  7;
	int ncolumns = 4;
	
	cout << ">> reating table..." << endl;
	TResultsTable t(nrows, ncolumns, true);

	cout << ">> Setting column titles..." << endl; 
	t.SetColumnTitle(0, "0 < /eta <= 1.0");
	t.SetColumnTitle(1, "1.0 < /eta <= 1.479");
	t.SetColumnTitle(2, "1.479 < /eta <= 2.0");
	t.SetColumnTitle(3, "2.0 < /eta <= 2.5");
	
	cout << ">> Setting row titles..." << endl;
	//t.SetRowTitle(0, "10 < pt <= 15");
	t.SetRowTitle(0, "15 < pt <= 20");
	t.SetRowTitle(1, "20 < pt <= 25"); 
	t.SetRowTitle(2, "25 < pt <= 30"); 
	t.SetRowTitle(3, "30 < pt <= 35");
	t.SetRowTitle(4, "35 < pt <= 40");
	t.SetRowTitle(5, "40 < pt <= 45");
	t.SetRowTitle(6, "45 < pt <= 50");
	

	cout << ">> Filling table..." << endl;
	
	for (unsigned int i = 0; i < nrows ; i++) {
	
	  for (unsigned int j = 0; j < ncolumns; j++) {

	    t[i][j] = Muon_FR_pT_eta_jet15->GetBinContent(i+1,j+1);
	    t[i][j].SetError(Muon_FR_pT_eta_jet15->GetBinError(i+1,j+1));       

	  }
	}

	cout << ">> Saving to somefile.tex..." << endl;
	t.SaveAs("MuFR_RunII_25ns_jet25_21Oct.tex");
	t.Print();
	
      } // END PRINTOUT
      
    }
  

    if (  lep == "Electron" ) {

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

	Ele_FR_pT_eta_jet15->Divide(Ele_signal_pT_eta_jet15, Ele_fake_pT_eta_jet15, 1., 1., "");

	cout << "-------- Results without MC corrections -------- " << endl;
	print2D_PosEta(Ele_FR_pT_eta_jet15);     

	Ele_FR_pT_eta_jet15->Write();
      }
    }

  }


  if (correction) {

    if (  lep == "Muon" ) {

     Muon_fake_pT_eta_corrW_jet15->Add(wjets_fake_pT_eta_jet15,-1); 
     Muon_signal_pT_eta_corrW_jet15->Add(wjets_signal_pT_eta_jet15,-1);

     Muon_fake_pT_eta_corrWZ_jet15->Add(wjets_fake_pT_eta_jet15,-1); 
     Muon_signal_pT_eta_corrWZ_jet15->Add(wjets_signal_pT_eta_jet15,-1);
  
     Muon_fake_pT_eta_corrWZ_jet15->Add(zjets_fake_pT_eta_jet15,-1); 
     Muon_signal_pT_eta_corrWZ_jet15->Add(zjets_signal_pT_eta_jet15,-1);


     TH1F *Muon_signal_eta_corrW_jet15 = (TH1F*)  Muon_signal_pT_eta_corrW_jet15->ProjectionY()->Clone();
     TH1F *Muon_fake_pT_corrW_jet15 =  (TH1F*) Muon_fake_pT_eta_corrW_jet15->ProjectionX()->Clone();
     TH1F *Muon_fake_eta_corrW_jet15 =  (TH1F*) Muon_fake_pT_eta_corrW_jet15->ProjectionY()->Clone();
     TH1F *Muon_signal_pT_corrW_jet15 = (TH1F*)  Muon_signal_pT_eta_corrW_jet15->ProjectionX()->Clone();
     TH1F *Muon_fake_pT_corrWZ_jet15 =  (TH1F*) Muon_fake_pT_eta_corrWZ_jet15->ProjectionX()->Clone();
     TH1F *Muon_fake_eta_corrWZ_jet15 =  (TH1F*) Muon_fake_pT_eta_corrWZ_jet15->ProjectionY()->Clone();
     TH1F *Muon_signal_pT_corrWZ_jet15 =  (TH1F*) Muon_signal_pT_eta_corrWZ_jet15->ProjectionX()->Clone();
     TH1F *Muon_signal_eta_corrWZ_jet15 = (TH1F*)  Muon_signal_pT_eta_corrWZ_jet15->ProjectionY()->Clone();
   
     Muon_FR_pT_corrW_jet15->Divide(Muon_signal_pT_corrW_jet15 , Muon_fake_pT_corrW_jet15, 1., 1., "");
     Muon_FR_pT_corrWZ_jet15->Divide(Muon_signal_pT_corrWZ_jet15 , Muon_fake_pT_corrWZ_jet15, 1., 1., "");

     Muon_FR_eta_corrW_jet15->Divide(Muon_signal_eta_corrW_jet15 , Muon_fake_eta_corrW_jet15, 1., 1., "");
     Muon_FR_eta_corrWZ_jet15->Divide(Muon_signal_eta_corrWZ_jet15 , Muon_fake_eta_corrWZ_jet15, 1., 1., "");

     Muon_FR_pT_eta_corrW_jet15->Divide(Muon_signal_pT_eta_corrW_jet15, Muon_fake_pT_eta_corrW_jet15, 1., 1., "");
     Muon_FR_pT_eta_corrWZ_jet15->Divide(Muon_signal_pT_eta_corrWZ_jet15, Muon_fake_pT_eta_corrWZ_jet15, 1., 1., "");

     


     TCanvas * NUM_Correction_pT = new TCanvas("NUM_Correction_pT", "NUM_Correction_pT", 750, 750);
     NUM_Correction_pT->cd(); 

     //Muon_FR_pT_eta_corrWZ_jet15->Draw("COLZTEXT");

     WZjets_fake_pT_jet15->Add(zjets_fake_pT_jet15, +1);
     WZjets_signal_pT_jet15->Add(zjets_signal_pT_jet15, +1);

     DrawTH1F(NUM_Correction_pT,Muon_signal_pT_jet15 ,"","","pT (GeV/c)","");
     LineOpt(Muon_signal_pT_jet15,kBlue,2,kSolid);
     MarkerOpt(Muon_signal_pT_jet15,kBlue,1,kFullCircle);

     DrawTH1F(NUM_Correction_pT, WZjets_signal_pT_jet15,"same","","pT (GeV/c)","");
     LineOpt(WZjets_signal_pT_jet15,kRed,2,kSolid);
     MarkerOpt(WZjets_signal_pT_jet15,kRed,1,kFullCircle);

     DrawLegend(0.25,0.7,Muon_signal_pT_jet15, "data", "LP",0.030,0.12, 0.10);
     DrawLegend(0.25,0.6,WZjets_signal_pT_jet15, "EWK contribution", "LP",0.030,0.12, 0.10);

     tex1->Draw();
     tex2->Draw();

   

     TCanvas * NUM_RelCorrection_pT = new TCanvas("NUM_RelCorrection_pT", "NUM_RelCorrection_pT", 750, 750);
     NUM_RelCorrection_pT->cd(); 

 
     Muon_RelCorrection_signal_pT_jet15->Divide(WZjets_signal_pT_jet15,Muon_signal_pT_jet15, 1., 1., "");

     DrawTH1F(NUM_RelCorrection_pT, Muon_RelCorrection_signal_pT_jet15 ,"","","pT (GeV/c)","");
     LineOpt( Muon_RelCorrection_signal_pT_jet15,kBlue,2,kSolid);
     MarkerOpt( Muon_RelCorrection_signal_pT_jet15,kBlue,1,kFullCircle);
    

     tex1->Draw();
     tex2->Draw();



     TCanvas * DEN_RelCorrection_pT = new TCanvas("DEN_RelCorrection_pT", "DEN_RelCorrection_pT", 750, 750);
     DEN_RelCorrection_pT->cd(); 

 
     Muon_RelCorrection_fake_pT_jet15->Divide(WZjets_fake_pT_jet15,Muon_fake_pT_jet15, 1., 1., "");

     DrawTH1F(DEN_RelCorrection_pT, Muon_RelCorrection_fake_pT_jet15 ,"","","pT (GeV/c)","");
     LineOpt( Muon_RelCorrection_fake_pT_jet15,kBlue,2,kSolid);
     MarkerOpt( Muon_RelCorrection_fake_pT_jet15,kBlue,1,kFullCircle);
    

     tex1->Draw();
     tex2->Draw();




     TCanvas * DEN_Correction_pT = new TCanvas("DEN_Correction_pT", "DEN_Correction_pT", 750, 750);
     DEN_Correction_pT->cd(); 


     DrawTH1F(DEN_Correction_pT,Muon_fake_pT_jet15 ,"","","pT (GeV/c)","");
     LineOpt(Muon_fake_pT_jet15,kBlue,2,kSolid);
     MarkerOpt(Muon_fake_pT_jet15,kBlue,1,kFullCircle);

     DrawTH1F(DEN_Correction_pT,WZjets_fake_pT_jet15,"same","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(WZjets_fake_pT_jet15,kRed,2,kSolid);
     MarkerOpt(WZjets_fake_pT_jet15,kRed,1,kFullCircle);

     DrawLegend(0.25,0.7,Muon_fake_pT_jet15, "data", "LP",0.030,0.12, 0.10);
     DrawLegend(0.25,0.6,WZjets_fake_pT_jet15, "EWK contribution", "LP",0.030,0.12, 0.10);

     tex1->Draw();
     tex2->Draw();


     TCanvas * fake_Correction_pT = new TCanvas("fake_Correction_pT", "fake_Correction_pT", 750, 750);
     fake_Correction_pT->cd(); 


     DrawTH1F(fake_Correction_pT,Muon_FR_pT_jet15 ,"","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(Muon_FR_pT_jet15,kBlue,2,kSolid);
     MarkerOpt(Muon_FR_pT_jet15,kBlue,1,kFullCircle);

     DrawTH1F(fake_Correction_pT,Muon_FR_pT_corrWZ_jet15 ,"same","","pT (GeV/c)","#Tight / #Loose");
     LineOpt(Muon_FR_pT_corrWZ_jet15,kRed,2,kSolid);
     MarkerOpt(Muon_FR_pT_corrWZ_jet15,kRed,1,kFullCircle);

     DrawLegend(0.25,0.7,Muon_FR_pT_jet15, "Before EWK removal", "LP",0.030,0.12, 0.10);
     DrawLegend(0.25,0.6,Muon_FR_pT_corrWZ_jet15, "After EWK removal", "LP",0.030,0.12, 0.10);

     tex1->Draw();
     tex2->Draw();


     TCanvas * fake_Correction_eta = new TCanvas("fake_Correction_eta", "fake_Correction_eta", 750, 750);
     fake_Correction_eta->cd(); 


     DrawTH1F(fake_Correction_eta,Muon_FR_eta_jet15 ,"","","jetEt (GeV/c)","#Tight / #Loose");
     LineOpt(Muon_FR_eta_jet15,kBlue,2,kSolid);
     MarkerOpt(Muon_FR_eta_jet15,kBlue,1,kFullCircle);

     DrawTH1F(fake_Correction_eta,Muon_FR_eta_corrWZ_jet15 ,"same","","","#Tight / #Loose");
     LineOpt(Muon_FR_eta_corrWZ_jet15,kRed,2,kSolid);
     MarkerOpt(Muon_FR_eta_corrWZ_jet15,kRed,1,kFullCircle);

     DrawLegend(0.25,0.7,Muon_FR_eta_jet15, "Before EWK removal", "LP",0.030,0.12, 0.10);
     DrawLegend(0.25,0.6,Muon_FR_eta_corrWZ_jet15, "After EWK removal", "LP",0.030,0.12, 0.10);


     tex1->Draw();
     tex2->Draw();




     if (printout) {

       MuFR_Data_file->cd();

       cout << "-------- Results for jet 15, corrected by MC  -------- " << endl;
       print2D_PosEta(Muon_FR_pT_eta_corrWZ_jet15);

       Muon_FR_pT_eta_jet15->Write("FR_pT_eta");
       Muon_FR_pT_eta_corrWZ_jet15->Write("FR_pT_eta_EWKcorr");

       MuFR_Data_file->Close();

       int nrows =  7;
       int ncolumns = 4;

       cout << ">> Creating table..." << endl;
       TResultsTable t(nrows, ncolumns, true);

       cout << ">> Setting column titles..." << endl; 
       t.SetColumnTitle(0, "0 < /eta <= 1.0");
       t.SetColumnTitle(1, "1.0 < /eta <= 1.479");
       t.SetColumnTitle(2, "1.479 < /eta <= 2.0");
       t.SetColumnTitle(3, "2.0 < /eta <= 2.5");

       cout << ">> Setting row titles..." << endl;
       //t.SetRowTitle(0, "10 < pt <= 15");
       t.SetRowTitle(0, "15 < pt <= 20");
       t.SetRowTitle(1, "20 < pt <= 25"); 
       t.SetRowTitle(2, "25 < pt <= 30"); 
       t.SetRowTitle(3, "30 < pt <= 35");
       t.SetRowTitle(4, "35 < pt <= 40");
       t.SetRowTitle(5, "40 < pt <= 45");
       t.SetRowTitle(6, "45 < pt <= 50");

       cout << ">> Filling table..." << endl;

       for (unsigned int i = 0; i < nrows ; i++) {
	
	 for (unsigned int j = 0; j < ncolumns; j++) {

	   t[i][j] = Muon_FR_pT_eta_corrWZ_jet15->GetBinContent(i+1,j+1);
	   t[i][j].SetError(Muon_FR_pT_eta_corrWZ_jet15->GetBinError(i+1,j+1));
       

	 }
       }

       cout << ">> Saving to somefile.tex..." << endl;
       t.SaveAs("MuFR_RunII_25ns_jet25_21Oct_EWKcorr.tex");
       t.Print();
  
     }

   
    }
  }



} // END
