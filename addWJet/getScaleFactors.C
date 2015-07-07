//------------------------------------------------------------------------------
//
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/MuFR_Moriond13_jet15_EWKcorr.root             ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/MuFR_Moriond13_jet20_EWKcorr.root             ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/MuFR_Moriond13_jet30_EWKcorr.root             ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/MuFR_Moriond13_jet50_EWKcorr.root             ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/fakerates_trigger_Moriond13_ewksub.root       ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/fakerates_trigger_Moriond13_ewksub_jet15.root ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/work/c/calderon/public/LatinoTrees/leptonEffAndFakesHWW/fakerates_trigger_Moriond13_ewksub_jet50.root ~/www/AuxiliaryFilesWZXS8TeV/.
//
// scp lxplus.cern.ch:/afs/cern.ch/user/h/hbrun/www/promptRateAll2012/prompt_rateMuons.root     ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/user/h/hbrun/www/promptRateAll2012/prompt_rateElectrons.root ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/user/h/hbrun/www/SFbyPeriod/allMuonsSF.root                  ~/www/AuxiliaryFilesWZXS8TeV/.
// scp lxplus.cern.ch:/afs/cern.ch/user/h/hbrun/www/SFbyPeriod/allElectronsSF.root              ~/www/AuxiliaryFilesWZXS8TeV/.
//
//------------------------------------------------------------------------------


const TString path = "../auxiliaryFiles8TeV/";


const Int_t ptBinsSizeElFake = 5;

Float_t ptBinsElFake[ptBinsSizeElFake+1] = {10, 15, 20, 25, 30, 35};

const Int_t etaBinsSizeElFake = 4;

Float_t etaBinsElFake[etaBinsSizeElFake+1] = {0.0, 1.0, 1.479, 2.0, 2.5};


void getScaleFactors()
{
  MakeElectronFR("fakerates_trigger_Moriond13_ewksub_jet15");
  MakeElectronFR("fakerates_trigger_Moriond13_ewksub_jet35");
  MakeElectronFR("fakerates_trigger_Moriond13_ewksub_jet50");

  Invert("prompt_rateMuons_Moriond13",     "effDATA_prompt_rate",      "MuPR");
  Invert("prompt_rateElectrons_Moriond13", "effDATA_All_selec",        "ElePR");
  //Invert("allMuonsSF",           "muonsDATAMCratio_all",     "MuSF");
  //Invert("allElectronsSF",       "electronsDATAMCratio_all", "EleSF");
}


//------------------------------------------------------------------------------
// MakeElectronFR
//------------------------------------------------------------------------------
void MakeElectronFR(TString filename)
{
  TString jetpt = "jet35";
 
  if (filename.Contains("jet15")) jetpt = "jet15";
  if (filename.Contains("jet50")) jetpt = "jet50";
 
  TH2F* fakeElH2 = new TH2F("fakeElH2_" + jetpt,
			    "fakeElH2_" + jetpt,
			    ptBinsSizeElFake, ptBinsElFake,
			    etaBinsSizeElFake, etaBinsElFake);

  TFile* input = new TFile(path + filename + ".root", "read");

  TH1F* histoFR_B1 = (TH1F*)input->Get("TrgElenewWPHWWPtBarrel1");
  TH1F* histoFR_B2 = (TH1F*)input->Get("TrgElenewWPHWWPtBarrel2");
  TH1F* histoFR_E1 = (TH1F*)input->Get("TrgElenewWPHWWPtEndcap1");
  TH1F* histoFR_E2 = (TH1F*)input->Get("TrgElenewWPHWWPtEndcap2");

  for (UInt_t ptBin=4; ptBin<4+ptBinsSizeElFake; ptBin++) {

    Float_t valueB1 = histoFR_B1->GetBinContent(ptBin);
    Float_t valueB2 = histoFR_B2->GetBinContent(ptBin);
    Float_t valueE1 = histoFR_E1->GetBinContent(ptBin);
    Float_t valueE2 = histoFR_E2->GetBinContent(ptBin);  
    
    fakeElH2->SetBinContent(ptBin-3, 1, valueB1);
    fakeElH2->SetBinContent(ptBin-3, 2, valueB2);
    fakeElH2->SetBinContent(ptBin-3, 3, valueE1);
    fakeElH2->SetBinContent(ptBin-3, 4, valueE2);
  }

  TFile* output = new TFile(path + "EleFR_Moriond13_" + jetpt + "_EWKcorr.root", "recreate");

  fakeElH2->Write("fakeElH2");

  output->Close();

  input->Close();
}


//------------------------------------------------------------------------------
// InvertAxis
//------------------------------------------------------------------------------
TH2F* InvertAxis(TH2* h2, TString hname)
{
  TArrayD* xarray = (TArrayD*)h2->GetXaxis()->GetXbins();
  TArrayD* yarray = (TArrayD*)h2->GetYaxis()->GetXbins();

  const UInt_t xsize = xarray->GetSize();
  const UInt_t ysize = yarray->GetSize();

  Double_t* xbins = new Double_t[xsize];
  Double_t* ybins = new Double_t[ysize];
  
  for (UInt_t i=0; i<xsize; i++) xbins[i] = xarray->GetAt(i);
  for (UInt_t i=0; i<ysize; i++) ybins[i] = yarray->GetAt(i);

  TH2F* h2inverted = new TH2F(hname, hname,
			      h2->GetNbinsY(), ybins,
			      h2->GetNbinsX(), xbins);

  for (UInt_t i=1; i<=h2inverted->GetNbinsX(); i++) {
    for (UInt_t j=1; j<=h2inverted->GetNbinsY(); j++) {

      h2inverted->SetBinContent(i, j, h2->GetBinContent(j, i));
      h2inverted->SetBinError  (i, j, h2->GetBinError  (j, i));
    }
  }

  delete [] xbins;
  delete [] ybins;

  return h2inverted;
}


//------------------------------------------------------------------------------
// Invert
//------------------------------------------------------------------------------
void Invert(TString filename, TString hname, TString prefix)
{
  TFile* input = new TFile(path + filename + ".root", "read");

  TH2F* h2 = (TH2F*)input->Get(hname);

  TH2F* h2inverted = InvertAxis(h2, "h2inverted");

  TFile* output = new TFile(path + prefix + "_2012.root", "recreate");

  h2inverted->Write();

  output->Close();

  input->Close();
}
