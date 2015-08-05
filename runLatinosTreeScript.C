void runLatinosTreeScript(Double_t luminosity    = 12.0,
			  Int_t    jetChannel    = 1,
			  TString  flavorChannel = "MuE",
			  TString  theSample     = "DYtautau",
			  Bool_t   verbose       = false)
{
  gInterpreter->LoadMacro("LatinosTreeScript.C+");

  LatinosTreeScript(luminosity, jetChannel, flavorChannel, theSample, verbose);
}
