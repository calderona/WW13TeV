void runLatinosTreeScriptLL(Double_t luminosity    = 12.103,
			  Int_t    jetChannel    = 1,
			  TString  flavorChannel = "MuE",
			  TString  theSample     = "DYtautau",
			  Bool_t   verbose       = false)
{
  gInterpreter->LoadMacro("LatinosTreeScriptLL.C+");

  LatinosTreeScriptLL(luminosity, jetChannel, flavorChannel, theSample, verbose);
}
