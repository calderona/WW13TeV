void runComputeFR(TString theSample)
{
  gInterpreter->LoadMacro("computeFR.C+");

  computeFR(theSample);
}
