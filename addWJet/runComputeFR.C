void runComputeFR (TString theSample, bool doGen = false)
{

  //TProof *p = TProof::Open("");  

  gInterpreter->LoadMacro("computeFR.C+");

  //p->Exec(".L computeFR.C+"); 

  //p->Load("computeFR.C+");
  
  //p->SetParameter("theSample", theSample)
  //p->SetParameter("doGen", doGen)

  //p->Exec("computeFR(\"Top\", 1)");

  //p->GetOutputList();

  computeFR(theSample, doGen); 

}
