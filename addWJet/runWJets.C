void runWJets(string input="latino_ll_LP_test.root",
	      float luminosity = 1.0,
	      string indir = "input/",
	      string outdir = "output/",
	      bool addWeight=false)
{
	//-- Reset ROOT environment
	gROOT->Reset();
	gROOT->SetBatch();
	
	cout << "----------> Starting WJets <----------" << endl << endl;

	//-- Compile and run
	//if( !gSystem->CompileMacro( "addWJetsWeights.C"       ,"gk" ) ) return;
	gSystem->Load("addWJetsWeightsv8_C.so");
	addWJetsWeights(input,1,indir,outdir,true) ;
	
}
