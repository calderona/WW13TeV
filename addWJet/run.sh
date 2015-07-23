SAMPLES="QCDEM     \
QCD                \  
WJets              \
ZJets"


for SAMPLE in $SAMPLES; do 
        
    root -l -b -q "runComputeFR.C(\"$SAMPLE\")"
  
done


