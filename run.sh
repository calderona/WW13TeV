if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
    echo "  "
    exit -1
fi


LUMINOSITY=0.04003

NJETS=$1

CHANNELS="OF"

SAMPLES="data"          
#WWTo2L2Nu_pow          \
#WZ                     \
#ZZ                     \
#Top                    \
#WJets                  \
#DY "

#DataRun2012_Total  \
#ggWWto2L           \
#WWTo2L2Nu_mcnlo    \
#GamGamWW           \
#WZ                 \
#ZZ                 \
#TTbar              \
#TW                 \
#DY                 \
#WgammaStar         \
#WgammaNoStar       \
#Zgamma             \
#VVV                \
#WJets              \
#"

#rm -rf rootfiles/${NJETS}jet

# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	root -l -b -q "runLatinosTreeScript.C($LUMINOSITY,$NJETS,\"$CHANNEL\",\"$SAMPLE\")"
  
    done

    OUTPATH=rootfiles/${NJETS}jet/${CHANNEL}
    
done


