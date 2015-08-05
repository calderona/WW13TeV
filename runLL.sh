if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
    echo "  "
    exit -1
fi


LUMINOSITY=5.0

NJETS=$1

CHANNELS="MuMu"

SAMPLES="WJetsFakes_Total"

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
	
	root -l -b -q "runLatinosTreeScriptLL.C($LUMINOSITY,$NJETS,\"$CHANNEL\",\"$SAMPLE\")"
  
    done

    OUTPATH=rootfiles/${NJETS}jet/${CHANNEL}

    # Merge some MC files
    hadd ${OUTPATH}/WW.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu.root 
    hadd ${OUTPATH}/WW_pow.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_pow.root 
    hadd ${OUTPATH}/WW_mc.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu_mcnlo.root
    hadd ${OUTPATH}/Wgamma.root ${OUTPATH}/WgammaNoStar.root ${OUTPATH}/WgammaStar.root ${OUTPATH}/GamGamWW.root 

done


