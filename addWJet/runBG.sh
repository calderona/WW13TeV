#
# This is a commented example GridEngine job for the IFCA cluster.
#

###################
# FLAG definition #
###################

# Request the Bourne Shell
#$ -S /bin/bash

# Change the job name to "hello_parallel"
#$ -N hxx_9

# We are using the "l.gaes" project
#$ -P l.gaes

#if you want to receive notifications when your job starts/ends
#$ -M calderon@ifca.unican.es

# Submit a parallel job of 16 slots
#$ -pe mpi 1

#$ -cwd
#$ -o output.out
#$ -e output.err


# Make a reservation. Since the job is paralell and it might be that all the
# resources are used, ensure that our resources are reserved.
#$ -R y


#################
# Actual script #
#################

#source run.sh

#root -l -b -q "runComputeFR.C()"

root -l -b -q "runComputeFR.C(\"WJets\", 1)" 

#root -l -b -q "runComputeFR.C(\"Top\", 1)"

#root -l -b -q "runComputeFR.C(\"ZJets\")"       

root -l -b -q "runComputeFR.C(\"QCD\")"

#root -l -b -q "runComputeFR.C(\"QCDEM\")"  