#!/bin/sh

if [ $# -ne 2 ]; then
    echo "Usage : " $0 " [start run#] [run_num in one job]"
    exit 1
fi

###########################################################

bsub_opt="-q s"			# for a short job
# bsub_opt="-q l"			# for a long job

###########################################################

dirWork="$HOME/HYPS_Geant4/CFT_LEPS_pub/Analyzer"

pathBin=$dirWork
pathLog="log"

###########################################################
bin="$pathBin/FilterSubmitter.sh"
log="$pathLog/run0$1.log"

echo "$> bsub -o $log $bin $1 $2 2>&1"
bsub $bsub_opt -o $log $bin $1 $2 2>&1
