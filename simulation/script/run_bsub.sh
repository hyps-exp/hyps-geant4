#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Usage : " $0 " [start run#]"
    exit 1
fi

###########################################################

bsub_opt="-q s"			# for a short job
# bsub_opt="-q l"			# for a long job

###########################################################

dirWork="$HOME/HYPS_Geant4/CFT_LEPS_pub/Geant4"

pathBin=$dirWork
pathLog="log"

###########################################################
bin="$pathBin/FilterAllRun.sh"
log="$pathLog/run0$1.log"

echo "$> bsub -o $log $bin $1 2>&1"
bsub $bsub_opt -o $log $bin $1 2>&1
