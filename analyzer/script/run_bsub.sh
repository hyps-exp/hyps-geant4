#!/bin/sh

if [ $# -ne 2 ]; then
    echo "Usage : " $0 " [start run#] [run_num in one job]"
    exit 1
fi

###########################################################

bsub_opt="-q s"			# for a short job
# bsub_opt="-q l"			# for a long job

###########################################################

dirWork=$(cd $(dirname $0); pwd)/..

pathBin=$dirWork/script
pathLog="log"

###########################################################
bin="$pathBin/FilterSubmitter.sh"
log="$dirWork/$pathLog/run0$1.log"

echo "$> bsub -o $log $bin $1 $2 2>&1"
bsub $bsub_opt -o $log $bin $1 $2 2>&1
