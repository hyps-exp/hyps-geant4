#!/bin/sh

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
    echo "Usage : FilterSubmitter.sh [START RUN#] [RUN_NUM_IN_ONE_JOB]"
    exit 1
fi

RUNNUM=$1
JOB_NUM=1

if [ $# = 2 ]; then
    JOB_NUM=$2
fi

i=1
while [ $i -le $JOB_NUM ]
do
    ./RunLambdaPForCondor.sh $RUNNUM
    RUNNUM=$(($RUNNUM+1))
    i=$(($i+1))
done
