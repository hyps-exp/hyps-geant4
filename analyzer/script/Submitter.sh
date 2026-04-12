#!/bin/sh

RUN_START=92
RUN_END=93
RUN_NUM_1JOB=1

i=$RUN_START
while [ $i -le $RUN_END ]
do
    ./run_bsub.sh $i $RUN_NUM_1JOB
    i=$(($i+$RUN_NUM_1JOB))
#    i=$(($i+40))
done
