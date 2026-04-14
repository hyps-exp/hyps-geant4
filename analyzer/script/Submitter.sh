#!/bin/sh

RUN_START=1
RUN_END=10
RUN_NUM_1JOB=1

cd $(dirname $0)

i=$RUN_START
while [ $i -le $RUN_END ]
do
    ./run_bsub.sh $i $RUN_NUM_1JOB
    i=$(($i+$RUN_NUM_1JOB))
#    i=$(($i+40))
done
