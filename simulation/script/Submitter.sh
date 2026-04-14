#!/bin/sh

RUN_START=1
RUN_END=10

cd $(dirname $0)

i=$RUN_START
while [ $i -le $RUN_END ]
do
    ./run_bsub.sh $i
    i=$(($i+1))
done
