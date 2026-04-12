#!/bin/sh

RUN_START=306
RUN_END=306

# GEANT_INSTALL=/home/had/miwaq/cern/geant4.10.05/install/bin
# source $GEANT_INSTALL/geant4.sh

i=$RUN_START
while [ $i -le $RUN_END ]
do
    ./run_bsub.sh $i
    i=$(($i+1))
done
