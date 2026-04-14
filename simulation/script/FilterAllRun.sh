#!/bin/sh

if [ $# -ne 1 ]
    then
    echo Usage: FilterAllRun.sh [Run#1]
    exit 0
fi

cd $(dirname $0)/..

./bin/G4Hyps ../param_sim/conf/geant4_hyps.conf root/run$1.root data/run$1.dat run.mac 

gzip -f data/run$1.dat

