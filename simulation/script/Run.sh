#!/bin/sh

if [ $# -ne 2 ]
    then
    echo Usage: Run.sh [output root file] [output dat file]
    exit 0
fi

source ~/cern/geant4/install/bin/geant4.sh

SksMinus=./CFT
MacroFile=run.mac
ConfFile=conf/sksm.conf.default
HistFile=$1
DataFile=$2

$SksMinus $MacroFile $ConfFile $HistFile $DataFile

