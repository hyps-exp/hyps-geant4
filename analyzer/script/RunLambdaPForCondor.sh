#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Usage : RunLambdaPForCondor.sh [RUN#]"
    exit 1
fi

RUNNUM=$1
# RUNNUM=101

ANADIR=$HOME/HYPS_Geant4/CFT_LEPS_pub/Analyzer
DATADIR=$HOME/HYPS_Geant4/CFT_LEPS_pub/Analyzer/data
ROOTDIR=$HOME/HYPS_Geant4/CFT_LEPS_pub/Analyzer/rootfile

cd $ANADIR

IN_FILE=run$RUNNUM.dat
echo cp $DATADIR/$IN_FILE\.gz $DATADIR/tmp
cp -f $DATADIR/$IN_FILE\.gz $DATADIR/tmp

BIN=analysLambdaProtonScat

#if [ -f conf/analyzer.conf.CFT.$RUNNUM ] 
#then
#    CONF_FILE=conf/analyzer.conf.CFT.$RUNNUM
#else
#    CONF_FILE=conf/analyzer.conf.CFT
#fi
CONF_FILE=conf/analyzer.conf.CFT
ROOT_FILE=run$RUNNUM\_ana.root
# ROOT_FILE=run$RUNNUM\_ana\_100.root

echo gzip -d $DATADIR/tmp/$IN_FILE\.gz
gzip -d $DATADIR/tmp/$IN_FILE\.gz

echo ./bin/$BIN  $CONF_FILE  $DATADIR/tmp/$IN_FILE $ROOTDIR/$ROOT_FILE  
./bin/$BIN  $CONF_FILE  $DATADIR/tmp/$IN_FILE $ROOTDIR/$ROOT_FILE  

echo rm -f $DATADIR/tmp/$IN_FILE
rm -f $DATADIR/tmp/$IN_FILE
