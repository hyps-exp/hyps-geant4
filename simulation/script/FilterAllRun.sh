#!/bin/sh

if [ $# -ne 1 ]
    then
    echo Usage: FilterAllRun.sh [Run#1]
    exit 0
fi

./CFT run.mac conf/hyps.conf.default root/run$1.root data/run$1.dat
#./CFT_K0Sks run.mac.36 conf/sksm.conf.default.36 root/run$1.root data_dhcp26ad/run$1.dat

gzip -f data/run$1.dat

