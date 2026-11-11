#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Usage : " $0 " [run#]"
    exit 1
fi

cd $(dirname $0)

./run_bsub_analysis-tools.sh $1 1
