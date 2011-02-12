#!/bin/bash
INPUT=$1
PATTERN=$2
OUTPUT=$3

echo 'zcat '$INPUT'| ./foliage sfilter' '--pt' $PATTERN '| gzip -c >' $OUTPUT |bash

