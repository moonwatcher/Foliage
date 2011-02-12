#!/bin/bash
ARGML="/software/jdk1.6.0_01/bin/java -Djava.awt.headless=true -Xmx1000M -Xmn1000M -jar /lustre/scratch1/sanger/lg8/argml/trunk/argmllib.jar"

echo 'zcat '$1'/'$3'.perm.stats.gz |' $ARGML 'msd' '--I' $1 '--O' $2 '--tl 2000 --pt "[0-9]+_A1B8"| gzip -c > ' $2'/'$3'.msd.gz'|bash

