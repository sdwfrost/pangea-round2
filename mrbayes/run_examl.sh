#!/usr/bin/env bash
/usr/local/bin/parse-examl -s $1.phy -m DNA -q $2 -n $1
/usr/local/bin/raxmlHPC-PTHREADS-SSE3 -y -m GTRCAT -p 12345 -s $1.phy -n $1_startingtree
mpirun -np 2 /usr/local/bin/examl -D -t RAxML_parsimonyTree.$1_startingtree -m PSR -s $1.binary -n $1
