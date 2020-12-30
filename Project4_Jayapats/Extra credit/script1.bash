#!/bin/bash
# number of threads:
for t in 1 2 3 4 6
do
echo NUMT = $t
# number of subdivisions:
for s in 1024	2048	4096	8192	16384	32768	65536	131072	262144	524288	1048576	2097152
do
echo ArraySize = $s
g++ -DArraySize=$s -DNUMT=$t project4.cpp -o project4 -lm -fopenmp
./project4
done
done
