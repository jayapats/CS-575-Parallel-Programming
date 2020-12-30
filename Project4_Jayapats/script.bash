#!/bin/bash
# number of subdivisions:
for s in 1024	2048	4096	8192	16384	32768	65536	131072	262144	524288	1048576	2097152
do
echo ArraySize = $s
g++ -DArraySize=$s Project41.cpp -o Project41 -lm -fopenmp
./Project41
done

