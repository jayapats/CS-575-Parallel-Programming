#!/bin/bash
#SBATCH -J montecarlo
#SBATCH -A cs475-575
#SBATCH -p class
#SBATCH --gres=gpu:1
#SBATCH -o montecarlo.out
#SBATCH -e montecarlo.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jayapats@oregonstate.edu
for n in 16384 32768 64512 128000 256000 512000 1000448
do
 for b in 16 32 64 128
 do
	/usr/local/apps/cuda/cuda-10.1/bin/nvcc  -DNUMTRIALS=$n -DBLOCKSIZE=$b -o montecarlo montecarlo.cu
	./montecarlo
 done	
done