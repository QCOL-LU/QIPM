#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_qlsa_01) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 3 4 5 6 7 8 9 10)
CONDITION=(64 128 256 512 1024 2048) #1 2 4 8 16 32
NORMA=(1 10 20 50 100)


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${CONDITION[@]}; do #Change this
        	DIRC="/home/quantum/QIPM/results/QLSA_cond_A" #Change this
        	mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_condA${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
