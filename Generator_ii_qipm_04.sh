#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_II_QIPM_04) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 4 8 16 32)
LO_Precision=(1 2 4 6 8 10)


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${NORMRHS[@]}; do #Change this
        	DIRC="/home/quantum/interior_point_methods/results/II_QIPM_norm_b" #Change this
        	mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_normb${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
