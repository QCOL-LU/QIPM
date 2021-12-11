#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_II_IPM_01) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 4 6 8 10)
CONDITION=(1 2 4 8 16)  #


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${CONDITION[@]}; do #Change this
        	DIRC="/home/quantum/interior_point_methods/results/II_IPM_A_cond"
        	mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_cond${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
