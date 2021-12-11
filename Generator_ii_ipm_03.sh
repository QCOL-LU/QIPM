#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_II_IPM_03) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 4 8 16 32)
LO_Precision=(1e-1 1e-2 1e-3 1e-4 1e-5)

for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${LO_Precision[@]}; do #Change this
        	DIRC="/home/quantum/interior_point_methods/results/II_IPM_precision" #Change this
        	mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_precision${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
