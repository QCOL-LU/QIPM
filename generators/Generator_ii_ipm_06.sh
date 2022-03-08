#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_II_IPM_06) #Change this
NORMA=(1 2 4 6 8 10)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 4 8 16 32)
LO_Precision=(1 2 4 6 8 10)
SEEDS=(11 23 38 39 81 91 9)
OMEGA=(5 10 50 100 500 1000)


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${OMEGA[@]}; do #Change this
        	DIRC="/home/quantum/interior_point_methods/results/II_IPM_omega" #Change this
        	mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_omega${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
