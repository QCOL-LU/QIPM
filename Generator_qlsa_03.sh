#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_qlsa_04) #Change this
NORMRHS=(1 2 4 8 16 32)
PRECISION=(1 5 10 50 100 500 1000)
DIMENSIONS=(2 3 4 5 6 7 8 9 10)
CONDITION=(1 2 4 8 16 32)
NORMA=(1 10 20 50 100)


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${PRECISION[@]}; do #Change this
            DIRC="/home/quantum/interior_point_methods/results/QLSA_precision" #Change this
            mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_precision${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done