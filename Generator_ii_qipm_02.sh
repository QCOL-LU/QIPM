#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_II_QIPM_02) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(10 50 100 200 500 1000)
CONDITION=(1 2 4 8 16 32)


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${DIMENSIONS[@]}; do #Change this
        	DIRC="/home/quantum/interior_point_methods/results/II_QIPM_var_num" #Change this
        	mkdir -p ${DIRC}
            COMMAND="python3.9 ${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_var${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
