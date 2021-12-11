#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_IR_II_QIPM_01) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 4 6 8 10)
CONDITION=(1 2 4 8 16 32)
METHODS=(0 )


for seed in {1..100}; do
    for runner in ${RUNNERS[@]}; do
        for method in ${METHODS[@]}; do #Change this
        	DIRC="/home/quantum/interior_point_methods/results/IR_II_QIPM_precision"
        	mkdir -p ${DIRC}
            if [[ "$method" -eq 1 ]]; then
                str_method="II_QIPM"
            else
                str_method="IR_II_QIPM"
            fi
            COMMAND="python3.9 ${runner}.py ${seed} ${method} > ${DIRC}/${runner}_${str_method}_seed${seed}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
