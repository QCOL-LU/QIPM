#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_ir_ls_01) #Change this
NORMRHS=(1 2 4 8 16 32)
ANCILLAES=(1 2 3 4 5 6 7 8 9)
DIMENSIONS=(2 3 4 5 6 7 8 9 10)
CONDITION=(64 128 256 512 1024 2048) #1 2 4 8 16 32
NORMA=(1 10 20 50 100)
EPSILON=(1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10)

DIRC="/home/quantum/QIPM/results/IR_LS_epsilon" #Change this
# mkdir -p ${DIRC}

for seed in {1..10}; do
    for runner in ${RUNNERS[@]}; do
        for param in ${EPSILON[@]}; do #Change this
    
            COMMAND="python3.9 ../runners/${runner}.py ${seed} ${param} > ${DIRC}/${runner}_seed${seed}_epsilon${param}.txt" #Change this folder and name
            eval $COMMAND
        done
    done
done
