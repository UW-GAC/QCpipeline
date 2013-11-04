#! /bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -o $JOB_NAME.o$JOB_ID

args=("$@") # all arguments
unset args[0] # remove first argument (R script name)

len=${#args[@]} # number of arguments

export R_LIBS=/projects/geneva/gcc-fs2/R_packages/library

if [ $len -gt 0 ]; then
    R -q --vanilla --args ${args[@]} < $1
else
    R -q --vanilla < $1
fi
