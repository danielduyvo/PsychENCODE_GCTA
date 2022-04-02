#!/usr/bin/env bash

# qsub options
#$ -w e
#$ -l h_rt=12:00:00,highp
#$ -pe shared 1
#$ -cwd
#$ -V
#$ -m a
#$ -M danieldu

julia scripts/run_cis_qsub.jl ${SGE_TASK_ID} $1
