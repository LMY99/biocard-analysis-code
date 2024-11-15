#!usr/bin/env bash

sbatch -J "real_data_20" --time=7-00:00:00 --ntasks 1 --cpus-per-task 11 run_real.sh
