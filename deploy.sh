#!/usr/bin/env bash

for cv in {1..5}; do
  if ! test -f "S_cv${cv}.rda"; then
    sbatch --export="cv=${cv}" -J "SCV${cv}"  --nodes=1 --ntasks=1 --cpus-per-task=11 --time=7-00:00:00 run_S.sh
  fi
  if ! test -f "para_cv${cv}.rda"; then
    sbatch --export="cv=${cv}" -J "paraCV${cv}"  --nodes=1 --ntasks=1 --cpus-per-task=11 --time=7-00:00:00 run_para.sh
  fi
  if ! test -f "flex_cv${cv}.rda"; then
    sbatch --export="cv=${cv}" -J "flexCV${cv}"  --nodes=1 --ntasks=1 --cpus-per-task=1 --time=7-00:00:00 run_flex.sh
  fi
done