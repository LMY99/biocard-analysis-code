#!/usr/bin/env bash

for cv in {1..5}; do
  if ! test -f "S_CV${cv}.rda"; then
    sbatch --export="cv=${cv}" -J "SCV{cv}" --time=7-00:00:00 run_S.sh
  fi
done