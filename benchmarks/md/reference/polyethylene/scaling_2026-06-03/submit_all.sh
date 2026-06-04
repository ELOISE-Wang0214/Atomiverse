#!/bin/bash
ref=/home/jicli594/work/molcrafts/Atomiverse/benchmarks/md/reference/polyethylene
base=$ref/scaling_2026-06-03
for N in 64; do
  d=$base/lammps_${N}r
  mkdir -p "$d"
  cp "$ref/lammps.in" "$ref/msu_pe.data" "$d"/
  sbatch --chdir="$d" --ntasks=$N --cpus-per-task=1 --job-name=lmp_${N}r "$base/job_lammps.sbatch"

  d=$base/atomiverse_${N}r
  mkdir -p "$d"
  sbatch --chdir="$d" --ntasks=$N --cpus-per-task=1 --job-name=av_${N}r "$base/job_atomiverse.sbatch"
done
