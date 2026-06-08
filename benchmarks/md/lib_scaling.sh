# Shared helpers for bench_scaling.sbatch. Sourced, not executed.
# Keep ALL parsing/measurement logic here so every run uses identical code.

cbrt(){ awk -v n="$1" 'BEGIN{printf "%d", int(n^(1/3)+0.5)}'; }

# peak total VRAM [GB] over the run of "$@" — uniform metric for every cell
# (sum of nvidia-smi memory.used across visible GPUs, sampled at 0.3 s).
pv(){ local f; f=$(mktemp); ( while :; do nvidia-smi --query-gpu=memory.used \
      --format=csv,noheader,nounits 2>/dev/null|awk '{s+=$1}END{print s}'>>"$f"; \
      sleep 0.3; done ) & local p=$!; "$@" >/dev/null 2>&1; kill "$p" 2>/dev/null; \
      awk 'BEGIN{m=0}{if($1>m)m=$1}END{printf "%.2f",m/1024}' "$f"; rm -f "$f"; }

# Parse the PRODUCTION "Loop time" from a LAMMPS log (the LAST one — the deck
# does an untimed warmup `run` first, then the timed `run`). Echoes
# "matom_step_s atoms" or nothing on failure.
lmp_throughput(){ local log="$1"; local ln; ln=$(grep "Loop time of" "$log"|tail -1); \
  [ -z "$ln" ] && return 1; \
  awk -v ln="$ln" 'BEGIN{n=split(ln,a," "); t=a[4]; s=a[9]; at=a[12]; \
    if(t>0) printf "%.1f %d", at*s/t/1e6, at}'; }
