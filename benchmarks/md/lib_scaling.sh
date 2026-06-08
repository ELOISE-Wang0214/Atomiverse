# Shared helpers for bench_scaling.sbatch. Sourced, not executed.
# Keep ALL parsing/measurement logic here so every run uses identical code.

cbrt(){ awk -v n="$1" 'BEGIN{printf "%d", int(n^(1/3)+0.5)}'; }

# peak total VRAM [GB] over the run of "$@" — uniform metric for every cell
# (sum of nvidia-smi memory.used across visible GPUs, sampled at 0.3 s).
pv(){ local f; f=$(mktemp); ( while :; do nvidia-smi --query-gpu=memory.used \
      --format=csv,noheader,nounits 2>/dev/null|awk '{s+=$1}END{print s}'>>"$f"; \
      sleep 0.3; done ) & local p=$!; "$@" >/dev/null 2>&1; kill "$p" 2>/dev/null; \
      awk 'BEGIN{m=0}{if($1>m)m=$1}END{printf "%.2f",m/1024}' "$f"; rm -f "$f"; }

# peak job RSS [GB] over the run of "$@" — node used-memory delta (the CPU
# analog of pv; run on an --exclusive node so the delta is this job).
# DEPRECATED for AV-vs-LAMMPS comparison: the node-delta under-counts LAMMPS to
# ~0 (it reuses/pre-touches differently), so AV's true RSS got compared against
# LAMMPS's spurious near-zero → a bogus ~12x. Use rss_hwm_gb instead.
pv_cpu(){ local f base; f=$(mktemp); base=$(free -m|awk '/^Mem:/{print $3}'); \
  ( while :; do free -m|awk '/^Mem:/{print $3}'>>"$f"; sleep 0.5; done ) & local p=$!; \
  "$@" >/dev/null 2>&1; kill "$p" 2>/dev/null; \
  awk -v b="$base" 'BEGIN{m=0}{if($1>m)m=$1}END{printf "%.2f",(m-b)/1024}' "$f"; rm -f "$f"; }

# Sum of per-rank peak VmHWM [GB] for an srun job whose ranks are EACH wrapped in
# `/usr/bin/time -v` (so every rank prints "Maximum resident set size (kbytes)").
# This is the SAME metric AV's `scaling` binary self-reports as peak_rss_total_kb
# (MPI_Reduce SUM of /proc/self/status VmHWM), so it makes AV vs LAMMPS an
# apples-to-apples total-resident-footprint comparison. Usage:
#   rss=$(rss_hwm_gb srun --mpi=pmi2 -n $p -c 1 /usr/bin/time -v $BIN ...args)
# (throughput is parsed separately from the LAMMPS log / AV out-dir, so routing
# the binary's stdout to /dev/null here is fine.)
rss_hwm_gb(){ local ef; ef=$(mktemp); "$@" >/dev/null 2>"$ef"; \
  awk -F: '/Maximum resident set size/{v=$NF; gsub(/[^0-9]/,"",v); s+=v} \
           END{printf "%.2f", s/1048576}' "$ef"; rm -f "$ef"; }

# Parse the PRODUCTION "Loop time" from a LAMMPS log (the LAST one — the deck
# does an untimed warmup `run` first, then the timed `run`). Echoes
# "matom_step_s atoms" or nothing on failure.
lmp_throughput(){ local log="$1"; local ln; ln=$(grep "^Loop time of" "$log"|tail -1); \
  [ -z "$ln" ] && return 1; \
  awk -v ln="$ln" 'BEGIN{n=split(ln,a," "); t=a[4]; s=a[9]; at=a[12]; \
    if(t>0) printf "%.1f %d", at*s/t/1e6, at}'; }
