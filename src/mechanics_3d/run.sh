#!/usr/bin/env bash
# Run one simulation in its own folder (the solver writes its output to the working directory).
# Usage: ./run.sh input.xml [final_time_seconds=1800] [output_dir=runs/<input>_<T>s]
# Set BEND_MODEL to use an executable other than build/bend-model.
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
xml=$(realpath -- "${1:?usage: ./run.sh input.xml [final_time_seconds] [output_dir]}")
secs=${2:-1800}
out=${3:-$here/runs/$(basename -- "$xml" .xml)_${secs}s}
bin=$(realpath -- "${BEND_MODEL:-$here/build/bend-model}")
[[ -x $bin ]] || { echo "No executable at $bin (build first, or set BEND_MODEL)" >&2; exit 2; }
mkdir -p -- "$out"
[[ -z $(ls -A -- "$out") ]] || { echo "Output folder is not empty: $out" >&2; exit 3; }
cp -- "$xml" "$out/input.xml"
cd -- "$out"
echo "Running $secs simulated seconds (dt=0.001) -> $out"
status=0
MICHAEL_FINAL_TIME_SECONDS=$secs "$bin" -dt=0.001 input.xml > execution.log 2>&1 || status=$?
echo "$status" > exit_code.txt
echo "Exit code $status. Log: $out/execution.log"
exit "$status"
