#!/bin/bash
# Benchmark the current package on ST-colon3, one run at a time (run on an idle machine).
# Usage: bash scripts/bench.sh [prefix]   -> runs/<prefix>_<pipeline>_<env>/
set -u
prefix=${1:-bench}
specs=(
  "1a v5 run_tutorial_1a.R"
  "1a v4 run_tutorial_1a.R"
  "gene v5 run_celltype_gene.R"
  "gene v4 run_celltype_gene.R"
  "1b v5 run_tutorial_1b.R"
  "1b v4 run_tutorial_1b.R"
)
for spec in "${specs[@]}"; do
  read -r name env script <<< "$spec"
  run="runs/${prefix}_${name}_${env}"
  rm -rf "$run"; mkdir -p "$run"
  python3 scripts/monitor.py "$run/monitor.tsv" -- \
    pixi run -e "$env" Rscript "scripts/$script" pkg "$run" > "$run/log.txt" 2>&1
  echo "$run: $(grep -h -E '^TOTAL|^MONITOR' "$run/log.txt" | tr '\n' ' ')"
done
