# dev_test: reproduction and regression environment

Not part of the package (excluded via `.Rbuildignore`). Data, environments,
the baseline snapshot and run outputs are git-ignored.

## Setup

```bash
cd dev_test
pixi install --all                    # envs: v4 (R 4.3 + Seurat 4.4) and v5 (R 4.4 + Seurat 5.5)
pixi run -e v4 install-extra          # Bioconductor data packages, then IntNMF, add2ggplot (CRAN),
pixi run -e v5 install-extra          # quickSparseM and scSeqComm (GitLab, ~400 MB + 150 MB downloads)
mkdir -p baseline && git -C .. archive <commit> R NAMESPACE DESCRIPTION data | tar -x -C baseline
```

Reference outputs published by the authors (ST-colon3, Baidu share in
`tutorial/2a.multisample_ecotype.md`) go to `data/baidu/ST-colon3/`:

```bash
BaiduPCS-Go transfer https://pan.baidu.com/s/1l7g5skCZkaz6MIpzvH8F4w 7bsc
BaiduPCS-Go download /<saved dir>/CRC/output/ST-colon3 --saveto data/baidu
```

## Scripts (run from `dev_test/`)

| command | what it checks |
|---|---|
| `pixi run -e v5 probe-baseline` | every Seurat-object access of the original code, one by one |
| `pixi run -e v5 probe` | the same for the current package |
| `pixi run -e v5 test` | package unit tests (new code vs original dense implementation) |
| `python3 scripts/monitor.py runs/X/monitor.tsv -- pixi run -e v5 Rscript scripts/run_tutorial_1a.R pkg runs/X [n_threads=N]` | tutorial 1a end to end, with thread count / busy cores sampled every second |
| `... scripts/run_tutorial_1b.R pkg runs/X [n_threads=N num_core=M]` | tutorial 1b (`spaniche_nmf_comm`, needs scSeqComm) |
| `pixi run -e v5 Rscript scripts/smoke_downstream.R runs/X` | plotting / downstream functions with nothing but SpaNiche loaded |
| `bash scripts/bench.sh [prefix]` | benchmark: tutorial 1a, celltype_and_gene, tutorial 1b under v5 and v4, one at a time |
| `python3 scripts/stage_times.py runs/X ...` | per-stage wall time and busy cores of finished runs |
| `pixi run -e v5 Rscript scripts/compare_runs.R runs/X <runs/Y or author>` | numerical comparison of two runs, or against the authors' published output |

`code_dir` is `baseline` (original code, sourced) or `pkg` (`devtools::load_all("..")`).
Set `SPANICHE_PROFILE=1` to print a per-stage profile.
