# Stochastic Knapsack Problem (SKP) – Code-to-Paper Reproducibility Guide

This repository contains Java and Python code to reproduce the computational elements described in the work:

R. Rossi, S. Prestwich, S. A. Tarim, "Mixed-Integer Linear Programming Approximations for the Stochastic Knapsack Problem," working paper, 2025

This guide maps paper components to code, shows how to build/run experiments, and points to curated CSV outputs under results/.

## What's Implemented
- Normal SKP (static): Piecewise-linear MILP, Lazy-Cuts warm-start, SAA/SAA-LD, SDP baseline
	- Batch driver: [src/main/java/skp/batch/SKPNormalBatch.java](src/main/java/skp/batch/SKPNormalBatch.java)
	- MILP PWLA (piecewise first-order loss): outputs solved_normal_instances_MILP_{partitions}.csv
	- MILP Lazy Cuts (LC) and LC warm-start: outputs solved_normal_instances_LC.csv and solved_normal_instances_LC_warm_start_25.csv
	- SAA and SAA_LD: outputs solved_normal_instances_SAA.csv and solved_normal_instances_SAA_LD.csv
	- Dynamic SKP (SDP baseline): outputs solved_normal_instances_DSKP.csv
	- Curated results: [results/normal_results](results/normal_results)

- Multivariate Normal SKP (correlated): PWLA MILP, Lazy-Cuts, SAA/SAA-LD, SDP (special covariance structure)
	- Batch driver: [src/main/java/skp/batch/SKPMultinormalBatch.java](src/main/java/skp/batch/SKPMultinormalBatch.java)
	- MILP PWLA: outputs solved_multinormal_instances_MILP_{partitions}.csv
	- MILP Lazy Cuts (LC): outputs solved_multinormal_instances_LC.csv (with optional warm start in the same batch driver)
	- SAA and SAA_LD: outputs solved_multinormal_instances_SAA.csv and solved_multinormal_instances_SAA_LD.csv
	- Dynamic SKP (SDP) when covariance has ρ^{|i−j|}σ_iσ_j structure: outputs solved_multinormal_instances_DSKP.csv
	- Curated results: [results/mvn_large_R_results](results/mvn_large_R_results), [results/mvn_small_R_results](results/mvn_small_R_results)

- Generic distributions (Lognormal, Gamma): MILP with Lazy Cuts (exact and normal-approximation), SAA/SAA-LD, SDP baseline
	- Batch driver: [src/main/java/skp/batch/SKPGenericDistributionBatch.java](src/main/java/skp/batch/SKPGenericDistributionBatch.java)
	- Methods: LAZY_CUTS, LAZY_CUTS_NORMAL_APPROXIMATION, SAA, SAA_LD, and DSKP (for small N)
	- Curated results: [results/lognormal_results](results/lognormal_results), [results/gamma_results](results/gamma_results)

- Merzifonluoğlu et al. (2012) exact B&B and heuristic baselines (normal sizes)
	- Heuristic (two-stage rounding): [src/main/python/binary_heuristic.py](src/main/python/binary_heuristic.py)
	- Branch-and-Bound (exact) + heuristic batcher: [src/main/python/skp_batch.py](src/main/python/skp_batch.py)
	- Produces: solved_normal_instances_Merzifonluoglu.csv, solved_normal_instances_Merzifonluoglu_heuristic.csv
	- Curated results: [results/normal_results](results/normal_results)

## Core Models and Components
- MILP PWLA (Piecewise First-Order Loss) formulations
	- Normal: [src/main/java/skp/milp/SKPNormalMILP.java](src/main/java/skp/milp/SKPNormalMILP.java)
	- Multivariate Normal: [src/main/java/skp/milp/SKPMultinormalMILP.java](src/main/java/skp/milp/SKPMultinormalMILP.java)
	- First-order loss utilities: [src/main/java/skp/folf](src/main/java/skp/folf)
	- Optional OPL models and data: [src/main/resources/opl_models](src/main/resources/opl_models); Java batchers also package .dat archives for OPL.

- MILP with Lazy Cuts (LC)
	- Generic distributions: [src/main/java/skp/milp/SKPGenericDistributionLazyCuts.java](src/main/java/skp/milp/SKPGenericDistributionLazyCuts.java)
	- Multivariate Normal: [src/main/java/skp/milp/SKPMultinormalLazyCuts.java](src/main/java/skp/milp/SKPMultinormalLazyCuts.java)

- Sample Average Approximation (SAA) and SAA_LD (Kleywegt et al. bound-guided)
	- Generic distributions: [src/main/java/skp/saa/SKPGenericDistributionSAA.java](src/main/java/skp/saa/SKPGenericDistributionSAA.java), [src/main/java/skp/saa/SKPGenericDistributionSAA_LD.java](src/main/java/skp/saa/SKPGenericDistributionSAA_LD.java)
	- Multivariate Normal: [src/main/java/skp/saa/SKPMultinormalSAA.java](src/main/java/skp/saa/SKPMultinormalSAA.java), [src/main/java/skp/saa/SKPMultinormalSAA_LD.java](src/main/java/skp/saa/SKPMultinormalSAA_LD.java)

- Dynamic SKP (SDP) baselines
	- Normal: [src/main/java/skp/sdp/DSKPNormal.java](src/main/java/skp/sdp/DSKPNormal.java)
	- Multivariate Normal (special structure): [src/main/java/skp/sdp/DSKPMultinormal.java](src/main/java/skp/sdp/DSKPMultinormal.java)
	- Generic distributions: [src/main/java/skp/sdp/DSKPGenericDistribution.java](src/main/java/skp/sdp/DSKPGenericDistribution.java)

- Instance generators (Pisinger-style classes, sizes/correlations)
	- Normal: generated in [SKPNormalBatch](src/main/java/skp/batch/SKPNormalBatch.java)
	- Multivariate Normal (ρ controls correlation): [SKPMultinormalBatch](src/main/java/skp/batch/SKPMultinormalBatch.java)
	- Lognormal/Gamma: [SKPGenericDistributionBatch](src/main/java/skp/batch/SKPGenericDistributionBatch.java)

## Prerequisites
- macOS or Linux; Java (JDK 11+), Maven, and Python 3.10+
- IBM ILOG CPLEX Optimization Studio (OPL/Concert) with a valid license. The Java MILP code uses the OPL/Concert APIs (`ilog.opl.*`, `ilog.concert.*`, `IloCplex`) to build/solve models and callbacks. On macOS you may need `DYLD_LIBRARY_PATH` and `-Djava.library.path` as shown in batch headers.
- Python baselines (B&B and heuristic) do not require CPLEX.

## Build (Java)
Use Maven to build the project and resolve dependencies:

```bash
mvn -q -DskipTests package
```

## Running Experiments (Java)
You can run batch drivers from your IDE (recommended) or from CLI with a full runtime classpath. Each batch writes instance JSON and solution CSVs under `batch/<type>/<size>/<cv>[/<rho>]` and prints progress to stdout.

- Normal SKP (all methods; PWLA, LC, LC warm-start, SAA_LD, DSKP for N=25)
	- Main class: skp.batch.SKPNormalBatch
	- Output examples placed next to their inputs; curated copies under [results/normal_results](results/normal_results)

- Multivariate Normal SKP (PWLA, LC, SAA_LD, DSKP for N=25 when structure holds)
	- Main class: skp.batch.SKPMultinormalBatch
	- Curated copies under [results/mvn_large_R_results](results/mvn_large_R_results) and [results/mvn_small_R_results](results/mvn_small_R_results)

- Lognormal/Gamma (Generic) SKP (LC/LC-normal-approx, SAA, SAA_LD, DSKP for N=25)
	- Main class: skp.batch.SKPGenericDistributionBatch
	- To switch distribution (GAMMA vs LOGNORMAL), edit Dist in main()
	- Curated copies under [results/lognormal_results](results/lognormal_results) and [results/gamma_results](results/gamma_results)

## Notes on CLI execution
- If running from CLI, you need the runtime classpath (target/classes plus dependencies). A simple route is to run from your IDE. If you prefer CLI, you can print the runtime classpath using your Maven setup and pass it to java -cp.

## Baselines (Python) – Merzifonluoğlu et al. (2012)
- Install requirements:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy
```

- Run the heuristic/B&B over generated normal instances:
	1) Generate normal instances via [SKPNormalBatch](src/main/java/skp/batch/SKPNormalBatch.java) (it creates normal_instances.json under batch/...)
	2) From repo root or src/main/python, run one of:

```bash
python src/main/python/skp_batch.py               # solve normal_instances.json in CWD
python src/main/python/skp_batch.py               # or call recursive_solve() inside the script to traverse a root directory
```

- Outputs
	- solved_normal_instances_Merzifonluoglu.csv (B&B exact + heuristic, N≲50)
	- solved_normal_instances_Merzifonluoglu_heuristic.csv (heuristic only, for larger N)
	- Curated examples are in [results/normal_results](results/normal_results)

## OPL Models
- OPL .mod files: [src/main/resources/opl_models](src/main/resources/opl_models)
- Batch drivers export zipped .dat sets alongside instance JSONs (see `storeBatchAsOPLDataFiles` in `SKPNormalBatch`/`SKPMultinormalBatch`) for reproducing MILP PWLA in IBM ILOG OPL. This is aligned with the Java runs, which already use the OPL/Concert Java APIs.

## Reproducing Tables/Figures from the Paper
- Normal SKP tables/figures: reproduced by SKPNormalBatch outputs (MILP PWLA, LC/LC warm-start, SAA/SAA_LD, DSKP). See [results/normal_results](results/normal_results) for curated CSVs.
- Multivariate Normal SKP: reproduced by SKPMultinormalBatch outputs (PWLA, LC, SAA_LD, DSKP where applicable). See [results/mvn_large_R_results](results/mvn_large_R_results) and [results/mvn_small_R_results](results/mvn_small_R_results).
- Lognormal/Gamma SKP: reproduced by SKPGenericDistributionBatch (LC vs LC-normal-approx, SAA/SAA_LD, DSKP). See [results/lognormal_results](results/lognormal_results) and [results/gamma_results](results/gamma_results).
- Merzifonluoğlu baselines used for comparison in the normal case: use the Python drivers; curated outputs are in [results/normal_results](results/normal_results).

## Instance Classes and Settings
- Instance types follow Pisinger’s families (e.g., UNCORRELATED, WEAKLY/STRONGLY CORRELATED, SUBSET SUM, etc.) parameterized by item counts (e.g., 25/50/100/200/500), coefficient of variation (cv ∈ {0.1, 0.2}), and for MVN, correlation ρ ∈ {0.75, 0.95}.
- Random seed is fixed in [src/main/java/skp/batch/SKPBatch.java](src/main/java/skp/batch/SKPBatch.java) for reproducibility.

## Where outputs are written
- Batch runs create per-experiment folders under `batch/<instance_type>/<N>/<cv>[/<rho>]`, with:
	- normal_instances.json or multinormal_instances.json (inputs)
	- solved_*.json and solved_*.csv (solutions)
- Curated CSVs used in the paper are mirrored under [results/](results).

## Troubleshooting
- For OPL-only reproduction, configure DYLD_LIBRARY_PATH and -Djava.library.path as noted in batch headers (macOS).
- For large N, prefer the Python heuristic exporter; the exact B&B is intended for moderate sizes.
