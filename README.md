# Stochastic Knapsack Problem (SKP) – Code-to-Paper Reproducibility Guide

This repository contains Java and Python code to reproduce the computational elements described in the work:

R. Rossi, S. Prestwich, S. A. Tarim, "Mixed-Integer Linear Programming Approximations for the Stochastic Knapsack Problem," [arXiv:2512.14912](https://arxiv.org/abs/2512.14912), 2025

## What's Implemented

### Static

- Normal: Piecewise-linear MILP, Lazy-Cuts warm-start, SAA/SAA-LD, SDP baseline, Merzifonluoğlu et al. (2012) exact B&B and heuristic baselines,
Range et al. (2018) heuristic baseline,
	- Batch driver: [src/main/java/skp/batch/SKPNormalBatch.java](src/main/java/skp/batch/SKPNormalBatch.java)
	- Batch driver: [src/main/python/sskp_merzifonluoglu/skp_batch.py](src/main/python/sskp_merzifonluoglu/skp_batch.py)
	- Batch driver: [src/main/python/sskp_range/skp_batch.py](src/main/python/sskp_range/skp_batch.py)
	- Curated results: [results/normal_results](results/normal_results)

- Multivariate Normal (correlated): PWLA MILP, Lazy-Cuts, SAA/SAA-LD, SDP (special covariance structure)
	- Batch driver: [src/main/java/skp/batch/SKPMultinormalBatch.java](src/main/java/skp/batch/SKPMultinormalBatch.java)
	- Curated results: [results/mvn_large_R_results](results/mvn_large_R_results), [results/mvn_small_R_results](results/mvn_small_R_results)

- Generic distributions (Lognormal, Gamma): MILP with Lazy Cuts (exact and normal-approximation), SAA/SAA-LD, SDP baseline
	- Batch driver: [src/main/java/skp/batch/SKPGenericDistributionBatch.java](src/main/java/skp/batch/SKPGenericDistributionBatch.java)
	- Curated results: [results/lognormal_results](results/lognormal_results), [results/gamma_results](results/gamma_results)

### Dynamic

- Normal
	- Batch driver: [src/main/java/skp/batch/SKPNormalRecedingBatch.java](src/main/java/skp/batch/SKPNormalRecedingBatch.java)

- Multivariate Normal (correlated with special covariance structure)
	- Batch driver: [src/main/java/skp/batch/SKPMultinormalRecedingBatch.java](src/main/java/skp/batch/SKPMultinormalRecedingBatch.java)

- Generic distributions (Lognormal, Gamma)
	- Batch driver: [src/main/java/skp/batch/SKPGenericDistributionRecedingBatch.java](src/main/java/skp/batch/SKPGenericDistributionRecedingBatch.java)

## Instances

- Instance generators (Pisinger-style classes, sizes/correlations)
	- Normal: generated in [SKPNormalBatch](src/main/java/skp/batch/SKPNormalBatch.java)
	- Multivariate Normal (ρ controls correlation): [SKPMultinormalBatch](src/main/java/skp/batch/SKPMultinormalBatch.java)
	- Lognormal/Gamma: [SKPGenericDistributionBatch](src/main/java/skp/batch/SKPGenericDistributionBatch.java)

## Prerequisites
- Java (JDK 11+), Maven, and Python 3.10+
- IBM ILOG CPLEX Optimization Studio (OPL/Concert) with a valid license. The Java MILP code uses the OPL/Concert APIs (`ilog.opl.*`, `ilog.concert.*`, `IloCplex`) to build/solve models and callbacks. 

## Running Experiments (Java)
You can run batch drivers from your IDE (recommended) or from CLI with a full runtime classpath. Each batch writes instance JSON and solution CSVs under `batch/<type>/<size>/<cv>[/<rho>]` and prints progress to stdout.

## Baselines (Python)
- Run the heuristic/B&B over generated normal instances:
	1) Generate normal instances via [SKPNormalBatch](src/main/java/skp/batch/SKPNormalBatch.java) (it creates normal_instances.json under batch/...)
	2) From repo root or src/main/python/..., run one of:

```bash
python src/main/python/.../skp_batch.py    # solve normal_instances.json in the current directory
							               # or call recursive_solve() inside the script to traverse a root directory
```

## Where outputs are written
- Batch runs create per-experiment folders under `batch/<instance_type>/<N>/<cv>[/<rho>]`, with:
	- normal_instances.json or multinormal_instances.json (inputs)
	- solved_*.json and solved_*.csv (solutions)
- Curated CSVs used in the paper are mirrored under [results/](results).


