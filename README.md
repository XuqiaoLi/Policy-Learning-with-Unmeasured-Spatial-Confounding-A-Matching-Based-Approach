# The R code for Policy Learning with Unmeasured Spatial Confounding: A Matching-Based Approach

[`simulation.R`](simulation.R) implements the simulation studies in Section 3 of our work *Policy Learning with Unmeasured Spatial Confounding: A Matching-Based Approach*.

The script compares the proposed spatial matching-based methods with nonspatial matching-based methods, Policy Tree, and Q-learning in four scenarios: Tree-X, Non-tree-X, Tree-XU, and Non-tree-XU.

Each scenario is replicated 100 times, with 800 training observations and 2,000 test observations per replication. The script generates the data and computes empirical policy values, accuracy, and regret. It also produces boxplots of empirical policy values.

## Main functions

| Function | Description |
| --- | --- |
| `generate_dataset` | Generates data with unmeasured spatial confounding. |
| `generate_train_test_data` | Constructs training and test samples. |
| `find_optimal_weight` | Selects the spatial weighting parameter based on ASDM. |
| `MBlearning` | Implements the proposed spatial matching-based methods, with or without bias correction. |
| `MBlearning_non` | Implements the nonspatial matching-based methods using observed covariates. |
| `evaluate_policy_sim` | Computes empirical policy value, accuracy relative to the oracle policy, and regret. |
| `compareAllmethod_sim` | Fits and evaluates all methods in one simulation replication. |
| `run_simulation_reps` | Repeats the simulation and collects results across replications. |

## Run

Install the R packages listed at the beginning of the script, then run from this directory:

```bash
Rscript simulation.R
```
