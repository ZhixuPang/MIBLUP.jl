# MIBLUP.jl
![Fig 1 MIBLUP](https://github.com/ZhixuPang/MIBLUP.jl/assets/137134445/6b5e056c-1209-4b4b-920d-c2c19c7979f9)

**MIBLUP.jl** implements Mutual Information-based Best Linear Unbiased Prediction (MIBLUP), a genomic prediction method that integrates mutual information theory with BLUP to improve prediction accuracy, especially when feature selection or marker weighting is beneficial.

## Features

* Mutual information-based marker selection (`mRMR`)
* Kinship matrix calculation with SNP weighting
* Variance component estimation via HE, AI-REML, EMMA, or hybrid methods
* Flexible phenotype support: continuous or discrete
* BLUP using the MME solver
* Cross-validation and multiple accuracy metrics: correlation, AUC, RMSE

## Installation

Here is the revised **Installation** section for your `README.md`:

---

## Installation

You can install this package from source as follows:

### Step 1: Download the source code

```bash
wget https://github.com/ZhixuPang/MIBLUP.jl/archive/refs/heads/master.zip
unzip master.zip
```

### Step 2: Launch Julia and include the package

```julia
include("MIBLUP.jl-master/src/MIBLUP.jl")
```


## Quick Start

```julia
using MIBLUP

res = runMIBLUP(
    phe,               # Phenotype vector
    geno;              # Genotype matrix (individuals × SNPs)
    covariance = cov,  # Optional covariate matrix
    kinship = nothing,
    ids = sample_ids,
    snp_names = snp_ids,
    cov_names = cov_names,
    vc_method = "he+ai",
    phe_type = "continuous",
    mi_target = "phe",
    nbins = 10,
    acc_type = "cor",
    kinship_weight = true,
    marker_selection = true,
    n_cv = 5,
    n_sels = 100,
    n_pre_sels = 1000,
    verbose = true
)
```

## Arguments

| Argument           | Type     | Description                                                 |
| ------------------ | -------- | ----------------------------------------------------------- |
| `phe`              | `Vector` | Phenotypic values (with `missing` allowed)                  |
| `geno`             | `Matrix` | Genotype matrix (individuals × SNPs)                        |
| `covariance`       | `Matrix` | Covariates (default: intercept-only)                        |
| `kinship`          | `Matrix` | Precomputed kinship matrix                                  |
| `vc_method`        | `String` | Variance estimation: `"he"`, `"emma"`, `"ai"`, or `"he+ai"` |
| `mi_target`        | `String` | `"phe"` or `"gebv"` for mutual information target           |
| `phe_type`         | `String` | `"continuous"` or `"discrete"` phenotype type               |
| `nbins`            | `Int`    | Number of bins for discretization                           |
| `n_cv`             | `Int`    | Number of CV folds                                          |
| `acc_type`         | `String` | `"cor"`, `"auc"`, or `"rmse"`                               |
| `kinship_weight`   | `Bool`   | Whether to use SNP-weighted kinship                         |
| `marker_selection` | `Bool`   | Whether to perform mRMR marker selection                    |
| `n_sels`           | `Int`    | Number of markers selected                                  |
| `n_pre_sels`       | `Int`    | Number of markers pre-filtered                              |

## Outputs

`runMIBLUP` returns a dictionary with:

* `res`: `DataFrame` with individual `ID`, `GEBV`, and `yHat`
* `cov_effects`: Covariate effect estimates
* `mutual_information`: MI values for SNPs
* `selected_marker_effects`: Effects for selected SNPs
* `vc`: Estimated variance components
* `kinship`: Final kinship matrix used


## Usage

### 1. Load Phenotype Data

The phenotype file should be in `.csv` format with at least two columns: `ID` and `Phenotype`, for example:

#### `phenotype.csv`

```csv
ID,Phenotype
Ind1,12.5
Ind2,15.3
Ind3,14.8
Ind4,NA
Ind5,13.7
```

In your Julia script, use the following code to load the phenotype values:

```julia
using CSV, DataFrames

y_df = CSV.read("phenotype.csv", DataFrame;missingstring = "NA")
y = y_df.Phenotype
```

> **Note:** The order of IDs in `phenotype.csv` should match the order in the PLINK `.fam` file used in `MIBLUP.read_plink`. The `ID` column corresponds to `Individual_ID` in `data.fam`.



### 2. Load PLINK Genotype Data

```julia
bfile = "data"  # Base name of the PLINK files (i.e., data.bed, data.bim, data.fam)
data = MIBLUP.read_plink(bfile)
geno = data.geno
```

---

## Run MIBLUP

```julia
res = MIBLUP.runMIBLUP(y, geno;
    ids = data.fam[:, "Individual_ID"],
    snp_names = data.bim[:, "SNP_ID"],
    n_pre_sels = 20
)
```

---

## Run GBLUP (as a baseline)

```julia
res = MIBLUP.runMIBLUP(yNa, geno;
    ids = data.fam[:, "Individual_ID"],
    snp_names = data.bim[:, "SNP_ID"],
    kinship_weight = false,
    marker_selection = false
)
```

---

## Output Explanation

### Estimated Genetic Parameters

```julia
vc = DataFrame(vg = res.vc.vg, ve = res.vc.ve)
```

### Estimated GEBVs and Predicted Phenotypes

```julia
res.res.GEBV  # Genomic Estimated Breeding Values
res.res.yHat  # Predicted phenotype values
```

### Estimated Mutual Information

```julia
res.mutual_information
```

### Effects of Selected SNPs

```julia
res.selected_marker_effects
```

---

## 🧠 Method Description

**MIBLUP (Mutual Information-based Best Linear Unbiased Prediction)** is a two-stage genomic prediction method that combines mutual information-based SNP selection with a weighted GBLUP model. It improves prediction accuracy by identifying and emphasizing informative SNPs prior to model fitting.

### Model Formulation

The MIBLUP model is defined as:

```
y = Xb + Mₛf + Zu + e                      (1)
```

* `y`: vector of phenotypes
* `X`: design matrix for fixed effects `b`
* `Mₛ`: matrix of selected SNP markers
* `f`: fixed effects of selected SNPs
* `Z`: design matrix for random effect `u`
* `u ~ N(0, G_w * σ²_u)`: random genetic effects modeled using the mutual information-weighted GRM `G_w`
* `e ~ N(0, I * σ²_e)`: residual error

The genomic estimated breeding value (GEBV) is computed as:

```
g = Mₛf + u = Mₛf + G_w Z' V⁻¹ (y - Mₛf - Xb)       (2)
```

where `V = ZG_w Z'σ²_g + Iσ²_e`.

### SNP Selection Procedure

MIBLUP uses a two-stage marker selection strategy:

1. **Pre-selection with mRMR (Minimum Redundancy Maximum Relevance):**
   The top `n_pre_sels` SNPs are selected by ranking their mutual information with the phenotype, while minimizing redundancy among selected markers. This step ensures the initial subset of SNPs is both relevant and non-redundant.

2. **Refinement via Incremental Feature Selection and Cross-validation:**
   The pre-selected SNPs are added sequentially to the model. For each SNP, the improvement in prediction accuracy is evaluated using N-fold cross-validation. A SNP is retained if it consistently improves accuracy (p < 0.05) across folds. The process terminates if a preset number `k` of consecutive SNPs fail to enhance performance.

### Mutual Information-weighted GRM

The genomic relationship matrix `G_w` is constructed as:

```
G_w = W D W' / m                         (3)
```

* `W`: standardized genotype matrix
* `D`: diagonal matrix of SNP weights `d_k`

The diagonal elements `d_k` are defined as:

```
d_k = w_k / (2p_k(1 - p_k))             (4)
```

where:

* `w_k = [I(m_k, y) × m] / Σ I(m_k, y)`
* `I(m_k, y)`: mutual information between marker `k` and phenotype `y`
* `p_k`: allele frequency at SNP `k`
* `w_k`: reflects the informativeness of SNP `k`

