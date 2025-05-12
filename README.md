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
    marker_selection = false,
    n_sels = 20,
    n_pre_sels = round(Int, 0.1 * size(data.bim, 1))
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

