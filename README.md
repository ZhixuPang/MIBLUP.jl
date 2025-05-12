# MIBLUP.jl
![Fig 1 MIBLUP](https://github.com/ZhixuPang/MIBLUP.jl/assets/137134445/6b5e056c-1209-4b4b-920d-c2c19c7979f9)

## Usage
```
include("./MIBLUP.jl")
```
### 1.Read data
```
bfile = "data" # bfile path
data = MIBLUP.read_plink(bfile) #to load PLINK binary files into memory
geno = data.bed
y = Vector(data.fam[:, "Trait_1"]) # phenotype data
covar = readdlm("covar.txt", Float64) # covariates data
```

### 2.Usage
```
res = MIBLUP.runMIBLUP(y, geno; covariance = covar,
                       ids = data.fam[:, "Individual_ID"],
                       snp_names = data.bim[:, "SNP_ID"])
```

以下是为您的 `MIBLUP` 函数编写的 GitHub `README.md` 文件草稿，适合用于项目主页展示其功能、安装和使用方法：

---

# MIBLUP.jl

**MIBLUP.jl** implements Mutual Information-based Best Linear Unbiased Prediction (MIBLUP), a genomic prediction method that integrates mutual information theory with BLUP to improve prediction accuracy, especially when feature selection or marker weighting is beneficial.

## Features

* Mutual information-based marker selection (`mRMR`)
* Kinship matrix calculation with SNP weighting
* Variance component estimation via HE, AI-REML, EMMA, or hybrid methods
* Flexible phenotype support: continuous or discrete
* BLUP using the MME solver
* Cross-validation and multiple accuracy metrics: correlation, AUC, RMSE

## Installation

You can install this package from source via:

```julia
using Pkg
include("./src/MIBLUP.jl")
```

Replace `YOUR_USERNAME` with your GitHub handle.

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
