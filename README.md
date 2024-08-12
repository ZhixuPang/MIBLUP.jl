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

