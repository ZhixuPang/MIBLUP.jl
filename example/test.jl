# 引入所需的包
using Pkg
# using ArgParse
using CSV
using DataFrames
using DelimitedFiles
using Base.Threads

include("../src/MIBLUP.jl")
bfile = "data"
println("Start read plink")
data = MIBLUP.read_plink(bfile)
#covar = readdlm("covar.txt", Float64)
# 进行基因型插补
#geno = MIBLUP.geno_imputation(data.bed)
geno = data.bed
# 获取y向量
pheno = CSV.read("pheno.txt", DataFrame)
y = Vector(pheno[:, 1])

res = MIBLUP.runMIBLUP(y, geno;
                       ids = data.fam[:, "Individual_ID"],
                       snp_names = data.bim[:, "SNP_ID"],
                       n_pre_sels = 20)

println("MIBLUP END.")
