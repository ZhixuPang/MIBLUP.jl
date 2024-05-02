# 引入所需的包
using Pkg
# using ArgParse
using CSV
using DataFrames
using DelimitedFiles
using Base.Threads

# Get the current number of threads
current_threads = Threads.nthreads()
println("Current number of threads: $current_threads")

# Set the desired number of threads
desired_threads = 16
Threads.nthreads() = desired_threads
println("Number of threads set to $desired_threads")

# Check the updated number of threads
updated_threads = Threads.nthreads()
println("Updated number of threads: $updated_threads")
println("Start activate MIBLUP")
include("/share/nas2/data/pangzx/Software/julia_dev/MIBLUP/src/MIBLUP.jl")
output_dir = "MIBLUP"
bfile = "data_qc"
mkdir(output_dir)

println("Start read plink")
data = MIBLUP.read_plink(bfile)
#covar = readdlm("covar.txt", Float64)
# 进行基因型插补
#geno = MIBLUP.geno_imputation(data.bed)
geno = data.bed
# 获取y向量
y = Vector(data.fam[:, "Trait_1"])

test_index = Matrix{Int}(readdlm("test_index.txt"))

yHat = zeros(size(data.fam, 1), 20)

for i in 1:20
    fold = "fold_$(i)_"
    whichNa = test_index[:, i]
    yNa = Vector{Union{Missing, typeof(y[1])}}(copy(y))
    yNa[whichNa] .= missing
    res = MIBLUP.runMIBLUP(yNa, geno;
                           ids = data.fam[:, "Individual_ID"],
                           snp_names = data.bim[:, "SNP_ID"],
                           n_pre_sels = 20)
    vc = DataFrame(vg = res.vc.vg, ve = res.vc.ve)
    CSV.write(joinpath(output_dir, fold*"solve.txt"), res.res)
    CSV.write(joinpath(output_dir, fold*"covar_effevts.txt"), res.cov_effects)
    CSV.write(joinpath(output_dir, fold*"mutual_information.txt"), res.mutual_information)
    CSV.write(joinpath(output_dir, fold*"selected_marker_effects.txt"), res.selected_marker_effects)
    CSV.write(joinpath(output_dir, fold*"vc.txt"), vc)
    yHat[:, i] .= res.res.yHat
end

CSV.write("MIBLUP_results.txt", DataFrame(yHat, :auto);
          delim = '\t', header=false)

println("MIBLUP分析已完成")
