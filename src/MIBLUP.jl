module MIBLUP
include("VarianceCompoent.jl")
include("read_file.jl")

using .VarianceCompoent
using StatsBase, Statistics, LinearAlgebra, DataFrames, Base.Threads, CSV, Random, HypothesisTests
#using MKL
export cal_vc, runMIBLUP, read_geno_file, compute_accuracy, discretize, cal_kinship, cal_entropy, cal_mi, mrmr, read_plink, geno_imputation

BLAS.set_num_threads(Threads.nthreads())

"""
    runMIBLUP(phe::Vector, geno::Matrix;
                  covariance = nothing,
                  kinship = nothing,
                  ids = nothing,
                  snp_names = nothing,
                  cov_names = nothing,
                  vc = (vg = nothing, ve = nothing),
                  kinship_weight::Bool = true,
                  marker_selection::Bool = true,
                  cal_mode = "speed",
                  vc_method = "he+ai",
                  phe_type = "continuous",
                  mi_target = "phe",
                  nbins = 10,
                  ngrids=100,
                  llim=log(0.01),
                  ulim=log(100),
                  esp=1e-10,
                  init = (vg = 0.2, ve = 0.8),
                  max_iter = 30,
                  cc = 1.0e-6,
                  n_sels::Int = 20,
                  n_pre_sels::Int = 1000,
                  verbose::Bool = true)

Run Mutual Information BLUP analysis.

# Arguments
- `phe::Vector`: Phenotypic values.
- `geno::Matrix`: Genotypic data matrix.
- `covariance=nothing`: Covariate matrix.
- `kinship=nothing`: Kinship matrix.
- `ids=nothing`: Individual IDs.
- `snp_names=nothing`: SNP names.
- `cov_names=nothing`: Covariate names.
- `vc=(vg=nothing, ve=nothing)`: Variance component.
- `kinship_weight::Bool=true`: Use kinship weighting.
- `marker_selection::Bool=true`: Perform marker selection.
- `cal_mode="speed"`: Kinship calculation mode.
- `vc_method="he+ai"`: Variance component estimation method ("he", "emma", "ai", "he+ai").
- `phe_type="continuous"`: Type of phenotype (continuous or discrete).
- `mi_target="phe"`: Target variable for mutual information and mrmr (gebv or pheno).\n
                     The main purpose of using `mi_target="gebv"` is to remove the effect of covariance. \n
                     If gebv is used, it will use `solve_mme` to calculate gebv, which will increase the computate time.
- `nbins=10`: Number of bins for discretizing phe or gebv. \n
              This parameter is ignored if the `phe_type=="discrete" && mi_target==“phe”`. \n
              If`mi_target==“gebv” && phe_type=="discrete"`, then it is better to set nbins equal to the number of levels.
-`acc_type`::String (optional): Type of accuracy or correlation to compute. Valid options are "cor", "auc", "rmse". Default is "cor".
- `ngrids=100`:(`vc_method = "emma"`) Number of grid points for log delta (default: 100).
- `llim=log(0.01)`: (`vc_method = "emma"`)Lower limit for log delta grid (default: log(0.01)).
- `ulim=log(100)`: (`vc_method = "emma"`)Upper limit for log delta grid (default: log(100)).
- `esp=1e-10`: (`vc_method = "emma"`)Tolerance for checking convergence (default: 1e-10).
- `init=(vg=0.2, ve=0.8)`: (`vc_method = "ai"`)Initial values for variance components in AI-REML (default: (vg = 0.2, ve = 0.8)).
- `max_iter=30`: (`vc_method = "ai" or "he+ai"`)Maximum number of iterations for AI-REML (default: 30).
- `cc=1.0e-6`: (`vc_method = "ai" or "he+ai"`) Convergence criterion for AI-REML (default: 1.0e-6).
- `n_sels::Int=20`: Number of marker selections.
- `n_pre_sels::Int=1000`: Number of pre-selections.
- `verbose::Bool=true`: Verbose output.

# Returns
- `res::DataFrame`: Results including ID, GEBV, and yHat.
- `cov_effects::DataFrame`: Covariate effects.
- `mutual_information::DataFrame`: Mutual information for SNPs.
- `selected_marker_effects::DataFrame`: Selected marker effects.
- `vc`: Variance components.
- `kinship`: Kinship matrix.
"""
function runMIBLUP(
    phe::Vector,
    geno::Matrix;
    covariance = nothing,
    kinship = nothing, # the kinship matrix if it is defined
    ids = nothing,
    snp_names = nothing,
    cov_names = nothing,
    vc = (vg = nothing, ve = nothing), # variance component
    kinship_weight::Bool = true,
    marker_selection::Bool = true,
    cal_mode = "speed", # calculation mode of kinship matrix
    vc_method = "he+ai",
    phe_type = "continuous", # continuous or discrete
    mi_target = "phe", # the target variable for mutual information, which can be either “gebv” or “pheno”.
    nbins = 10, # the number of bins for discretizing a continuous variable or a gebv variable when calculating mutual information.
    n_cv = 5,
    random_seed = 123,
    acc_type = "cor", #::String (optional): Type of accuracy or correlation to compute. Valid options are "cor", "auc", "rmse". Default is "cor".
    acc_incre = 0.001,
    ngrids=100, llim=log(0.01), ulim=log(100), esp=1e-10, init = (vg = 0.2, ve = 0.8), max_iter = 30, cc = 1.0e-6, # Argument of cal_vc
    n_pca::Int = 10, n_k::Int = 5, n_sels::Int = 20, n_pre_sels::Int = 1000, verbose::Bool = true, # Argument of mRMR
)
    # t1 = time()
    inf_index = ismissing.(phe) # 检查表型缺失值
    ref_index = .!inf_index     # 生成非缺失值的索引
    y = phe[ref_index]          # 获取非缺失值的表型数据
    N = length(phe)             # 总个体数量
    n = length(y)               # 非缺失值的个体数量
    m = size(geno, 2)           # SNP数量

    if isnothing(covariance)
        covariance = ones(N, 1) # 如果没有指定协变量矩阵，则将其设置为全1矩阵
    end
    n_cov = size(covariance, 2) # 协变量的数量
    X = covariance
    K = kinship
    @assert size(covariance, 1) == N
    if !isnothing(kinship)
        @assert size(kinship, 1) == size(kinship, 2) == N
    end

    if isnothing(ids)
        ids = 1:N               # 如果没有指定个体ID，则设置默认值为1到N
    end
    if isnothing(snp_names)
        snp_names = 1:m         # 如果没有指定SNP名称，则设置默认值为1到m
    end
    if isnothing(cov_names)
        cov_names = 1:n_cov     # 如果没有指定协变量名称，则设置默认值为1到n_cov
    end

    if verbose
        println("
Mutual Information BLUP\n
--------------------------------------------------------------------------------
")
        println("   Number of individuals: ", N)
        println("   Number of individuals with phenotype: ", n)
        println("   Number of SNPs: ", m)
        println("   Number of covariance: ", n_cov)
        println("")
    end

    # Calculate Variance Component

    if isnothing(vc.vg) || isnothing(vc.ve)
        # 计算 kinship 矩阵
        if isnothing(K)
            K = cal_kinship(geno; cal_mode = cal_mode, verbose = verbose)
        end
        # K = cal_kinship(geno; cal_mode = cal_mode, verbose = verbose)
        # 计算方差组分，得到 lambda 值
        vc = cal_vc(phe, covariance, K; method = vc_method,
                     ngrids=ngrids, llim=llim, ulim=ulim, esp=esp,
                     init = init, max_iter = max_iter, cc = cc, verbose = verbose)
        # vc = cal_vc(phe, covariance, K)
    end
    lambda = vc.ve / vc.vg

    # Calculate Mutual Information
    mi = nothing
    if kinship_weight || marker_selection
        if mi_target == "phe"
            if phe_type == "discrete"
                dis_y = Vector{Int}(y)
            elseif phe_type == "continuous"
                dis_y = discretize(y, nbins)# 将连续表型离散化
            else
                throw(DomainError(phe_type, "The phe_type=$(phe_type) argument is error!"))
            end
        elseif mi_target == "gebv"
            if isnothing(K)
                K = cal_kinship(geno; cal_mode = cal_mode, verbose = verbose)
            end
            if verbose
                println("To estimate the Best Linear Unbiased Predictions (BLUP) using the MME method.\n")
            end
            vals, vecs = eigen(K)
            
            beta, gebv = solve_mme(phe, hcat(X, vecs[:,sortperm(vals, rev=true)[1:10]]), K, lambda) # 估计GEBV
            dis_y = discretize(gebv, nbins)[ref_index] # 将GEBV离散化
        else
            throw(DomainError(mi_target, "The mi_target=$(mi_target) argument is error!"))
        end
    end

    # 计算权重（snp weights）

    if kinship_weight && isnothing(kinship)
        # 计算权重
        mi = cal_mi(geno[ref_index, :], dis_y) # calculate mi
        snp_weights = mi ./ sum(mi) .* m       # calculate weights
        Kw = cal_kinship(geno; weights = snp_weights, type = "scale", cal_mode = cal_mode, verbose = verbose) # calculate weighted kinship
    end

    if isnothing(K)
        K = cal_kinship(geno; cal_mode = cal_mode, verbose = verbose)
    end

    # Marker Selection

    if marker_selection
        # 使用 mrmr 选择 SNP
        if isnothing(mi)
            mrmr_selected_snps = mrmr(geno[ref_index, :], dis_y; n_sels = n_sels, n_pre_sels = n_pre_sels, verbose = verbose)
        else
            mrmr_selected_snps = mrmr(geno[ref_index, :], dis_y; n_sels = n_sels, n_pre_sels = n_pre_sels, mis_geno_y = mi, verbose = verbose)
        end
        # X = hcat(X, geno[:, selected_snps]) # 将选择的SNP加入协变量矩阵
    end

    function cv(y, covariance, kinship;lambda=1.0,n_cv = 5,acc_type = "cor", random_seed = 123)
        # 设置随机数种子
        Random.seed!(random_seed)
        # 确保输入数据的维度匹配
        n_samples = length(y)
        n_covar = size(covariance, 2)
        indices = collect(1:n_samples)
        shuffle!(indices)
        fold_size = Int(div(n_samples, n_cv))
        cv_accs = []
        for fold in 1:n_cv
            start_index = (fold - 1) * fold_size + 1
            end_index = fold * fold_size
            # 在当前折中选择索引
            if fold == n_cv && end_index < n_samples
                current_indices = indices[start_index:end]
            else
                current_indices = indices[start_index:end_index]
            end
            train_indices = setdiff(indices, current_indices)
            y_train = Vector{Union{Missing, typeof(y[1])}}(y)
            y_train[current_indices] .= missing
            beta, u = solve_mme(y_train, covariance, kinship, lambda)
            yHat = u .+ covariance * beta
            acc = compute_accuracy(y[current_indices], yHat[current_indices], acc_type)
            push!(cv_accs, acc)
        end
        # println(indices[1:5])
        # return mean(cv_accs)
        return cv_accs
    end

    function snp_selection(y, geno, covariance, kinship;set_acc = [], snps = [], lambda=1.0,n_cv = 5,acc_type = "cor", acc_incre = acc_incre, verbose = true, random_seed = 123)
        # set_acc = cv(y, covariance, kinship; lambda=lambda,n_cv = n_cv,acc_type = acc_type, random_seed = random_seed)
        X = copy(covariance)
        selected_snps = []
        k = 0
        for snp in snps
            acc = cv(y, hcat(X, geno[:, snp]), kinship; lambda=lambda,n_cv = n_cv,acc_type = acc_type, random_seed = random_seed)
            p = pvalue(OneSampleTTest(acc .- set_acc))
            improved_acc = mean(acc) - mean(set_acc)
            if improved_acc >= acc_incre && p < 0.05
                push!(selected_snps, snp)
                X = hcat(X, geno[:, snp])
                if verbose
                    println("    Selecting the $(snp)th marker as a covariate improves the accuracy by $(round(improved_acc, digits = 4)), with a p-value of $(p).\n")
                end
                set_acc = acc
                k = 0
            else
                if verbose
                    println("    Rejecting the selection of the $(snp)th marker as a covariate, the accuracy is $(round(improved_acc, digits = 4)), the p-value is $(p).\n")
                end
                k += 1
                if k >= n_k
                    break
                end
            end
        end
        if verbose
            println("Selection completed. $(length(selected_snps)) SNPs were selected as covariates.\n")
        end
        return selected_snps
    end

    # Select kinship
    kinship_type = "varanden kinship"
    if kinship_weight && isnothing(kinship)
        set_acc = cv(y, X[ref_index, :], Symmetric(K[ref_index, ref_index]); lambda = lambda,n_cv = n_cv,acc_type = acc_type, random_seed = random_seed)
        acc = cv(y, X[ref_index, :], Symmetric(Kw[ref_index, ref_index]); lambda=lambda,n_cv = n_cv,acc_type = acc_type, random_seed = random_seed)
        p = pvalue(OneSampleTTest(acc .- set_acc))
        improved_acc = mean(acc) - mean(set_acc)
        if improved_acc >= acc_incre && p < 0.05
            K = Kw
            kinship_type = "weighted kinship"
            if verbose
                println("\nUsing weighted kinship significantly improves the accuracy by $(improved_acc), with a p-value of $(p)\n")
            end
            set_acc .= acc
        end
    end

    if verbose
        println("Use $(kinship_type) as genetic relationship matrix.\n")
    end

    if marker_selection
        selected_snps = snp_selection(y, geno[ref_index, :], X[ref_index, :], Symmetric(K[ref_index, ref_index]);
                                      set_acc = set_acc,
                                      snps = mrmr_selected_snps,
                                      lambda = lambda,
                                      n_cv = n_cv,
                                      acc_type = acc_type,
                                      acc_incre = acc_incre,
                                      verbose = verbose,
                                      random_seed = random_seed)
        X = hcat(covariance, geno[:, selected_snps]) # 将选择的SNP加入协变量矩阵
    end

    if verbose
        println("To estimate the Best Linear Unbiased Predictions (BLUP) using the MME method.\n")
    end
    # solve MME
    beta, u = solve_mme(phe, X, K, lambda)
    yHat = u .+ X * beta # y的估计值
    cov_effects = DataFrame(;covariance = cov_names, values = beta[1:n_cov]) # 构建协变量的DataFrame
    gebv = u
    mutual_information = DataFrame(;snp_name = snp_names, mutual_information = mi) # 构建互信息的DataFrame
    if marker_selection
        selected_marker_effects = DataFrame(;snp_name = snp_names[selected_snps], values = beta[n_cov + 1: length(beta)]) # 构建选择的SNP的DataFrame
        gebv += X[:, n_cov + 1 : length(beta)] * beta[n_cov + 1: length(beta)] # 更新GEBV，加上选择的SNP效应
    else
        selected_marker_effects = nothing
    end
    res = DataFrame(;ID = ids, GEBV = gebv, yHat = yHat)
    if verbose
        println("End of MIBLUP.\n")
    end
    return (res = res,
            cov_effects = cov_effects,
            mutual_information = mutual_information,
            selected_marker_effects = selected_marker_effects,
            vc = vc, kinship = K, kinship_type = kinship_type)
end

# MIBLUP(y, geno)
# cal_vc(phe, covariance, K)

"""
solve_mme(phe::Vector, covariance::Matrix, K::Symmetric, lambda::Float64 = 1.0)

To estimate the Best Linear Unbiased Predictions (BLUP) using the MME method.

# Arguments
-`phe`: The phenotypes vector
-`covariance`: covariance, design matrix(n * x) for the fixed effects(covariance must contain a column of 1's)
-`K`: Kinship matrix for all individuals
-`lambda`: ve/vg, ve is the variance of residual, vg is the variance of additive effect (default: 1.0)

# Returns
-`beta`: The estimated fixed covariance effects coefficients
-`u`: The Best Linear Unbiased Predictions (BLUP) for each individual
"""
function solve_mme(phe::Vector, covariance::Matrix, K::Symmetric, lambda::Float64 = 1.0)
    inf_index = ismissing.(phe)
    ref_index = .!inf_index
    refphe = phe[ref_index]
    y = Vector{Float64}(refphe)
    X = Matrix{Float64}(covariance[ref_index, :])

    n = sum(ref_index)
    m = size(K, 1)

    # Compute Cholesky factorization of K[ref_index, ref_index] + lambda * I
    L = cholesky(K[ref_index, ref_index] + I * lambda)
    # L = inv(L.U)
    # Calculate y^T * U and X^T * U
    yt = L.L \ y
    Xt = L.L \ X
    # Compute X^T * X
    # XtX = Xt' * Xt
    # Calculate beta using ldiv!
    # beta = XtX \ (Xt' * yt)
    beta = Xt \ yt
    # Dt = L.L \ (y .- X * beta)
    S = L \ (y .- X * beta)
    # Calculate u
    u = Vector{Float64}(undef, length(phe))
    mul!(u, K[:, ref_index], S)
    # Calculate yHat = u + X * beta
    # yHat = u .+ covariance * beta
    return (beta = beta, u = u)
end


"""
cal_kinship(geno::Matrix; weights::Union{Vector, Nothing} = nothing, type = "vanraden", cal_mode = "speed", verbose::Bool = true)

To calculate the kinship matrix based on genotype data.

# Arguments
-`geno`: The genotype matrix (rows represent individuals, columns represent markers)
-`weights`: Weights for each marker (default: nothing)
-`type`: Type of kinship calculation ("vanraden" or "scale", default: "vanraden")

# Returns
-`K`: The kinship matrix
"""
function cal_kinship(geno::Matrix; weights::Union{Vector, Nothing} = nothing, type = "vanraden", cal_mode = "speed", verbose::Bool = true)
    if verbose
        println("Use $(cal_mode) mode to construct the kinship matrix.\n")
    end
    M = Matrix{Float64}(geno)
    n, m = size(M)
    p = mean(M, dims=1) ./ 2
    M .-= 2 .* p
    if type == "scale"
        sd = sqrt.(2 * p .* (1 .- p))
        M ./= sd
    end
    if verbose
        if !isnothing(weights)
            println("Weighting the kinship matrix.\n")
        end
    end
    SUM = ifelse(type == "vanraden", 2 * sum(p .* (1 .- p)), m)
    if cal_mode == "speed"
        if !isnothing(weights)
            tMw = M'
            @views tMw .*= sqrt.(weights)
            M = tMw'
        end
        K = M * M' ./ SUM
    elseif cal_mode == "memory"
        block_size = 1024
        K = zeros(n, n)
        if !isnothing(weights)
            @threads for j in 1:m
                M[:, j] .*= sqrt(weights[j])
            end
        end
        if verbose
            progress_total = div(n, block_size) + (n % block_size > 0 ? 1 : 0)
            progress_current = 0
            println("Calculating kinship: ")
        end
        @threads for i_start in 1:block_size:n
            i_end = min(i_start + block_size - 1, n)
            @threads for j_start in i_start:block_size:n
                j_end = min(j_start + block_size - 1, n)
                block = M[i_start:i_end, :] * M[j_start:j_end, :]'
                K[i_start:i_end, j_start:j_end] += block
                if i_start != j_start
                    K[j_start:j_end, i_start:i_end] += block'
                end
            end
            if verbose
                progress_current += 1
                progress_percent = round(progress_current / progress_total * 100, digits=1)
                println("[$(repeat("=", round(Int, progress_percent * 0.5)))| $(progress_percent)%] ")
                flush(stdout)
            end
        end
        K /= SUM
    else
        error("Invalid calculation method specified.")
    end
    if verbose
        println("The kinship matrix has been constructed.\n")
    end
    return Symmetric(K)
end

function cal_entropy(x)
    entropy(counts(x) / length(x), 2)
end

function cal_entropy(X::Matrix)
    m = size(X, 2)
    entropies = Vector{Float64}(undef, m)
    @threads for j in 1:m
        entropies[j] = cal_entropy(view(X, :, j))
    end
    return entropies
end

function cal_mi(vec::Vector, y::Vector)
    # counts_y = counts(y)
    counts_map = countmap(y)
    counts_y = values(counts_map)
    p_y = counts_y ./ length(y)
    entropy_x = cal_entropy(vec)
    cond_entropy_xy = 0.0
    # cond_entropy_xy = [cal_entropy(vec[y .== v]) for v in StatsBase.levels(y)]' * p_y
    cond_entropy_xy = [cal_entropy(vec[y .== v]) for v in keys(counts_map)]' * p_y
    mi = entropy_x - cond_entropy_xy
    return mi
end

function cal_mi(X::Matrix, y::Vector)
    n, m = size(X)
    # counts_y = counts(y)
    counts_map = countmap(y)
    counts_y = values(counts_map)
    keys_y = keys(counts_map)
    p_y = counts_y ./ n
    entropy_y = entropy(p_y, 2)

    mis = Vector{Float64}(undef, m)
    # Threads.@threads for j in 1:m
    @inbounds @simd for j in 1:m
        vec = view(X, :, j)
        entropy_x = cal_entropy(vec)
        cond_entropy_xy = [cal_entropy(vec[y .== v]) for v in keys(counts_map)]' * p_y
        # cond_entropy = Vector{Float64}(undef, length(keys_y))
        # @threads for i in 1:length(keys_y)
        #     v = keys_y[i]
        #     cond_entropy[i] = cal_entropy(vec[y .== v])
        # end
        # cond_entropy_xy = cond_entropy' * p_y
        # print(mi)
        mis[j] = entropy_x - cond_entropy_xy
    end

    return mis
end


"""
mrmr(geno::Matrix, y::Vector; n_sels::Int = 20, n_pre_sels::Int = 1000, mis_geno_y = nothing, verbose::Bool = true)

The mrmr function performs feature selection using the Minimum Redundancy Maximum Relevance (mRMR) algorithm.

# Arguments
-`geno`::Matrix: A matrix representing the genotype data. Each row corresponds to a sample, and each column corresponds to a genetic feature.
-`y`::Vector: A vector representing the target variable or class labels. It should have the same length as the number of samples in geno.
-`n_sels`::Int = 20 (optional): The number of features to select. Default is 20.
-`n_pre_sels`::Int = 1000 (optional): The number of pre-selected features. Default is 1000.
-`verbose`::Bool = true (optional): Whether to display the progress and intermediate results. Default is true.

# Returns
-`sels`: The selected markers.
"""
function mrmr(geno::Matrix, y::Vector; n_sels::Int = 20, n_pre_sels::Int = 1000, mis_geno_y = nothing, verbose::Bool = true)
    t1 = time()
    if verbose
        println("Minimum Redundancy Maximum Relevance (mRMR)")
        println("
--------------------------------------------------------------------------------
")
    end
    n, m = size(geno)
    t = zeros(m)
    if verbose
        println("Pre-select the top $(n_pre_sels) SNPs.\n")
        println("Select $(n_sels) SNPs using mRMR algorithm.\n")
    end
    if isnothing(mis_geno_y)
        mis_geno_y = cal_mi(geno, y)
    end
    sort_idxs = sortperm(mis_geno_y, rev = true)
    fea_base = sort_idxs[1:n_sels]
    pre_sels = sort_idxs[1:n_pre_sels]
    idxleft = sort_idxs[2:n_pre_sels]
    k = 1
    if verbose
        progress_total = n_sels - 1
        progress_current = 0
        println("Selecting SNPs: ")
        # println("Currently selected 1 remaining 19.\n")
        # println("k=1 cost_time=$(time()-t1) cur_fea=$(sort_idxs[1]) n_left=$(length(idxleft))")
    end
    sels = Vector{Int}()
    push!(sels, sort_idxs[1])
    for k in 2:n_sels
        n_left = length(idxleft)
        n_cur_sels = length(sels)
        mis_left_y = mis_geno_y[idxleft]
        redundancy = similar(mis_left_y)
        for i in 1:n_left
            mis_i_sels = mean(cal_mi(geno[:, sels], geno[:, idxleft[i]]))
            redundancy[i] = mis_i_sels
        end
        tmp_idx = argmax(mis_left_y - redundancy)
        push!(sels, idxleft[tmp_idx])
        deleteat!(idxleft, tmp_idx)
        if verbose
            progress_current += 1
            if progress_current % (div(progress_total, 4)) == 0
                progress_percent = round(progress_current / progress_total * 100, digits=1)
                println("[$(repeat("=", round(Int, progress_percent * 0.5)))| $(progress_percent)%]")
                flush(stdout)
            end
            # println("Currently selected $(k) remaining $(n_sels - k).\n")
            # println("k=$(k) cost_time=$(time()-t1) cur_fea=$(sels[k]) #left_cand=$(length(idxleft))")
        end
    end
    # close(pb)
    return sels
end

"""
discretize(x, nbins::Int64 = 10)

Discretizes a continuous variable into nbins bins.
"""
function discretize(x, nbins::Int64 = 10)
    min_x = minimum(x)
    max_x = maximum(x)
    bin_width = (max_x - min_x) / nbins
    bin_edges = Vector{Float64}(undef, nbins+1)
    bins = Vector{Int64}(undef, length(x))

    for i in 0:nbins
        bin_edges[i+1] = min_x + i * bin_width
    end

    for (i, val) in enumerate(x)
        if val == max_x
            bins[i] = nbins
        else
            bin = findfirst(edge -> edge > val, bin_edges) - 1
            bins[i] = bin
        end
    end

    return bins
end


"""
Function to compute accuracy or correlation of predictions.

# Arguments
-`x1`: Vector of true values.
-`x2`: Vector of predicted values.
-`type`::String (optional): Type of accuracy or correlation to compute. Valid options are "cor", "auc", "rmse". Default is "cor".
Returns:
-`res`: Computed accuracy or correlation.
"""
function compute_accuracy(x1, x2, type::String = "cor")
    if type == "cor"
        res = cor(x1, x2)
    elseif type == "auc"
        if !(0 in unique(x1) && 1 in unique(x1))
            error("Only two levels (case/1, control/0) are allowed!")
        end
        # println("use auc")
        X = hcat(x1, x2)
        X = X[sortperm(X[:, 2]), :]
        N_case = sum(X[:, 1] .== 1)
        N_cont = sum(X[:, 1] .== 0)
        case_index = findall(X[:, 1] .== 1)
        case_index_mean = mean(case_index)
        res = (1 / N_cont) * (case_index_mean - 0.5 * N_case - 0.5)
    elseif type == "rmse"
        res = sqrt(sum((x2 - x1).^2) / length(x2))
    else
        error("Invalid type. Available options: 'cor', 'auc', 'rmse'.")
    end
    return res
end

end
