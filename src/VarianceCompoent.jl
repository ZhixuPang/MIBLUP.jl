module VarianceCompoent

using LinearAlgebra, Optim, Statistics
# using MKL

export cal_vc

"""
To estimate variance component using Haseman-Elston regression.
Input:
y: The phenotypes vector
X: The fixed effect(X must contain a column of 1's)
K: Kinship for all individuals
"""
function he(y, X, K; verbose::Bool = true)
    if verbose
        println("Genetic parameters were estimated by HE (Haseman–Elston) regression.\n")
    end
    n = size(K, 1)
    # center the Kinship
    w = ones(n)
    Gw = sum(K, dims=2)
    alpha = -1.0 / n
    K = K + alpha * (w * transpose(Gw) + Gw * transpose(w))
    function calc_v_che(y, X, K)
        n = size(K, 1)
        r = n / (n - size(X, 2))
        v_traceG = mean(diag(K))
        # center and scale K by X
        WtW = X' * X
        WtWi = inv(WtW)
        WtWiWt = WtWi * X'
        GW = K * X
        Gtmp = GW * WtWiWt
        K = K - Gtmp
        Gtmp = transpose(Gtmp)
        K = K - Gtmp
        WtGW = X' * GW
        GW = WtWiWt' * WtGW
        Gtmp = GW * WtWiWt
        K = K + Gtmp
        d = mean(diag(K))
        traceG_new = d
        if d != 0
            K = K / d
        end
        # center y by X, and standardize it to have variance 1 (transpose(y) * y / n = 1)
        Wty = X' * y
        WtWiWty = WtWi * Wty
        y_scale = y - X * WtWiWty
        function vector_var(x)
            m = mean(x)
            m2 = sum(x .^ 2) / length(x)
            return m2 - m * m
        end
        function standardize_vector(x)
            m = mean(x)
            v = sum((x .- m) .^ 2) / length(x)
            v = v - m * m
            x = (x .- m) ./ sqrt(v)
            return x
        end
        var_y = vector_var(y)
        var_y_new = vector_var(y_scale)
        y_scale = standardize_vector(y_scale)
        # compute Kry, which is used for confidence interval; also compute q_vec (*n^2)
        Kry = K * y_scale - r * y_scale
        q_vec = Kry' * y_scale
        # compute yKrKKry, which is used later for confidence interval
        KKry = K * Kry
        yKrKKry = Kry' * KKry
        d = Kry' * Kry
        yKrKKry = vcat(yKrKKry, d)
        # compute Sij (*n^2)
        tr = sum(K .^ 2)
        S_mat = tr - r * n
        # compute S^{-1}q
        Si_mat = 1 / S_mat
        # compute pve (on the transformed scale)
        pve = Si_mat * q_vec
        # compute q_var (*n^4)
        s = 1
        qvar_mat = yKrKKry[1] * pve
        s = s - pve
        tmp_mat = yKrKKry[2] * s
        qvar_mat = (qvar_mat + tmp_mat) * 2.0
        # compute S^{-1}var_qS^{-1}
        Var_mat = Si_mat * qvar_mat * Si_mat
        # transform pve back to the original scale and save data
        s = 1
        vyNtgN = var_y_new / traceG_new
        vtgvy = v_traceG / var_y
        v_sigma2 = pve * vyNtgN
        v_pve = pve * vyNtgN * vtgvy
        s = s - pve
        pve_total = pve * vyNtgN * vtgvy
        d = sqrt(Var_mat)
        v_se_sigma2 = d * vyNtgN
        v_se_pve = d * vyNtgN * vtgvy
        se_pve_total = Var_mat * vyNtgN * vtgvy * vyNtgN * vtgvy
        v_sigma2 = vcat(v_sigma2, s * r * var_y_new)
        v_se_sigma2 = vcat(v_se_sigma2, sqrt(Var_mat) * r * var_y_new)
        se_pve_total = sqrt(se_pve_total)
        return v_sigma2
    end
    # initialize sigma2/log_sigma2
    v_sigma2 = calc_v_che(y, X, K)
    log_sigma2 = Vector{Union{Float64, Missing}}(missing, length(v_sigma2))
    for i in 1:length(v_sigma2)
        if v_sigma2[i] <= 0
            log_sigma2[i] = log(0.1)
        else
            log_sigma2[i] = log(v_sigma2[i])
        end
    end

    function log_rl_dev1(log_sigma2, parms)
        y = parms[:y]
        X = parms[:X]
        K = parms[:K]
        n = size(K, 1)

        P = K * exp(log_sigma2[1]) + Matrix{Float64}(I, n, n) * exp(log_sigma2[2])

        # calculate H^{-1}
        P = inv(P)

        # calculate P = H^{-1} - H^{-1}X(X^TH^{-1}X)^{-1}X^TH^{-1}
        HiW = P * X
        WtHiW = X' * HiW
        WtHiWi = inv(WtHiW)
        WtHiWiWtHi = WtHiWi * HiW'
        P = P - HiW * WtHiWiWtHi

        # calculate Py, KPy, PKPy
        Py = P * y
        KPy = K * Py

        # calculate dev1 = -0.5 * trace(PK_i) + 0.5 * yPKPy
        para1 = (-0.5 * sum(P .* K) + 0.5 * transpose(Py) * KPy) * exp(log_sigma2[1])
        para2 = (-0.5 * sum(diag(P)) + 0.5 * transpose(Py) * ones(n)) * exp(log_sigma2[2])

        return (para1=para1, para2=para2)
    end

    parms = Dict(:y => y, :X => X, :K => K)
    res = log_rl_dev1(log_sigma2, parms)

    vg = exp(log_sigma2[1])
    ve = exp(log_sigma2[2])
    delta = ve / vg

    return (vg=vg, ve=ve)
end

"""
To estimate variance component using optimized EMMA.
Input:
y: The phenotypes vector
X: The fixed effect(X must contain a column of 1's)
K: Kinship for all individuals
ngrids: Number of grid points for delta (default: 100)
llim: Lower limit for log delta grid (default: 0)
ulim: Upper limit for log delta grid (default: 1)
esp: Tolerance for checking convergence (default: 1e-10)
"""
function emma_reml(y, X, K; ngrids=100, llim=log(0.01), ulim=log(100), esp=1e-10, verbose::Bool = true)
    if verbose
        println("Genetic parameters were estimated by approximating the REML using the EMMA algorithm.\n")
    end

    if !isa(y, Vector{Float64})
        y = Vector{Float64}(y)
    end
    function eigen_R_wo_Z(K, X)
        n = size(X, 1)
        q = size(X, 2)
        cXX = X' * X
        iXX = inv(cXX)
        # iXX = LinearAlgebra.svdfact!(cXX).Vt
        multiply(X, Y) = sum(X .* Y, dims=2)
        SS1 = X * iXX
        SS2 = SS1 * X'
        S = Matrix{Float64}(I, n, n) - SS2
        eig = eigen(S * (K + Matrix{Float64}(I, n, n)) * S)
        sorted_indices = sortperm(eig.values, rev=true)
        eigvalues = eig.values[sorted_indices]
        eigvectors = eig.vectors[:, sorted_indices]
        return (values=eigvalues[1:(n-q)] .- 1, vectors=eigvectors[:, 1:(n-q)])
    end
    function delta_REML_dLL_wo_Z(logdelta, lambda, etas)
        nq = length(etas)
        delta = exp(logdelta)
        etasq = etas .* etas
        ldelta = lambda .+ delta
        return 0.5 * (nq * sum(etasq ./ (ldelta .* ldelta)) / sum(etasq ./ ldelta) - sum(1 ./ ldelta))
    end
    function delta_REML_LL_wo_Z(logdelta, lambda, etas)
        nq = length(etas)
        delta = exp(logdelta)
        return 0.5 * (nq * (log(nq / (2 * pi)) - 1 - log(sum(etas .* etas ./ (lambda .+ delta)))) - sum(log.(lambda .+ delta)))
    end
    function multi(i, logdelta, dLL, eig_values, etas)
        # if (dLL[i] * dLL[i + 1] < 0) && (dLL[i] > 0) && (dLL[i + 1] < 0)
        #     println(i)
        #     println((dLL[i] * dLL[i + 1] < 0) && (dLL[i] > 0) && (dLL[i + 1] < 0))
        # end
        if (dLL[i] * dLL[i + 1] < 0) && (dLL[i] > 0) && (dLL[i + 1] < 0)
            r = optimize(logdelta -> delta_REML_dLL_wo_Z(logdelta, eig_values, etas), logdelta[i], logdelta[i + 1])
            optlogdelta_mult = r.minimum
            optLL_mult = delta_REML_LL_wo_Z(r.minimum, eig_values, etas)
            return (optlogdelta_mult=optlogdelta_mult, optLL_mult=optLL_mult)
        else
            return nothing
        end
    end
    eig_R = eigen_R_wo_Z(K, X)
    n = length(y)
    t = size(K, 1)
    q = size(X, 2)

    if size(K, 2) != t || size(X, 1) != n
        error("Incorrect dimensions of K and X")
    end

    cXX = X' * X
    if det(cXX) ≈ 0
        @warn("X is singular")
        return (REML=0, delta=0, ve=0, vg=0)
    end

    etas = eig_R.vectors' * y
    logdelta = range(llim, stop=ulim, length=ngrids+1)
    m = length(logdelta)
    delta = exp.(logdelta)

    Lambdas = eig_R.values .+ delta'
    Etasq = (etas .^ 2) * ones(1, m)
    LL = 0.5 * ((n - q) * (log((n - q) / (2 * pi)) .- 1 .- log.(sum(Etasq ./ Lambdas, dims=1))) - sum(log.(Lambdas), dims=1))
    dLL = 0.5 * delta .* ((n - q) * sum(Etasq ./ (Lambdas .* Lambdas), dims=1) ./ sum(Etasq ./ Lambdas, dims=1) .- sum(1 ./ Lambdas, dims=1))'

    optlogdelta = Float64[]
    optLL = Float64[]

    if dLL[1] < esp
        push!(optlogdelta, llim)
        push!(optLL, delta_REML_LL_wo_Z(llim, eig_R.values, etas))
    end

    if dLL[m-1] > 0 - esp
        push!(optlogdelta, ulim)
        push!(optLL, delta_REML_LL_wo_Z(ulim, eig_R.values, etas))
    end

    multi_res = Vector{Float64}[]
    for i in 1:(m-1)
        # println(i)
        # println(logdelta)
        res = multi(i, logdelta, dLL, eig_R.values, etas)
        # println(res)
        if res != nothing
            push!(multi_res, [res.optlogdelta_mult, res.optLL_mult])
        end
    end

    if !isempty(multi_res)
        multi_res = hcat(multi_res...)
        optlogdelta = vcat(optlogdelta, multi_res[1, :])
        optLL = vcat(optLL, multi_res[2, :])
    end

    maxdelta = exp(optlogdelta[argmax(optLL)])
    optLL = replace(x -> isnan(x) ? minimum(optLL[.!isnan.(optLL)]) : x, optLL)
    maxLL = maximum(optLL)
    maxva = sum(etas .* etas ./ (eig_R.values .+ maxdelta)) / (n - q)
    maxve = maxva * maxdelta

    return (ve=maxve, vg=maxva)
end

"""
To estimate variance component using average information REML (AI-REML).
Input:
y: The phenotypes vector
X: The fixed effect (X must contain a column of 1's)
K: Kinship for all individuals
init: Initial values for variance components (default: (vg = 0.2, ve = 0.8))
max_iter: Maximum number of iterations (default: 30)
cc: Convergence criterion (default: 1.0e-6)
"""
function ai_reml(y, X, K; init = (vg = 0.2, ve = 0.8), max_iter = 30, cc = 1.0e-6, verbose::Bool = true)
    if verbose
        println("Genetic parameters were estimated by approximating the REML using the AI algorithm.\n")
    end
    function trace(vec)
        return sum(diag(vec))
    end
    t1 = time()
    n = length(y)
    init = [init.vg, init.ve]
    vars = reshape(init, (2, 1)) # 方差组分
    if verbose
        println("   Initial variance: genetic variance is $(round(vars[1, 1]; digits=3)) and residual variance is $(round(vars[2, 1]; digits=3)).\n")
    end
    AI = zeros(Float64, (2, 2)) # Hessian矩阵
    s = zeros(Float64, (2, 1)) # 一阶偏导向量
    # eig_values, eig_vectors = eigen(K)
    for i in 1:max_iter
        if verbose
            print("   Iteration $(i): ")
        end
        G = K * vars[1, 1]
        R = diagm(ones(n)) * vars[2, 1]
        R = I(n) * vars[2, 1]
        V = G + R
        L = cholesky(V)
        Vinv = inv(L)
        Vinv_X = L \ X
        X_Vinv_X = X' * Vinv_X
        X_Vinv_X_inv = inv(X_Vinv_X)
        P = Vinv - Vinv_X * X_Vinv_X_inv * Vinv_X'
        PK = P * K
        yPK = y' * PK
        Py = P * y
        yPKP = yPK * P
        s[1, 1] = sum(0.5 .* (-trace(PK) .+ yPK * Py))
        s[2, 1] = sum(0.5 .* (-trace(P) .+ Py' * Py))
        AI[1, 1] = 0.5 * sum(yPKP * yPK')
        AI[2, 1] = 0.5 * sum(yPKP * Py)
        AI[1, 2] = AI[2, 1]
        AI[2, 2] = 0.5 * sum(Py' * P * Py)
        delta = AI \ s
        new_vars = vars + delta
        if new_vars[1, 1] < 0
            new_vars[1, 1] = 0.01
        end
        if new_vars[2, 1] < 0
            new_vars[2, 1] = 0.01
        end
        if verbose
            println("     genetic variance is $(round(new_vars[1, 1]; digits=3)) and residual variance is $(round(new_vars[2, 1]; digits=3)).\n")
        end
        # println("Time: ", time()-t1)
        vars = copy(new_vars)
        diff1 = abs(delta[1, 1])
        diff2 = abs(delta[2, 1])
        if diff1 < cc && diff2 < cc
            break
        end
    end
    return (vg=vars[1], ve=vars[2])
end


"""
To estimate variance component.

Input:

y: The phenotypes vector

covariance: The fixed effect(covariance must contain a column of 1's)

K: Kinship for all individuals

method: Method for estimating variance components ("he", "emma", "ai", "he+ai")

ngrids: Number of grid points for log delta (default: 100)

llim: Lower limit for log delta grid (default: log(0.01))

ulim: Upper limit for log delta grid (default: log(100))

esp: Tolerance for checking convergence (default: 1e-10)

init: Initial values for variance components in AI-REML (default: (vg = 0.2, ve = 0.8))

max_iter: Maximum number of iterations for AI-REML (default: 30)

cc: Convergence criterion for AI-REML (default: 1.0e-6)

Output:

vg: Genetic variance

ve: Residual variance
"""

function cal_vc(phe::Vector, covariance::Matrix, K::Symmetric;
    method = "he", ngrids=100, llim=log(0.01), ulim=log(100), esp=1e-10,
    init = (vg = 0.2, ve = 0.8), max_iter = 30, cc = 1.0e-6, verbose::Bool = true)
    inf_index = ismissing.(phe)
    ref_index = .!inf_index
    refphe = phe[ref_index]
    y = Float64.(refphe)
    N = length(phe)
    n = length(refphe)
    X = covariance[ref_index, :]
    # print(size(K))
    K = K[ref_index, ref_index]
    # test = X
    if verbose
        println("
Estimate Variance Component

--------------------------------------------------------------------------------
        ")
    end
    if method == "he"
        vcs = he(y, X, K; verbose = verbose)
    elseif method == "emma"
        vcs = emma_reml(y, X, K; ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, verbose = verbose)
    elseif method == "ai"
        vcs = ai_reml(y, X, K; init = init, max_iter = max_iter, cc = cc, verbose = verbose)
    elseif method == "he+ai"
        vg, ve = he(y, X, K; verbose = verbose)
        # print(vg)
        vcs = ai_reml(y, X, K; init = (vg = vg, ve = ve), max_iter = max_iter, cc = cc, verbose = verbose)
    else
        throw(DomainError(method, "The method argument is error!"))
    end
    if verbose
        println("
Estimated Result:

    Genetic variance    : $(round(vcs.vg; digits=3))
    Rresidual variance  : $(round(vcs.ve; digits=3))
    Heritability        : $(round(vcs.vg / (vcs.vg + vcs.ve); digits=3))
                 ")
    end
    return vcs
end

end  # module VarianceCompoent

# VarianceCompoent.cal_vc(phe, X, K; method = "ai")
# VarianceCompoent.cal_vc(phe, X, K; method = "emma", ngrids= 500, llim = log(0.01), ulim = log(100))
# VarianceCompoent.cal_vc(phe, X, K; method = "he")
# vcs = VarianceCompoent.cal_vc(phe, X, K; method = "he+ai")
