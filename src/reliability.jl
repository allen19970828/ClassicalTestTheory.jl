function cronbach_alpha(X::AbstractMatrix)
    J = size(X, 2)
    J <= 1 && return 0.0
    C = pairwise_cov(X)
    item_vars = [C[i,i] for i in 1:J]
    total_var = sum(C)
    total_var == 0 && return 0.0
    return (J / (J - 1)) * (1.0 - sum(item_vars) / total_var)
end

"""
計算指定分量表的 Cronbach's Alpha
"""
function cronbach_alpha(X::AbstractMatrix, scales::Dict{String, Vector{Int}})
    alphas = Dict{String, Float64}()
    for (scale_name, item_indices) in scales
        alphas[scale_name] = cronbach_alpha(X[:, item_indices])
    end
    return alphas
end

"""
    split_half(X::AbstractMatrix{<:Real})

計算折半信度 (Split-half reliability)。
預設將題目分為奇數題與偶數題兩半，計算兩半總分的相關後，
再使用 Spearman-Brown 公式進行校正。
"""
function split_half(X::AbstractMatrix)
    J = size(X, 2)
    J < 2 && return 0.0
    
    odd_items = X[:, 1:2:end]
    even_items = X[:, 2:2:end]
    
    # 漏答當作沒分 (補 0) 計算總分
    odd_scores = vec(sum(coalesce.(odd_items, 0), dims=2))
    even_scores = vec(sum(coalesce.(even_items, 0), dims=2))
    
    # 避免標準差為 0 的情況
    if std(odd_scores) == 0 || std(even_scores) == 0
        return 0.0
    end
    
    r = cor(odd_scores, even_scores)
    # 使用 Spearman-Brown 公式校正 (n=2)
    # 由於我們沒有在 reliability.jl 引入 scoring.jl 的 spearman_brown，
    # 這裡直接實作或是透過 ClassicalTestTheory 模組調用。為避免循環依賴，直接寫公式：
    return (2 * r) / (1 + r)
end

"""
    guttman_lambda(X::AbstractMatrix{<:Real})

計算 Guttman's lambda 信度 (Lambda 1 到 6)。
回傳一個包含各 Lambda 值的 Dict。
其中 Lambda 3 等同於 Cronbach's Alpha。
"""
function guttman_lambda(X::AbstractMatrix)
    n_subj, n_items = size(X)
    if n_items < 2 || n_subj < 2
        return Dict("L1"=>0.0, "L2"=>0.0, "L3"=>0.0, "L4"=>0.0, "L5"=>0.0, "L6"=>0.0)
    end
    
    C = pairwise_cov(X)
    Vt = sum(C)
    if Vt == 0
        return Dict("L1"=>0.0, "L2"=>0.0, "L3"=>0.0, "L4"=>0.0, "L5"=>0.0, "L6"=>0.0)
    end
    
    Vj = [C[i,i] for i in 1:n_items]
    sum_Vj = sum(Vj)
    
    L1 = 1.0 - sum_Vj / Vt
    
    C_offdiag = C - diagm(Vj)
    sum_sq_off = sum(C_offdiag.^2)
    L2 = L1 + sqrt((n_items / (n_items - 1)) * sum_sq_off) / Vt
    
    L3 = (n_items / (n_items - 1)) * L1
    
    L4 = NaN # Max split-half is computationally heavy to implement simply
    
    max_cov_sq = maximum(C_offdiag.^2)
    L5 = L1 + 2 * sqrt(max_cov_sq) / Vt
    
    R = pairwise_cor(X)
    smc = zeros(n_items)
    try
        R_inv = inv(R)
        for i in 1:n_items
            smc[i] = 1.0 - 1.0 / R_inv[i, i]
        end
    catch
        smc .= 0.0
    end
    L6 = 1.0 - sum((1 .- smc) .* Vj) / Vt
    
    return Dict("L1" => L1, "L2" => L2, "L3" => L3, "L4" => L4, "L5" => L5, "L6" => L6)
end

"""
    mcdonald_omega(X::AbstractMatrix{<:Real})

計算 McDonald's Omega 信度。
此函數利用主軸因素法 (PAF) 擷取第一個潛在因素的載荷量，
以作為單一因素模型的近似，計算出 Omega 信度，對齊 psych 套件的嚴謹算法。
"""
function mcdonald_omega(X::AbstractMatrix)
    n_subj, n_items = size(X)
    if n_items < 2 || n_subj < 2
        return 0.0
    end
    
    R = pairwise_cor(X)
    # 若有變異數為 0 導致 NaN
    if any(isnan, R)
        return 0.0
    end
    
    # 使用主軸因素法 (PAF) 抽取 1 個因素
    loadings, _, _ = paf(R, 1)
    loadings = vec(loadings)
    
    # 確保載荷量總和為正 (雖然 paf 內已有做符號統一，但再確認一次)
    if sum(loadings) < 0
        loadings = -loadings
    end
    
    sum_loadings = sum(loadings)
    sum_uniqueness = sum(1.0 .- loadings.^2)
    
    omega = (sum_loadings^2) / (sum_loadings^2 + sum_uniqueness)
    return omega
end


