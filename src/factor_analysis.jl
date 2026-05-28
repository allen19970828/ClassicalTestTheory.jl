using LinearAlgebra
using Statistics
using FactorRotations

export paf, efa, EFAResult

"""
    EFAResult

探索性因素分析的結果型別。
- `loadings`: 轉軸後的載荷量矩陣
- `communality`: 各題的共同性 (h²)
- `eigenvalues`: 初始的縮減相關矩陣特徵值
- `phi`: 因素間的相關矩陣 (斜交轉軸適用，直交轉軸則為單位矩陣)
"""
struct EFAResult
    loadings::Matrix{Float64}
    communality::Vector{Float64}
    eigenvalues::Vector{Float64}
    phi::Matrix{Float64}
end

"""
    paf(R::AbstractMatrix, n_factors::Int; max_iter=100, tol=1e-5)

主軸因素法 (Principal Axis Factoring, PAF) 因素萃取。
輸入 `R` 必須為相關係數矩陣。
回傳 Tuple: `(未轉軸的載荷量矩陣, 最終共同性, 特徵值)`。
"""
function paf(R::AbstractMatrix, n_factors::Int; max_iter=100, tol=1e-5)
    p = size(R, 1)
    
    # 初始共同性估計：SMC (Squared Multiple Correlation)
    # 若矩陣不可逆，則使用矩陣中的最大絕對相關係數
    smc = zeros(p)
    try
        R_inv = inv(R)
        for i in 1:p
            smc[i] = 1.0 - 1.0 / R_inv[i, i]
        end
    catch
        for i in 1:p
            smc[i] = maximum(abs.(R[i, filter(x -> x != i, 1:p)]))
        end
    end
    
    # 避免 SMC 超過 1 或小於 0
    smc = clamp.(smc, 0.0, 1.0)
    
    R_reduced = copy(R)
    for i in 1:p
        R_reduced[i, i] = smc[i]
    end
    
    h2 = copy(smc)
    loadings = zeros(p, n_factors)
    eigen_vals = zeros(p)
    
    for iter in 1:max_iter
        vals, vecs = eigen(Symmetric(R_reduced))
        
        # 排序特徵值
        idx = sortperm(vals, rev=true)
        vals = vals[idx]
        vecs = vecs[:, idx]
        
        if iter == 1
            eigen_vals = vals
        end
        
        # 計算載荷量
        for j in 1:n_factors
            if vals[j] > 0
                loadings[:, j] = vecs[:, j] .* sqrt(vals[j])
            else
                loadings[:, j] .= 0.0
            end
        end
        
        # 計算新的共同性
        new_h2 = vec(sum(loadings.^2, dims=2))
        new_h2 = clamp.(new_h2, 0.0, 1.0) # 避免 Heywood case 導致發散
        
        # 檢查收斂
        diff = maximum(abs.(new_h2 .- h2))
        if diff < tol
            h2 = new_h2
            break
        end
        
        h2 = new_h2
        for i in 1:p
            R_reduced[i, i] = h2[i]
        end
    end
    
    # 統一符號：確保每個因素中絕對值最大的載荷量為正
    for j in 1:n_factors
        max_idx = argmax(abs.(loadings[:, j]))
        if loadings[max_idx, j] < 0
            loadings[:, j] = -loadings[:, j]
        end
    end
    
    return loadings, h2, eigen_vals
end

"""
    efa(X::AbstractMatrix, n_factors::Int; cor_type=:pearson, rotation=Varimax(), smooth=true)

探索性因素分析 (Exploratory Factor Analysis) 完整流程。
支援的 `cor_type`：`:pearson` 或 `:polychoric`。
支援的 `rotation`：來自 FactorRotations.jl 的轉軸方法，如 `Varimax()`, `Promax()`, `Oblimin()` 等。
若不轉軸，請傳入 `nothing`。
預設 `smooth=true` 會自動將非正定矩陣平滑化。
回傳 `EFAResult`。
"""
function efa(X::AbstractMatrix, n_factors::Int; cor_type=:pearson, rotation=Varimax(), smooth=true)
    if cor_type == :polychoric
        R = polychoric_matrix(X)
    else
        R = pairwise_cor(X)
    end
    
    if smooth
        R = cor_smooth(R)
    end
    
    L_unrotated, h2, evs = paf(R, n_factors)
    
    if rotation === nothing
        phi = Matrix{Float64}(I, n_factors, n_factors)
        return EFAResult(L_unrotated, h2, evs, phi)
    end
    
    rot_res = rotate(L_unrotated, rotation)
    L_rot = rot_res.L
    
    # 斜交轉軸會有 phi (因素相關矩陣)，直交轉軸則沒有這個 field
    if hasproperty(rot_res, :phi)
        Phi = rot_res.phi
    else
        Phi = Matrix{Float64}(I, n_factors, n_factors)
    end
    
    return EFAResult(L_rot, h2, evs, Phi)
end

"""
    parallel_analysis(X::AbstractMatrix; n_iter=20, cor_type=:pearson, centile=95.0)

執行平行分析 (Parallel Analysis)，以決定應抽取之因素數量。
此方法透過蒙地卡羅模擬產生與原始資料相同維度的隨機常態矩陣，
計算其相關矩陣的特徵值，並與原始資料的特徵值進行比較。
若原始特徵值大於隨機特徵值的特定百分位數 (預設 95%)，則建議保留該因素。

回傳建議的因素數量 `n_factors`。
"""
function parallel_analysis(X::AbstractMatrix; n_iter=20, cor_type=:pearson, centile=95.0)
    n_subj, n_vars = size(X)
    
    # 原始資料特徵值
    if cor_type == :polychoric
        R_obs = polychoric_matrix(X)
    else
        R_obs = pairwise_cor(X)
    end
    # 平滑化以策安全
    R_obs = cor_smooth(R_obs)
    
    vals_obs = reverse(sort(real.(eigvals(Symmetric(R_obs)))))
    
    # 隨機資料特徵值矩陣
    random_vals = zeros(n_iter, n_vars)
    
    for i in 1:n_iter
        # 產生隨機常態資料
        X_rand = randn(n_subj, n_vars)
        R_rand = cor(X_rand)
        vals_rand = reverse(sort(real.(eigvals(Symmetric(R_rand)))))
        random_vals[i, :] = vals_rand
    end
    
    # 計算百分位數 (手動寫一個簡單的 percentile，避免依賴不同版本的 Statistics quantile，
    # 雖然 Statistics 有 quantile，但自己算 sort 後取 index 也很穩)
    target_vals = zeros(n_vars)
    for j in 1:n_vars
        col = sort(random_vals[:, j])
        idx = max(1, min(n_iter, round(Int, n_iter * (centile / 100.0))))
        target_vals[j] = col[idx]
    end
    
    # 決定因素數
    n_factors = 0
    for j in 1:n_vars
        if vals_obs[j] > target_vals[j]
            n_factors += 1
        else
            break
        end
    end
    
    # 至少為 1
    return max(1, n_factors)
end
