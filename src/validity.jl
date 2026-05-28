"""
    criterion_validity(X::AbstractMatrix{<:Real}, Y::AbstractVector{<:Real})

計算測驗總分與外部效標 Y 之間的效標關聯效度 (Criterion-related Validity)。
使用 Pearson 相關係數。
"""
function criterion_validity(X::AbstractMatrix, Y::AbstractVector)
    total_scores = vec(sum(coalesce.(X, 0), dims=2))
    valid_idx = .!(ismissing.(total_scores) .| ismissing.(Y))
    ts_v = Float64.(total_scores[valid_idx])
    y_v = Float64.(Y[valid_idx])
    if std(ts_v) == 0 || std(y_v) == 0
        return NaN
    end
    return cor(ts_v, y_v)
end

"""
    scale_cor(X::AbstractMatrix, scales::Dict{String, Vector{Int}})

計算多個子量表之間的相關矩陣。
回傳 `(scale_names, cor_matrix)`。
"""
function scale_cor(X::AbstractMatrix, scales::Dict{String, Vector{Int}})
    scale_names = collect(keys(scales))
    n_scales = length(scale_names)
    scores = zeros(Float64, size(X, 1), n_scales)
    for (i, name) in enumerate(scale_names)
        scores[:, i] = sum(coalesce.(X[:, scales[name]], 0), dims=2)
    end
    corr_matrix = pairwise_cor(scores)
    return scale_names, corr_matrix
end

"""
    disattenuated_cor(r_xy::Real, r_xx::Real, r_yy::Real)

計算校正衰減相關 (Disattenuated Correlation)。
用以估計當兩個測量工具皆沒有測量誤差時，兩者之間的真實相關係數。
`r_xy` 為兩變數間的觀察相關係數。
`r_xx` 為變數 X 的信度。
`r_yy` 為變數 Y 的信度。
"""
function disattenuated_cor(r_xy::Real, r_xx::Real, r_yy::Real)
    if r_xx <= 0 || r_yy <= 0
        return NaN
    end
    return r_xy / sqrt(r_xx * r_yy)
end

"""
    pca_loading(X::AbstractMatrix, n_components::Int=1)

透過主成份分析 (PCA) 取得題目在各主成份上的載荷量 (Loadings)。
此函數可用來初步檢查量表的單向度 (Unidimensionality) 或建構效度。
回傳 `(loadings_matrix, explained_variance_ratio)`。
"""
function pca_loading(X::AbstractMatrix, n_components::Int=1)
    R = pairwise_cor(X)
    vals, vecs = eigen(R)
    
    real_vals = real.(vals)
    sorted_idx = sortperm(real_vals, rev=true)
    
    vals_sorted = real_vals[sorted_idx]
    vecs_sorted = real.(vecs[:, sorted_idx])
    
    total_var = sum(vals_sorted)
    k = min(n_components, size(R, 1))
    
    loadings = zeros(Float64, size(R, 1), k)
    var_explained = zeros(Float64, k)
    
    for i in 1:k
        v = vals_sorted[i]
        if v > 0
            loadings[:, i] = vecs_sorted[:, i] .* sqrt(v)
            var_explained[i] = v / total_var
        end
    end
    
    for i in 1:k
        if sum(loadings[:, i]) < 0
            loadings[:, i] = -loadings[:, i]
        end
    end
    
    return loadings, var_explained
end
