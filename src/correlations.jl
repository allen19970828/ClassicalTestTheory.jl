using Distributions
using QuadGK
using Optim

"""
    biserial_cor(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

計算二系列相關 (Biserial Correlation)。
`x` 為二元變數 (例如 0, 1)，`y` 為連續變數 (例如總分)。
"""
function biserial_cor(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("x 與 y 長度必須相同"))
    end
    
    valid_idx = .!(ismissing.(x) .| ismissing.(y))
    x_v = x[valid_idx]
    y_v = Float64.(y[valid_idx])
    
    unique_x = unique(x_v)
    if length(unique_x) != 2
        throw(ArgumentError("x 必須為二元變數 (僅包含兩種數值)"))
    end
    
    val1, val0 = maximum(unique_x), minimum(unique_x)
    idx1 = x_v .== val1
    idx0 = x_v .== val0
    
    mean1 = mean(y_v[idx1])
    mean0 = mean(y_v[idx0])
    p = sum(idx1) / length(x_v)
    q = 1.0 - p
    
    sy = std(y_v)
    if sy == 0
        return NaN
    end
    
    # 標準常態分配中面積為 q 對應的 z 值
    z = quantile(Normal(), q)
    # 在 z 處的常態分配高度
    u = pdf(Normal(), z)
    
    if u == 0
        return NaN
    end
    
    rb = ((mean1 - mean0) / sy) * (p * q / u)
    return max(min(rb, 1.0), -1.0)
end

"""
    polyserial_cor(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

計算多分系列相關 (Polyserial Correlation) 的兩步估計值。
`x` 為次序變數 (ordinal)，`y` 為連續變數。
"""
function polyserial_cor(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("x 與 y 長度必須相同"))
    end
    
    valid_idx = .!(ismissing.(x) .| ismissing.(y))
    x_v = x[valid_idx]
    y_v = Float64.(y[valid_idx])
    
    n = length(x_v)
    unique_x = sort(unique(x_v))
    k = length(unique_x)
    
    if k < 2
        return NaN
    elseif k == 2
        return biserial_cor(x, y)
    end
    
    # 計算各類別累積比例與閾值 (thresholds)
    freqs = [sum(x_v .== val) for val in unique_x]
    props = freqs ./ n
    cum_props = cumsum(props)[1:end-1]
    
    taus = quantile.(Normal(), cum_props)
    
    # 計算各閾值處的常態分配高度，然後加總
    sum_phi = sum(pdf.(Normal(), taus))
    
    r_xy = cor(x_v, y_v)
    sx = std(x_v)
    
    rps = r_xy * sqrt((n - 1) / n) * (sx / sum_phi)
    
    # 限制在 [-1, 1] 之間
    rps = max(min(rps, 1.0), -1.0)
    return rps
end

"""
    bivnor(h::Real, k::Real, r::Real)

計算標準二元常態分配在 (-Inf, h) 與 (-Inf, k) 區間內的累積機率。
`r` 為相關係數。
"""
function bivnor(h::Real, k::Real, r::Real)
    if h == -Inf || k == -Inf
        return 0.0
    elseif h == Inf
        return cdf(Normal(), k)
    elseif k == Inf
        return cdf(Normal(), h)
    end
    
    if r == 1.0
        return cdf(Normal(), min(h, k))
    elseif r == -1.0
        return max(0.0, cdf(Normal(), h) + cdf(Normal(), k) - 1.0)
    end
    
    f(x) = pdf(Normal(), x) * cdf(Normal(), (k - r * x) / sqrt(1.0 - r^2))
    res, _ = quadgk(f, -Inf, h)
    return max(0.0, min(1.0, res))
end

"""
    polychoric_cor(x::AbstractVector, y::AbstractVector)

使用兩步最大概似估計法 (Two-step MLE) 計算兩個次序變數之間的多分相關 (Polychoric Correlation)。
"""
function polychoric_cor(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        throw(ArgumentError("x 與 y 長度必須相同"))
    end
    
    valid_idx = .!(ismissing.(x) .| ismissing.(y))
    x_v = x[valid_idx]
    y_v = y[valid_idx]
    n = length(x_v)
    
    unique_x = sort(unique(x_v))
    unique_y = sort(unique(y_v))
    r_levels = length(unique_x)
    c_levels = length(unique_y)
    
    if r_levels < 2 || c_levels < 2
        return NaN
    end
    
    # 建立列聯表
    freqs = zeros(Int, r_levels, c_levels)
    for i in 1:n
        r_idx = findfirst(==(x_v[i]), unique_x)
        c_idx = findfirst(==(y_v[i]), unique_y)
        freqs[r_idx, c_idx] += 1
    end
    
    # 邊際機率與閾值
    row_props = sum(freqs, dims=2) ./ n
    col_props = sum(freqs, dims=1) ./ n
    
    tau_x = [-Inf; quantile.(Normal(), cumsum(vec(row_props))[1:end-1]); Inf]
    tau_y = [-Inf; quantile.(Normal(), cumsum(vec(col_props))[1:end-1]); Inf]
    
    # 對數概似函數 (Log-Likelihood)
    function neg_loglik(rho)
        ll = 0.0
        for i in 1:r_levels
            for j in 1:c_levels
                if freqs[i, j] > 0
                    p = bivnor(tau_x[i+1], tau_y[j+1], rho) - 
                        bivnor(tau_x[i], tau_y[j+1], rho) - 
                        bivnor(tau_x[i+1], tau_y[j], rho) + 
                        bivnor(tau_x[i], tau_y[j], rho)
                    if p <= 0
                        ll -= freqs[i, j] * -1e10 # penalty for invalid prob
                    else
                        ll += freqs[i, j] * log(p)
                    end
                end
            end
        end
        return -ll
    end
    
    res = optimize(neg_loglik, -0.999, 0.999, Brent())
    return Optim.minimizer(res)
end

"""
    tetrachoric_cor(x::AbstractVector, y::AbstractVector)

計算四分相關 (Tetrachoric Correlation)，這是多分相關在二元變數上的特例。
"""
function tetrachoric_cor(x::AbstractVector, y::AbstractVector)
    return polychoric_cor(x, y)
end

"""
    polychoric_matrix(X::AbstractMatrix)

計算矩陣內變數兩兩之間的多分相關，回傳相關矩陣。
"""
function polychoric_matrix(X::AbstractMatrix)
    J = size(X, 2)
    R = zeros(Float64, J, J)
    for i in 1:J
        for j in i:J
            if i == j
                R[i, i] = 1.0
            else
                r_val = polychoric_cor(X[:, i], X[:, j])
                R[i, j] = R[j, i] = r_val
            end
        end
    end
    return R
end
