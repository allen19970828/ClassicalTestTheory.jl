"""
    spearman_brown(r::Real, n::Real)

使用 Spearman-Brown 預言公式計算改變測驗長度後的信度。
`r` 為原始測驗信度，`n` 為測驗長度的倍數（例如 `n = 2` 代表測驗長度加倍）。
"""
function spearman_brown(r::Real, n::Real)
    return (n * r) / (1 + (n - 1) * r)
end

"""
    zscore(X::AbstractVecOrMat{<:Real})

計算標準分數 (Z-scores)，使平均數為 0，標準差為 1。
如果 X 是一個矩陣，將會按列 (column) 進行標準化。
"""
function zscore(X::AbstractVector{<:Real})
    return (X .- mean(X)) ./ std(X)
end

function zscore(X::AbstractMatrix{<:Real})
    return (X .- mean(X, dims=1)) ./ std(X, dims=1)
end

"""
    tscore(X::AbstractVecOrMat{<:Real})

計算 T 分數 (T-scores)，使平均數為 50，標準差為 10。
"""
function tscore(X::AbstractVecOrMat{<:Real})
    return 50 .+ 10 .* zscore(X)
end
