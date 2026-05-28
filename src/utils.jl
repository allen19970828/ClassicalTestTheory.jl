using Statistics

"""
    pairwise_cov(X::AbstractMatrix)

計算成對刪除 (Pairwise Deletion) 後的共變異數矩陣。
會忽略含有 `missing` 的配對資料。
"""
function pairwise_cov(X::AbstractMatrix)
    n, p = size(X)
    C = zeros(Float64, p, p)
    for i in 1:p
        for j in i:p
            valid_rows = .!(ismissing.(X[:, i]) .| ismissing.(X[:, j]))
            if sum(valid_rows) > 1
                c_val = cov(X[valid_rows, i], X[valid_rows, j])
                C[i, j] = C[j, i] = c_val
            else
                C[i, j] = C[j, i] = NaN
            end
        end
    end
    return C
end

"""
    pairwise_cor(X::AbstractMatrix)

計算成對刪除 (Pairwise Deletion) 後的相關係數矩陣。
會忽略含有 `missing` 的配對資料。
"""
function pairwise_cor(X::AbstractMatrix)
    n, p = size(X)
    R = zeros(Float64, p, p)
    for i in 1:p
        for j in i:p
            if i == j
                R[i, i] = 1.0
            else
                valid_rows = .!(ismissing.(X[:, i]) .| ismissing.(X[:, j]))
                if sum(valid_rows) > 1
                    r_val = cor(X[valid_rows, i], X[valid_rows, j])
                    R[i, j] = R[j, i] = r_val
                else
                    R[i, j] = R[j, i] = NaN
                end
            end
        end
    end
    return R
end

"""
    cor_smooth(R::AbstractMatrix; tol=1e-12)

平滑化相關矩陣，確保矩陣為正定 (Positive Definite)。
此函數對標 R 的 `psych::cor.smooth`。當成對刪除 (Pairwise deletion) 
產生非正定矩陣時，此函數會將非正特徵值抹平，避免後續因素分析出錯。
"""
function cor_smooth(R::AbstractMatrix; tol=1e-12)
    p = size(R, 1)
    vals, vecs = eigen(Symmetric(R))
    
    if all(vals .>= tol)
        return R
    end
    
    # 將小於 tol 的特徵值替換為微小正數
    vals_new = copy(vals)
    vals_new[vals .< tol] .= 100 * tol
    
    # 重新縮放特徵值，使其總和維持為 p (對角線總和)
    vals_new .= vals_new .* (p / sum(vals_new))
    
    # 重建矩陣
    R_new = vecs * diagm(vals_new) * vecs'
    
    # 將共變異數矩陣轉換回相關矩陣 (cov2cor)
    D = sqrt.(abs.(diag(R_new)))
    R_smooth = R_new ./ (D * D')
    
    # 確保對角線精確為 1.0
    for i in 1:p
        R_smooth[i, i] = 1.0
    end
    
    return R_smooth
end

"""
    kmo(R::AbstractMatrix)

計算 Kaiser-Meyer-Olkin (KMO) 抽樣適切性檢定指標。
用來評估相關矩陣是否適合進行因素分析。
回傳 (總體 KMO 值, 各變數 KMO 值)。
"""
function kmo(R::AbstractMatrix)
    p = size(R, 1)
    
    # 計算反矩陣
    Q = try
        inv(R)
    catch
        pinv(R)
    end
    
    # Anti-image 相關矩陣的平方
    D = sqrt.(abs.(diag(Q)))
    Q_norm = Q ./ (D * D')
    q_sq = Q_norm.^2
    
    r_sq = R.^2
    
    # 將對角線歸零
    for i in 1:p
        q_sq[i, i] = 0.0
        r_sq[i, i] = 0.0
    end
    
    sum_r_sq = sum(r_sq)
    sum_q_sq = sum(q_sq)
    
    if sum_r_sq + sum_q_sq == 0
        return NaN, fill(NaN, p)
    end
    
    overall_kmo = sum_r_sq / (sum_r_sq + sum_q_sq)
    
    col_r_sq = vec(sum(r_sq, dims=1))
    col_q_sq = vec(sum(q_sq, dims=1))
    
    item_kmo = col_r_sq ./ (col_r_sq .+ col_q_sq)
    
    return overall_kmo, item_kmo
end
