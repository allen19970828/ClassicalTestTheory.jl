"""
    item_analysis(X::AbstractMatrix{<:Real}; corrected=true)

進行項目分析，回傳 `ItemAnalysisResult`。
若 `corrected` 為 true，則計算點二系列相關與二系列相關時會扣除該題分數 (Item-rest correlation)。
"""
function item_analysis(X::AbstractMatrix; corrected=true)
    n_subj, J = size(X)
    
    # 計算單題難度 (通過率或平均數) 與總分
    difficulties = zeros(Float64, J)
    total_scores = zeros(Float64, n_subj)
    
    for i in 1:n_subj
        # 遺漏值當作 0 分 (或不計分)
        total_scores[i] = sum(skipmissing(X[i, :]))
    end
    
    for j in 1:J
        valid_items = collect(skipmissing(X[:, j]))
        difficulties[j] = isempty(valid_items) ? NaN : mean(valid_items)
    end
    pbis = zeros(Float64, J)
    bis = zeros(Float64, J)
    alphas_if_deleted = zeros(Float64, J)
    
    for j in 1:J
        item_scores = X[:, j]
        valid_idx = .!ismissing.(item_scores)
        
        if sum(valid_idx) < 2
            pbis[j] = NaN
            bis[j] = NaN
        else
            item_scores_v = Float64.(item_scores[valid_idx])
            total_scores_v = Float64.(total_scores[valid_idx])
            
            if std(item_scores_v) == 0 || std(total_scores_v) == 0
                pbis[j] = NaN
                bis[j] = NaN
            else
                if corrected
                    total_scores_sub = total_scores_v .- item_scores_v
                    if std(total_scores_sub) == 0
                        pbis[j] = NaN
                        bis[j] = NaN
                    else
                        pbis[j] = cor(item_scores_v, total_scores_sub)
                        bis[j] = polyserial_cor(item_scores_v, total_scores_sub)
                    end
                else
                    pbis[j] = cor(item_scores_v, total_scores_v)
                    bis[j] = polyserial_cor(item_scores_v, total_scores_v)
                end
            end
        end
        
        # 刪除該題後的 Cronbach's Alpha
        X_sub = X[:, filter(idx -> idx != j, 1:J)]
        alphas_if_deleted[j] = cronbach_alpha(X_sub)
    end
    
    return ItemAnalysisResult(difficulties, pbis, bis, alphas_if_deleted)
end

"""
    item_analysis(X::AbstractMatrix{<:Real}, scales::Dict{String, Vector{Int}}; corrected=true)

針對多個子量表進行項目分析，回傳 Dict。
"""
function item_analysis(X::AbstractMatrix, scales::Dict{String, Vector{Int}}; corrected=true)
    results = Dict{String, ItemAnalysisResult}()
    for (scale_name, items) in scales
        results[scale_name] = item_analysis(X[:, items]; corrected=corrected)
    end
    return results
end

"""
    distractor_analysis(responses::AbstractMatrix, key::AbstractVector)

進行誘答力分析 (Distractor Analysis)。
`responses` 為 N x J 的矩陣，代表 N 個受試者在 J 個題目上的作答反應 (例如選項 A, B, C, D)。
`key` 為長度 J 的向量，代表正確答案。
回傳長度為 J 的 `DistractorAnalysisResult` 陣列。
"""
function distractor_analysis(responses::AbstractMatrix, key::AbstractVector)
    N, J = size(responses)
    if length(key) != J
        throw(ArgumentError("長度不一致：key 的長度必須與 responses 的欄數相同"))
    end
    
    # 計算總分：先轉換 responses 為二元計分 (答對為 1，答錯為 0，漏答為 0)
    binary_scores = zeros(Int, N, J)
    for j in 1:J
        binary_scores[:, j] .= coalesce.(responses[:, j] .== key[j], false)
    end
    total_scores = vec(sum(binary_scores, dims=2))
    
    results = Vector{DistractorAnalysisResult}(undef, J)
    for j in 1:J
        item_responses = string.(responses[:, j])
        unique_options = unique(item_responses)
        
        freqs = Dict{String, Float64}()
        pbis_dict = Dict{String, Float64}()
        
        for opt in unique_options
            # 選擇此選項的人 (1 vs 0)
            chosen = (item_responses .== opt)
            freqs[opt] = sum(chosen) / N
            
            # 計算點二系列相關
            if std(chosen) == 0 || std(total_scores) == 0
                pbis_dict[opt] = NaN
            else
                pbis_dict[opt] = cor(chosen, total_scores)
            end
        end
        results[j] = DistractorAnalysisResult(freqs, pbis_dict)
    end
    
    return results
end
