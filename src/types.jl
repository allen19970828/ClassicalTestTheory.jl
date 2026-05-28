"""
    ItemAnalysisResult

儲存項目分析結果的結構體。包含 `difficulty` (難度)、`pbis` (點二系列相關)、`bis` (二系列相關) 以及 `alpha_if_deleted`。
"""
struct ItemAnalysisResult
    difficulty::Vector{Float64}
    pbis::Vector{Float64}
    bis::Vector{Float64}
    alpha_if_deleted::Vector{Float64}
end

"""
    DistractorAnalysisResult

儲存選項分析結果的結構體。包含 `frequencies` (選項頻率) 與 `pbis` (選項點二系列相關)。
"""
struct DistractorAnalysisResult
    frequencies::Dict{String, Float64}
    pbis::Dict{String, Float64}
end


