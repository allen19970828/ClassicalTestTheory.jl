module CTT

using Statistics
using LinearAlgebra

include("utils.jl")
include("types.jl")
include("reliability.jl")
include("analysis.jl")
include("validity.jl")
include("scoring.jl")
include("correlations.jl")
include("factor_analysis.jl")

export ItemAnalysisResult, DistractorAnalysisResult
export score_items, item_analysis, distractor_analysis
export cronbach_alpha, guttman_lambda, mcdonald_omega, split_half
export criterion_validity, scale_cor, disattenuated_cor, pca_loading
export spearman_brown, zscore, tscore
export biserial_cor, polyserial_cor, polychoric_cor, tetrachoric_cor, polychoric_matrix
export paf, efa, EFAResult, parallel_analysis
export cor_smooth, kmo

end
