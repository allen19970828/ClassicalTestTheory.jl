using Test
using CTT
using Statistics

@testset "CTT.jl" begin
    @testset "Reliability" begin
        X = [1 1 0;
             1 0 1;
             0 1 1;
             1 1 1;
             0 0 0]
             
        # Manually calculated alpha
        # Total var = 1.2
        # Sum of item vars = 0.9
        # Alpha = (3/2) * (1 - 0.9/1.2) = 1.5 * 0.25 = 0.375
        alpha = cronbach_alpha(X)
        @test isapprox(alpha, 0.375, atol=1e-5)
    end
    
    @testset "Item Analysis" begin
        X = [1 1 0;
             1 0 1;
             0 1 1;
             1 1 1;
             0 0 0]
             
        res = item_analysis(X)
        
        # Difficulties
        @test isapprox(res.difficulty[1], 0.6, atol=1e-5)
        @test isapprox(res.difficulty[2], 0.6, atol=1e-5)
        @test isapprox(res.difficulty[3], 0.6, atol=1e-5)
        
        # Test alpha if deleted
        # If we delete item 1, X_sub = [1 0; 0 1; 1 1; 1 1; 0 0]
        # Total scores sub = [1, 1, 2, 2, 0] -> mean = 1.2, var = 0.7
        # Item vars = 0.3, 0.3 -> sum = 0.6
        # Alpha if item 1 deleted = (2/1) * (1 - 0.6/0.7) = 2 * (1 - 0.857142) = 0.285714
        @test isapprox(res.alpha_if_deleted[1], 0.285714, atol=1e-5)
        
        @test length(res.bis) == 3
    end
    
    @testset "Subscales" begin
        X = [1 1 0 2 3;
             1 0 1 1 2;
             0 1 1 2 2;
             1 1 1 3 4;
             0 0 0 1 1]
        
        scales = Dict("Scale1" => [1, 2, 3], "Scale2" => [4, 5])
        
        # Test reliability
        alphas = cronbach_alpha(X, scales)
        @test haskey(alphas, "Scale1")
        @test haskey(alphas, "Scale2")
        @test isapprox(alphas["Scale1"], 0.375, atol=1e-5)
        
        # Test item analysis
        res = item_analysis(X, scales)
        @test haskey(res, "Scale1")
        @test haskey(res, "Scale2")
        @test length(res["Scale1"].difficulty) == 3
        @test length(res["Scale1"].bis) == 3
        @test length(res["Scale2"].difficulty) == 2
        @test length(res["Scale2"].bis) == 2
    end
    
    @testset "Validity" begin
        X = [1 1 0;
             1 0 1;
             0 1 1;
             1 1 1;
             0 0 0]
        Y = [2.5, 2.0, 1.5, 3.0, 0.5] # criterion
        
        val = criterion_validity(X, Y)
        @test -1.0 <= val <= 1.0
        
        scales = Dict("S1" => [1, 2], "S2" => [2, 3])
        names, corr = scale_cor(X, scales)
        @test length(names) == 2
        @test size(corr) == (2, 2)
        @test isapprox(corr[1, 1], 1.0)
    end
    @testset "Scoring" begin
        # spearman_brown
        r = 0.5
        @test isapprox(spearman_brown(r, 2), 0.66666, atol=1e-4)
        
        # zscore
        X_z = [1 2 3; 3 4 5; 5 6 7] # means: 3, 4, 5. stds: 2, 2, 2
        Z = zscore(X_z)
        @test isapprox(mean(Z, dims=1)[1], 0.0, atol=1e-7)
        @test isapprox(std(Z, dims=1)[1], 1.0, atol=1e-7)
        
        # tscore
        T = tscore(X_z)
        @test isapprox(mean(T, dims=1)[1], 50.0, atol=1e-7)
        @test isapprox(std(T, dims=1)[1], 10.0, atol=1e-7)
    end
    
    @testset "Phase 1 Extensions" begin
        # split_half
        X_sh = [1 1 0 1; 1 0 1 0; 0 1 1 1; 1 1 1 1; 0 0 0 0]
        sh_rel = split_half(X_sh)
        @test sh_rel >= -1.0 && sh_rel <= 1.0
        
        # distractor_analysis
        responses = ["A" "B" "C"; "A" "C" "C"; "B" "B" "A"; "A" "B" "C"]
        key = ["A", "B", "C"]
        res = distractor_analysis(responses, key)
        @test length(res) == 3
        @test res[1].frequencies["A"] == 0.75
        @test haskey(res[1].pbis, "B")
        
        # disattenuated_cor
        r_xy = 0.5
        r_xx = 0.8
        r_yy = 0.9
        @test isapprox(disattenuated_cor(r_xy, r_xx, r_yy), 0.58925, atol=1e-4)
    end
    
    @testset "Phase 2: Reliability & Correlations" begin
        # guttman_lambda
        X = [1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0]
        gl = guttman_lambda(X)
        @test haskey(gl, "L1")
        @test haskey(gl, "L6")
        @test isapprox(gl["L3"], cronbach_alpha(X), atol=1e-5)
        
        # mcdonald_omega
        # For a simple matrix, omega should be computable and > 0
        omega = mcdonald_omega(X)
        @test omega >= 0.0 && omega <= 1.0
        
        # biserial_cor
        x_bin = [1, 1, 0, 0, 1, 0, 1, 1, 0, 0]
        y_cont = [5.5, 6.2, 3.1, 4.0, 5.8, 2.9, 7.1, 6.0, 4.2, 3.8]
        rb = biserial_cor(x_bin, y_cont)
        @test rb > 0.0 && rb <= 1.0
        
        # polyserial_cor
        x_poly = [3, 2, 1, 1, 3, 1, 2, 3, 2, 1]
        rps = polyserial_cor(x_poly, y_cont)
        @test rps >= -1.0 && rps <= 1.0
        
        # polychoric_cor and tetrachoric_cor
        y_poly = [1, 2, 3, 3, 1, 3, 2, 1, 2, 3]
        rpc = polychoric_cor(x_poly, y_poly)
        @test rpc >= -1.0 && rpc <= 1.0
        
        x_bin = [1, 0, 1, 0, 1, 1, 0, 0, 1, 0]
        y_bin = [1, 1, 0, 0, 1, 0, 1, 0, 1, 0]
        rtc = tetrachoric_cor(x_bin, y_bin)
        @test rtc >= -1.0 && rtc <= 1.0
        
        # polychoric_matrix
        X_poly = [x_poly y_poly x_bin]
        R_poly = polychoric_matrix(X_poly)
        @test size(R_poly) == (3, 3)
        @test isapprox(R_poly[1, 1], 1.0, atol=1e-5)
    end
    
    @testset "Phase 3: Construct Validity" begin
        X = [1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0]
        loadings, var_explained = pca_loading(X, 2)
        @test size(loadings, 1) == 3
        @test size(loadings, 2) <= 2
        @test length(var_explained) <= 2
        @test sum(var_explained) > 0.0
        
        # Test EFA with FactorRotations
        X_efa = rand(100, 5) # 100 samples, 5 items
        res_efa = efa(X_efa, 2; cor_type=:pearson, rotation=nothing)
        @test size(res_efa.loadings) == (5, 2)
        @test length(res_efa.communality) == 5
        
        # Test EFA with Varimax
        using FactorRotations
        res_varimax = efa(X_efa, 2; rotation=Varimax())
        @test size(res_varimax.loadings) == (5, 2)
        @test size(res_varimax.phi) == (2, 2)
        
        # Test EFA with Polychoric
        X_poly_efa = rand(1:4, 100, 3)
        res_poly = efa(X_poly_efa, 1; cor_type=:polychoric, rotation=nothing)
        @test size(res_poly.loadings) == (3, 1)
    end
    
    @testset "Missing Data Handling" begin
        # Create a matrix with some missing values
        X_miss = [1 1 missing; 1 0 1; 0 1 1; 1 missing 1; 0 0 0]
        
        # Test item analyses
        res = item_analysis(X_miss)
        @test isapprox(res.difficulty[3], 0.75) # (1+1+1+0)/4 = 3/4 = 0.75
        @test !isnan(res.pbis[1])
        
        # Test reliability
        alpha = cronbach_alpha(X_miss)
        @test alpha > -1.0 && alpha < 1.0
        
        sh = split_half(X_miss)
        @test sh > -1.0 && sh < 1.0
        
        gl = guttman_lambda(X_miss)
        @test haskey(gl, "L1")
        
        omega = mcdonald_omega(X_miss)
        @test omega >= 0.0 && omega <= 1.0
        
        # Test construct validity
        loadings, var_exp = pca_loading(X_miss)
        @test size(loadings, 1) == 3
        
        # Test correlations
        y_cont = [5.0, 4.0, 3.5, 4.5, 2.0]
        cv = criterion_validity(X_miss, y_cont)
        @test cv >= -1.0 && cv <= 1.0
        
        x_bin = [1, missing, 0, 1, 0]
        rb = biserial_cor(x_bin, y_cont)
        @test rb >= -1.0 && rb <= 1.0
    end
end

