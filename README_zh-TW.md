# ClassicalTestTheory.jl

這是一個基於 Julia 語言開發，專為**古典測驗理論 (Classical Test Theory, CTT)** 與**心理計量學 (Psychometrics)** 設計的高效能、現代化套件。

受到 R 語言中廣受歡迎的 `CTT` 與 `psych` 套件啟發，`ClassicalTestTheory.jl` 提供了一套完整且嚴謹的分析工具，涵蓋了項目分析、信度與效度檢定、進階相關係數 (如 Polychoric)，以及探索性因素分析 (EFA)。本套件特別針對真實世界中不完美的資料進行了優化，內建成對刪除 (Pairwise Deletion) 與相關矩陣平滑化 (Matrix Smoothing) 等安全防護機制。

## 安裝方式

您可以在 Julia 的 REPL 中安裝 `ClassicalTestTheory.jl`。按下 `]` 進入 Pkg 模式後輸入：

```julia
pkg> add ClassicalTestTheory
```

*(註：若套件尚未註冊至官方，您可能需要透過 GitHub 網址來安裝。)*

## 核心功能

- **項目分析 (`item_analysis`)**
  - 難度 (Difficulty) / 通過率
  - 點二系列相關 (Point-Biserial Correlation)
  - 二系列相關 (Biserial Correlation)
  - 刪題後 Alpha (Alpha if Item Deleted)
- **誘答力分析 (`distractor_analysis`)**
  - 選項選擇頻率
  - 錯誤選項之點二系列相關
- **信度分析 (Reliability)**
  - **Cronbach's Alpha (`cronbach_alpha`)**
  - **McDonald's Omega (`mcdonald_omega`)** (基於 PAF 萃取，精準度對標 `psych` 套件)
  - **Guttman's Lambda (`guttman_lambda`)** (1 至 6)
- **效度與縮放 (Validity & Scaling)**
  - **效標關聯效度 (`criterion_validity`)**
  - **子量表相關 (`scale_cor`)**
  - **衰減校正 (`disattenuated_cor`)**
  - **標準分數轉換**：`zscore`, `tscore`
- **進階相關係數 (Advanced Correlations)**：純 Julia 實作的二系列 (Biserial)、多分系列 (Polyserial)、四分 (Tetrachoric) 與多分相關 (Polychoric) 係數 (採用穩健的兩步最大概似估計法)。
- **探索性因素分析 (EFA)**：
  - 內建主軸因素法 (Principal Axis Factoring, PAF)，有效避免最大概似法 (ML) 常見的 Heywood Cases 崩潰問題。
  - 無縫整合 `FactorRotations.jl`，輕鬆支援 Varimax、Promax 等數十種直交與斜交轉軸。
  - 蒙地卡羅平行分析 (Parallel Analysis)，客觀建議應保留的因素數量。
  - KMO (Kaiser-Meyer-Olkin) 抽樣適切性檢定。
- **遺漏值與髒資料防護**：全域支援遺漏值 (`missing`)，內建自動化的成對刪除，以及在因素分析前會自動進行的「相關矩陣平滑化 (`cor_smooth`)」，確保非正定矩陣不再讓分析中斷。

## 快速上手

### 1. 基礎項目分析與信度

```julia
using ClassicalTestTheory

# 準備受試者反應矩陣 (列 = 受試者，欄 = 題目)
responses = [
    1 1 0 1;
    1 0 1 1;
    0 0 1 0
]

# 1. 進行項目分析
item_res = item_analysis(responses)
println("各題難度: ", item_res.difficulty)
println("各題點二系列相關: ", item_res.pbis)

# 2. 計算 Cronbach's Alpha 與 McDonald's Omega
alpha = cronbach_alpha(responses)
omega = mcdonald_omega(responses)
println("Cronbach's Alpha: ", alpha)
println("McDonald's Omega: ", omega)

# 3. 探索性因素分析 (EFA) - 使用主軸因素法 (PAF)
efa_res = efa(responses, 1, rotation="varimax")
println("因素載荷量:\n", efa_res.loadings)
```

### 2. 進階相關係數 (Polychoric)

處理李克特量表 (Likert scale) 這類次序尺度資料時，傳統皮爾森相關常會低估真實關聯。您可以輕鬆計算多分相關矩陣：

```julia
# 運算成對刪除的多分相關 (Polychoric) 矩陣
R_poly = polychoric_matrix(X)
```

### 3. 探索性因素分析 (EFA)

不再依靠通靈或陡坡圖盲猜因素數量。使用平行分析搭配嚴謹的 EFA 流程：

```julia
using FactorRotations

# 1. 透過平行分析 (Parallel Analysis) 決定因素數量
n_factors = parallel_analysis(X; cor_type=:polychoric)

# 2. 檢查抽樣適切性 (KMO)
overall_kmo, item_kmo = kmo(pairwise_cor(X))
println("總體 KMO: ", overall_kmo)

# 3. 執行 EFA (內部預設啟用 smooth=true 保護機制)
efa_res = efa(X, n_factors; cor_type=:polychoric, rotation=Varimax())

# 取得轉軸後的載荷矩陣與共同性
println(efa_res.loadings)
println(efa_res.communality)
```

## 型別支援

`ClassicalTestTheory.jl` 中的函數皆接受 `AbstractMatrix` 作為輸入，這意味著它能無縫接軌 `Matrix{Float64}`, `Matrix{Int}`，或者是含有遺漏值的 `Matrix{Union{Missing, Float64}}` 型別矩陣。

## 授權

本專案採用 MIT 授權條款 (MIT License)。
