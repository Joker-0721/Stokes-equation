# Stokes方程推导

**叙述**:对于两场形式（$\mathbf{u} $和$p$）强制将不可压缩约束与动量方程耦合，导致离散系统需同时满足两种矛盾的数值需求，则可能导致LBB条件不满足，压力场出现高频振荡等問題(如：相同阶数的连续速度$( \mathbf{u}) $ - 压力元$(p)$)。所以对于不可压缩牛顿流体的稳态流动可表述为 $( \mathbf{u}, \boldsymbol{\varepsilon}, \boldsymbol{\sigma}, p )$ 四场问题，四場分別對應

* **速度场**​$( \mathbf{u}) $：流体速度分布
* ​**应变率张量场** $(\boldsymbol{\varepsilon}) $：速度梯度的对称部分（$\boldsymbol{\varepsilon}  = \nabla^s \mathbf{u}$）
* ​**应力张量场**​ $( \boldsymbol{\sigma} )$：流体内部分布的张应力
* ​**压力场**​$(p)$：不可压缩约束下的拉格朗日乘子​

通过引入更多约束或自由度，缓解原始离散系统的过约束问题，并使得每个方程对应不同场变量的离散化需求。这种分离允许​**针对不同变量选择最合适的有限元空间**​，例如：

* 速度场$( \mathbf{u}) $需要高阶连续性以描述流动。
* 应力场 $( \boldsymbol{\sigma} )$和应变率场 $(\boldsymbol{\varepsilon}) $可选用低阶非连续空间。
* 压力场$(p)$需要与速度场匹配以避免数值振荡。
* 不可压缩条件 $\nabla \cdot \mathbf{u} = 0$和本构方程 $\boldsymbol{\sigma} = 2\mu\boldsymbol{\varepsilon}  - p\mathbf{1} $的离散化可分别处理，避免直接耦合导致矩阵病态。

---

针对不同变量选择最合适的有限元空间从而满足稳定性条件（如**Ladyzhenskaya–Babuška–Brezzi条件**​：简称LBB条件或inf-sup条件）

**边界条件**:

1. ​**速度边界条件**​（Dirichlet条件）：
   
   $$
   \mathbf{u} = 0 \quad \text{on } \Omega \quad
   $$
   
   其中 $ \Omega  $ 为速度已知的边界（如固定壁面或入口流速）。
2. ​**应力边界条件**​（Neumann条件）：
   
   $$
   \boldsymbol{\sigma} \cdot \mathbf{n} = \mathbf{t} \quad \text{on } \Gamma
   $$
   
   其中 $ \Gamma $ 为应力已知的边界，$ \mathbf{n} $ 为边界法向量，$ \mathbf{t} $ 为施加的表面力（如自由表面或出口条件）。

---

## 四场问题表述

$$
\begin{cases} 
\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = \mathbf{0} & \text{in } \Omega \quad \text{(平衡方程)} \\
\boldsymbol{\sigma} = 2\mu\boldsymbol{\varepsilon}  - p\mathbf{1} & \text{in } \Omega \quad \text{(本构方程)} \\
\boldsymbol{\varepsilon}  = \nabla^s \mathbf{u} & \text{in } \Omega \quad \text{(协调方程)} \\
\nabla \cdot \mathbf{u} = 0 & \text{in } \Omega \quad \text{(不可压缩条件)}
\end{cases}
$$

- $ \mathbf{u} $: 流体速度场
- $ \mathbf{b} $: 单位体积的分配体荷载
- $ \boldsymbol{\sigma} $: 应力张量
- $ \mu$：膨胀粘性系数
- $ \boldsymbol{\varepsilon} $: 速度梯度的对称部分（应变率张量）
- $ p $: 压力（拉格朗日乘子）/各向同性压力
- $\mathbf{1}：\delta_{ij}$
- $\nabla \cdot ：\frac{D\mathbf{}}{Dx}e_{i}+\frac{D\mathbf{}}{Dy}e_{j}+\frac{D\mathbf{}}{Dz}e_{z}$
- $\nabla^s \mathbf{u} ： \dfrac{1}{2}\Bigl(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\Bigr)$
- $\nabla \mathbf{u}：\Big(\frac{D\mathbf{}}{Dx},\frac{D\mathbf{}}{Dy},\frac{D\mathbf{}}{Dz}\Big)$

---

## 动量守恒方程推导

- 根据牛顿第二定律，对流体微团有

$$
\rho\,\frac{D\mathbf{u}}{Dt} = \nabla \cdot \boldsymbol{\sigma}  + \mathbf{b}.
$$

- 考虑稳态流动（$\frac{D\mathbf{u}}{Dt} = 0$）且低雷诺数，惯性项忽略：

$$
\nabla \cdot \boldsymbol{\sigma}  + \mathbf{b} = \mathbf{0}.
$$

---

## 牛顿流体本构方程推导

本构方程以及协调方程可以根据牛顿黏性定律（一維剪應力）推广到**應力張量** $\tau_{ij}$ 和**應變率張量** $S_{ij}$的關係，以下是推导过程

$$
\begin{cases}
F=\mu A  \frac{u}{y}\\
\tau =\frac{F}{A}
\end{cases}
$$

$$
\tau = \mu \frac{du}{dy} \tag{1}
$$

- ​$\tau$：剪應力
- ​$\mu$：流體動力黏性係數
- ​$\frac{du}{dy}$：速度梯度

在實際流體中，應力不僅包含剪應力，還包含法向應力（正應力）。因此需將標量方程推廣到**應力張量** $\tau_{ij}$ 和**應變率張量** $S_{ij}$的關係，下標 $i= 1,2,3$ 分別對應 $x,y,z$ 方向。

### 應力張量

$$
\tau_{ij} = \begin{bmatrix}
\tau_{xx} & \tau_{xy} & \tau_{xz} \\
\tau_{yx} & \tau_{yy} & \tau_{yz} \\
\tau_{zx} & \tau_{zy} & \tau_{zz}
\end{bmatrix}
$$

### 應變率張量

$$
S_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right)
$$

- ($\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}$）：表示剪切變形(張量對稱部分)

把 **應力張量** $\tau_{ij}$ 和**應變率張量** $S_{ij}$代入(1)，并根据笛卡尔座标中写成一般关系式

$$
\tau_{ij} = 2\mu S_{ij}\tag{2}
$$

$$
\tau_{ij} = \sum_k \sum_\ell \mu_{ijkl} \frac{\partial u_k}{\partial x_\ell}
$$

- $\mu_{ijkl}$：速度梯度張量 $\frac{\partial u_k}{\partial x_\ell}$ 映射到黏性應力張量 $\tau_{ij}$ 的黏度張量。
- 由於指數 $i, j, k, \ell$ 可以在 1~3 之間變化，因此總共有 ​**81 個**​「黏度係數」$\mu_{ijkl}$。

假設黏度張量是各向同性的，則將 81 個係數簡化為 ​**三個獨立參數** $\alpha$、$\beta$、$\gamma$：

$$
\mu_{ijkl} = \alpha \delta_{ij} \delta_{kl} + \beta \delta_{ik} \delta_{jl} + \gamma \delta_{il} \delta_{jk}
$$

- $\delta_{ij}$：Kronecker delta 符號（單位張量）。
- $\alpha, \beta, \gamma$：各向同性假設下的獨立黏度參數。

假設流體黏性與方向無關（各向同性），應力與應變率滿足線性關係：

$$
\tau_{ij} = \alpha S_{kk} \delta_{ij} + 2\mu S_{ij} \tag{3}
$$

- $\alpha$：與體積黏性相關的係數
- $\delta_{ij}$：Kronecker delta（單位張量）
- $S_{kk} = S_{11} + S_{22} + S_{33}$：應變率張量的跡（體積膨脹速率）

流體總應力包含靜壓力 $p$ 和黏性應力：

$$
\sigma_{ij} = -p \delta_{ij} + \tau_{ji} \tag{4}
$$

將 (3) 代入 (4)，得：

$$
\sigma_{ij} = -p \delta_{ij} + \alpha S_{kk} \delta_{ij} + 2\mu S_{ij}
$$

定義 ​**第一黏度** $\mu$（剪切黏度）和 ​**第二黏度** $\mu'$（體積黏度）：

$$
\mu' = \alpha + \frac{2}{3}\mu
$$

最終得到廣義牛頓黏性定律：

$$
\sigma_{ij} = -p \delta_{ij} + 2\mu \left( S_{ij} - \frac{1}{3} S_{kk} \delta_{ij} \right) + \mu' S_{kk} \delta_{ij} \tag{5}
$$

因為體積黏度(第二黏度)只有當流體被迅速壓縮或擴展時，如聲音和衝擊波，體積黏度顯得很重要，在流體動力學問題中通常不是必需的，不可壓縮流體滿足$\nabla \cdot \mathbf{u} = 0$ 因此不含 $\mu'$（體積黏度）。则可得到以下公式：

根据：

$$
\begin{cases}
\mathbf{1} = \delta_{ij}\\
{\varepsilon_m}=\frac{1}{3} S_{kk}\\
\boldsymbol{\varepsilon}={\varepsilon_{ij}} - {\varepsilon_m}\delta_{ij}\\
\end{cases}
$$

可以得到：

$$
\sigma_{ij} =2\mu \boldsymbol{\varepsilon} -p \mathbf{1}\\
$$

- 如果$(i \neq j)$則$\sigma_{ij}$可以寫成$\tau_{ij}$

---

## 应变率定义推导

1.位移場描述相鄰兩點的位移差為：

$$
du_i = \frac{\partial u_i}{\partial x_j} dx_j
$$

其中 $\frac{\partial u_i}{\partial x_j}$ 為位移梯度張量。

2.位移梯度分解將位移梯度分解為對稱（應變張量）和反對稱（轉動張量）部分：

$$
\frac{\partial u_i}{\partial x_j} = e_{ij} + \omega_{ij}
$$

其中：

$$
e_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right), \quad \omega_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial x_j} - \frac{\partial u_j}{\partial x_i} \right)
$$

3.引入速度場將位移場推廣為速度場 $\boldsymbol{v}(\boldsymbol{x}, t)$，其梯度為：

$$
\frac{\partial v_i}{\partial x_j}
$$

4.對稱化速度梯度的對稱部分即為應變率張量：

$$
D_{ij} = \frac{1}{2} \left( \frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i} \right)
$$

- 根据牛津黏性定律可以得到应变率张量：

$$
S_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right)
$$

$$
S_{ij} = \begin{bmatrix}
\frac{1}{2} \left( \frac{\partial u_1}{\partial x} + \frac{\partial u_1}{\partial x} \right) & 
\frac{1}{2} \left( \frac{\partial u_1}{\partial y} + \frac{\partial u_2}{\partial x} \right)&
 \frac{1}{2} \left( \frac{\partial u_1}{\partial z} + \frac{\partial u_3}{\partial x} \right)\\
\frac{1}{2} \left( \frac{\partial u_2}{\partial x} + \frac{\partial u_1}{\partial y} \right) &
\frac{1}{2} \left( \frac{\partial u_2}{\partial y} + \frac{\partial u_2}{\partial y} \right)&
\frac{1}{2} \left( \frac{\partial u_2}{\partial z} + \frac{\partial u_3}{\partial y} \right)\\
\frac{1}{2} \left( \frac{\partial u_3}{\partial x} + \frac{\partial u_1}{\partial z} \right) &
\frac{1}{2} \left( \frac{\partial u_3}{\partial y} + \frac{\partial u_2}{\partial z} \right)&
\frac{1}{2} \left( \frac{\partial u_3}{\partial z} + \frac{\partial u_3}{\partial z} \right)
\end{bmatrix}
$$

若写成矢量形式，有:

$$
{\frac{1}{2}\Bigl(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\Bigr)}={\nabla^s \mathbf{u}}
$$

- 黏性应力仅与对称部分相关：

$$
\boldsymbol{\varepsilon} = \nabla^s \mathbf{u}.
$$

---

## 方程简化：四场→两场

1. ​**代入本构方程到平衡方程**：

$$
\nabla \cdot (2\mu \nabla^s \mathbf{u} - p\mathbf{1}) + \mathbf{b} = 0.
$$

2. ​**展开散度运算**：

2.1 单位矩阵的定义：

$$
\mathbf{1} = \begin{bmatrix} 
1 & 0 & 0 \\ 
0 & 1 & 0 \\ 
0 & 0 & 1 
\end{bmatrix}
$$

2.2 散度运算的本质：

$$
\nabla \cdot \mathbf{1} = \sum_{i=1}^{3} \frac{\partial I_{ij}}{\partial x_i} \mathbf{e}_j
$$

- ​**数学解释**：
  - 对每一列 $ j $，计算其所有行方向的偏导数之和：
    $$
    \frac{\partial 1_{1j}}{\partial x_1} + \frac{\partial 1_{2j}}{\partial x_2} + \frac{\partial 1_{3j}}{\partial x_3}
    $$
  - 由于单位矩阵元素 $ 1_{ij} $ 均为常数（0或1），所有偏导数均为零。
  - ​**结果**：
    $$
    \nabla \cdot \mathbf{1} = \mathbf{0}
    $$

2.3 压力项简化：

压力与单位矩阵的乘积散度运算：

$$
\nabla \cdot (-p\mathbf{1}) = -p(\nabla \cdot \mathbf{1}) - \mathbf{1} \cdot \nabla p
$$

- ​**分步推导**：
  1. 根据散度运算的乘积法则：
     $$
     \nabla \cdot (a\mathbf{A}) = a(\nabla \cdot \mathbf{A}) + \mathbf{A} \cdot \nabla a
     $$
  2. 代入 $ a = -p $ 和 $ \mathbf{A} = \mathbf{I} $:
     $$
     \nabla \cdot (-p\mathbf{1}) = -p(\nabla \cdot \mathbf{1}) - \mathbf{1} \cdot \nabla p
     $$
  3. 利用 $ \nabla \cdot \mathbf{1} = 0 $，简化为：
     $$
     \nabla \cdot (-p\mathbf{1}) = -\nabla p
     $$

代入后得到：

$$
2\mu \, \nabla \cdot  (\nabla^s \mathbf{u}) - \nabla p + \mathbf{b} = 0.
$$

- $ \nabla \cdot \mathbf{A} $: 张量场的散度运算
- $ \mathbf{e}_j $: 笛卡尔坐标系的基向量
- $ \mathbf{0} $: 零向量

3. ​**利用对称梯度特性**
   对称梯度 $\nabla^s \mathbf{u} = \frac{1}{2}(\nabla \mathbf{u} + (\nabla \mathbf{u})^T)$ 的散度为：

$$
\nabla \cdot (\nabla^s \mathbf{u}) = \frac{1}{2} \Delta \mathbf{u}.
$$

代入后系数抵消，得到最终两场方程：

$$
\mu \Delta \mathbf{u} - \nabla p + \mathbf{b} = 0.
$$

**符号标注**

$$
\Delta = \nabla \cdot \nabla
$$

---

## 黏性流体的自由能泛函

**黏性耗散能的泛函形式**

$$
\Pi(\mathbf{u}) = \frac{\mu}{2} \int_{\Omega} \nabla \mathbf{u} : \nabla \mathbf{u} \, d\Omega - \int_{\Omega} \mathbf{b} \cdot \mathbf{u} \, d\Omega.
$$

- 第一項：黏性耗散能（與速度梯度相關）
- 第二項：外力做功（與體積力相關）

---

## 变分原理推导

約束泛函的推導依據基於**黏性流體的自由能泛函**和**不可壓縮約束條件**兩個核心公式的結合：

$$
\begin{cases}
\Pi (\mathbf{u}) = \frac{\mu}{2} \int_\Omega \nabla \mathbf{u} : \nabla \mathbf{u} \, d\Omega - \int_\Omega \mathbf{b} \cdot \mathbf{u} \, d\Omega \\
\nabla \cdot \mathbf{u} = 0 & \text{in } \Omega  \\
\boldsymbol{\sigma} \cdot \boldsymbol{r} = t & \text{on } \Gamma  \\
\end{cases}
$$

通過拉格朗日乘子法，將此約束嵌入能量泛函，構造混合泛函：

$$
L(\mathbf{u}, p) = \Pi(\mathbf{u}) - \int_\Omega p \, (\nabla \cdot \mathbf{u}) \, d\Omega
$$

- $ p $ 既是拉格朗日乘子，也是物理壓力場

將不可壓縮條件以拉格朗日乘子形式加入自由能泛函：

$$
L(\mathbf{u}, p) = \frac{\mu}{2} \int_\Omega \nabla \mathbf{u} : \nabla \mathbf{u} \, d\Omega - \int_\Omega \mathbf{b} \cdot \mathbf{u} \, d\Omega - \int_\Omega p \, (\nabla \cdot \mathbf{u}) \, d\Omega
$$

#### ​**對速度和壓力取變分**

1. ​**速度變分 $\delta \mathbf{u}$**
   對速度場施加微小擾動 $\mathbf{u} \to \mathbf{u} + \delta \mathbf{u}$，要求泛函的極值條件：
   
   $$
   \delta L = \mu \int_\Omega \nabla \delta \mathbf{u} : \nabla \mathbf{u} \, d\Omega - \int_\Omega \delta \mathbf{u} \cdot \mathbf{b} \, d\Omega - \int_\Omega p \, (\nabla \cdot \delta \mathbf{u}) \, d\Omega = 0
   $$
2. ​**壓力變分 $\delta p$**
   對壓力場施加微小擾動 $p \to p + \delta p$，得到約束條件：
   
   $$
   \delta L = -\int_\Omega \delta p \, (\nabla \cdot \mathbf{u}) \, d\Omega = 0
   $$

---

## 鞍点问题与离散化

**弱形式**：

$$
\begin{cases}
A(\mathbf{u},\mathbf{v}) + B(p,\mathbf{v}) = f(\mathbf{v})  \\
B(\mathbf{u},q) = 0
\end{cases}
$$

把流體速度張量用形函數方式表示，得出：

$$
A|_{ij} = 2\mu \int \left[ \nabla N^{u}_i : \nabla N^{u}_j \right] d\Omega
$$

$$
B|_{rj} = \int \left[ N^{p}_r \, \mathrm{div}(N^{u}_j) \right] d\Omega
$$

$$
f|_{i} = \int \left[ b \cdot N^{u}_i \right] d\Omega + \int_{\Gamma} \left[ N^{u}_i \cdot t \right] d\Gamma
$$

​簡化積分

$$
A|_{ij} = \iint 2\mu \nabla u \nabla v \, dxdy
$$

$$
B|_{rj} = \int p \, \mathrm{div}v \, dxdy
$$

$$
f|_{i} = \int b v \, dxdy + \nabla u t \, dxdy
$$


- $\mathbf{N_i}$：形函數向量 ($i=1,2,3 / x,y,z$)
- $\mathbf{t}$：邊界切向量
- $\Gamma$：積分邊界路徑/曲面

**离散化系统**：

$$
\begin{pmatrix}
k^{uu} & k^{up} \\
k^{pu} & k^{pp}
\end{pmatrix}
\begin{pmatrix}
\mathbf{u} \\
\mathbf{p}
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{f} \\
0
\end{pmatrix}.
$$

- $ k^{uu} = A|_{ij}$
- $ k^{up} = B|_{jr}$
- $ k^{pu} = B|_{rj}$
- $ k^{pp} = 0$

---

## 结论

1. ​**四场形式**：

$$
\begin{cases}
\nabla \cdot \,\boldsymbol{\sigma} + \mathbf{b} = 0, \\
\boldsymbol{\sigma} = 2\mu\,\boldsymbol{\varepsilon} - p\,\mathbf{1}, \\
\boldsymbol{\varepsilon} = \nabla^s \mathbf{u}, \\
\nabla \cdot \,\mathbf{u} = 0.
\end{cases}
$$

2. ​**两场形式**：

$$
\nabla \cdot \,(2\mu\,\nabla^s \mathbf{u} - p\,\mathbf{1}) + \mathbf{b} = 0, \quad \nabla \cdot \,\mathbf{u} = 0.
$$

3. ​**数值实现**：通过混合有限元方法求解鞍点矩阵系统。
