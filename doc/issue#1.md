# Stokes方程推导

**叙述**:对于两场形式（$\mathbf{u} $和$p$）强制将不可压缩约束与动量方程耦合，导致离散系统需同时满足两种矛盾的数值需求，则可能导致LBB条件不满足，压力场出现高频振荡等問題。所以对于不可压缩牛顿流体的稳态流动可表述为 $( \mathbf{u}, \boldsymbol{\varepsilon}, \boldsymbol{\sigma}, p )$ ，分別對應

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
   \mathbf{u} = \mathbf{0} \quad \text{on } \Omega \quad
   $$
   
   其中 $ \Omega  $ 为速度已知的边界（如固定壁面或入口流速）。
2. ​**应力边界条件**​（Neumann条件）：
   
   $$
   \boldsymbol{\sigma} \cdot \mathbf{n} = \mathbf{t} \quad \text{on } \Gamma
   $$
   
   其中 $ \Gamma $ 为应力已知的边界，$ \mathbf{n} $ 为边界法向量，$ \mathbf{t} $ 为施加的表面力（如自由表面或出口条件）。

---

## 问题表述

$$
\begin{cases} 
\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = \mathbf{0} & \text{in } \Omega \quad \text{(平衡方程)} \\
\boldsymbol{\sigma} = 2\mu\boldsymbol{\varepsilon}  - p\mathbf{1} & \text{in } \Omega \quad \text{(本构方程)} \\
\boldsymbol{\varepsilon}  = \nabla^s \mathbf{u} & \text{in } \Omega \quad \text{(幾何方程)} \\
\nabla \cdot \mathbf{u} = \mathbf{0} & \text{in } \Omega \quad \text{(不可压缩条件)}
\end{cases}
$$

- $ \mathbf{u} $: 流体速度场
- $ \mathbf{b} $: 单位体积的分配体荷载
- $ \boldsymbol{\sigma} $: 应力张量
- $ \mu$：膨胀粘性系数
- $ \boldsymbol{\varepsilon} $: 速度梯度的对称部分（应变率张量）
- $ p $: 压力（拉格朗日乘子）/各向同性压力
- $\mathbf{1}：\delta_{ij}$
- $\nabla ：\frac{D\mathbf{}}{Dx_i}e_{i}$
- $\nabla^s \mathbf{u} ： \dfrac{1}{2}(u_{i,j} +u_{j,i})$

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

本构方程以及协调方程可以根据牛顿黏性定律（一維剪應力）推广到**應力張量** $\sigma_{ij}$ 和**應變率張量** $\varepsilon_{ij}$的關係，以下是推导过程

$$
F=\mu A  \frac{u}{y}\\
\sigma_{ij} =2\mu \varepsilon_{ij} - p \delta_{ij}
$$

下標 $i= 1,2,3$ 分別對應 $x,y,z$ 方向。

### 應力張量

$$
\sigma_{ij} = \begin{bmatrix}
\sigma_{xx} & \sigma_{xy} & \sigma_{xz} \\
\sigma_{yx} & \sigma_{yy} & \sigma_{yz} \\
\sigma_{zx} & \sigma_{zy} & \sigma_{zz}
\end{bmatrix}
$$

### 應變率張量

$$
\varepsilon_{ij} = \frac{1}{2} \left(u_{i,j}+u_{j,i} \right)
$$

- $\frac{1}{2} \left(u_{i,j}+u_{j,i} \right)$：表示剪切變形(張量對稱部分)

---

## 应变率定义推导

1.位移場描述相鄰兩點的位移差為：

$$
du_i =  u_{i,j} dx_j
$$

其中 $u_{i,j}$ 為位移梯度張量。

2.位移梯度分解將位移梯度分解為對稱（應變張量）和反對稱（轉動張量）部分：

$$
\frac{\partial u_i}{\partial x_j} = \varepsilon_{ij} + \omega_{ij}
$$

其中：

$$
\varepsilon_{ij} = \frac{1}{2} \left( u_{i,j} + u_{j,i} \right), \quad \omega_{ij} = \frac{1}{2} \left(u_{i,j} - u_{j,i}  \right)
$$

$$
\varepsilon_{ij} = \begin{bmatrix}
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
{\frac{1}{2} (u_{i,j} + u_{j,i})}={\nabla^s \mathbf{u}}
$$

---

## 方程简化

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

## 离散化

**弱形式**：

$$
\begin{cases}
 \delta L =2 \mu \int_\Omega \nabla \delta \mathbf{u} : \nabla \mathbf{u} \, d\Omega \\
\delta L = -\int_\Omega \delta p \, (\nabla \cdot \mathbf{u}) \, d\Omega 
\end{cases}
$$

**速度离散**

$$
u^h_{i}(\boldsymbol{x}) = \sum_{I=1}^{n_u} N_{I}(\boldsymbol{x}) d_{iI}
$$

**应变率张量场**

$$
\begin{split}
\varepsilon^h_{ij} 
&= \frac{1}{2} \left( u^h_{i,j} + u^h_{j,i} \right) \\
&= \sum_{I=1}^{n_p} \frac{1}{2} \left( N_{I,j} d_{iI} + N_{I,i} d_{jI} \right) \\
&= \sum_{I=1}^{n_p} \frac{1}{2} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) d_{kI}
\end{split}
$$

**代入弱形式（动量方程）**

$$
\begin{split}
\int_{\Omega} 2\mu \, \delta \varepsilon^h_{ij} \varepsilon^h_{ij} \, \mathrm{d}\Omega 
&= \sum_{I,J}^{n_p} \int_{\Omega} 2\mu \cdot \frac{1}{4} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) \delta d_{kI} \left( N_{J,j} \delta_{il} + N_{J,i} \delta_{jl} \right) d_{lJ} \, \mathrm{d}\Omega \\
&= \sum_{I,J}^{n_p} \delta d_{kI} \int_{\Omega} 2\mu \cdot \frac{1}{4} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) \left( N_{J,j} \delta_{il} + N_{J,i} \delta_{jl} \right) \mathrm{d}\Omega \, d_{lJ}
\end{split}
$$

**壓力離散**

$$
\begin{gathered}
\delta p^h(\boldsymbol{x}) = \sum_{J=1}^{n_p} \Psi_J(\boldsymbol{x}) \delta p_J \\
u^h_i(\boldsymbol{x}) = \sum_{I=1}^{n_u} N_I(\boldsymbol{x}) d_{iI}
\end{gathered}
$$

**散度项表达式：**

$$
\nabla \cdot \mathbf{u}^h = u^h_{k,k} = \sum_{I=1}^{n_u} N_{I,k} d_{kI}
$$

**代入弱形式（连续性方程）**

$$
\int_\Omega \delta p^h (\nabla \cdot \mathbf{u}^h) \, \mathrm{d}\Omega \mathbf{K}^{up} = \left( \int_\Omega \Psi_J N_{I,k} \, \mathrm{d}\Omega \right)
= \sum_{I,J=1}^{n_p} \delta p_J \left( \int_\Omega N_J N_{I,k} \, \mathrm{d}\Omega \right) d_{kI}
$$

**离散化系统**：

$$
\mathbf{K}_{IJkl} = \int_{\Omega} 2\mu \cdot \frac{1}{4} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) \left( N_{J,j} \delta_{il} + N_{J,i} \delta_{jl} \right)
$$

$$
\mathbf{K}^{up} = \left( \int_\Omega \Psi_J N_{I,k} \, \mathrm{d}\Omega \right)
$$

$$
\mathbf{K}^{pp} = \mathbf{0}
$$

#### Component-wise Expansion:

**Case $(k,l)=(1,1)$**:

$$
\mathbf{K}_{IJ11} = \mathbf{K}[2I-1,2J-1] = \mu \int_{\Omega} \left( 2 N_{I,x} N_{J,x} + N_{I,y} N_{J,y} \right) d\Omega
$$

. ​**Case $(k,l)=(2,1)$**:

$$
\mathbf{K}_{IJ21} = \mathbf{K}[2I,2J-1] = \mu \int_{\Omega} N_{I,x} N_{J,y} \, d\Omega
$$

. ​**Case $(k,l)=(1,2)$**:

$$
\mathbf{K}_{IJ12} = \mathbf{K}[2I-1,2J] = \mu \int_{\Omega} N_{I,y} N_{J,x} \, d\Omega
$$

. ​**Case $(k,l)=(2,2)$**:

$$
\mathbf{K}_{IJ22} = \mathbf{K}[2I,2J]= \mu \int_{\Omega} \left( N_{I,x} N_{J,x} + 2 N_{I,y} N_{J,y} \right) d\Omega
$$

$$
\begin{bmatrix}
\delta d_{1I}\\
\delta d_{2I}
\end{bmatrix}

\begin{bmatrix}
k[2I-1,2J-1] & k[2I,2J-1] \\
k[2I-1,2J] & k[2I,2J]
\end{bmatrix}

\begin{bmatrix}
d_{1J}\\
d_{2J}
\end{bmatrix}
$$

. ​**Case $(k)=(1)$**:

$$
\mathbf{K}^{up} = \left( \int_\Omega N_{I,x} \Psi_J \, \mathrm{d}\Omega \right)
$$

. ​**Case $(k)=(2)$**:

$$
\mathbf{K}^{up} = \left( \int_\Omega N_{I,y}\Psi_J  \, \mathrm{d}\Omega \right)
$$

$$
\begin{bmatrix}
\delta d_{1I}\\
\delta d_{2I}
\end{bmatrix}

\begin{bmatrix}
k[2I-1,2J-1] \\
k[2I-1,2J]
\end{bmatrix}

\begin{bmatrix}
d_{1J}\\
d_{2J}
\end{bmatrix}
$$

$$
\begin{bmatrix}
\mathbf{K}^{uu} & \mathbf{K}^{up} \\
(\mathbf{K}^{up})^{T} & \mathbf{K}^{pp}
\end{bmatrix}
\left\{
\begin{matrix}
\mathbf{d}^u \\
\mathbf{d}^p
\end{matrix}
\right\}=
\left\{
\begin{matrix}
\mathbf{f} \\
\mathbf{0}
\end{matrix}
\right\}
$$



1. ​**速度變分 $\delta \mathbf{u}$**

$$
\delta L = \mu \int_\Omega \nabla \delta \mathbf{u} : \nabla \mathbf{u} \, d\Omega - \int_\Omega \delta\mathbf{u} \cdot \mathbf{b} \,d\Omega（1）
$$


2.​**壓力變分 $\delta p$**

$$
\delta L = -\int_\Omega \delta p \, (\nabla \cdot \mathbf{u}) \, d\Omega （2）
$$


**速度离散**

$$
\begin{split}
\delta u^h_{i}(\boldsymbol{x}) = \sum_{I=1}^{n_u} N_{I}(\boldsymbol{x}) \delta d_{iI}（3）\\
u^h_{i}(\boldsymbol{x}) = \sum_{J=1}^{n_u} N_{J}(\boldsymbol{x}) d_{iJ} （4）
\end{split}
$$


**应变率张量场**

$$
\begin{split}
\varepsilon^h_{ij}
&= \frac{1}{2} \left( u^h_{i,j} + u^h_{j,i} \right) \\
&= \sum_{I=1}^{n_p} \frac{1}{2} \left( N_{I,j} d_{iI} + N_{I,i} d_{jI} \right) （5）\\
&= \sum_{I=1}^{n_p} \frac{1}{2} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) d_{kI}（6）
\end{split}
$$


**（3）（4）（6）代入（1）**

$$
\begin{split}
\int_{\Omega} 2\mu \, \delta \varepsilon^h_{ij} \varepsilon^h_{ij} \, \mathrm{d}\Omega
&= \sum_{I,J}^{n_p} \int_{\Omega} 2\mu \frac{1}{4} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) \delta d_{kI} \left( N_{J,j} \delta_{il} + N_{J,i} \delta_{jl} \right) d_{lJ} \, \mathrm{d}\Omega \\
&= \sum_{I,J}^{n_p} \delta d_{kI} \int_{\Omega} 2\mu \frac{1}{4} \left( N_{I,j} \delta_{ik} + N_{I,i} \delta_{jk} \right) \left( N_{J,j} \delta_{il} + N_{J,i} \delta_{jl} \right) \mathrm{d}\Omega \, d_{lJ}\\
&=\sum_{I,J}^{n_p} \delta d_{kI} \int_{\Omega} \frac{1}{2} \mu \left( N_{I,j}N_{J,j}\delta_{kl} + 2N_{I,l}N_{J,k} + N_{I,i}N_{J,i}\delta_{kl} \right) \mathrm{d}\Omega \, d_{lJ}\\
&=\sum_{I,J}^{n_p} \delta d_{kI} \int_{\Omega} \frac{1}{2} \mu \left(2N_{I,i}N_{J,i}\delta_{kl} + 2N_{I,l}N_{J,k} \right) \mathrm{d}\Omega \, d_{lJ}\\
&=\sum_{I,J}^{n_p} \delta d_{kI} \mu \int_{\Omega} \left(N_{I,i}N_{J,i}\delta_{kl} + N_{I,l}N_{J,k} \right) \mathrm{d}\Omega \, d_{lJ}（7）
\end{split}
$$

$$
K_{IJkl} = \mu \int_{\Omega} (N_{I,i}N_{J,i}\delta_{kl} + N_{I,l}N_{J,k})d\Omega
$$

$$
K_{IJ} =\begin{bmatrix}\mu \int_{\Omega} 2(N_{I,x}N_{J,x}) + N_{I,x}N_{J,y}d\Omega &\mu \int_{\Omega} N_{I,x}N_{J,y}d\Omega\\\mu \int_{\Omega} N_{I,x}N_{J,y}d\Omega &\mu \int_{\Omega} N_{I,x}N_{J,x} + 2(N_{I,x}N_{J,y})d\Omega\end{bmatrix}
$$


