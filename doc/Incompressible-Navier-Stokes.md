# Incompressible Navier-Stokes Equation

## 有限体积法离散化
#### 守恒方程微分形式:
$$
\begin{aligned}
\frac{\partial \rho}{\partial t} + \frac{\partial}{\partial x_i} (\rho u_i) &= 0\\

\frac{\partial}{\partial t} (\rho u_i) + \frac{\partial}{\partial x_j}(\rho u_i u_j) &= -\frac{\partial}{\partial x_i} p + \frac{\partial}{\partial x_i} (\tau_{ij})\\

\tau_{ij} &= \mu \left[ \left( \frac{\partial}{\partial x_j} u_i + \frac{\partial}{\partial x_j} u_i^{T} \right) - \frac{2}{3} \frac{\partial}{\partial x_i} u_i I \right]
\end{aligned} \quad (1)
$$


- $\rho : 流体密度$
- $t : 时间 $
- $u_i:速度矢量$
- $p:压力$
- $\tau_{ij} : 粘性应力张量$
- $\mu:动力粘度$
- $I:单位张量$ 

以标量$\phi$的瞬态输运方程为例说明控制方程的离散化$u_i ->\phi$。对于任意控制体$\mathbf V$，以积分形式写成下面的形式：
$$
\int_{\mathbf V} \frac {\partial \rho \phi}{\partial t}d\mathbf V \int_f \rho \phi u_j \cdot d \mathbf A = \int_f \Gamma_\phi \frac{\partial}{\partial x_i} \phi \cdot d \mathbf A + \int_\mathbf V S_\phi d \mathbf V \quad(2)
$$

式中，$\rho $为密度； $u_j$为速度向量； $A$为面积向量；$\Gamma_\phi$为扩散系数；$\nabla \phi为\phi$的梯度；$S_\phi$为单位体积的$\phi$的源项。
<img src=/home/ddd/文件/stokes/Stokes-equation/20220114102032.png>

#### 离散方程为：
$$
\frac{\partial \rho \phi}{\partial t} V + \sum_{f}^{N_{\text{faces}}} \rho_{f} \mathbf{v}_{f} \phi_{f} \cdot \mathbf{A}_{f} = \sum_{f}^{N_{\text{faces}}} \Gamma_{\phi} \frac{\partial}{\partial x_i} \phi_{f} \cdot \mathbf{A}_{f} + S_{\phi} V \quad (3)
$$

式中，$N_{faces}$为封闭网格单元的面数；$\phi_f$为通过面$f$的对流量；$\rho_f \mathbf v_f \cdot a_f$为网格面f的面积向量；$\frac{\partial}{\partial x_i} \phi_f$为网格面$f$上变量$\phi$的梯度；$\mathbf V$为网格体积。

$$a_p \phi + \sum_{ub} a_{nb} \phi_{nb} = b_p \quad (4)$$

## Pressure - based

<img src=/home/ddd/文件/stokes/Stokes-equation/pressure-based.svg>

无散度约束转化为压力的椭圆型方程，并且假设空间导数已经完成离散化处理

$$\frac {\partial}{\partial x_i}\left(\frac{\partial p}{\partial x_i}\right) = \frac {\partial}{\partial x_i}\left[\frac{\partial \rho u_iu_j}{\partial x_j}\right] \quad(5)$$

$$\frac {\partial}{\partial x_i} \Rightarrow \frac {\delta}{\delta x_i} $$

#### 显式时间积分:
$$\frac {\phi ^{n+1}-\phi^n}{\Delta t} = \mathbf F(\phi^n)\quad (6)$$
代入（1）

$$
\frac{\partial\left(\rho u_{i}\right)}{\partial t} = -\frac{\delta\left(\rho u_{i}u_{j}\right)}{\delta x_{j}} + \frac{\delta\tau_{ij}}{\delta x_{i}} - \frac{\delta p}{\delta x_{i}} = H_{i} - \frac{\delta p}{\delta x_{i}}
$$

$$
\begin {aligned}
\left(\rho u_{i}\right)^{n+1} - \left(\rho u_{i}\right)^{n} &= \Delta t\left(H_{i}^{n} - \frac{\delta p}{\delta x_{i}}^{n}\right) \quad(7)\\
\frac {\delta (\rho u_i)^{n+1}}{\delta x_i} &\neq 0
\end {aligned}
$$
压力泊松方程:
$$
\frac{\delta}{\delta x_{i}}\left(\frac{\delta p^{n}}{\delta x_{i}}\right) = \frac{\delta}{\delta x_{i}}H_{i}^{n}\quad(8)
$$

#### 隐式时间积分:
$$\begin{aligned}
\frac {\phi ^{n+1}-\phi^n}{\Delta t} &= \mathbf F(\phi^{n+1})\\
\phi ^{n+1} &= \phi^n + \Delta t \mathbf F (\phi ^{n+1})
\end{aligned}\quad (9)$$
代入（1）

$$
\frac{\partial\left(\rho u_{i}\right)}{\partial t} = -\frac{\delta\left(\rho u_{i}u_{j}\right)}{\delta x_{j}} + \frac{\delta\tau_{ij}}{\delta x_{i}} - \frac{\delta p}{\delta x_{i}} = H_{i} - \frac{\delta p}{\delta x_{i}}
$$

$$
\left(\rho u_{i}\right)^{n+1} - \left(\rho u_{i}\right)^{n} = \Delta t\left(H_{i}^{n+1} - \frac{\delta p}{\delta x_{i}}^{n+1}\right)\quad(10)
$$
压力泊松方程:
$$
\frac{\delta}{\delta x_{i}}\left(\frac{\delta p^{n+1}}{\delta x_{i}}\right) = \frac{\delta}{\delta x_{i}}H_{i}^{n+1}\quad(11)
$$


### 投影法
$$\frac{\partial}{\partial t} u_i+ \frac{\partial}{\partial x_j}(\rho u_i u_j) = -\frac{\partial}{\partial x_i} p + \frac{1}{\rho}\frac{\partial}{\partial x_i} (\tau_{ij})\quad(1)$$

$u_i^{n+1}-u_i^n \Rightarrow (u_i^{n+1} - u_i^*)+(u_i^* - u_i^n)$

$$\frac {u_i^* - u_i^n}{\Delta t} = -\frac{\partial}{\partial x_j}(\rho u_i^t u_j^t) + \frac{\partial}{\partial x_i} (\tau_{ij}^*) \quad(12)$$

$$\frac {u_i^{n+1} - u_i^*}{\Delta t} = -\frac{\partial}{\partial x_i} p^* \quad(13)$$

$$
\begin{aligned}
\Rightarrow \frac {\delta}{\delta x_i}(\frac {u_i^{n+1}}{\Delta t} - \frac {u_i^*}{\Delta t}) &= \frac {\delta}{\delta x_i}(-\frac{\delta}{\delta x_i} p^*)\\
\frac {\delta}{\delta x_i}u_i^{n+1} &= 0\\

\frac{1}{\Delta t}\frac{\delta u_i^*}{\delta x_i} &= \frac {\delta ^2 p^*}{\delta^2 x_i}\quad(14)\\
\frac{u_i^{n+1}-u_i^*}{\Delta t} &= \frac{1}{\rho}\frac{\delta p^*}{\delta x_i}\quad(15)
\end{aligned}$$

### Uzawa演算法

$$
\left[\begin{array}{cc}
A & \nabla \\
\nabla \cdot  & 0
\end{array}\right]
\left[\begin{array}{c}
u^{n+1} \\
p^{n+1}
\end{array}\right]
=
\left[\begin{array}{c}
RHS \\
0
\end{array}\right]
$$

#### LU分解方法

$$
\left[\begin{array}{cc}
A & 0 \\
\nabla \cdot  & -\nabla \cdot A^{-1}\nabla
\end{array}\right]
\left[\begin{array}{cc}
I & A^{-1}\nabla \\
0 & I
\end{array}\right]
\left[\begin{array}{c}
u^{n+1} \\
p^{n+1}
\end{array}\right]
=
\left[\begin{array}{c}
RHS \\
0
\end{array}\right]
$$

#### 精确分裂方法（Exact Splitting）
$$
\left\{
\begin{array}{c}
\left[\begin{array}{cc}
A & 0 \\
\nabla \cdot  & -\nabla \cdot A^{-1}\nabla
\end{array}\right]
\left[\begin{array}{c}
u^{*} \\
p^{*}
\end{array}\right]
=
\left[\begin{array}{c}
RHS \\
0
\end{array}\right] \\[2em]
\left[\begin{array}{cc}
I & A^{-1}\nabla \\
0 & I
\end{array}\right]
\left[\begin{array}{c}
u^{n+1} \\
p^{n+1}
\end{array}\right]
=
\left[\begin{array}{c}
u^{*} \\
p^{*}
\end{array}\right]
\end{array}
\right.
$$

### 分步求解过程

#### 第一步：动量方程预测
$$ A u^{*} = RHS $$

#### 第二步：压力泊松方程
$$ \nabla \cdot  u^{*} - \nabla \cdot  A^{-1} \nabla p^{*} = 0 $$

#### 第三步：速度修正
$$ u^{n+1} + A^{-1} \nabla p^{n+1} = u^{*} $$

#### 第四步：压力更新
$$ p^{n+1} = p^{*} $$

- $\mathbf A :  (\mathbf{u} \cdot \nabla)  - \nu \nabla^2 $ 
- $u^{*} : 预测速度场$ 
- $p^{*}: 中间压力场 $
- $RHS :(\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \nabla^2 \mathbf{u}$

