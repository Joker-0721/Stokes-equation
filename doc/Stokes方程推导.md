# Stokes方程推导

**叙述**: 对于两场形式(u和 p)强制将不可压缩约束与动量方程耦合,导致离散系统需同时满足两种矛盾的数值需求，则可能导致LBB条件不满足，压力场出现高频振荡等問題。所以对于不可压缩牛顿流体的稳态流动可表述为 $(u,\varepsilon,\sigma, p)$ ,分别對應速度场 $(u)$ ,应变率张量场 $(\varepsilon)$ ,应力张量场 $(\sigma)$ ,压力场 $(p)$。

## 边界条件:

### 1. 速度边界:

$$
u=0\quad\text{ on}\Omega
$$

### 2. 自然边界:

$$
\sigma\cdot n=t\quad\text{ on}\Gamma_t
$$

## 1.1 体积不可压材料有限元分析

考虑如图 2.1所示维度为 $n_d$ 且具有边界 $\Gamma$ 的体积不可压材料弹性体 $\Gamma_t$ 和 $\Gamma_g$ 。分别表示其自然边界和本质边界,并满足 $\Gamma_t\cup\Gamma_g=\Gamma$ , $\Gamma_t\cap\Gamma_g=\varnothing$ 。该问题相应的强形式为:

$$
\left\{\begin{array}{ll}\nabla\cdot\sigma+b=0&\text{ in}\Omega\\ \sigma\cdot n=t&\text{ on}\Gamma_t\\ \nabla\cdot u=0&\text{ in}\Omega\end{array}\right.
$$

其中 $\sigma$ 为应力张量,对于各向同性线弹性材料,其本构关系表示为:

$$
\sigma=2\mu\varepsilon-p 1\quad(1.2)
$$

$1=\delta_{i j}$ ,是二阶恒等张量。式中 $\varepsilon$ 为应变率张量:

$$
\varepsilon=\frac{1}{2}\left(u_{i, j}+u_{j, i}\right)
$$

弹性力学方程组(1.1)相对应的能量泛函$\Pi$具有如下形式:

$$
L(u, p)=\frac{\mu}{2}\int_{\Omega}\varepsilon:\varepsilon d\Omega-\int_{\Omega} u\cdot b d\Omega-\int_{\Omega} u\cdot t d\Gamma-\int_{\Omega} p(\nabla\cdot u) d\Omega
$$

对上式进行变分:

$$
\text{ findu}\in V, a(\delta u, u)+a^1(\delta p, u)=f(\delta u)
$$

具有以下形式:

$$
\begin{align*}& a(\delta u, u)=\mu\int_{\Omega}\delta\varepsilon(u):\varepsilon(u) d\Omega\\ & a^1(\delta p, u)=-\int_{\Omega}\delta p(\nabla\cdot u) d\Omega\end{align*}
$$

$$
f(\delta u)=\int_{\Omega}\delta u\cdot b d\Omega-\int_{\Gamma_t}\delta u\cdot t d\Gamma\quad(1.7)
$$

## 速度离散

$$
\begin{align*}\delta u_i^h(x)&=\sum_{I=1}^{n_u} N_I(x)\delta d_{i I}, u_i^h(x)=\sum_{J=1}^{n_u} N_J(x) d_{i J}\end{align*}\qquad(1.8)
$$

### (1.8)代入 $\varepsilon$

$$
\begin{align*}\varepsilon_{i j}^h&=\frac{1}{2}\left(u_{i, j}^h+u_{j, i}^h\right)\\ &=\sum_{I=1}^{n_p}\frac{1}{2}\left(N_{I, j} d_{i I}+N_{I, i} d_{j I}\right)\\ &=\sum_{I=1}^{n_p}\frac{1}{2}\left(N_{I, j}\delta_{i k}+N_{I, i}\delta_{j k}\right) d_{k I}\end{align*}
$$

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

- Case $(k, l)=(1,1)$:

$$
K_{I J 11}=K[2 I-1,2 J-1]=\mu\int_{\Omega}\left(2 N_{I, x} N_{J, x}+N_{I, y} N_{J, y}\right) d\Omega
$$

- Case $(k, l)=(1,2)$:

$$
K_{I J 12}=K[2 I-1,2 J]=\mu\int_{\Omega} N_{I, y} N_{J, x} d\Omega
$$

- Case $(k, l)=(2,1)$:

$$
K_{I J 21}=K[2 I, 2 J-1]=\mu\int_{\Omega} N_{I, x} N_{J, y} d\Omega
$$

- Case $(k, l)=(2,2)$:

$$
\begin{gathered} K_{I J 22}=K[2 I, 2 J]=\mu\int_{\Omega}\left(N_{I, x} N_{J, x}+2 N_{I, y} N_{J, y}\right) d\Omega\\ {\left[\begin{array}{l}\delta d_{1 I}\\ \delta d_{2 I}\end{array}\right]\left[\begin{array}{cc}k[2 I-1,2 J-1]& k[2 I, 2 J-1]\\ k[2 I-1,2 J]& k[2 I, 2 J]\end{array}\right]\left[\begin{array}{l}d_{1 J}\\ d_{2 J}\end{array}\right]}\end{gathered}
$$

## 壓力離散和散度项表达式:

$$
\delta p^h(x)=\sum_{J=1}^{n_p}\Psi_J(x)\delta p_J,\nabla\cdot u^h=u_{k, k}^h=\sum_{I=1}^{n_u} N_{I, k} d_{k I}\quad(1.10)
$$

代入弱形式(连续性方程)

$$
\begin{gathered}\int_{\Omega}\delta p^h\left(\nabla\cdot u^h\right) d\Omega\\ K_{I J k}=\int_{\Omega}\Psi_J N_{I, k}~d\Omega\\ =\sum_{I, J=1}^{n_p}\delta p_J\left(\int_{\Omega} N_J N_{I, k}~d\Omega\right) d_{k I}\end{gathered}
$$

## 离散化系统:

$$
\begin{gathered} K_{I J k l}=\int_{\Omega} 2\mu\cdot\frac{1}{4}\left(N_{I, j}\delta_{i k}+N_{I, i}\delta_{j k}\right)\left(N_{J, j}\delta_{i l}+N_{J, i}\delta_{j l}\right)\\ K_{I J k}=\int_{\Omega}\Psi_J N_{I, k}~d\Omega\\ K^{p p}=0\\ {\left[\begin{array}{cc}K^{u u}& K^{u p}\\ \left(K^{u p}\right)^T& K^{p p}\end{array}\right]\left\{\begin{array}{l}d^u\\ d^p\end{array}\right\}=\left\{\begin{array}{l}f\\ 0\end{array}\right\}}\end{gathered}
$$


