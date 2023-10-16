# README

本仓库代码用于解决 RoboMaster2022 赛季至今的大能量机关运动方程解算问题。

## 仓库结构

```
.
├─ Filter.m                       # 测试程序，使用卡尔曼滤波
├─ Monte.m                        # 测试程序，使用蒙特卡洛法分析
├─ Sim.m                          # 测试程序，模拟真实运行环境分析
├─ log
└─ lib
    ├─ getData.m                  # 获取数据序列
    ├─ meanFilter.m               # 均值滤波
    ├─ RMSE.m                     # 均方根误差
    ├─ OLS.m                      # 最小二乘法
    ├─ KF.m                       # 卡尔曼滤波
    ├─ GN.m                       # 高斯牛顿法
    └─ FFT.m                      # 带窗傅里叶变换法
```

*注：由于 Matlab 默认使用行向量，而数学上一般倾向于使用列向量，故请在理解代码和公式时注意转换关系*

## 问题描述

现有一半径为 700mm 的叶片绕定轴旋转，其旋转速度 $spd$ 与时间 $t$ 的关系式为

$$
spd(t) = a \sin(\omega t) + b \tag{1} \\
\text{s.t.} \ 
\begin{cases}
a \in \[0.780, 1.045\] \\
\omega \in \[1.884, 2.000\] \\
b = 2.090 - a
\end{cases}
$$

已知过去 $n$ 个时刻的叶片角度测量值和时间序列 $(\alpha_i, t_i), \ t_i = i \Delta t, \ i \in \{0, \cdots, n-1\}$ ，求任意时刻 $t$ 时叶片的角度，并尽可能使叶片末端的位置误差小于 30mm。（以下简称为“大符问题”）

## 方案分析

由于直接观测量为角度，首先将速度方程积分得角度方程，可得多种形式

$$
\begin{align}
f(t) &= -\frac{a}{\omega}\cos(\omega t + \phi) + (2.090-a) t + c \\
&= A \sin(\omega t + \phi) + b t + c \\
&= B \sin(\omega t) + C \cos(\omega t) + b t + c
\end{align}
$$

对 $t$ 离散化后得到两种下文中主要使用到的表达形式

$$
f(i) = -\frac{a}{\omega}\cos(i \omega \Delta t + \phi) + (2.090-a) i \Delta t + c \tag{2.1}
$$

$$
f(i) = B \sin(i \omega \Delta t) + C \cos(i \omega \Delta t) + b i \Delta t + c, \\
\text{s.t.} \ \omega\sqrt{B^2 + C^2} + b = 2.090 \tag{2.2}
$$

可以认为该问题中共含有 4 个未知参数，若能求解所有参数，即可预测任意时刻下叶片的旋转角度。从式 $(2.2)$ 中可以发现，若已知 $\omega$ ， $f(\cdot)$ 可以表示为其余 4 个参数和已知量的线性表达式，可通过最小二乘法直接求解，于是问题的关键转化为对 $\omega$ 的求解。

第一种直观的方案，通过带窗的快速傅里叶变换求解速度方程的信号主频率 $\omega$ ；

第二种方案，通过高斯牛顿法求解式 $(2.1)$ 的最优参数；

第三种方案可以认为是第二种的改版，设计如下交替迭代法求解式 $(2.2)$ 的优化参数：

1. 以上一时刻计算的 $\omega$ 值为当前时刻的初值，若当前为初始时刻，则根据 $\omega$ 的取值范围设定初值为中间数 1.942；
2. 固定 $\omega$，根据最小二乘法求取 4 个线性参数；
3. 固定 4 个线性参数，使用高斯牛顿法迭代优化 $\omega$ ；
4. 重复步骤 2~3，直到 $\Delta \omega < \eta$ 。

此外，为了解决噪声对拟合算法的影响，我们在进入优化之前需要对角度数据进行滤波处理，实验的测试中主要使用了均值滤波，更好的选择是卡尔曼滤波；为了降低优化的复杂度，提升运行速度，我们可以对滤波后的数据进行降采样。

## 核心算法理论概述

### 最小二乘法

最小二乘法用于求解如下最小二乘问题：现有线性方程 $y = f(\boldsymbol{x}) = \boldsymbol{\theta}^T \boldsymbol{x}, \ \boldsymbol{x} \in \mathbb{R}^m$ 和 $n$ 个样本 $(\boldsymbol{x}\_{i}, y_{i}), \ i \in \{0, \cdots, n-1\}$ ，求使得残差平方和 $J = \sum\limits_{i=0}^{n-1} (y_{i} - f(\boldsymbol{x}_{i}))^2$ 最小的方程参数 $\boldsymbol{\theta}$。

令 $X = \[\boldsymbol{x}\_0, \cdots, \boldsymbol{x}\_{n-1}\]^T, \ Y = \[y_0, \cdots, y_{n-1}\]^T$ ，由于

$$
\begin{align}
J(\boldsymbol{\theta})
&= \sum_{i=1}^n (y_i - \boldsymbol{\theta}^T \boldsymbol{x}_i))^2 \\
&= (Y - X \boldsymbol{\theta})^T (Y - X \boldsymbol{\theta})
\end{align}
$$

求梯度，化简得 $\boldsymbol{\theta}$ 的显式解。

$$
\frac{\partial J(\boldsymbol{\theta})}{\partial \boldsymbol{\theta}} = X^T(Y - X^T \boldsymbol{\theta}) = 0 \\
\boldsymbol{\theta} = (X^T X)^{-1} X^T Y
$$

对于大符问题，我们令 $\boldsymbol{x}_i = \[\sin(i \omega \Delta t), \cos(i \omega \Delta t), i \Delta t, 1\]^T$ ，构造

$$
X = \begin{pmatrix}
0 & 1 & 0 & 1 \\
\sin(\omega \Delta t) & \cos(\omega \Delta t) & \Delta t & 1 \\
\sin(2\omega \Delta t) & \cos(2\omega \Delta t) & 2 \Delta t & 1 \\
\vdots & \vdots & \vdots & \vdots \\
\sin((n-1)\omega \Delta t) & \cos((n-1)\omega \Delta t) & (n-1) \Delta t & 1
\end{pmatrix}, \ 
Y = \begin{pmatrix}
\alpha_0 \\
\alpha_1 \\
\alpha_2 \\
\vdots \\
\alpha_{n-1}
\end{pmatrix}, \ 
\boldsymbol{\theta} = \begin{pmatrix}
B \\ C \\ b \\ c
\end{pmatrix} \\
$$

在固定 $\omega$ 的情况下，直接求解即可。

### 高斯牛顿法

高斯牛顿法适用于求解如下优化问题：设误差函数为 $e(\boldsymbol{\omega})$ ，优化目标为

$$
\min F(\boldsymbol{\omega}) = \|e(\boldsymbol{\omega})\|^2
$$

设 $\boldsymbol{\omega}$ 在 $\boldsymbol{\omega}^{(k)}$ 处具有增量 $\Delta \boldsymbol{\omega}^{(k)}$

$$
\begin{align}
F(\boldsymbol{\omega}^{(k)} + \Delta \boldsymbol{\omega}^{(k)})
&= \|e(\boldsymbol{\omega}^{(k)} + \Delta \boldsymbol{\omega}^{(k)})\|^2 \\
&\simeq \|e(\boldsymbol{\omega}^{(k)}) + J_e(\mathbf{\omega}^{(k)})^T \Delta \boldsymbol{\omega}^{(k)}\|^2
\end{align}
$$

其中 $J_e(\boldsymbol{\omega})$ 为 $e(\boldsymbol{\omega})$ 的一阶导数，称为雅各比矩阵。对其求关于 $\Delta \boldsymbol{\omega}^{(k)}$ 的偏导得

$$
\begin{align}
& \frac{\partial \|e(\boldsymbol{\omega}^{(k)}) + J_e(\boldsymbol{\omega}^{(k)})^T \Delta \boldsymbol{\omega}^{(k)}\|^2}{\partial \Delta \boldsymbol{\omega}^{(k)}} = 0 \\
\Rightarrow \; & \Delta \boldsymbol{\omega}^{(k)} = - (J_e(\boldsymbol{\omega}^{(k)}) J_e(\boldsymbol{\omega}^{(k)})^T)^{-1} J_e(\boldsymbol{\omega}^{(k)}) e(\boldsymbol{\omega}^{(k)})
\end{align}
$$

算法迭代步骤如下：

1. 给定初值 $\boldsymbol{\omega}^{(0)}$ ；
2. 在第 k 次迭代中，根据雅各比矩阵 $J_e(\boldsymbol{\omega}^{(k)})$ 和 $e(\boldsymbol{\omega}^{(k)})$ 求解增量 $\Delta \boldsymbol{\omega}^{(k)}$ ；
3. 令 $\boldsymbol{\omega}^{(k+1)} = \boldsymbol{\omega}^{(k)} + \beta \Delta \boldsymbol{\omega}^{(k)}$ ，其中 $\beta$ 为步长倍率；
4. 重复步骤 2~3，直到 $\|\Delta \boldsymbol\omega\|^2 < \eta$ 。

对于大符问题，我们令 $e(\omega) = \[e_0, \cdots, e_{n-1}\]^T, \ e_i = \alpha_i - \boldsymbol{\theta}^T \boldsymbol{x}$ ，从而有

$$
J_e(\omega) = \begin{pmatrix}
\cdots &  - B i \Delta t \cos(i \omega \Delta t) + C i \Delta t \sin(i \omega \Delta t) & \cdots
\end{pmatrix}
$$

在固定 $\boldsymbol{\theta}$ 的情况下，按上述步骤迭代求解即可；若对所有参数同时进行优化，则需要求解 $\hat{e}(\boldsymbol\omega) = \alpha - (-\frac{a}{\omega}\cos(\omega t) + (2.090-a) t + c)$ 关于四个原始参数 $\boldsymbol\omega = [a, \omega, \phi, c]^T$ 的偏导数，且终止条件需改为 $\|\Delta \boldsymbol\omega\|^2 < \eta$ ，其余同理。

$$
J_{\hat{e}}(\boldsymbol{\omega}) = \begin{pmatrix}
\cdots & \frac{1}{\omega} \cos(i \omega \Delta t + \phi) + i \Delta t & \cdots \\
\cdots & -\frac{a i \Delta t}{\omega^2} \cos(i \omega \Delta t + \phi) - \frac{a i \Delta t}{\omega} \sin(i \omega \Delta t + \phi) & \cdots \\
\cdots & -\frac{1}{\omega} \sin(i \omega \Delta t + \phi) & \cdots \\
\cdots & 1 & \cdots \\
\end{pmatrix}
$$

### 卡尔曼滤波

卡尔曼滤波的原理不再赘述，此处仅记录其五步公式

$$
\begin{align}
\bar{\boldsymbol{x}\_{k+1}} &= F \boldsymbol{x}\_k \\
\bar{P_{k+1}} &= F P_k F^T + Q \\
K_{k+1} &= \frac{\bar{P_{k+1}} H^T}{H \bar{P_{k+1}} H^T + R} \\
\boldsymbol{x}\_{k+1} &= \bar{\boldsymbol{x}\_{k+1}} + K_{k+1} (\boldsymbol{z}\_{k+1} - H \bar{\boldsymbol{x}\_{k+1}}) \\
P_{k+1} &= (I - K_{k+1} H) \bar{P_{k+1}}
\end{align}
$$

对于大符问题，为了拟合原方程中的正弦分量，联想到简谐振动 $m\ddot{x} + kx = 0$ 的解为 $x(t) = A \sin(\omega t + \phi), \ \omega = \sqrt{\frac{k}{m}}$ ，构造状态空间方程如下

$$
\boldsymbol{x}_k = \begin{pmatrix}
\mathtt{x} \\ \dot{\mathtt{x}}
\end{pmatrix}, \ 
\boldsymbol{z}_k = \alpha_i - b t_i - c, \\
F = \begin{pmatrix}
1 & \Delta t \\
-\omega^2 \Delta t & 1
\end{pmatrix}, \ 
H = \begin{pmatrix}
1 & 0
\end{pmatrix}
$$

## 实验效果

### 带窗傅里叶变换法

由于角度数据经过微分后误差放大，且数据中包含周期数极少，使用 600 帧数据通过该算法求解 $\omega$ 的均方根误差约为 0.015；结合最小二乘法求解角度方程后，在 0.3s 后预测坐标的均方根误差达到 15mm 左右，预测误差在可接受范围内，优点在于运算速度较快、受噪声影响相对小。

若采用不带窗的形式，则求解误差极大，结果基本无效。

### 最小二乘法

效果非常稳定； $\omega$ 越精准，效果越好。该算法可以认为是大符问题的基石。

### 高斯牛顿法

#### 直接优化法

直接对 4 个原始参数进行优化效果喜人，求解噪声较小时（高斯噪声， $\sigma < 0.02$ ）的运动方程具有相当良好的精度，使用 600 帧数据可以使得 $\omega$ 的均方根误差小于 0.0015，在 0.3s 后预测坐标的均方根误差达到 2mm 水平，散布保持在直径 10mm 以内。

#### 交替迭代法

使用交替迭代法直接求解噪声较小时（高斯噪声， $\sigma < 0.02$ ）的运动方程具有相当良好的精度，使用 600 帧数据可以使得 $\omega$ 的均方根误差小于 0.002，在 0.3s 后预测坐标的均方根误差达到 2mm 水平，散布保持在直径 10mm 以内。相比于直接优化法，  $a+b=2.090$ 的约束并不严格满足，且预测精度略低，但是求解速度略快（Matlab 测试并不代表 C 语言运行速度）。

上述两种方案就实验结果而言相距甚小，且可以通过终止条件进一步权衡精度与计算速度，可根据实际测试效果灵活选用。高斯牛顿法的共性问题是在噪声较大的情况下预测精度将持续下降，甚至经常出现反向优化的情况；且 $\omega$ 的优化易呈现与时间相关的有偏估计，导致一定时间段内的预测误差呈现有偏分布。

关于这两种方案还有两个现象。其一是直接优化法更容易在连续过程中计算出错误的 $\omega$ 的下降方向，这可能是由于梯度计算较为复杂导致的，但奇怪的是交替迭代法则更容易在随机测试中出现该问题，总体而言发生这类情况的概率都不大。其二是交替迭代法一般需要通过设置步长倍率加快下降速度，一般取 $\beta\in(2,5)$ ，而直接优化法则无须设置。以上现象暂时未知其原因。

### 卡尔曼滤波

使用 200 帧数据通过最小二乘法求取初始值，并逐步更新参数，可有效对 $\sigma < 0.1$ 高斯噪声与 $s < 0.02, \ \sigma < 0.5$ 稀疏噪声的融合噪声进行滤波（可通过调节滤波器参数优化性能，此处实验数据仅供参考）。另外，由于初值求取的条件一般较为简陋，初值一定范围的领域内滤波器输出将有较大偏差（依据噪声决定）。

理论上来说，线性分量亦可加入状态量中，但实际测试效果不佳，故不予采纳。若 $\omega$ 取值不精准，滤波器将在数个周期内产生明显可见的偏差，故不可摒弃对 $\omega$ 的优化求解。

