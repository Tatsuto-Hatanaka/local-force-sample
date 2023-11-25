# Local-force-sample

## Overview

This is the sample code to calculate the effective Heisenberg interaction on each atom by the local force method, originally developed by Liechtenstein[1]. Here, we consider the simple model, namely the single orbital Hubbard model on an infinite-dimensional Bethe lattice.

## Usage

Run the code below. Figures will be output in the `figures` directory.

```python
python local_force_sample.py
```

## Note

We consider the single orbital Hubbard model, and the Hamiltonian is defined below.

```math
\begin{align}
\mathcal{H} = -\sum_{ij}\sum_{\sigma}t_{ij}c^{\dagger}_{i\sigma}c_{j\sigma} + U\sum_{i}n^{\uparrow}_{i}n^{\downarrow}_{i}\nonumber
\end{align}
```

This Hamiltonian is transformed to the form below under the mean-field approximation of the second term $`n^{\uparrow}_{i}n^{\downarrow}_{i}`$, and this Hamiltonian corresponds to the single particle eigenstates problem (= band calculation).

```math
\begin{align}
\mathcal{H} = -\sum_{ij}\sum_{\sigma}t_{ij}c^{\dagger}_{i\sigma}c_{j\sigma} &+ \sum_{i}V_{i}n_{i} - \sum_{i}\boldsymbol{B}_{i}\cdot\boldsymbol{m}_{i}+\frac{U}{4}\left(\langle\boldsymbol{m}_{i}\rangle^2-\langle n_{i}\rangle^2\right)\nonumber\\
\boldsymbol{m}_{i} &= \sum_{\sigma\sigma'}c^{\dagger}_{i\sigma}\boldsymbol{\sigma}_{\sigma\sigma'}c_{j\sigma'}\nonumber\\
V_i &= \frac{U}{2}\langle n_{i}\rangle\nonumber\\
\boldsymbol{B}_{i} &= \frac{U}{2}\langle\boldsymbol{m}_{i}\rangle\nonumber
\end{align}
```

The magnetization $m_0\equiv\langle\boldsymbol{m}_{i}\rangle$ can be obtained from the self-consistent equation under the mean-field approximation. Here, we assumed the ferromagnetic state polarized along the $z$-axis for simplicity, and hence, the Green's function is diagonal in the spin space.

```math
\begin{align}
G^{\sigma\sigma}_{ij}(\omega) &= \frac{1}{N}\sum_{k}\frac{e^{i\boldsymbol{k\cdot(\boldsymbol{R}_i-\boldsymbol{R}_j)}}}{\omega-\epsilon_{\boldsymbol{k}}+\frac{U}{2}m_{0}\sigma}\nonumber\\
m_{0} &= -\frac{1}{\pi}\mathrm{Im}\int d\epsilon f(\epsilon)\left(G^{\uparrow\uparrow}_{00}(\epsilon)-G^{\downarrow\downarrow}_{00}(\epsilon)\right)\nonumber
\end{align}
```

Then, we can calculate the effective Heisenberg interaction on each atom (orbital) using self-consistently determined $m_{0}$.

```math
\begin{align}
J^{(i)}_{0}&\equiv \sum_{j\neq i}J_{ij}\\
J_{0} = -\frac{1}{4\pi}\mathrm{Im}\int d\epsilon f(\epsilon)&\left[Um_{0}\left(G^{\uparrow\uparrow}_{00}(\epsilon)-G^{\downarrow\downarrow}_{00}(\epsilon)\right)+(Um_{0})^{2}G^{\uparrow\uparrow}_{00}(\epsilon)G^{\downarrow\downarrow}_{00}(\epsilon)\right]\nonumber
\end{align}
```

Here, we consider an infinite-dimensional Bethe lattice with a hopping term between only nearest neighbors. The Green's functions can be calculated below for the problem on this lattice[3].

```math
\begin{align}
G_{00}(\epsilon)=\frac{2}{D}\left[\frac{\epsilon}{D}-\mathrm{sgn}(\epsilon)\sqrt{\left(\frac{\epsilon}{D}\right)^{2}-1}\right]\nonumber
\end{align}
```

The Green's functions for up/down spin can be obtained below.

```math
\begin{align}
G^{\sigma\sigma}_{00} (\epsilon)&=G_{00}\left(\epsilon+\frac{Um_{0}}{2}\sigma\right)\nonumber
\end{align}
```

## Reference

1. A. I. Liechtenstein,M. I. Katsnelson, andV. A. Gubanov, Exchange interactions and spin-wave stiffness in ferromagnetic metals, [J.Phys. F14, L125 (1984).](https://doi.org/10.1088/0305-4608/14/7/007)
2. A. Sakuma, Theoretical study on the exchange constants of the transition metal systems, [IEEE Trans. Magn. 35, 3349 (1999).](https://doi.org/10.1109/20.800521)
3. Green’s Functions for Tight-Binding Hamiltonians, pages 77–110. [Springer Berlin Heidelberg, Berlin,
Heidelberg, 2006.](https://doi.org/10.1007/3-540-28841-4_5)

## Author

Tatsuto Hatanaka

## License

MIT license
