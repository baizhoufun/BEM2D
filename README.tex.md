# BEM2D = Boundary Element Method in Two Dimensions

## Problem statement

Solve 2D Laplace equation

$$
\nabla^2\phi = 0\quad\textrm{in}\quad \omega~,
$$

subject to boundary conditions imposed on domain boundary $\gamma$ ,

$$\left\{
\begin{aligned}
\phi &= \varphi\quad\textrm{on}\quad \gamma_\mathrm{D}~,\\
\bm{n}\cdot\nabla\phi &= q\quad\textrm{on}\quad \gamma_\mathrm{N}~,
\end{aligned}\right.
$$

where boundary $\gamma=\gamma_\mathrm{D}\cup\gamma_\mathrm{N}$ 
is split into a Dirichlet part $\gamma_\mathrm{D}$ and a Neumann part $\gamma_\mathrm{N}$ .