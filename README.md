# BEM2D = Boundary Element Method in Two Dimensions

## Problem statement

Solve 2D Axisymmetric Laplace equation

<p align="center"><img src="/tex/dde9b62e32b758022822abed45e1d0d4.svg?invert_in_darkmode&sanitize=true" align=middle width=163.8661695pt height=18.312383099999998pt/></p>

subject to boundary conditions imposed on domain boundary <img src="/tex/11c596de17c342edeed29f489aa4b274.svg?invert_in_darkmode&sanitize=true" align=middle width=9.423880949999988pt height=14.15524440000002pt/> ,

<p align="center"><img src="/tex/40ee62c6197af3fa9ba6be91668c099b.svg?invert_in_darkmode&sanitize=true" align=middle width=170.64913635pt height=49.315569599999996pt/></p>

where boundary <img src="/tex/6a134afa087f6f4f0518d21f6324e905.svg?invert_in_darkmode&sanitize=true" align=middle width=87.04995749999999pt height=18.264896099999987pt/> 
is split into a Dirichlet part <img src="/tex/c6756650b41a2a1f7f4ab6b70d2f124d.svg?invert_in_darkmode&sanitize=true" align=middle width=18.407965949999994pt height=14.15524440000002pt/> and a Neumann part <img src="/tex/3ec459e527f7d6875940e3190e2c199d.svg?invert_in_darkmode&sanitize=true" align=middle width=18.21390614999999pt height=14.15524440000002pt/> .

## Compile with CMAKE

```bash
mkdir build/
cd build
cmake ..
make -j4
```
