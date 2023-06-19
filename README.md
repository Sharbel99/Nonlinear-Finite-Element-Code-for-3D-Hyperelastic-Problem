# Nonlinear-Finite-Element-Code-for-3D-Hyperelastic-Problem
A program showcasing the process of solving a 3D hyperelastic problem in an 8-node cubic element using the Newton-Raphson method
Two cases are presented: tension and shear. Various stress, strain, tangent components are plotted in addition to the residual log.

## Description of the problem
The code solves the problem having the following energy function: 
```math
 U \left( J^e \right) = \dfrac{1}{2} \kappa \left(\dfrac{1}{2}\left({J^e}^2 -1 \right) - ln J^e \right)
```
```math
\bar{W}\left(\bar{b}^e\right) = \mu_1 \left(tr \left[\bar{b}^e\right] -3 \right) +\dfrac{1}{2} \mu_2 \left(tr \left[\bar{b}^e\right] -3 \right)^2
```

## Instructions
First run the input file (Input_tension.m or Input_shear.m) then run FEA_Program.m

To edit the code to run a different energy functional, edit the parameters in the input file, and edit the residual and tangent in Elast3d_Elem.m

To change the plots, edit FEA_Program.m
