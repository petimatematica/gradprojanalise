# PG1.jl 

This file contains the implementation of the Projected Gradient Method (PG) equipped with Armijo's linear search along the projection arc.

## Function method1:
Implements an optimization algorithm with projection. It takes the following parameters:

- x0 (Vector): The initial point.
- f (Function): The objective function to be minimized.
- ∇f (Function): The gradient of the objective function.
- ε (Float64): The convergence tolerance.
- max_iter (Int): The maximum number of iterations allowed.
- PG1 (Function): The step size and projection update function.
- projection (Function): The projection function to ensure feasibility.

## Function PG1:
Armijo linesearch along the projection arc.

- x (Vector): The current point in the optimization process.
- f (Function): The objective function to be minimized.
- ∇f (Function): The gradient of the objective function.
- projection (Function): The projection function to ensure the point stays within the feasible region.
- η (Float64): The sufficient decrease parameter used in the line search condition.
- min_step (Float64): The minimum allowable step size before triggering an error.
- β_start (Float64): The initial step size guess.

### Returns:
- β (Float64): The computed step size.
- ierror (Int): Error flag, set to 1 if the step size becomes too small.
- evalf (Int): The number of function evaluations performed.
- projeval (Int): The number of projections performed.

# PG2.jl 

This file contains the implementation of the Projected Gradient (PG) Method equipped with Armijo's linesearch along feasible directions.

## Function method2:
Implements an optimization algorithm with projection. It takes the following parameters:

- x (Vector): The initial point for the optimization process.
- f (Function): The objective function to be minimized.
- ∇f (Function): The gradient of the objective function.
- ε (Float64): The convergence tolerance.
- max_iter (Int): The maximum number of iterations allowed.
- PG2 (Function): The function that computes the step sizes γ and β, as well as handles projection updates.
- projection (Function): The function that projects the updated point onto a feasible set.

## Function PG2:
Armijo linesearch along feasible directions.

- x (Vector): The current point in the optimization process.
- f (Function): The objective function to be minimized.
- ∇f (Function): The gradient of the objective function.
- projection (Function): The function that projects the updated point onto a feasible set.
- η (Float64): The sufficient decrease parameter used in the line search condition.
- min_step (Float64): The minimum allowable step size before triggering an error.
- γ_start (Float64): The initial guess for the step size γ.
- β2 (Float64): The upper bound for the step size β.

### Returns:
- γ (Float64): The computed step size for the projection update.
- β (Float64): The updated step size after evaluating the quadratic approximation.
- ierror (Int): Error flag, set to 1 if the step size γ becomes too small.
- evalf (Int): The number of function evaluations performed during the step size adjustment.

# PG1.jl 
This file contains the implementation of the Spectral Projected Gradient Method (SPG).

## Function spg:
Implements the Spectral Projected Gradient (SPG) Method. It performs optimization using projected gradients and adaptive step sizes. It takes the following parameters:

- x0 (Vector): The initial point for the optimization.
- f (Function): The objective function to be minimized.
- ∇f (Function): The gradient of the objective function.
- proj (Function): The projection function to ensure feasibility.
- ε (Float64): The convergence tolerance.
- max_iter (Int): The maximum number of iterations allowed.
- lambda_min (Float64): The minimum value for the step size adjustment parameter.
- lambda_max (Float64): The maximum value for the step size adjustment parameter.
- M (Int): History size for backtracking.
- sigma1 (Float64): Lower bound for backtracking step size adjustment.
- sigma2 (Float64): Upper bound for backtracking step size adjustment.
- η (Float64): The sufficient decrease parameter.
- linesearch (Function): The line search function to adjust step sizes.

## Function spg1:
Nonmonotone linesearch proposed by Grippo, Lampariello and Lucidi.

- k (Int): Current iteration number.
- lambda_k (Float64): Current step size parameter.
- x_k (Vector): Current point in the optimization process.
- gradf_x_k (Vector): Gradient at the current point.
- f_hist (Vector): History of function values.
- M (Int): History size for backtracking.
- sigma1 (Float64): Lower bound for step size adjustment.
- sigma2 (Float64): Upper bound for step size adjustment.
- η (Float64): Sufficient decrease parameter.

### Returns:
- x_plus (Vector): Updated point after step size adjustment.
- gradf_x_kp1 (Vector): Gradient at the updated point.
- s_k (Vector): Step direction.
- y_k (Vector): Difference in gradients.
- f_hist (Vector): Updated history of function values.
- alpha (Float64): Final step size.
- et (Float64): Elapsed time for the line search.
- evalf (Int): Number of function evaluations.
- evalproj (Int): Number of projection evaluations.

## Function spg2:
Nonmonotone linesearch proposed by Birgin, Raydan and Martínez.

- k (Int): Current iteration number.
- lambda_k (Float64): Current step size parameter.
- x_k (Vector): Current point in the optimization process.
- gradf_x_k (Vector): Gradient at the current point.
- f_hist (Vector): History of function values.
- M (Int): History size for backtracking.
- sigma1 (Float64): Lower bound for step size adjustment.
- sigma2 (Float64): Upper bound for step size adjustment.
- η (Float64): Sufficient decrease parameter.

### Returns:
- x_plus (Vector): Updated point after step size adjustment.
- gradf_x_kp1 (Vector): Gradient at the updated point.
- s_k (Vector): Step direction.
- y_k (Vector): Difference in gradients.
- f_hist (Vector): Updated history of function values.
- alpha (Float64): Final step size.
- et (Float64): Elapsed time for the line search.
- evalf (Int): Number of function evaluations.
- evalproj (Int): Number of projection evaluations.

# performances.jl 
This code compares different optimization methods (PG1, PG2, SPG1, SPG2) on various test problems using the CUTEst library. It measures the performance in terms of iterations, CPU time, function evaluations, and projection evaluations. The results are stored and performance profiles are generated.

## Requirements
Ensure the following Julia packages are installed:
- CUTEst
- NLPModels
- LinearAlgebra
- DataFrames
- Random
- Printf
- Plots
- BenchmarkProfiles
- JLD2

## Parameters

- η: 1.e-4 (The sufficient decrease parameter).
- ε: 1.e-5 (Convergence criterion).
- β_start: 1.0 (Initial value for the β parameter in the projected gradient method).
- β1: 0.01 (Minimum value for the β parameter).
- β2: 0.9 (Maximum value for the β parameter).
- γ_start: 1.0  (Initial step size γ in the linesearch).
- min_step: 1.e-5 (Minimum allowed step size γ).
- max_iter: 40000 (Maximum number of iterations).
- lambda_min: 1.e-30 (Minimum value for the spectral step parameter in the SPG method).
- lambda_max: 1.e+30 (Maximum value for the spectral step parameter in the SPG method).
- M: 10 (Number that determines the decrease of the function in nonmonotone linesearch).
- sigma1: 0.1 (Lower bound for step size update in the linesearch).
- sigma2: 0.9 (Upper bound for step size update in the linesearch).
