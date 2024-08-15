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
