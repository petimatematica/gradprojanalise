## Projected Gradient Method with GPA1 ##

function PG2(x, f, ∇f, projection, η, min_step, γ_start, β2)
    ierror = 0
    β = β2
    γ = γ_start
    j = 0
    zk = projection(x - β * ∇f(x))
    #println("zk =", zk)
    if x == x0
        evalf = 1
        else
        evalf = 0
    end

    while true
     evalf += 1 
     x_plus = x + 2.0^(-j) * (zk - x)   
     stptest = f(x_plus) - f(x) - η * 2.0^(-j) * dot(∇f(x), zk - x) 
        
        if stptest > 0.0   
           j += 1 
        else
           γ = 2.0^(-j)
            if γ < min_step
               ierror = 1
               println("Step length too small!")
               break
            end
            break  
        end       
    end
    
    β = (-dot(zk - x, ∇f(x)) * β^2) / (2.0 * (f(x + β * (zk - x)) - f(x) - β * dot(zk - x, ∇f(x))))
    if β < β1
       β = β1
       elseif β > β2
       β = β2
    end
    return (γ, β, ierror, evalf)
end

function method2(x, f, ∇f, ε, max_iter, GPA1, projection)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    stepsizes_γ = Float64[]
    evalf_γ = Float64[]
    projeval = Float64[]
    iteration_time = Float64[]
    
    # Initialization
    ierror = 0
    iter = 0
    projeval = 0
    seqx = x
    t0 = time()

    if norm(x - projection(x - ∇f(x))) < ε
        info = 0
        et = time - t0
        evalf_γ = 0
        println("x0 is a stationary point!")
        return (x, f(x), info, et, ierror, seqx, evalf_γ)
    end

    while true
        xk = copy(x)
        it0 = time()
        ∇fx = ∇f(x)
        fx = f(x) 
        gradnorm = norm(∇fx)
        (γ, β, ierror, evalf) = PG2(x, f, ∇f, projection, η, min_step, γ_start, β2)
        
        if ierror == 1
            break
        end 

        z = projection(x - β * ∇f(x)) 
        projeval += 1
        x = x + γ * (z - x)
        seqx = [seqx x] 
        it = time() - it0
        
        push!(iteration_time, it)
        push!(stepsizes_β, β)
        push!(stepsizes_γ, γ)
        push!(evalf_γ, evalf)
        push!(fvals, fx) 
        push!(gradnorms, gradnorm)
        
        # First stopping condition
        if norm(x - xk) < ε
           println("The solution has found!")
           break
        end

        # Update iteration
        iter += 1

        # Second stopping condition
        if iter > max_iter
            println("Maximum of iterations was achieved! Stopping...")
            ierror = 2
            break
        end
    end
    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, stepsizes_γ = stepsizes_γ, evalf_γ = evalf_γ, iteration_time = iteration_time)
    et = time() - t0
    println("iter = $iter  f(x) = $(f(x)) tempo = $et evalf = $(sum(evalf_γ)) evalproj = $projeval") 
    return (x, info, et, ierror, seqx, sum(evalf_γ), projeval)
end