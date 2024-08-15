## Projected Gradient Method with PG1 ##

function PG1(x, f, ∇f, projection, η, min_step, β_start)
    ierror = 0
    β = β_start
    j = 0

    if x == x0
        evalf = 1
        else
        evalf = 0
    end

    projeval = 0
 
    while true
     evalf += 1 
     projeval += 1   
     zkj = projection(x - β * 2.0^(-j) * ∇f(x)) 
     stptest = f(zkj) - f(x) + η * dot(∇f(x), x - zkj) 
 
        if stptest > 0.0 
         j += 1 
        else
            β = β_start * 2.0^(-j)
            if β < min_step
               ierror = 1
               println("Step length too small!")
               break
            end
            break
        end
    end 
    return (β, ierror, evalf, projeval)
end

function method1(x0, f, ∇f, ε, max_iter, PG1, projection)

    fvals = Float64[]
    gradnorms = Float64[]
    stepsizes_β = Float64[]
    evalf_β = Float64[]
    iteration_time = Float64[]
    projevals = Float64[]
    
    # Initialization
    ierror = 0
    iter = 0
    x = x0
    seqx = x
    t0 = time()

    if norm(x - projection(x - ∇f(x))) < ε
        info = 0
        et = time() - t0
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
        (β, ierror, evalf, projeval) = PG1(x, f, ∇f, projection, η, min_step, β_start)

        if ierror == 1
            break
        end  

        x = projection(x - β * ∇fx)
        seqx = [seqx x]
        it = time() - it0
        push!(fvals, fx)
        push!(gradnorms, gradnorm)
        push!(iteration_time, it)
        push!(stepsizes_β, β) 
        push!(evalf_β, evalf)
        push!(projevals, projeval)
        
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

    info = DataFrame(fvals = fvals, gradnorms = gradnorms, stepsizes_β = stepsizes_β, evalf_β = evalf_β, iteration_time = iteration_time)
    et = time() - t0
    println("iter = $iter  f(x) = $(f(x)) tempo = $et evalf = $(sum(evalf_β)) projeval = $(sum(projevals))") 
    return (x, info, et, ierror, seqx, sum(evalf_β), sum(projevals))
end
