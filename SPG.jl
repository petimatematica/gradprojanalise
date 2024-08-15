#
# Spectral Projected Gradient Method (SPG)
#
function spg(x0, f, ∇f, proj, ε, max_iter, lambda_min, lambda_max, M, sigma1, sigma2, η, linesearch)

    stplen = Float64[];
    fvals = Float64[];
    Gnorm = Float64[];

    info = DataFrame()

    f_hist = Float64[]
    feval = Int64[]
    projeval = Int64[]

    x = copy(x0)
    x = proj(x)
    seqx=x
    
    push!(f_hist,f(x))

    iter = 0

    push!(stplen,NaN)

    gradf_x = ∇f(x)
    
    lambda = 1.0 / norm(proj(x-gradf_x)-x,Inf)

    fbest = f(x)
    xbest= x

    t0 = time()

    while true
        P_x_m_grad = proj(x-gradf_x)
        norm_P = norm(P_x_m_grad - x,2)
        norm_G1_inf = norm(P_x_m_grad - x,Inf)
        #println(norm_G1_inf)
        
        fx = f(x)
        push!(fvals,fx)
        if fbest > fx 
            fbest = fx
            xbest = x
        end

        normgradfx = norm(gradf_x,2)
        push!(Gnorm,normgradfx)

        if norm_G1_inf < ε

            info.stplen = stplen 
            info.fvals = fvals
            info.gradnorm = Gnorm

            error=0

            et = time() - t0
            
            evalf_k = sum(feval) + 1
            evalproj_k = sum(projeval)

            println("iter = $iter  norm_G1_inf = $norm_G1_inf  f(x) = $(f(x)) tempo = $et evalf = $evalf_k evalproj = $evalproj_k")
            println("Solutions has found!")

            return(xbest, info, et, error, seqx, evalf_k, evalproj_k)
        end

        # Update number of iterations
        iter += 1

        if iter == 5000
            println("5000 iteradas...")
        end
        
        if iter == 10000
            println("10000 iteradas...")
        end        

        if iter > max_iter
            println("Maximum of iterations was achieved! Stoping...")

            info.stplen = stplen 
            info.fvals = fvals
            info.gradnorm = Gnorm

            error=1 

            et = time() - t0
            
            evalf_k = sum(feval) + 1
            evalproj_k = sum(projeval)

            return(xbest, info, et, error, seqx, evalf_k, evalproj_k)
        end

        # Backtrackin routine
        (x, gradf_x, s, y, f_hist, alpha, et, evalf, evalproj) = linesearch(iter, lambda, x, gradf_x, f_hist, M, sigma1, sigma2, η)
        push!(stplen,alpha)
        push!(feval,evalf)
        push!(projeval,evalproj)

        b = dot(s,y)
        if b > 0.0
            a = dot(s,s)
            lambda = max(lambda_min,min(lambda_max,a/b))
        else
            lambda = lambda_max               
        end

        seqx=[seqx x];
    end
end

#
# Backtrackin routine SPG1
#
function spg1(k,lambda_k,x_k,gradf_x_k,f_hist, M, sigma1, sigma2, η)
    alpha = copy(lambda_k)
    t0 = time()
    evalf = 0
    evalproj = 0

    while true
        x_plus = proj(x_k - alpha * gradf_x_k)
        evalproj += 1  
        d_k = x_plus - x_k
        m_k = min(k,M-1)
        f_max = maximum(f_hist[end-m_k+1:end])
        #m_k = min(k-1,M-1)
        #f_max = maximum(f_hist[end-m_k:end])
        f_x_plus = f(x_plus)
        evalf += 1  
        test = f_x_plus > f_max + η * dot(x_plus - x_k,gradf_x_k)
        if ~test
            s_k = x_plus - x_k
            gradf_x_kp1 = ∇f(x_plus)
            y_k = gradf_x_kp1 - gradf_x_k
            push!(f_hist,f_x_plus)
            et = time() - t0
            #evalf += 1
            return (x_plus,gradf_x_kp1,s_k,y_k,f_hist,alpha,et,evalf,evalproj)
        else
            fxx = f(x_k + alpha * d_k)
            evalf += 1 
            if alpha <= 0.1
                alpha = alpha / 2.0
            else
                atemp = (- dot(d_k,gradf_x_k) * alpha^2) / (2.0 * (fxx - f_hist[end] - alpha * dot(d_k,gradf_x_k)))
                if atemp < sigma1 || atemp > sigma2 * alpha
                    atemp = alpha / 2.0
                end
                alpha = atemp
            end
        end 
    end
end

#
# Backtrackin routine SPG2
#
function spg2(k,lambda_k,x_k,gradf_x_k,f_hist, M, sigma1, sigma2, η)
    lambda = copy(lambda_k)
    alpha = 1
    t0 = time()
    evalf = 0
    d_k = proj(x_k - lambda * gradf_x_k) - x_k
    evalproj = 1
    
    while true
        x_plus = x_k + alpha * d_k
        m_k = min(k-1,M-1)
        f_max = maximum(f_hist[end-m_k:end])
        f_x_plus = f(x_plus)
        evalf += 1
        test = f_x_plus > f_max + η * alpha * dot(d_k,gradf_x_k)
        if ~test
            s_k = x_plus - x_k
            gradf_x_kp1 = ∇f(x_plus)
            y_k = gradf_x_kp1 - gradf_x_k
            push!(f_hist,f_x_plus)
            #evalf += 1
            et = time() - t0
            return (x_plus,gradf_x_kp1,s_k,y_k,f_hist,alpha,et,evalf,evalproj)
        else
            fxx = f(x_k + alpha * d_k)
            evalf += 1
            if alpha <= 0.1
                alpha = alpha / 2.0
            else
                atemp = (- dot(d_k,gradf_x_k) * alpha^2) / (2.0 * (fxx - f_hist[end] - alpha * dot(d_k,gradf_x_k)))
                if atemp < sigma1 || atemp > sigma2 * alpha
                    atemp = alpha / 2.0
                end
                alpha = atemp
            end
        end 
    end
end