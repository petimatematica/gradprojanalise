# Performances 

using CUTEst, NLPModels, LinearAlgebra, DataFrames, Random, Printf, Plots, BenchmarkProfiles, JLD2

include("SPG.jl")
include("PG1.jl")
include("PG2.jl")

V = Float64[]
T = Float64[]
S = Float64[]
G = Float64[]

# Parameters

η = 1.e-4 
ε = 1.e-5 
β_start = 1.0 
β1 = 0.01
β2 = 0.9
γ_start = 1.0
min_step = 1.e-5
max_iter = 40000
lambda_min = 1.e-30
lambda_max = 1.e+30
M = 10
sigma1 = 0.1
sigma2 = 0.9

problems = ["TORSION1", "TORSION2", "TORSION3", "TORSION4", "TORSION5", "TORSION6", "TORSIONA", "TORSIONB", "TORSIONC", "TORSIOND", "TORSIONE", "TORSIONF"] 
dim_parameter = ["2", "5"]
strategies = ["SPG1", "SPG2", "PG1", "PG2"]

t = time()

for strategy in strategies
    for problem in problems
        for D in dim_parameter
            nlp = CUTEstModel(problem, "-param", "Q=$D")
            
            # Initial guess from CUTEst
            x0 = nlp.meta.x0
            println("x0 =", x0)
            global x0
        
            # Objective functions from CUTEst
            global function f(x)
                return obj(nlp, x) 
            end
            
            # Gradient of Objective function from CUTEst
            global function ∇f(x)
                return grad(nlp, x)
            end

            # Upper and lower bounds setting
            l = Array{Float64}(undef,size(x0))
            u = Array{Float64}(undef,size(x0))
            for i in 1 : size(x0,1)
                global l[i] = -100.0
                global u[i] = 50.0
            end

            # Orthogonal projection
            global function proj(x)
                n = size(x,1)
                z = Array{Float64}(undef,size(x0))
            
                for i in 1:n
                    z[i] = max(l[i],min(x[i],u[i]))
                end
                return z
            end

            t0 = time()
            println("Running test with: strategy = $strategy, problem = $problem, dim_parameter = $D")

            if strategy == "SPG1"
            (xbest, info, et, ierror, seqx, evalsf, evalsproj) = spg(x0, f, ∇f, proj, ε, max_iter, lambda_min, lambda_max, M, sigma1, sigma2, η, spg1)
            elseif strategy == "SPG2"                 
            (xbest, info, et, ierror, seqx, evalsf, evalsproj) = spg(x0, f, ∇f, proj, ε, max_iter, lambda_min, lambda_max, M, sigma1, sigma2, η, spg2)
            elseif strategy == "PG1"
            (xbest, info, et, ierror, seqx, evalsf, evalsproj) = method1(x0, f, ∇f, ε, max_iter, PG1, proj)
            elseif strategy == "PG2"
            (xbest, info, et, ierror, seqx, evalsf, evalsproj) = method2(x0, f, ∇f, ε, max_iter, PG2, proj)
            end

            filename = "echo/" * D * problem * strategy * string("ierror", ierror) * string("min", f(xbest)) * ".jld2"
            @save filename info 

            if ierror > 0
                push!(V, Inf)
                push!(T, Inf)
                push!(S, Inf)
                push!(G, Inf)
            else
                iters = size(seqx, 2)
                push!(V, iters)
                push!(T, et)
                push!(S, evalsf)
                push!(G, evalsproj)
            end   

            #ENV["LINES"] = 10000
            #println(info)
            # println("Minimum value of f: ", f(xbest))
            # println("Total time spent: ", et)
            #println("Function evaluations = ", evalsf)
            # println("Number of iterations: ", iters)
            # println("Ierror = ", ierror)

            finalize(nlp)
        end    
    end
end

t_final = time() - t
println("Total time spent = ", t_final/60)
println("Problems tested, generating performance profiles...")

h=length(problems) * length(dim_parameter);
W=[V[1:h] V[h+1:2h] V[2h+1:3h] V[3h+1:4h]]; #Matrix which stores iterations
Z=[T[1:h] T[h+1:2h] T[2h+1:3h] T[3h+1:4h]]; #Matrix which stores CPU time
R=[S[1:h] S[h+1:2h] S[2h+1:3h] S[3h+1:4h]]; #Matrix which stores function evaluation
E=[G[1:h] G[h+1:2h] G[2h+1:3h] G[3h+1:4h]]; #Matrix which stores projection evaluation

colors=[:blue, :green2, :red, :cyan4]

X = performance_profile(PlotsBackend(), W, ["SPG1", "SPG2", "PG1", "PG2"], xlabel="Number of iterations", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)
Y = performance_profile(PlotsBackend(), Z, ["SPG1", "SPG2", "PG1", "PG2"], xlabel="CPU time ratio", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)
Q = performance_profile(PlotsBackend(), R, ["SPG1", "SPG2", "PG1", "PG2"], xlabel="Function evaluation", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)
N = performance_profile(PlotsBackend(), E, ["SPG1", "SPG2", "PG1", "PG2"], xlabel="Projection evaluation", ylabel="Solved problems [%]", legend=:bottomright, palette=colors, lw=1.5, dpi=1000)

println("Performance profiles generated, saving figures...")

p = plot(X)
savefig(p, "performanceprofileiters")

q = plot(Y)
savefig(q,"performanceprofiletime")

r = plot(Q)
savefig(r,"performanceprofileevalf")

n = plot(N)
savefig(n,"performanceprofileevalproj")

println("Figures saved, the code has finished running.")