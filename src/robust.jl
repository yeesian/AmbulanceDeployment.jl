type Gamma
    _single::Vector{Int}
    _local::Vector{Int}
    _regional::Vector{Int}
    _global::Int
end

type Qrobust
    m::JuMP.Model
    I::UnitRange{Int}
    J::UnitRange{Int}
    d::Vector{JuMP.Variable}
    q::Vector{JuMP.Variable}
    γ::Gamma
end

type RobustDeployment <: DeploymentModel
    m::JuMP.Model
    Q::Qrobust

    I::UnitRange{Int}
    J::UnitRange{Int}

    x::Vector{JuMP.Variable}
    y::Vector{Matrix{JuMP.Variable}}
    z::Vector{Vector{JuMP.Variable}}
    η::JuMP.Variable

    scenarios::Vector{Vector{Int}}
    upperbounds::Vector{Float64}
    lowerbounds::Vector{Float64}
    deployment::Vector{Vector{Int}}
end
deployment(m::RobustDeployment) = m.deployment[end]

function Gamma(p::DeploymentProblem; α=params.α)
    demand = p.demand[p.train,:]
    γ_single = vec(maximum(demand,1) + 1*(maximum(demand,1) .== 0))
    γ_local = [quantile(Poisson(mean(sum(demand[:,vec(p.adjacency[i,:])], 2))),1-α) for i=1:p.nregions]
    γ_regional = [quantile(Poisson(mean(sum(demand[:,p.coverage[:,i]],2))),1-α) for i in 1:p.nlocations]
    γ_global = quantile(Poisson(mean(sum(demand,2))),1-α)
    Gamma(γ_single,γ_local,γ_regional,γ_global)
end

function Qrobust(problem::DeploymentProblem; α=params.α, verbose=false,
    solver=GurobiSolver(OutputFlag=0)) #, MIPGapAbs=0.9)) #, TimeLimit=30))
    if verbose
        solver=GurobiSolver(OutputFlag=1) #, MIPGapAbs=0.9)
    end
    γ = Gamma(problem, α=α)
    upp_bound = maximum(γ._single)
    I = 1:problem.nlocations
    J = 1:problem.nregions

    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, d[1:problem.nregions]>=0, Int)
    JuMP.@defVar(m, p[1:problem.nregions], Bin)
    JuMP.@defVar(m, q[1:problem.nlocations], Bin)

    for i in I, j in J
        problem.coverage[j, i] && JuMP.@addConstraint(m, p[j] <= q[i])
    end

    # Uncertainty
    for j in J
        JuMP.@addConstraint(m, d[j] <= γ._single[j]*p[j])
        adjacent_regions = filter(k->problem.adjacency[k,j],J)
        JuMP.@addConstraint(m, sum{d[k], k in adjacent_regions} <= γ._local[j])
    end
    for i in I
        covered_regions = filter(j->problem.coverage[j,i],J)
        JuMP.@addConstraint(m, sum{d[j], j in covered_regions} <= γ._regional[i])
    end
    JuMP.@addConstraint(m, sum{d[j], j in J} <= γ._global)

    Qrobust(m, I, J, d, q, γ)
end

function evaluate{T <: Real}(Q::Qrobust, x::Vector{T})
    JuMP.@setObjective(Q.m, Max, sum{Q.d[j], j in Q.J} - sum{x[i]*Q.q[i], i in Q.I})
    status = JuMP.solve(Q.m)
    JuMP.getObjectiveValue(Q.m), Int[round(Int,d) for d in JuMP.getValue(Q.d)]
end

function evaluate_objvalue{T <: Real}(Q::Qrobust, x::Vector{T})
    JuMP.@setObjective(Q.m, Max, sum{Q.d[j], j in Q.J} - sum{x[i]*Q.q[i], i in Q.I})
    status = JuMP.solve(Q.m)
    JuMP.getObjectiveValue(Q.m)
end

function RobustDeployment(p::DeploymentProblem; α=params.α, eps=params.ε, tol=params.δ,
    solver=GurobiSolver(OutputFlag=0, MIPGapAbs=0.9), verbose=false, master_verbose=false)
    if master_verbose
        solver=GurobiSolver(OutputFlag=1, MIPGapAbs=0.9)
    end
    I = 1:p.nlocations
    J = 1:p.nregions

    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, x[1:p.nlocations] >= 0, Int)
    JuMP.@defVar(m, η >= 0)
    y = Vector{Matrix{JuMP.Variable}}()
    z = Vector{JuMP.Variable}()

    # Initial Restricted Master Problem
    JuMP.@setObjective(m, Min, η)
    JuMP.@addConstraint(m, sum{x[i], i=I} <= p.nambulances)
    for j in J # coverage over all regions
        JuMP.@addConstraint(m, sum{x[i], i in filter(i->p.coverage[j,i], I)} >= 1)
    end

    RobustDeployment(m, Qrobust(p, α=α, verbose=verbose), I, J, x, y, z, η,
                     Vector{Vector{Int}}(), Vector{Float64}(), Vector{Float64}(), Vector{Int}[naive_solution(p)])
end

function add_scenario{T <: Real}(model::RobustDeployment, p::DeploymentProblem, scenario::Vector{T}; tol=params.δ)
    # Create variables yˡ
    push!(model.y, Array(JuMP.Variable, (p.nlocations,p.nregions)))
    l = length(model.y)
    for i in model.I, j in model.J
        model.y[l][i,j] = JuMP.Variable(model.m, 0, p.nambulances, :Int, UTF8String("y[$i,$j,$l]"))
    end
    push!(model.z, Array(JuMP.Variable, p.nregions))
    for j in model.J
        model.z[l][j] = JuMP.Variable(model.m, 0, Inf, :Int, UTF8String("z[$j,$l]"))
    end

    # (1) η >= 1ᵀ(dˡ + Bᴶyˡ)^+
    JuMP.@addConstraint(model.m, model.η >= sum{model.z[l][j], j=model.J} + tol*sum{model.y[l][i,j], i=model.I, j=model.J})
    for i in model.I # flow constraints at each station
        JuMP.@defExpr(outflow, sum{model.y[l][i,j], j in filter(j->p.coverage[j,i], model.J)})
        JuMP.@addConstraint(model.m, model.x[i] >= outflow)
    end
    # (2) yˡ ∈ Y(x)
    for j in model.J # shortfall from satisfying demand/calls
        JuMP.@defExpr(inflow, sum{model.y[l][i,j], i in filter(i->p.coverage[j,i], model.I)})
        JuMP.@addConstraint(model.m, model.z[l][j] >= scenario[j] - inflow)
    end
end

function solve(model::RobustDeployment, p::DeploymentProblem; verbose=false, maxiter=params.maxiter, eps=params.ε)
    LB = 0.0
    UB, scenario = evaluate(model.Q, model.deployment[end])
    push!(model.lowerbounds, LB)
    push!(model.upperbounds, UB)
    push!(model.scenarios, scenario)

    for k in 1:maxiter
        verbose && println("iteration $k: LB $LB, UB $UB")
        abs(UB - LB) < eps && break
        verbose && println("  solving Q with $(model.deployment)")

        add_scenario(model, p, scenario)
        status = JuMP.solve(model.m)
        @assert status == :Optimal

        LB = JuMP.getObjectiveValue(model.m)

        push!(model.deployment, [round(Int,x) for x in JuMP.getValue(model.x)])
        shortfall, scenario = evaluate(model.Q, model.deployment[end])
        UB = min(UB, shortfall)

        # for tracking convergence later
        push!(model.upperbounds, UB)
        push!(model.scenarios, scenario)
        push!(model.lowerbounds, LB)
    end
end