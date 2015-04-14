type Gamma
    _single::Vector{Int}
    _local::Vector{Int}
    _regional::Vector{Int}
    _global::Int
end

function Gamma(p::DeploymentProblem; α::Float64=0.01)
    γ_single = vec(maximum(p.demand,1) + 1*(maximum(p.demand,1) .== 0))
    γ_local = [quantile(Poisson(mean(sum(p.demand[:,vec(p.adjacency[i,:])], 2))),1-α) for i=1:p.nregions]
    γ_regional = [quantile(Poisson(mean(sum(p.demand[:,p.coverage[:,i]],2))),1-α) for i in 1:p.nlocations]
    γ_global = quantile(Poisson(mean(sum(p.demand,2))),1-α)
    Gamma(γ_single,γ_local,γ_regional,γ_global)
end

function Qrobust{T <: Real}(problem::DeploymentProblem, x::Vector{T}, LB::Float64; α::Float64=0.01,
    solver=GurobiSolver(PrePasses=3))#, TimeLimit=10))

    γ = Gamma(problem, α=α)
    upp_bound = maximum(γ._single)

    I = 1:problem.nlocations
    J = 1:problem.nregions
    K = 1:upp_bound

    m = Model(solver=solver)

    @defVar(m, 1 >= p[1:problem.nregions] >= 0)
    @defVar(m, q[1:problem.nlocations] >= 0)
    @defVar(m, τ[1:problem.nregions,1:upp_bound], Bin)
    @defVar(m, f[1:problem.nregions,1:upp_bound], Bin)

    for j in J, k in K
        @addConstraint(m, τ[j, k] <= p[j])
        @addConstraint(m, τ[j, k] <= f[j, k])
    end

    for i in I, j in J
        if problem.coverage[j, i]
            @addConstraint(m, p[j] <= q[i])
        end
    end

    # Uncertainty

    for region in J
        @addConstraint(m, sum{f[region, k], k in K} <= γ._single[region])
        adjacent_regions = filter(j->problem.adjacency[j,region],J)
        @addConstraint(m, sum{f[j, k], j in adjacent_regions, k in K} <= γ._local[region])
    end

    for location in I
        covered_regions = filter(j->problem.coverage[j,location],J)
        @addConstraint(m, sum{f[j, k], j in covered_regions, k in K} <= γ._regional[location])
    end

    @addConstraint(m, sum{f[j, k], j in J, k in K} <= γ._global)

    @setObjective(m, Max, sum{τ[j, k], j in J, k in K} - sum{x[i]*q[i], i in I})

    JuMP.solve(m)

    return (getObjectiveValue(m), vec(sum(getValue(τ)[:,:], 2)))
end

type RobustDeployment{T <: Real} <: DeploymentModel
    m::Model
    x::Vector{Int}
    scenarios::Vector{Vector{T}}
    upperbounds::Vector{Float64}
    lowerbounds::Vector{Float64}
end # module

function naive_solution(p::DeploymentProblem)
    x = zeros(Int, p.nlocations)
    for i in p.nambulances
        x[i % p.nlocations] += 1
    end
    x
end

function RobustDeployment(p::DeploymentProblem; eps=1, tol=1e-6,
    iteration_limit=100, solver=GurobiSolver(OutputFlag=0,PrePasses=3))

    initial_x = naive_solution(p)

    scenarios = Vector{Int}[]
    upperbounds = Float64[]
    lowerbounds = Float64[]

    I = 1:p.nlocations
    J = 1:p.nregions

    m = Model(solver=solver)
    @defVar(m, x[1:p.nlocations] >= 0, Int)
    @defVar(m, η >= 0)
    y = Matrix{Variable}[]
    z = Vector{Variable}[]

    function add_scenario{T <: Real}(scenario::Vector{T}, l)

        # Create variables yˡ
        push!(y, Array(Variable, (p.nlocations,p.nregions)))
        for i in I, j in J
            y[l][i,j] = Variable(m, 0, p.nambulances, :Int, "y[$i,$j,$l]")
        end
        push!(z, Array(Variable, p.nregions))
        for j in J
            z[l][j] = Variable(m, 0, Inf, :Int, "z[$j,$l]")
        end
        
        # Add the following Constraints

        # (1) η >= 1ᵀ(dˡ + Bᴶyˡ)^+
        @addConstraint(m, η >= sum{z[l][j], j=J} + tol*sum{y[l][i,j], i=I, j=J})
        for i in I # flow constraints at each station
            @defExpr(outflow, sum{y[l][i,j], j in filter(j->p.coverage[j,i], J)})
            @addConstraint(m, x[i] >= outflow)
        end

        # (2) yˡ ∈ Y(x)
        for j in J # shortfall from satisfying demand/calls
            @defExpr(inflow, sum{y[l][i,j], i in filter(i->p.coverage[j,i], I)})
            @addConstraint(m, z[l][j] >= scenario[j] - inflow)
        end
    end

    # Initial Restricted Master Problem
    k = 1
    UB, scenario = Qrobust(p, initial_x, 0.0) # generate an initial scenario

    @setObjective(m, Min, η)
    @addConstraint(m, sum{x[i], i=I} <= p.nambulances)
    add_scenario(scenario, k)
    # @assert length(y) == k
    status = JuMP.solve(m)
    @assert status == :Optimal
    LB = max(0, getObjectiveValue(m)) # Initial LB
    UB, scenario = Qrobust(p, initial_x, LB) # Initial UB

    while k < iteration_limit
        println("iteration $k: LB $LB, UB $UB")
        if abs(UB - LB) < eps
            k = iteration_limit # terminate
        else
            # Solve Q(xⁱ), and form d = τ1 + ... + τk
            # (1) scenario          := d
            # (2) worst_shortfall   := Q(xⁱ)
            k += 1
            println("  solving Q with $(map(int,getValue(x)[:]))")
            worst_shortfall, scenario = Qrobust(p, getValue(x)[:], LB)
            println("  Q solved")
            push!(scenarios, scenario)
            
            if worst_shortfall < UB # Update UB = min{UB, Q(xⁱ)}
                UB = min(UB, worst_shortfall)
                deployment = getValue(x)[:] # keep track of best deployment plan
            end
            push!(upperbounds, UB) # for tracking convergence later
            
            add_scenario(scenario, k)
            status = JuMP.solve(m)
            @assert status == :Optimal
            LB = getObjectiveValue(m)
            push!(lowerbounds, LB) # for tracking convergence later
        end
    end
    RobustDeployment(m, map(int,getValue(x)[:]), scenarios, upperbounds, lowerbounds)
end