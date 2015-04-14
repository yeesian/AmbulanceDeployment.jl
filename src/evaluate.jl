function within{T <: Real}(scenario::Vector{T}, gamma::Gamma, p::DeploymentProblem)
    # returns true if scenario ∈ gamma; false otherwise
    if any(scenario .> gamma._single) 
        return false
    elseif sum(scenario) > gamma._global
        return false
    else
        for j in 1:p.nregions
            if sum(scenario[vec(p.adjacency[j,:])]) > gamma._local[j]
                return false
            end
        end
        for i in 1:p.nlocations
            if sum(scenario[p.coverage[:,i]]) > gamma._regional[i]
                return false
            end
        end
    end
    true
end

type Evaluation
    objvalue::Int
    shortfall::Vector{Int}
    dispatch::Vector{Int}
end

function evaluate{T1, T2 <: Real}(
    x::Vector{T1},
    scenario::Vector{T2},
    p::DeploymentProblem;
    solver=GurobiSolver(OutputFlag=0,PrePasses=3))

    I = 1:p.nlocations
    J = 1:p.nregions

    m = Model(solver=solver)
    @defVar(m, y[1:p.nlocations,1:p.nregions] >= 0, Int)
    @defVar(m, z[1:p.nregions] >= 0, Int)

    @setObjective(m, Min, sum{z[j], j in J})

    # flow constraints at each station
    for i in I
        @defExpr(outflow, sum{y[i,j], j in filter(j->p.coverage[j,i], J)})
        @addConstraint(m, x[i] >= outflow)
    end

    # shortfall from satisfying demand/calls
    for j in J
        @defExpr(inflow, sum{y[i,j], i in filter(i->p.coverage[j,i], I)})
        @addConstraint(m, z[j,t] >= scenario[j] - inflow)
    end

    Evaluation(int(getObjectiveValue(m)),
               map(int,getValue(z)[:]),
               map(int,getValue(y)[:]))
end

evaluate{T1, T2 <: Real}(x::Vector{T1}, scenarios::Vector{Vector{T2}}) = Evaluation[evaluate(x,scene) for scene in scenarios]

