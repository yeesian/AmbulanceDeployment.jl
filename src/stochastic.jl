type StochasticDeployment <: DeploymentModel
    m::JuMP.Model
    x::Vector{JuMP.Variable}
end
deployment(m::StochasticDeployment) = [round(Int,x) for x in JuMP.getValue(m.x)]

function StochasticDeployment(p::DeploymentProblem; nperiods=params.nperiods, tol=params.δ,
    solver=GurobiSolver(OutputFlag=0))

    nperiods = min(length(p.train), nperiods)
    demand = p.demand[p.train,:]
    I = 1:p.nlocations
    J = 1:p.nregions
    T = 1:nperiods

    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, x[1:p.nlocations] >= 0, Int)
    JuMP.@defVar(m, y[1:p.nlocations,1:p.nregions,1:nperiods] >= 0, Int)
    JuMP.@defVar(m, z[1:p.nregions,1:nperiods] >= 0, Int)

    JuMP.@setObjective(m, Min, sum{z[j,t], j=J, t=T} + tol*sum{y[i,j,t], i=I, j=J, t=T})

    JuMP.@addConstraint(m, sum{x[i], i=I} <= p.nambulances)

    for j in J # coverage over all regions
        JuMP.@addConstraint(m, sum{x[i], i in filter(i->p.coverage[j,i], I)} >= 1)
    end

    # flow constraints at each station
    for i in I, t in T
        JuMP.@defExpr(outflow, sum{y[i,j,t], j in filter(j->p.coverage[j,i], J)})
        JuMP.@addConstraint(m, x[i] >= outflow)
    end

    # shortfall from satisfying demand/calls
    for j in J, t in T
        JuMP.@defExpr(inflow, sum{y[i,j,t], i in filter(i->p.coverage[j,i], I)})
        JuMP.@addConstraint(m, z[j,t] >= demand[t,j] - inflow)
    end

    StochasticDeployment(m, x)
end

solve(model::StochasticDeployment) = JuMP.solve(model.m)
solve(model::StochasticDeployment, p::DeploymentProblem) = JuMP.solve(model.m)