type StochasticDeployment <: DeploymentModel
    m::Model
    x::JuMP.JuMPArray{Variable,1}
end

function StochasticDeployment(  
    p::DeploymentProblem, nperiods::Int = 100;
    tol::Float64=1e-6,
    solver=GurobiSolver(OutputFlag=0,PrePasses=3))

    nperiods = min(size(p.demand, 1), nperiods)
    I = 1:p.nlocations
    J = 1:p.nregions
    T = 1:nperiods

    m = Model(solver=solver)
    @defVar(m, x[1:p.nlocations] >= 0, Int)
    @defVar(m, y[1:p.nlocations,1:p.nregions,1:nperiods] >= 0, Int)
    @defVar(m, z[1:p.nregions,1:nperiods] >= 0, Int)

    @setObjective(m, Min, sum{z[j,t], j=J, t=T} + tol*sum{y[i,j,t], i=I, j=J, t=T})

    @addConstraint(m, sum{x[i], i=I} <= p.nambulances)

    # flow constraints at each station
    for i in I, t in T
        @defExpr(outflow, sum{y[i,j,t], j in filter(j->p.coverage[j,i], J)})
        @addConstraint(m, x[i] >= outflow)
    end

    # shortfall from satisfying demand/calls
    for j in J, t in T
        @defExpr(inflow, sum{y[i,j,t], i in filter(i->p.coverage[j,i], I)})
        @addConstraint(m, z[j,t] >= p.demand[t,j] - inflow)
    end

    StochasticDeployment(m, x)
end