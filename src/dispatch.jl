type LPDispatchGreedy <: DispatchModel
    m::JuMP.Model
    candidates::Vector{Vector{Int}}
    location::Vector{JuMP.ConstraintRef}
    y::Array{JuMP.Variable,2}
    wait_queue::Vector{Int} # number of people waiitng for ambulances
end

type LPDispatchRandom <: DispatchModel
    m::JuMP.Model
    candidates::Vector{Vector{Int}}
    location::Vector{JuMP.ConstraintRef}
    y::Array{JuMP.Variable,2}
    wait_queue::Vector{Int} # number of people waiitng for ambulances
end

type LPDispatchBacklog <: DispatchModel
    m::JuMP.Model
    candidates::Vector{Vector{Int}}
    location::Vector{JuMP.ConstraintRef}
    y::Array{JuMP.Variable,2}
    yb::Array{JuMP.Variable,2}
end

function LPDispatchGreedy(p::DeploymentProblem,
                          available::Vector{Int},
                          solver=GurobiSolver(OutputFlag=0),
                          tol=params.δ)
    demand = vec(mean(p.demand[p.train,:],1))
    I = 1:p.nlocations ; J = 1:p.nregions
    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, x[1:p.nlocations] >= 0)
    JuMP.@defVar(m, y[1:p.nlocations,1:p.nregions] >= 0)
    JuMP.@defVar(m, z[1:p.nregions] >= 0)
    JuMP.@setObjective(m, Min, sum{z[j], j=J} + tol*sum{y[i,j], i=I, j=J})
    location = Array(JuMP.ConstraintRef, p.nlocations)
    for i in I # flow constraints at each station
        JuMP.@defExpr(outflow, sum{y[i,j], j in filter(j->p.coverage[j,i], J)})
        location[i] = JuMP.@addConstraint(m, outflow <= available[i])
    end
    for j in J # shortfall from satisfying demand/calls
        JuMP.@defExpr(inflow, sum{y[i,j], i in filter(i->p.coverage[j,i], I)})
        JuMP.@addConstraint(m, z[j] >= demand[j] - inflow)
    end
    candidates = Array(Vector{Int}, p.nregions)
    for region in 1:p.nregions
        candidates[region] = I[vec(p.coverage[region,:])]
    end
    status = JuMP.solve(m)
    @assert status == :Optimal
    LPDispatchGreedy(m, candidates, location, y, zeros(Int, p.nlocations))
end

function LPDispatchRandom(p::DeploymentProblem,
                          available::Vector{Int},
                          solver=GurobiSolver(OutputFlag=0),
                          tol=params.δ)
    demand = vec(mean(p.demand[p.train,:],1))
    I = 1:p.nlocations ; J = 1:p.nregions
    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, x[1:p.nlocations] >= 0)
    JuMP.@defVar(m, y[1:p.nlocations,1:p.nregions] >= 0)
    JuMP.@defVar(m, z[1:p.nregions] >= 0)
    JuMP.@setObjective(m, Min, sum{z[j], j=J} + tol*sum{y[i,j], i=I, j=J})
    location = Array(JuMP.ConstraintRef, p.nlocations)
    for i in I # flow constraints at each station
        JuMP.@defExpr(outflow, sum{y[i,j], j in filter(j->p.coverage[j,i], J)})
        location[i] = JuMP.@addConstraint(m, outflow <= available[i])
    end
    for j in J # shortfall from satisfying demand/calls
        JuMP.@defExpr(inflow, sum{y[i,j], i in filter(i->p.coverage[j,i], I)})
        JuMP.@addConstraint(m, z[j] >= demand[j] - inflow)
    end
    candidates = Array(Vector{Int}, p.nregions)
    for region in 1:p.nregions
        candidates[region] = I[vec(p.coverage[region,:])]
    end
    status = JuMP.solve(m)
    @assert status == :Optimal
    LPDispatchRandom(m, candidates, location, y, zeros(Int, p.nlocations))
end

function LPDispatchBacklog(p::DeploymentProblem,
                           available::Vector{Int},
                           solver=GurobiSolver(OutputFlag=0),
                           tol=params.δ)
    demand = vec(mean(p.demand[p.train,:],1))
    I = 1:p.nlocations ; J = 1:p.nregions
    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, x[1:p.nlocations] >= 0)
    JuMP.@defVar(m, y[1:p.nlocations,1:p.nregions] >= 0)
    JuMP.@defVar(m, yb[1:p.nlocations,1:p.nregions] >= 0)
    JuMP.@defVar(m, z[1:p.nregions] >= 0)
    JuMP.@setObjective(m, Min, sum{z[j], j=J} + tol*sum{y[i,j]+yb[i,j], i=I, j=J})
    location = Array(JuMP.ConstraintRef, p.nlocations)
    for i in I # flow constraints at each station
        JuMP.@defExpr(outflow, sum{y[i,j]-yb[i,j], j in filter(j->p.coverage[j,i], J)})
        location[i] = JuMP.@addConstraint(m, outflow <= available[i])
    end
    for j in J # shortfall from satisfying demand/calls
        JuMP.@defExpr(inflow, sum{y[i,j]-yb[i,j], i in filter(i->p.coverage[j,i], I)})
        JuMP.@addConstraint(m, z[j] >= demand[j] - inflow)
    end
    candidates = Array(Vector{Int}, p.nregions)
    for region in 1:p.nregions
        candidates[region] = I[vec(p.coverage[region,:])]
    end
    status = JuMP.solve(m)
    @assert status == :Optimal
    LPDispatchBacklog(m, candidates, location, y, yb)
end

function update_ambulances!(model::LPDispatchGreedy, i::Int, delta::Int)
    constr = model.location[i]
    prev = JuMP.rhs(JuMP.LinearConstraint(constr))
    @assert prev >= 0
    if prev + delta < 0
        model.wait_queue[i] -= delta
    elseif model.wait_queue[i] > 0 && delta > 0
        model.wait_queue[i] -= delta
        @assert model.wait_queue[i] >= 0
    else
        JuMP.chgConstrRHS(constr, prev + delta)
        # @assert abs(JuMP.rhs(constr) - prev - delta) < 1e-6
        #println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
        status = JuMP.solve(model.m)
        if status != :Optimal
            # println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
            @assert status == :Optimal
        end
    end
end

function update_ambulances!(model::LPDispatchRandom, i::Int, delta::Int)
    constr = model.location[i]
    prev = JuMP.rhs(JuMP.LinearConstraint(constr))
    @assert prev >= 0
    if prev + delta < 0
        model.wait_queue[i] -= delta
    elseif model.wait_queue[i] > 0 && delta > 0
        model.wait_queue[i] -= delta
        @assert model.wait_queue[i] >= 0
    else
        JuMP.chgConstrRHS(constr, prev + delta)
        # @assert abs(JuMP.rhs(constr) - prev - delta) < 1e-6
        #println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
        status = JuMP.solve(model.m)
        if status != :Optimal
            # println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
            @assert status == :Optimal
        end
    end
end

function update_ambulances!(model::LPDispatchBacklog, i::Int, delta::Int)
    constr = model.location[i]
    prev = JuMP.rhs(JuMP.LinearConstraint(constr))
    JuMP.chgConstrRHS(constr, prev + delta)
    # @assert abs(JuMP.rhs(constr) - prev - delta) < 1e-6
    #println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
    status = JuMP.solve(model.m)
    if status != :Optimal
        # println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
        @assert status == :Optimal
    end
end

function available_for(model::LPDispatchGreedy, j::Int, problem::DispatchProblem)
    # when an emergency call arrives from region j
    candidates = model.candidates[j]
    yvalues = JuMP.getValue(model.y[candidates,j])
    for (yi,xi) in enumerate(candidates)
        if problem.available[xi] > 0
            yvalues[yi] = yvalues[yi] / problem.available[xi]
        end
    end
    yindices = sortperm(yvalues) # order the candidates in decreasing "desirability"
    for i in yindices
        if problem.available[candidates[i]] > 0 # send the most "desirable" one
            return candidates[i]
        end
    end
    # if no ambulance is available
    return 0
end

function available_for(model::LPDispatchRandom, j::Int, problem::DispatchProblem)
    candidates = model.candidates[j]
    yvalues = JuMP.getValue(model.y[candidates,j])
    if sum(yvalues) < 1e-6 # if the differences are too small
        for i in candidates
            if problem.available[i] > 0
                return i
            end
        end
        return 0
    else
        for (yi,xi) in enumerate(candidates)
            if problem.available[xi] > 0
                yvalues[yi] = yvalues[yi] / problem.available[xi]
            end
        end
        i = Distributions.sample(StatsBase.WeightVec(yvalues))
        return candidates[i]
    end
end

function available_for(model::LPDispatchBacklog, j::Int, problem::DispatchProblem)
    # when an emergency call arrives from region j
    candidates = model.candidates[j]
    fwd_values = JuMP.getValue(model.y[model.candidates[j],j])
    bwd_values = JuMP.getValue(model.yb[model.candidates[j],j])

    yvalues = fwd_values - bwd_values
    for (yi,xi) in enumerate(candidates)
        if problem.available[xi] > 0
            yvalues[yi] = yvalues[yi] / problem.available[xi]
        end
    end

    yindices = sortperm(-yvalues)  # order the candidates in decreasing "desirability"
    for i in yindices
        if problem.available[candidate[i]] > 0
            return candidate[i]
        end
    end
    # if no ambulance is available
    return 0
end

type RobustDispatch <: DispatchModel
    Q::Qrobust
    candidates::Vector{Vector{Int}}
    available::Vector{Int}
end

function RobustDispatch(p::DeploymentProblem,
                        available::Vector{Int},
                        solver=GurobiSolver(OutputFlag=0),
                        tol=params.δ)
    I = 1:p.nlocations
    candidates = Array(Vector{Int}, p.nregions)
    for region in 1:p.nregions
        candidates[region] = I[vec(p.coverage[region,:])]
    end
    RobustDispatch(Qrobust(p), candidates, available)
end

function update_ambulances!(model::RobustDispatch, i::Int, delta::Int)
    model.available[i] += delta
end

function available_for(model::RobustDispatch, j::Int, problem::DispatchProblem)
    qvalues = Vector{Float64}()
    for i in model.candidates[j]
        model.available[i] -= 1
        push!(qvalues, evaluate_objvalue(model.Q, model.available))
        model.available[i] += 1
    end
    for i in sortperm(qvalues)
        location = model.candidates[j][i]
        if problem.available[location] > 0 # send the most "desirable" one
            return location
        end
    end
    return 0 # if no ambulance is available
end

type StochasticDispatch <: DispatchModel
    m::JuMP.Model
    candidates::Vector{Vector{Int}}
    location::Vector{JuMP.ConstraintRef}
end

function StochasticDispatch(p::DeploymentProblem,
                            available::Vector{Int},
                            nperiods = params.nperiods,
                            solver = GurobiSolver(OutputFlag=0),
                            tol = params.δ)
    @assert sum(available) <= p.nambulances
    nperiods = min(length(p.train), nperiods)
    demand = vec(mean(p.demand[p.train,:],1))
    I = 1:p.nlocations
    J = 1:p.nregions
    T = 1:nperiods

    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, y[1:p.nlocations,1:p.nregions,1:nperiods] >= 0, Int)
    JuMP.@defVar(m, z[1:p.nregions,1:nperiods] >= 0, Int)

    JuMP.@setObjective(m, Min, sum{z[j,t], j=J, t=T} + tol*sum{y[i,j,t], i=I, j=J, t=T})

    location = Array(JuMP.ConstraintRef, p.nlocations)
    for t in T
        for i in I # flow constraints at each station
            JuMP.@defExpr(outflow, sum{y[i,j,t], j in filter(j->p.coverage[j,i], J)})
            location[i] = JuMP.@addConstraint(m, outflow <= available[i])
        end
        for j in J # shortfall from satisfying demand/calls
            JuMP.@defExpr(inflow, sum{y[i,j,t], i in filter(i->p.coverage[j,i], I)})
            JuMP.@addConstraint(m, z[j,t] >= demand[t,j] - inflow)
        end
    end

    candidates = Array(Vector{Int}, p.nregions)
    for j in J
        candidates[region] = I[vec(p.coverage[j,:])]
    end
    # status = JuMP.solve(m)
    # @assert status == :Optimal
    StochasticDispatch(m, candidates, location)
end

function update_ambulances!(model::StochasticDispatch, i::Int, delta::Int)
    constr = model.location[i]
    prev = JuMP.rhs(JuMP.LinearConstraint(constr))
    JuMP.chgConstrRHS(constr, prev + delta)
    # @assert abs(JuMP.rhs(constr) - prev - delta) < 1e-6
    #println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
    # status = JuMP.solve(model.m)
    # if status != :Optimal
    #     println("available: $([iround(JuMP.rhs(model.m.linconstr[loc])) for loc in 1:35])")
    #     @assert status == :Optimal
    # end
end

function available_for(model::StochasticDispatch, j::Int, problem::DispatchProblem)
    qvalues = Vector{Float64}()
    for i in model.candidates[j]
        constr = model.location[i]
        prev = JuMP.rhs(JuMP.LinearConstraint(constr))
        if prev > 0.5 # consider sending an ambulance only if it's available
            JuMP.chgConstrRHS(constr, prev - 1)
            status = JuMP.solve(model.m)
            @assert status == :Optimal
            push!(qvalues, JuMP.getObjectiveValue(model.m))
            JuMP.chgConstrRHS(constr, prev)
        else
            push!(qvalues, Inf)
        end
    end
    for i in sortperm(qvalues)
        location = model.candidates[j][i]
        if problem.available[location] > 0 # send the most "desirable" one
            return location
        end
    end
    return 0 # if no ambulance is available
end

type MEXCLPDispatch <: DispatchModel
    x::Vector{Int}  # nlocation
    candidates::Vector{Vector{Int}}
    q::Float64

    function MEXCLPDispatch(p::DeploymentProblem, available::Vector{Int}, q::Float64)
        I = 1:p.nlocations
        candidates = Array(Vector{Int}, p.nregions)
        for region in 1:p.nregions
            candidates[region] = I[vec(p.coverage[region,:])]
        end
        new(available, candidates, q)
    end
end

function update_ambulances!(model::MEXCLPDispatch, i::Int, delta::Int)
    model.x[i] += delta
end

"the math says you'll return the location with the highest number of ambulances"
function available_for(model::MEXCLPDispatch, j::Int, problem::DispatchProblem)
    location = 0
    max_x = 0
    for i in model.candidates[j]
        if model.x[i] > max_x
            location = i
            max_x = model.x[i]
        end
    end
    location
end

type MALPDispatch{BM <: AbstractMatrix{Bool}} <: DispatchModel
    m::JuMP.Model
    available::Vector{Int}
    coverage::BM    # nregion x nlocation
    region::Vector{JuMP.ConstraintRef}
end

function MALPDispatch(p::DeploymentProblem,
                      available::Vector{Int},
                      q::Float64; # busy fraction
                      α::Float64 = 0.99, # reliability level
                      solver = GurobiSolver(OutputFlag=0))
    demand = vec(mean(p.demand[p.train,:],1))
    @assert length(demand) == p.nregions
    
    I = 1:p.nlocations
    J = 1:p.nregions
    b = ceil(Int, log(1-α)/log(q))

    m = JuMP.Model(solver=solver)
    JuMP.@defVar(m, x[1:p.nlocations] >= 0, Int)
    JuMP.@defVar(m, z[1:p.nregions, 1:b], Bin)

    JuMP.@setObjective(m, Max, sum{demand[j]*z[j,b], j in J})

    JuMP.@addConstraint(m, sum{x[i], i in I} <= p.nambulances)
    region = Array(JuMP.ConstraintRef, p.nregions)
    for j in J # coverage over all regions
        region_j = sum([available[i] for i in filter(i->p.coverage[j,i], I)])
        region[j] = JuMP.@addConstraint(m, sum{z[j,k], k in 1:b} <= region_j)
        for k in 2:b
            JuMP.@addConstraint(m, z[j,k] <= z[j,k-1])
        end
    end

    MALPDispatch(m, available, p.coverage, region)
end

function update_ambulances!(model::MALPDispatch, i::Int, delta::Int)
    model.available[i] += delta
    @assert model.available[i] >= 0
    for j in 1:size(model.coverage,1)
        if model.coverage[j,i]
            constr = model.region[j]
            prev = JuMP.rhs(JuMP.LinearConstraint(constr))
            JuMP.chgConstrRHS(constr, prev + delta)
            # @assert abs(JuMP.rhs(model.m.linconstr[i]) - prev - delta) < 1e-6
        end
    end
end

function available_for(model::MALPDispatch, j::Int, problem::DispatchProblem)
    location = 0
    max_q = 0.0
    for i in 1:size(model.coverage,2)
        if model.coverage[j,i] && model.available[i] > 0
            update_ambulances!(model, i, -1)
            status = JuMP.solve(model.m)
            @assert status == :Optimal
            qvalue = JuMP.getObjectiveValue(model.m)
            update_ambulances!(model, i, 1)
            if qvalue > max_q
                location = i
                max_q = qvalue
            end
        end
    end
    location
end

type ClosestDispatch <: DispatchModel
    drivetime::DataFrame
    available::Vector{Int}
    candidates::Vector{Vector{Int}}
end

function ClosestDispatch(p::DeploymentProblem,
                         drivetime::DataFrame,
                         available::Vector{Int})
    candidates = Array(Vector{Int}, p.nregions)
    I = 1:p.nlocations
    for region in 1:p.nregions
        candidates[region] = I[vec(p.coverage[region,:])]
    end
    ClosestDispatch(drivetime, deepcopy(available), candidates)
end

function update_ambulances!(model::ClosestDispatch, i::Int, delta::Int)
    model.available[i] += delta
end

function available_for(model::ClosestDispatch, id::Int, problem::DispatchProblem)
    location = 0
    min_time = typemax(Int)
    for i in model.candidates[problem.emergency_calls[id, :neighborhood]]
        if problem.available[i] > 0 && model.drivetime[id, i] < min_time
            location = i
            min_time = model.drivetime[id, i]
        end
    end
    location
end

# function LP_simulate_events(deployment::Vector{Int},
#                             problem::DeploymentProblem,
#                             emergency_calls::DataFrames.DataFrame;
#                             verbose=false)
#     ncalls = DataFrames.nrow(emergency_calls)
#     coverage = problem.coverage
#     nlocations = problem.nlocations
#     nregions = problem.nregions
#     locations = 1:problem.nlocations
#     regions = 1:problem.nregions
#     delays = (Int,Int,Int,Int)[]
#     available = deployment[:]
#     online_p = LPDispatch(problem, deployment)
#     next_in_location = Array(Vector{(Int,Int)}, nlocations)
    
#     for i in locations
#         next_in_location[i] = Vector{(Int,Int)}()
#     end
#     events_queue = form_queue(emergency_calls)
#     while !isempty(events_queue)
#         (id, event, t, value) = dequeue!(events_queue)
#         if event == :call
#             j = value # incoming emergency call from region j
#             verbose && println("incoming emergency call from region $j at time $t")
#             delay = 0
#             i = LP_ambulance_for(j, available, online_p)
#             if i == 0 # no ambulances available at the moment
#                 verbose && println("  no ambulances available for $j at time $(t)!")
#                 (i,tnew) = =next_ambulance_for(j, next_in_location, coverage)

#                 if i==0 # no coverage for region j
#                     println("unable to serve one call from region $j")
#                     continue
#                 end

#                 @assert i != 0
#                 verbose && println("  dispatching ambulance from $i, from queue $(next_in_location[i])")
#                 shift!(next_in_location[i])

#                 @assert tnew >= t
#                 delay = tnew - t
#                 push!(delays, (id, delay, i, j))
#             else
#                 verbose && println("  dispatching ambulance from $i")
#                 update_ambulances(online_p, i, -1)
#                 available[i] -= 1
#             end
#             tend = t + delay + 3600 # 3600seconds = 1 hour
#             enqueue!(events_queue, (id+ncalls, :done, tend, i), tend)
#             push!(next_in_location[i], (id+ncalls,tend))
#             verbose && println("  with $(available[i]) remaining at $i, and a queue of $(next_in_location[i]).")
#         else
#             @assert event == :done
#             i = value
#             @assert i != 0
#             for (n,(index, tend)) in enumerate(next_in_location[i])
#                 @assert t <= tend
#                 if t == tend
#                     if id == index
#                         splice!(next_in_location[i], n)
#                         update_ambulances(online_p, i, 1)
#                         available[i] += 1
#                         verbose && println("1 ambulance returned to location $i at time $t")
#                         break
#                     end
#                 else
#                     @assert t < tend
#                     break
#                 end
#             end
#         end
#     end
#     delays
# end

