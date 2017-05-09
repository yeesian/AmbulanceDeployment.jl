type AssignmentModel <: RedeployModel
    model::Gurobi.Model
    lambda::Float64

    hosp2stn::Matrix{Float64}
    stn2stn::Matrix{Float64}

    assignment::Vector{Int} # which location the ambulance is assigned to
    ambulances::Vector{Vector{Int}} # list of ambulances assigned to each location
    status::Vector{Symbol} # the current status of the ambulance
        # possible statuses: :available, :responding, :atscene, :conveying, :returning
    fromtime::Vector{Int} # the time it started the new status
    hospital::Vector{Int} # the hospital the ambulance is at (0 otherwise)

    soln::Vector{Float64} # buffer for storing dynamic assignment solutions
end

function AssignmentModel(
        p::DeploymentProblem,
        available::Vector{Int},
        # utilization::Vector{Float64},
        hospitals::DataFrame,
        stations::DataFrame;
        lambda::Float64 = 100.0
    )

    nambulances = sum(available)
    nlocations = length(available)
    assignment = zeros(Int, nambulances)

    hosp2stn = convert(Matrix{Float64}, hcat([hospitals[Symbol("stn$(i)_min")] for i in 1:nlocations]...))
    stn2stn =  convert(Matrix{Float64}, hcat([stations[Symbol("stn$(i)_min")] for i in 1:nlocations]...))

    k = 1
    ambulances = [Int[] for i in 1:nlocations]
    for i in eachindex(available), j in 1:available[i]
        assignment[k] = i
        push!(ambulances[i], k)
        k += 1
    end
    @assert k == nambulances + 1

    @assert sum(length(a) for a in ambulances) == nambulances
    status = fill(:available, nambulances)
    fromtime = zeros(Int, nambulances)
    hospital = zeros(Int, nambulances)

    m = Gurobi.Model(Gurobi.Env(), "redeploy", :minimize)
    Gurobi.setparam!(m, "OutputFlag", 0)
    for a in 1:nambulances, i in 1:nlocations # w variables
        Gurobi.add_bvar!(m, 0.)
    end
    for i in 1:nlocations # eta1
        Gurobi.add_cvar!(m, 1., 0., Inf)
    end
    for a in 1:nambulances, i in 1:nlocations # eta2 variables
        Gurobi.add_cvar!(m, lambda, 0., Inf)
    end

    # η >= available[i] - sum(w[a,i] for a in 1:nambulances)
    #     reformulated to
    # sum(w[a,i] for a in 1:nambulances) + η >= available[i]
    for i in 1:nlocations
        Gurobi.add_constr!(m,
            [((1:nambulances)-1)*nlocations+i; nambulances*nlocations + i], # inds
            ones(nambulances + 1), # coeffs
            '>', Float64(available[i]))
    end
    # sum(w[a,i] for i in 1:nlocations) == 1       [a=1:nambulances]
    for a in 1:nambulances
        Gurobi.add_constr!(m,
            collect((a-1)*nlocations + (1:nlocations)), # inds
            ones(nlocations), # coeffs
            '=', 1.)
    end

    # eta2[a,i] >= |w[a,i] - (assignment[a] == i)|   [a=1:nambulances, i=1:nlocations]
    for a in 1:nambulances, i in 1:nlocations
        offset = (a-1)*nlocations + i
        inds = [(nambulances+1)*nlocations + offset, offset]
        Gurobi.add_constr!(m, inds, [1., -1.], '>', - Float64(assignment[a] == i))
        Gurobi.add_constr!(m, inds, [1., 1.], '>', Float64(assignment[a] == i))
    end

    AssignmentModel(m, lambda, hosp2stn, stn2stn, assignment, ambulances,
                    status, fromtime, hospital, zeros(nambulances*nlocations))
end

function reassign_ambulances!(ems, problem::DispatchProblem, redeploy::AssignmentModel, t::Int)
    nlocations = length(redeploy.ambulances)
    nambulances = length(redeploy.assignment)
    # DEBUG
    # for a in 1:nambulances
    #     @assert t >= redeploy.fromtime[a]
    # end
    cost(a,i) = if redeploy.status[a] == :available
            0 + redeploy.stn2stn[redeploy.assignment[a],i]
        elseif redeploy.status[a] == :responding
            55 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i] # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :atscene
            45 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i] # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :conveying
            30 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i] # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :returning
            15 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i]
        elseif redeploy.status[a] == :redeploying
            10 - (t - redeploy.fromtime[a])/60 + redeploy.stn2stn[redeploy.assignment[a],i]
        end
    # (1) Optimize Dynamic Assignment Problem
    # (1i) change RHS to account for wait_queue
    nwait = [length(problem.wait_queue[nbhd]) for nbhd in 1:size(problem.coverage,1)]
    Gurobi.set_dblattrarray!(redeploy.model, "RHS", 1, nlocations, [
        Float64(problem.deployment[i] + sum(
            nwait[nbhd] for nbhd in 1:size(problem.coverage,1)
            if problem.coverage[nbhd,i]
        ))
        for i in 1:nlocations
    ])
    # (1ii) change coeffs to account for traveling time
    let con = Cint[0], ind = Cint[0], val = Float64[0.0]
        for i in 1:nlocations, a in 1:nambulances
            # JuMP.setRHS(s, problem.deployment[loc] + length(problem.wait_queue[loc])) # include backlog from wait_queue
            con[1] = i
            ind[1] = (a-1)*nlocations + i
            val[1] = 1-max(min(60,cost(a,i)),0)/60
            Gurobi.chg_coeffs!(redeploy.model, con, ind, val)
        end
    end
    Gurobi.update_model!(redeploy.model)
    Gurobi.optimize(redeploy.model)
    # @show Gurobi.get_objval(redeploy.model)

    # (2) Reassign ambulances based on JuMP.getvalue(redeploy.w)
    Gurobi.get_dblattrarray!(redeploy.soln, redeploy.model, "X", 1)
    # @show [(1-max(min(60,cost(a,i)),0)/60)*soln[redeploy.w[a,i].col] for a in 1:size(redeploy.w,1), i in 1:size(redeploy.w,2)]
    for a in 1:nambulances
        if redeploy.status[a] == :available
            for i in 1:nlocations
                let stn = redeploy.assignment[a]
                    # redeploy an existing ambulance
                    if stn != i && redeploy.soln[(a-1)*nlocations + i] > 0.5
                        @assert problem.available[redeploy.assignment[a]] > 0
                        problem.available[stn] -= 1
                        t_end = t + ceil(Int, 0*60*redeploy.stn2stn[stn,i])
                        redeploying_to!(redeploy, a, i, t)
                        enqueue!(ems.eventqueue, (:done, 0, t_end, a), t_end)
                    end
                end
            end
        end
    end
end
