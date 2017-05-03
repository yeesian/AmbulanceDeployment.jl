abstract DeployModel
    # Interface
    # =========
    # assignment::Vector{Int} # which location the ambulance is assigned to
    # status::Vector{Int} # the current status of the ambulance
    # fromtime::Vector{Int} # the time it started the new status

type NoRedeployModel <: DeployModel
    model::JuMP.Model
    w::Matrix{JuMP.Variable}
    eta1::Vector{JuMP.Variable}
    eta2::Matrix{JuMP.Variable}
    shortfall::Vector{JuMP.ConstraintRef}
    lambda::Float64

    hosp2stn::Matrix{Float64}
    stn2stn::Matrix{Float64}

    assignment::Vector{Int} # which location the ambulance is assigned to
    ambulances::Vector{Vector{Int}} # list of ambulances assigned to each location
    status::Vector{Symbol} # the current status of the ambulance
        # possible statuses: :available, :responding, :atscene, :conveying, :returning
    fromtime::Vector{Int} # the time it started the new status
    hospital::Vector{Int} # the hospital the ambulance is at (0 otherwise)
end

function NoRedeployModel(
        p::DeploymentProblem,
        available::Vector{Int},
        hospitals::DataFrame,
        stations::DataFrame;
        lambda::Float64 = 100.0,
        solver=GurobiSolver(OutputFlag=0)
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

    m = JuMP.Model(solver=solver)
    JuMP.@variable(m, w[1:nambulances, 1:nlocations], Bin)
    JuMP.@variable(m, eta1[1:nlocations] >= 0)
    JuMP.@variable(m, eta2[1:nambulances, 1:nlocations] >= 0)
    
    JuMP.@constraint(m, [a=1:nambulances], sum(w[a,i] for i in 1:nlocations) == 1)
    JuMP.@constraint(m, shortfall[i=1:nlocations], eta1[i] >= available[i] - sum(w[a,i] for a in 1:nambulances))
    JuMP.@constraint(m, [a=1:nambulances, i=1:nlocations], eta2[a,i] >= w[a,i] - (assignment[a] == i))
    JuMP.@constraint(m, [a=1:nambulances, i=1:nlocations], eta2[a,i] >= - (w[a,i] - (assignment[a] == i)))
    
    JuMP.@objective(m, Min,
        sum(eta1[i] for i in 1:nlocations) +
        lambda * sum(eta2[a,i] for a in 1:nambulances, i in 1:nlocations)
    )
    JuMP.build(m)

    NoRedeployModel(m, w, eta1, eta2, shortfall, lambda, hosp2stn, stn2stn, assignment, ambulances, status, fromtime, hospital)
end

function reassign_ambulances!(ems, problem::DispatchProblem, redeploy::DeployModel, t::Int)
    # (1) Optimize Dynamic Assignment Problem
    cost(a::Int,i::Int) =
        if redeploy.status[a] == :available
            redeploy.stn2stn[redeploy.assignment[a],i]
        elseif redeploy.status[a] == :responding
            50 # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :atscene
            40 # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :conveying
            30 # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :returning
            15 - redeploy.fromtime[a]
        elseif redeploy.status[a] == :redeploying
            t - redeploy.fromtime[a]
        end
    
    for s in redeploy.shortfall, a in 1:size(redeploy.w,1), i in 1:size(redeploy.w,2)
        chg_coeffs!(redeploy.model.internalModel.inner, [s.idx], [redeploy.w[a,i].col], [cost(a,i)/60-1])
    end
    JuMP.solve(redeploy.model)

    # (2) Reassign ambulances based on JuMP.getvalue(redeploy.w)
    for a in 1:size(redeploy.w,1), i in 1:size(redeploy.w,2)
        if JuMP.getvalue(redeploy.w[a,i]) > 0.5
            if redeploy.assignment[a] != i # redeploy an existing ambulance
                if redeploy.status[a] == :available
                    @assert problem.available[redeploy.assignment[a]] > 0
                    problem.available[redeploy.assignment[a]] -= 1
                    t_end = t + ceil(Int, 0*60*redeploy.stn2stn[redeploy.assignment[a],i])
                    redeploying_to!(redeploy, a, i, t)
                    enqueue!(ems.eventqueue, (:done, 0, t_end, a), t_end)
                end
            end
        end
    end
end

function respond_to!(redeploy::DeployModel, i::Int, t::Int)
    @assert length(redeploy.ambulances[i]) > 0 "$i: $(redeploy.ambulances[i])"
    amb = shift!(redeploy.ambulances[i])
    # @assert redeploy.hospital[amb] == 0
    @assert amb != 0
    @assert redeploy.status[amb] == :available "$amb: $(redeploy.status[amb])"
    redeploy.status[amb] = :responding
    redeploy.fromtime[amb] = t
    amb
end

function arriveatscene!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :responding "$amb: $(redeploy.status[amb])"
    @assert redeploy.hospital[amb] == 0
    redeploy.status[amb] = :atscene
    redeploy.fromtime[amb] = t
end

function conveying!(redeploy::DeployModel, amb::Int, hosp::Int, t::Int)
    @assert redeploy.status[amb] == :atscene "$amb: $(redeploy.status[amb])"
    @assert redeploy.hospital[amb] != 0
    redeploy.status[amb] = :conveying
    redeploy.fromtime[amb] = t
    redeploy.hospital[amb] = hosp
end

function returning_to!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :conveying "$amb: $(redeploy.status[amb])"
    @assert redeploy.hospital[amb] != 0
    redeploy.status[amb] = :returning
    redeploy.fromtime[amb] = t
    redeploy.assignment[amb]
end

function redeploying_to!(redeploy::DeployModel, amb::Int, i::Int, t::Int)
    ambulances = redeploy.ambulances[redeploy.assignment[amb]]
    @assert !in(amb, redeploy.ambulances[i])
    @assert in(amb, ambulances)
    deleteat!(ambulances, findfirst(ambulances, amb))
    @assert !in(amb, ambulances)
    redeploy.assignment[amb] = i
    redeploy.status[amb] = :redeploying
    redeploy.fromtime[amb] = t
end

function returned_to!(redeploy::DeployModel, amb::Int, t::Int)
    @assert in(redeploy.status[amb], (:returning, :redeploying)) "$amb: $(redeploy.status[amb])"
    redeploy.hospital[amb] = 0
    redeploy.status[amb] = :available
    redeploy.fromtime[amb] = t
    @assert !in(amb, redeploy.ambulances[redeploy.assignment[amb]])
    push!(redeploy.ambulances[redeploy.assignment[amb]], amb)
end

function redirected!(redeploy::DeployModel, amb::Int, t::Int)
    @assert in(redeploy.status[amb], (:returning, :redeploying)) "$amb: $(redeploy.status[amb])"
    redeploy.hospital[amb] = 0
    redeploy.status[amb] = :responding
    redeploy.fromtime[amb] = t
end
