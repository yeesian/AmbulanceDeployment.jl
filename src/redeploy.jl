abstract DeployModel
    # Interface
    # =========
    # assignment::Vector{Int} # which location the ambulance is assigned to
    # status::Vector{Int} # the current status of the ambulance
    # fromtime::Vector{Int} # the time it started the new status

type NoRedeployModel <: DeployModel
    model::JuMP.Model
    w::Matrix{JuMP.Variable}
    eta::Matrix{JuMP.Variable}
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
    JuMP.@variable(m, eta[1:nambulances, 1:nlocations] >= 0)
    JuMP.@constraint(m, [a=1:nambulances], sum(w[a,i] for i in 1:nlocations) == 1)
    JuMP.@constraint(m, [a=1:nambulances, i=1:nlocations], eta[a,i] >= w[a,i] - (assignment[a] == i))
    JuMP.@constraint(m, [a=1:nambulances, i=1:nlocations], eta[a,i] >= - (w[a,i] - (assignment[a] == i)))

    NoRedeployModel(m, w, eta, lambda, hosp2stn, stn2stn, assignment, ambulances, status, fromtime, hospital)
end

function reassign_ambulances!(problem::DispatchProblem, redeploy::DeployModel, t::Int)
    # (1) Optimize Dynamic Assignment Problem
    cost(a::Int,i::Int) =
        if redeploy.status[a] == :available
            100 + redeploy.stn2stn[redeploy.assignment[a],i]
        elseif redeploy.status[a] == :responding
            50 # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :atscene
            40 # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :conveying
            30 # + redeploy.hosp2stn[redeploy.hospital[a],i]
        elseif redeploy.status[a] == :returning
            10
        end
    JuMP.setobjective(redeploy.model, :Min, sum(
        cost(a,i)*redeploy.w[a,i] + redeploy.lambda*redeploy.eta[a,i]
        for a in eachindex(redeploy.assignment), i in eachindex(redeploy.ambulances)
    ))
    JuMP.solve(redeploy.model)

    # (2) Reassign ambulances based on JuMP.getvalue(redeploy.w)
    for a in 1:size(redeploy.w,1), i in 1:size(redeploy.w,2)
        if JuMP.getvalue(redeploy.w[a,i]) > 0.5
            if redeploy.assignment[a] != i # redeploy an existing ambulance
                if redeploy.status[a] == :available
                    @assert problem.available[redeploy.assignment[a]] > 0
                    problem.available[i] += 1; problem.available[redeploy.assignment[a]] -= 1
                    t_end = t + ceil(Int, 60*redeploy.stn2stn[redeploy.assignment[a],i])
                    enqueue!(ems.eventqueue, (:done, 0, t_end, a), t_end)
                    redeploying_to!(redeploy, a, t)
                end
            end
            redeploy.assignment[a] = i
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

function redeploying_to!(redeploy::DeployModel, amb::Int, t::Int)
    @assert in(amb, redeploy.ambulances)
    deleteat!(redeploy.ambulances, findfirst(redeploy.ambulances, amb))
    @assert !in(amb, redeploy.ambulances)
    redeploy.status[amb] = :redeploying
    redeploy.fromtime[amb] = t
end

function returned_to!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :returning "$amb: $(redeploy.status[amb])"
    redeploy.hospital[amb] = 0
    redeploy.status[amb] = :available
    redeploy.fromtime[amb] = t
    @assert !in(amb, redeploy.ambulances[redeploy.assignment[amb]])
    push!(redeploy.ambulances[redeploy.assignment[amb]], amb)
end

function redirected!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :returning "$amb: $(redeploy.status[amb])"
    redeploy.hospital[amb] = 0
    redeploy.status[amb] = :responding
    redeploy.fromtime[amb] = t
end
