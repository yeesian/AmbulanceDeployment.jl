abstract DeployModel
    # Interface
    # =========
    # assignment::Vector{Int} # which location the ambulance is assigned to
    # status::Vector{Int} # the current status of the ambulance
    # fromtime::Vector{Int} # the time it started the new status

type NoRedeployModel <: DeployModel
    assignment::Vector{Int} # which location the ambulance is assigned to
    ambulances::Vector{Vector{Int}} # list of ambulances assigned to each location
    status::Vector{Symbol} # the current status of the ambulance
        # possible statuses: :available, :responding, :atscene, :conveying, :returning
    fromtime::Vector{Int} # the time it started the new status
    hospital::Vector{Int} # the hospital the ambulance is at (0 otherwise)
end

function NoRedeployModel(p::DeploymentProblem, available::Vector{Int})
    nambulances = sum(available)
    assignment = zeros(Int, nambulances)
    k = 1
    ambulances = [Int[] for i in 1:nambulances]
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
    NoRedeployModel(assignment, ambulances, status, hospital)
end

reassign_ambulances!(redeploy::DeployModel) = nothing

function respond_to!(redeploy::DeployModel, i::Int, t::Int)
    @assert length(redeploy.ambulances[i]) > 0
    amb = shift!(redeploy.ambulances[i])
    redeploy.status[amb] = :responding
    redeploy.fromtime[amb] = t
    amb
end

function arriveatscene!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :responding
    redeploy.status[amb] = :atscene
    redeploy.fromtime[amb] = t
end

function conveying!(redeploy::DeployModel, amb::Int, hosp::Int, t::Int)
    @assert redeploy.status[amb] == :atscene
    redeploy.status[amb] = :conveying
    redeploy.fromtime[amb] = t
    redeploy.hospital[amb] = hosp
end

function returning_to!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :conveying
    redeploy.status[amb] = :returning
    redeploy.fromtime[amb] = t
    redeploy.assignment[amb]
    redeploy.hospital[amb] = 0
end

function returned_to!(redeploy::DeployModel, amb::Int, t::Int)
    @assert redeploy.status[amb] == :returning
    @assert redeploy.hospital[amb] == 0
    redeploy.status[amb] = :available
    redeploy.fromtime[amb] = t
    push!(redeploy.ambulances[redeploy.assignment[amb]], amb)
end