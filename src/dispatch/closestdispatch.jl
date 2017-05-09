type ClosestDispatch <: DispatchModel
    drivetime::DataFrame
    candidates::Vector{Vector{Int}}
end

function ambulance_for(model::ClosestDispatch,
                       id::Int,
                       problem::DispatchProblem; verbose=false)
    i = available_for(model, id, problem)
    if i == 0
        verbose && println(problem.available)
        verbose && println(model.available)
        return waiting_for(problem.emergency_calls[id, :neighborhood],
                           problem)
    else
        update_ambulances!(model, i, -1)
        return i
    end
end

function ClosestDispatch(p::DeploymentProblem, drivetime::DataFrame)
    candidates = Array(Vector{Int}, p.nregions)
    I = 1:p.nlocations
    for region in 1:p.nregions
        candidates[region] = I[vec(p.coverage[region,:])]
    end
    ClosestDispatch(drivetime, candidates)
end

update_ambulances!(model::ClosestDispatch, i::Int, delta::Int) = nothing

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
