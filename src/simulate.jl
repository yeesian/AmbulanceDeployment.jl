
function returned_to!(problem::DispatchProblem, location::Int, t::Int)
    @assert problem.deployment[location] > 0
    @assert problem.available[location] >= 0
    problem.available[location] += 1
end

function form_queue(emergency_calls::DataFrame)
    pq = PriorityQueue{Tuple{Symbol,Int,Int,Int},Int,Base.Order.ForwardOrdering}()
    for i in 1:nrow(emergency_calls)
        enqueue!(pq, (:call,
                      i,
                      emergency_calls[i, :arrival_seconds],
                      emergency_calls[i, :neighborhood]),
                 emergency_calls[i, :arrival_seconds])
    end
    pq
end
form_queue(problem::DispatchProblem) = form_queue(problem.emergency_calls)

function call_event!(
        events::PriorityQueue,
        problem::DispatchProblem,
        dispatch::DispatchModel,
        redeploy::DeployModel,
        id::Int, # the id of the emergency call
        t::Int, # the time of the emergency call
        nbhd::Int # the neighborhood the call is from
    )
    if sum(problem.deployment[problem.coverage[nbhd,:]]) == 0
        verbose && println("no ambulance reachable for call at $nbhd")
    elseif sum(problem.available[problem.coverage[nbhd,:]]) > 0
        i = ambulance_for(dispatch, id, problem) # assume valid i
        update_ambulances!(dispatch, i, -1)
        problem.emergency_calls[id, :dispatch_from] = i
        problem.available[i] -= 1

        travel_time = ceil(Int, 60*problem.emergency_calls[id, Symbol("stn$(i)_min")])
        @assert travel_time >= 0
        calls[id, :delay] = travel_time / 60 # minutes
        
        amb = respond_to!(redeploy, i, t)
        enqueue!(events, (:arrive, id, t + travel_time, amb), t + travel_time)
    else
        push!(problem.wait_queue[nbhd], id) # queue the emergency call
    end
end

function arrive_event!(
        events::PriorityQueue,
        redeploy::DeployModel,
        id::Int, # the id of the emergency call
        t::Int, # the time of the emergency call
        amb::Int,
        turnaround
    )
    arriveatscene!(redeploy, amb, t)
    scene_time = ceil(Int,15*rand(turnaround)) # time the ambulance spends at the scene
    @assert scene_time > 0
    enqueue!(events, (:convey, id, t + scene_time, amb), t + scene_time)
end

"determine the hospital to convey the patient to (currently it's based on the closest hospital)"
function convey_event!(
        events::PriorityQueue,
        problem::DispatchProblem,
        redeploy::DeployModel,
        id::Int, # the id of the emergency call
        t::Int, # the time of the emergency call
        amb::Int
    )
    h = let mintime = Inf, minindex = 0
        for h in 1:nrow(problem.hospitals)
            traveltime = problem.emergency_calls[id, Symbol("hosp$(h)_min")]
            if !isna(traveltime) && traveltime < mintime
                @assert traveltime >= 0
                minindex = h; mintime = traveltime
            end
        end
        minindex
    end
    redeploy.hospital[amb] = problem.emergency_calls[id, :hospital] = h
    conveytime = 15 + ceil(Int, 60*problem.emergency_calls[id, Symbol("hosp$(h)_min")])
    conveying!(redeploy, amb, h, t)
    enqueue!(events, (:return, id, t+conveytime, amb), t+conveytime)
end

function return_event!(
        events::PriorityQueue,
        problem::DispatchProblem,
        redeploy::DeployModel,
        id::Int,
        t::Int,
        amb::Int
    )
    stn = returning_to!(redeploy, amb, t)
    h = redeploy.hospital[amb]
    returntime = ceil(Int,60*problem.hospitals[h, Symbol("stn$(stn)_min")])
    t_end = t + returntime
    enqueue!(events, (:done, id, t_end, amb), t_end)
end

function done_event!(
        events::PriorityQueue,
        problem::DispatchProblem,
        model::DispatchModel,
        redeploy::DeployModel,
        id::Int,
        t::Int,
        amb::Int
    )
    stn = redeploy.assignment[amb]
    if sum(length(wq) for wq in problem.wait_queue[problem.coverage[:,stn]]) > 0
        # people are waiting in a queue
        @assert redeploy.status[amb] == :returning
        redeploy.fromtime[amb] = t
        redeploy.status[amb] = :responding

        # determine the person who has waited the longest
        minindex = 0; mintime = Inf
        for nbhd in 1:size(problem.coverage,1)
            if problem.coverage[nbhd,stn] && length(problem.wait_queue[nbhd]) > 0
                arrivaltime = problem.emergency_calls[problem.wait_queue[nbhd][1], :arrival_seconds]
                if arrivaltime < mintime
                    mintime = arrivaltime
                    minindex = nbhd
                end
            end
        end
        @assert minindex != 0 && t - mintime >= 0 && mintime < Inf

        # respond to the person
        let id = shift!(problem.wait_queue[minindex])
            problem.emergency_calls[id, :dispatch_from] = stn
            travel_time = ceil(Int,60*problem.emergency_calls[id, Symbol("stn$(stn)_min")])
            @assert t - mintime >= 0 && travel_time > 0 && t > 0
            total_delay = (t - mintime) + travel_time
            problem.emergency_calls[id, delay_col] = total_delay / 60 # minutes
            enqueue!(events, (:arrive, id, t + total_delay, amb), t + total_delay)
        end
    else # returned to base location
        returned_to!(redeploy, amb, t)
        returned_to!(problem, stn, t)
        update_ambulances!(model, stn, 1)
    end
end

function simulate_events!(problem::DispatchProblem,
                          model::DispatchModel,
                          redeploy::DeployModel,
                          turnaround::Distributions.LogNormal,
                          model_name::String="model",
                          verbose::Bool=false; mini_verbose=false)
    events = form_queue(problem)
    while !isempty(events)
        (event, id, t, value) = dequeue!(events)
        if event == :call
            call_event!(events, problem, model, redeploy, id, t, value)
        elseif event == :arrive
            arrive_event!(events, redeploy, id, t, value, turnaround)
        elseif event == :convey
            convey_event!(events, problem, redeploy, id, t, value)
        elseif event == :return
            reassign_ambulances!(redeploy)
            return_event!(events, problem, redeploy, id, t, value)
        else
            @assert event == :done
            done_event!(events, problem, model, redeploy, id, t, value)
        end
    end
    @assert all(available .== problem.deployment)
    @assert all(calls[dispatch_col] .> 0)
end
