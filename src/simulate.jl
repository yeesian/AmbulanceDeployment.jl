
type EMSEngine{T}
    eventlog::DataFrame
    eventqueue::PriorityQueue{T,Int,Base.Order.ForwardOrdering}
end

function EMSEngine(problem::DispatchProblem)
    ncalls = nrow(problem.emergency_calls)
    eventlog = DataFrame(
        id = 1:ncalls,
        dispatch_from = zeros(Int, ncalls),
        delay = fill(Inf, ncalls),
        hospital = zeros(Int, ncalls)
    )
    eventqueue = PriorityQueue{Tuple{Symbol,Int,Int,Int},Int,Base.Order.ForwardOrdering}()
    for i in 1:nrow(problem.emergency_calls)
        t = problem.emergency_calls[i, :arrival_seconds]
        enqueue!(eventqueue, (:call, i, t, problem.emergency_calls[i, :neighborhood]), t)
    end
    EMSEngine{Tuple{Symbol,Int,Int,Int}}(eventlog, eventqueue)
end

function call_event!(
        ems::EMSEngine,
        problem::DispatchProblem,
        dispatch::DispatchModel,
        redeploy::DeployModel,
        id::Int, # the id of the emergency call
        t::Int, # the time of the emergency call
        nbhd::Int; # the neighborhood the call is from
        verbose::Bool = false
    )
    if sum(problem.deployment[problem.coverage[nbhd,:]]) == 0
        verbose && println("no ambulance reachable for call at $nbhd")
    elseif sum(problem.available[problem.coverage[nbhd,:]]) > 0
        i = ambulance_for(dispatch, id, problem)
        @assert i > 0 "dispatch from $i" # assume valid i (enforced by <if> condition) 
        update_ambulances!(dispatch, i, -1)
        ems.eventlog[id, :dispatch_from] = i
        @assert problem.available[i] > 0
        problem.available[i] -= 1

        travel_time = ceil(Int, 60*problem.emergency_calls[id, Symbol("stn$(i)_min")])
        @assert travel_time >= 0
        ems.eventlog[id, :delay] = travel_time / 60 # minutes
        
        amb = respond_to!(redeploy, i, t)
        enqueue!(ems.eventqueue, (:arrive, id, t + travel_time, amb), t + travel_time)
    else
        push!(problem.wait_queue[nbhd], id) # queue the emergency call
    end
end

function arrive_event!(
        ems::EMSEngine,
        problem::DispatchProblem,
        redeploy::DeployModel,
        id::Int, # the id of the emergency call
        t::Int, # the time of the emergency call
        amb::Int
    )
    arriveatscene!(redeploy, amb, t)
    # time the ambulance spends at the scene
    scene_time = ceil(Int,60*0.4*rand(problem.turnaround)) # 60sec*0.4*mean(40minutes) ~ 15minutes
    @assert scene_time > 0
    enqueue!(ems.eventqueue, (:convey, id, t + scene_time, amb), t + scene_time)
end

"determine the hospital to convey the patient to (currently it's based on the closest hospital)"
function convey_event!(
        ems::EMSEngine,
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
    @assert h != 0
    redeploy.hospital[amb] = ems.eventlog[id, :hospital] = h
    conveytime = 60*15 + ceil(Int, 60*problem.emergency_calls[id, Symbol("hosp$(h)_min")]) # ~20minutes
    @assert conveytime >= 0 conveytime
    conveying!(redeploy, amb, h, t)
    enqueue!(ems.eventqueue, (:return, id, t+conveytime, amb), t+conveytime)
end

function return_event!(
        ems::EMSEngine,
        problem::DispatchProblem,
        redeploy::DeployModel,
        id::Int,
        t::Int,
        amb::Int
    )
    stn = returning_to!(redeploy, amb, t)
    h = redeploy.hospital[amb]
    returntime = ceil(Int,60*problem.hospitals[h, Symbol("stn$(stn)_min")]) # ~ 10minutes
    @assert returntime >= 0 returntime
    t_end = t + returntime
    enqueue!(ems.eventqueue, (:done, id, t_end, amb), t_end)
end

function done_event!(
        ems::EMSEngine,
        problem::DispatchProblem,
        dispatch::DispatchModel,
        redeploy::DeployModel,
        id::Int,
        t::Int,
        amb::Int
    )
    stn = redeploy.assignment[amb]; @assert stn > 0
    if sum(length(wq) for wq in problem.wait_queue[problem.coverage[:,stn]]) > 0
        # people are waiting in a queue
        redirected!(redeploy, amb, t)
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
        @assert minindex != 0
        @assert t - mintime >= 0
        @assert 0 <= mintime < Inf
        # respond to the person
        let id = shift!(problem.wait_queue[minindex])
            ems.eventlog[id, :dispatch_from] = stn
            travel_time = ceil(Int,60*problem.emergency_calls[id, Symbol("stn$(stn)_min")])
            @assert travel_time >= 0
            total_delay = (t - mintime) + travel_time; @assert total_delay >= 0
            tarrive = t + total_delay; @assert t + total_delay >= 0 "$t, $total_delay"
            ems.eventlog[id, :delay] = total_delay / 60 # minutes
            enqueue!(ems.eventqueue, (:arrive, id, tarrive, amb), tarrive)
        end
    else # returned to base location
        returned_to!(redeploy, amb, t)
        @assert problem.available[stn] >= 0
        returned_to!(problem, stn, t)
        update_ambulances!(dispatch, stn, 1)
    end
end

function simulate_events!(
        problem::DispatchProblem,
        dispatch::DispatchModel,
        redeploy::DeployModel;
        verbose::Bool=false
    )
    ems = EMSEngine(problem)
    # @show problem.available
    # @show redeploy.ambulances
    k = 0
    while !isempty(ems.eventqueue)
        if k > 10_000
            break
        else
            k += 1
        end
        (event, id, t, value) = dequeue!(ems.eventqueue)
        # @show (event, id, t, value)
        @assert t >= 0 # in case of integer overflow (when calls > ambulances)
        if event == :call
            call_event!(ems, problem, dispatch, redeploy, id, t, value, verbose=verbose)
        elseif event == :arrive
            arrive_event!(ems, problem, redeploy, id, t, value)
        elseif event == :convey
            convey_event!(ems, problem, redeploy, id, t, value)
        elseif event == :return
            reassign_ambulances!(problem, redeploy, t)
            return_event!(ems, problem, redeploy, id, t, value)
        else
            @assert event == :done
            done_event!(ems, problem, dispatch, redeploy, id, t, value)
        end
        # @show problem.available
        # @show redeploy.ambulances
        for i in eachindex(problem.available)
            @assert problem.available[i] == length(redeploy.ambulances[i])
        end
    end
    # @assert all(problem.available .== problem.deployment)
    @assert all(ems.eventlog[:dispatch_from] .>= 0)
    ems.eventlog
end
