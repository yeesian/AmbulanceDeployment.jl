function request_for!(problem::DispatchProblem, location::Int, t::Int)
    @assert problem.deployment[location] > 0
    if problem.available[location] == 0 # no ambulances available
        problem.wait_queue[location] += 1
        @assert length(problem.amb_queue[location]) >= 1
        delay = problem.amb_queue[location][1] - t
    else # ambulance available
        problem.available[location] -= 1
        delay = 0
    end
    delay
end

function returned_to!(problem::DispatchProblem, location::Int, t::Int)
    @assert problem.deployment[location] > 0
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

function simulate_events!(problem::DispatchProblem,
                          model::DispatchModel,
                          turnaround::Distributions.LogNormal,
                          model_name::String="model",
                          verbose::Bool=false; mini_verbose=false)
    
    amb_allocation = problem.deployment # for checking of invariance properties
    events = form_queue(problem.emergency_calls)
    ncalls = nrow(problem.emergency_calls)
    dispatch_col = Symbol("$(model_name)_dispatch")
    delay_col = Symbol("$(model_name)_delay")
    hospital_col = Symbol("$(model_name)_hospital")
    problem.emergency_calls[dispatch_col] = 0
    problem.emergency_calls[delay_col] = 1000.0 # Inf; should be filled with smaller values after
    problem.emergency_calls[hospital_col] = 0

    while !isempty(events)
        (event, id, t, value) = dequeue!(events)
        @assert t > 0 "$((event, id, t, value))"
        if event == :call
            if sum(problem.deployment[problem.coverage[value,:]]) == 0
                verbose && println("no ambulance reachable for call at $value")
            elseif sum(problem.available[problem.coverage[value,:]]) > 0
                i = ambulance_for(model, id, problem, verbose=mini_verbose) # assume valid i
                problem.emergency_calls[id, dispatch_col] = i
                update_ambulances!(model, i, -1)
                problem.available[i] -= 1
                travel_time = ceil(Int,60*problem.emergency_calls[id, Symbol("stn$(i)_min")])
                @assert travel_time >= 0 "$id, $i"
                problem.emergency_calls[id, delay_col] = travel_time / 60 # minutes
                enqueue!(events, (:arrive, id, t + travel_time, i), t + travel_time)
            else
                push!(problem.wait_queue[value], id) # queue the emergency call
            end
        elseif event == :arrive
            scene_time = ceil(Int,15*rand(turnaround)) # time the ambulance ends service (back at base)
            @assert scene_time > 0
            enqueue!(events, (:convey, id, t + scene_time, value), t + scene_time)
        elseif event == :convey
            h = let mintime = Inf, minindex = 0
                for i in 1:nrow(problem.hospitals)
                    if problem.emergency_calls[id, Symbol("hosp$(i)_min")] < mintime
                        minindex = i
                        mintime = problem.emergency_calls[id, Symbol("hosp$(i)_min")]
                    end
                end
                minindex
            end
            problem.emergency_calls[id, hospital_col] = h
            conveytime = 15 + ceil(Int, 60*problem.emergency_calls[id, Symbol("hosp$(h)_min")])
            returntime = ceil(Int,60*problem.hospitals[h, Symbol("stn$(value)_min")])
            t_end = t + conveytime + returntime
            enqueue!(events, (:done, id, t_end, value), t_end)
        else
            @assert event == :done
            if sum(length(wq) for wq in problem.wait_queue[problem.coverage[:,value]]) > 0
                minindex = 0; mintime = Inf
                for nbhd in 1:size(problem.coverage,1)
                    if problem.coverage[nbhd,value] && length(problem.wait_queue[nbhd]) > 0
                        arrivaltime = problem.emergency_calls[problem.wait_queue[nbhd][1], :arrival_seconds]
                        if arrivaltime < mintime
                            mintime = arrivaltime
                            minindex = nbhd
                        end
                    end
                end
                @assert minindex != 0 && t - mintime >= 0 && mintime < Inf "$minindex, $t, $mintime"

                let id = shift!(problem.wait_queue[minindex])
                    problem.emergency_calls[id, dispatch_col] = value
                    travel_time = ceil(Int,60*problem.emergency_calls[id, Symbol("stn$(i)_min")])
                    @assert t - mintime > 0 && travel_time > 0 && t > 0
                    total_delay = (t - mintime) + travel_time
                    problem.emergency_calls[id, delay_col] = total_delay / 60 # minutes
                    enqueue!(events, (:arrive, id, t + total_delay, value), t + total_delay)
                end
            else
                returned_to!(problem, value, t)
                update_ambulances!(model, value, 1)
            end
        end
    end

    @assert all(problem.available .== amb_allocation)
    @assert all(problem.emergency_calls[dispatch_col] .> 0)
end
