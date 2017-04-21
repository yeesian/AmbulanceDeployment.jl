
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

function simulate_events!(problem::DispatchProblem,
                          model::DispatchModel,
                          redeploy::DeployModel,
                          turnaround::Distributions.LogNormal,
                          model_name::String="model",
                          verbose::Bool=false; mini_verbose=false)
    
    calls = problem.emergency_calls; ncalls = nrow(calls)
    coverage = problem.coverage; available = problem.available
    events = form_queue(calls); wait_queue = problem.wait_queue

    dispatch_col = Symbol("$(model_name)_dispatch");  calls[dispatch_col] = 0
    delay_col = Symbol("$(model_name)_delay");        calls[delay_col] = Inf
    hospital_col = Symbol("$(model_name)_hospital");  calls[hospital_col] = 0

    while !isempty(events)
        (event, id, t, value) = dequeue!(events)
        if event == :call
            if sum(problem.deployment[coverage[value,:]]) == 0
                verbose && println("no ambulance reachable for call at $value")
            elseif sum(available[coverage[value,:]]) > 0
                i = ambulance_for(model, id, problem, verbose=mini_verbose) # assume valid i
                a = respond_to!(redeploy, i, t)

                calls[id, dispatch_col] = i
                available[i] -= 1; update_ambulances!(model, i, -1)

                travel_time = ceil(Int,60*calls[id, Symbol("stn$(i)_min")])
                @assert travel_time >= 0
                calls[id, delay_col] = travel_time / 60 # minutes
                
                enqueue!(events, (:arrive, id, t + travel_time, a), t + travel_time)
            else
                push!(wait_queue[value], id) # queue the emergency call
            end
        elseif event == :arrive
            arriveatscene!(redeploy, value, t)
            scene_time = ceil(Int,15*rand(turnaround)) # time the ambulance spends at the scene
            @assert scene_time > 0
            enqueue!(events, (:convey, id, t + scene_time, value), t + scene_time)
        elseif event == :convey
            # determine the hospital to convey the patient to
            # (currently it's based on the closest hospital)
            h = let mintime = Inf, minindex = 0
                for i in 1:nrow(problem.hospitals)
                    traveltime = calls[id, Symbol("hosp$(i)_min")]
                    if !isna(traveltime) && traveltime < mintime
                        @assert traveltime >= 0
                        minindex = i
                        mintime = calls[id, Symbol("hosp$(i)_min")]
                    end
                end
                minindex
            end
            calls[id, hospital_col] = h
            conveytime = 15 + ceil(Int, 60*calls[id, Symbol("hosp$(h)_min")])

            conveying!(redeploy, value, t)
            enqueue!(events, (:return, id, t+conveytime, value), t+conveytime)
        elseif event == :return
            ### add redeployment model here
            reassign_ambulances!(redeploy)
            ### end redeployment model here
            stn = returning_to!(redeploy, value, t)
            returntime = ceil(Int,60*problem.hospitals[h, Symbol("stn$(stn)_min")])
            t_end = t + conveytime + returntime
            enqueue!(events, (:done, id, t_end, value), t_end)
        else
            @assert event == :done
            stn = redeploy.assignment[value]
            if sum(length(wq) for wq in wait_queue[coverage[:,stn]]) > 0
                # people are waiting in a queue
                @assert redeploy.status[value] == :returning
                redeploy.fromtime[value] = t
                redeploy.status[value] = :responding

                # determine the person who has waited the longest
                minindex = 0; mintime = Inf
                for nbhd in 1:size(coverage,1)
                    if coverage[nbhd,stn] && length(wait_queue[nbhd]) > 0
                        arrivaltime = calls[wait_queue[nbhd][1], :arrival_seconds]
                        if arrivaltime < mintime
                            mintime = arrivaltime
                            minindex = nbhd
                        end
                    end
                end
                @assert minindex != 0 && t - mintime >= 0 && mintime < Inf

                # respond to the person
                let id = shift!(wait_queue[minindex])
                    calls[id, dispatch_col] = stn
                    travel_time = ceil(Int,60*calls[id, Symbol("stn$(stn)_min")])
                    @assert t - mintime >= 0 && travel_time > 0 && t > 0
                    total_delay = (t - mintime) + travel_time
                    calls[id, delay_col] = total_delay / 60 # minutes
                    enqueue!(events, (:arrive, id, t + total_delay, value), t + total_delay)
                end
            else # returned to base location
                returned_to!(redeploy, value, t)
                returned_to!(problem, stn, t)
                update_ambulances!(model, stn, 1)
            end
        end
    end

    @assert all(available .== problem.deployment)
    @assert all(calls[dispatch_col] .> 0)
end
