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
    problem.emergency_calls[dispatch_col] = 0
    problem.emergency_calls[delay_col] = 1000.0 # Inf; should be filled with smaller values after

    while !isempty(events)
        (event, id, t, value) = dequeue!(events)
        if event == :call
            if sum(problem.deployment[problem.coverage[value,:]]) == 0
                verbose && println("no ambulance reachable for call at $value")
            elseif sum(problem.available[problem.coverage[value,:]]) > 0
                i = ambulance_for(model, id, problem, verbose=mini_verbose) # assume valid i
                problem.emergency_calls[id, dispatch_col] = i
                update_ambulances!(model, i, -1)
                travel_time = ceil(Int,60*problem.emergency_calls[id, Symbol("stn$(i)_min")])
                problem.emergency_calls[id, delay_col] = travel_time / 60 # minutes
                enqueue!(events, (:arrive, id, t + travel_time, i), t + travel_time)
            else
                sort!(push!(problem.wait_queue[value], t)) # queue the emergency call
            end
        elseif event == :arrive
            t_end = t + ceil(Int,45*rand(turnaround)) # time the ambulance ends service (back at base)
            enqueue!(events, (:done, id, t_end, i), t_end)
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
                @assert minindex != 0 && t - mintime >= 0 && mintime < Inf

                let id = problem.wait_queue[nbhd][1]
                    problem.emergency_calls[id, dispatch_col] = value
                    travel_time = ceil(Int,60*problem.emergency_calls[id, Symbol("stn$(i)_min")])
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

    # @assert sum(problem.wait_queue) == 0
    # @assert all(problem.available .== amb_allocation)
    # @assert all(problem.emergency_calls[dispatch_col] .> 0)
end

# function form_queue(emergency_calls::DataFrame)
#     pq = PriorityQueue{(Int,Symbol,Int,Int),Int}()
#     for i in 1:nrow(emergency_calls)
#         enqueue!(pq, (i, :call, emergency_calls[i, :arrival_seconds], emergency_calls[i, :neighborhood]),
#                  emergency_calls[i, :arrival_seconds])
#     end
#     pq
# end

# function random_ambulance_for(j::Int, available::Vector{Int}, coverage::BitArray{2})
#     nlocations = size(coverage,2)
#     available_locations = [1:nlocations][vec(coverage[j,:])]
#     for loc in sample(available_locations,length(available_locations),replace=false)
#         if available[loc] > 0
#             return loc
#         end
#     end
#     return 0
# end

# function ambulance_for(model::DispatchModel,
#                        region::Int,
#                        problem::DispatchProblem,
#                        coverage::BitArray{2},
#                        verbose::Bool=false)
#     i = available_for(region, problem.available, model)
#     if i == 0
#         return waiting_for(region, problem.amb_queue, problem.coverage)
#     else
#         update_ambulances(model, i, -1)
#         return i
#     end
# end

# function waiting_for(j::Int, next_in_location::Vector{Vector{Int}}, coverage::BitArray{2})
#     nlocations = length(next_in_location)
#     earliest_time = Inf
#     earliest_location = 0
#     for i in 1:nlocations
#         if !isempty(next_in_location[i]) && coverage[j,i]
#             if next_in_location[i][1] < earliest_time
#                 earliest_location = i
#                 earliest_time = next_in_location[i][1]
#             end
#         end
#     end
#     @assert earliest_location != 0
#     (earliest_location, earliest_time)
# end

# function simulate_events(deployment::Vector{Int},
#                          coverage::BitArray{2},
#                          emergency_calls::DataFrame;
#                          verbose=false)
#     ncalls = nrow(emergency_calls)
#     (nregions, nlocations) = size(coverage)
#     locations = 1:nlocations
#     regions = 1:nregions
#     delays = (Int,Int,Int,Int)[]
#     available = deployment[:]
#     next_in_location = Array(Vector{Int}, nlocations)

#     for i in locations
#         next_in_location[i] = Vector{Int}()
#     end
#     events_queue = form_queue(emergency_calls)
#     while !isempty(events_queue)
#         (id, event, t, value) = dequeue!(events_queue)
#         if event == :call
#             j = value # incoming emergency call from region j
#             verbose && println("incoming emergency call from region $j at time $t")
#             delay = 0
#             i = ambulance_for(j, available, coverage)
#             if i == 0 # no ambulances available at the moment
#                 verbose && println("  no ambulances available for $j at time $(t)!")
#                 (i,tnew) = next_ambulance_for(j, next_in_location, coverage)
#                 verbose && println("  dispatching ambulance from $i, from queue $(next_in_location[i])")
#                 shift!(next_in_location[i])
#                 delay = tnew - t
#                 push!(delays, (id, delay, i, j))
#             else
#                 verbose && println("  dispatching ambulance from $i")
#                 available[i] -= 1
#             end
#             tend = t + delay + 3600
#             enqueue!(events_queue, (id+ncalls, :done, tend, i), tend)
#             push!(next_in_location[i], tend)
#             verbose && println("  with $(available[i]) remaining at $i, and a queue of $(next_in_location[i]).")
#         else
#             @assert event == :done
#             i = value
#             available[i] += 1
#             verbose && println("1 ambulance returned to location $i at time $t")
#             if t == next_in_location[i][1]
#                 shift!(next_in_location[i])
#             end
#         end
#     end
#     delays
# end
