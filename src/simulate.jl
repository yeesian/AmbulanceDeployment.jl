function request_for!(problem::DispatchProblem, location::Int, t::Int)
    @assert problem.deployment[location] > 0
    if problem.available[location] == 0
        problem.wait_queue[location] += 1
        @assert length(problem.amb_queue[location]) >= 1
        delay = problem.amb_queue[location][1] - t
    else
        problem.available[location] -= 1
        delay = 0
    end
    delay
end

function returned_to!(problem::DispatchProblem, location::Int, t::Int)
    @assert problem.deployment[location] > 0
    @assert !isempty(problem.amb_queue[location])
    shift!(problem.amb_queue[location])
    problem.available[location] += (problem.wait_queue[location] == 0)
    problem.wait_queue[location] -= (problem.wait_queue[location] > 0)
end

function form_queue(emergency_calls::DataFrame)
    pq = pq = PriorityQueue{Tuple{Int,Int,Int},Int,Base.Order.ForwardOrdering}()
    for i in 1:nrow(emergency_calls)
        enqueue!(pq, (i, emergency_calls[i, :arrival_seconds],
                         emergency_calls[i, :neighborhood]),
                 emergency_calls[i, :arrival_seconds])
    end
    pq
end

function simulate_events!(problem::DispatchProblem,
                          model::DispatchModel,
                          turnaround::Distributions.LogNormal,
                          model_name::ASCIIString="model",
                          verbose::Bool=false; mini_verbose=false)
    
    amb_allocation = problem.deployment # for checking of invariance properties
    events = form_queue(problem.emergency_calls)
    ncalls = nrow(problem.emergency_calls)
    dispatch_col = symbol("$(model_name)_dispatch")
    delay_col = symbol("$(model_name)_delay")
    problem.emergency_calls[dispatch_col] = 0
    problem.emergency_calls[delay_col] = 0

    while !isempty(events)
        (id, t, value) = dequeue!(events)
        if id <= ncalls # event == :call
            i = ambulance_for(model, id, problem, verbose=mini_verbose)
            if i == 0
                region = problem.emergency_calls[id, :neighborhood]
                println(amb_allocation)
                println(vec(problem.coverage[region,:]))
                println(amb_allocation[vec(problem.coverage[region,:])])
                println(problem.amb_queue)
                @assert sum(amb_allocation[vec(problem.coverage[region,:])]) == 0
                println("       out-of-range: $((id,t,value))")
            else
                problem.emergency_calls[id, dispatch_col] = i
                verbose && println("time $t: incoming emergency call from region $value")
                verbose && println("time $t:  requesting for ambulance from location $i")

                update_ambulances!(model, i, -1)
                delay = request_for!(problem, i, t) # queue wait time if busy
                if delay < 0
                    println(problem.amb_queue[i])
                    println(t)
                    @assert delay >= 0
                end
                delay += ceil(Int,60*problem.emergency_calls[id, i+2]) # include road travel time
                problem.emergency_calls[id, delay_col] = delay
                verbose && (delay  > 0) && println("time $t:   send from $(problem.amb_queue[i])")
                verbose && (delay == 0) && println("time $t:   ambulance dispatched from $i")

                t_end = t + delay + ceil(Int,60*rand(turnaround)) # time the ambualnce ends service (back at base)
                push!(problem.amb_queue[i], t_end)
                sort!(problem.amb_queue[i]) # may not be in order
                enqueue!(events, (id+ncalls, t_end, i), t_end)
            end
        else # event == :done
            verbose && println("time $t: 1 ambulance returned to location $value")
            returned_to!(problem, value, t)
            update_ambulances!(model, value, 1)
        end
        # Sanity check
        for i in eachindex(problem.amb_queue)
            @assert amb_allocation[i] == problem.available[i] + length(problem.amb_queue[i]) - problem.wait_queue[i]
        end
    end

    for i in eachindex(problem.amb_queue)
        @assert isempty(problem.amb_queue[i])
    end
    @assert sum(problem.wait_queue) == 0
    @assert all(problem.available .== amb_allocation)
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
