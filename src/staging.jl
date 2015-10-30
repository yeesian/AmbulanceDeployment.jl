

function request_for!(problem::DispatchProblem, location::Int,
                     t::Int, verbose::Bool=true)
    if problem.available[i] == 0
        problem.wait_queue[i] += 1
        delay = problem.amb_queue[i][1] - t
    else
        problem.available[i] -= 1
        delay = 0
    end
    delay
end

function returned_to!(problem::DispatchProblem, location::Int, t::Int)
    shift!(problem.amb_queue[value])
    problem.available[location] += (problem.wait_queue[location] == 0)
    problem.wait_queue[location] -= (problem.wait_queue[location] > 0)
end

function ambulance_for(model::DispatchModel,
                       region::Int,
                       problem::DispatchProblem,
                       coverage::BitArray{2},
                       verbose::Bool=false)
    i = available_for(region, problem.available, model)
    if i == 0
        return waiting_for(region, problem.amb_queue, problem.coverage)
    else
        update_ambulances!(model, i, -1)
        return i
    end
end

function waiting_for(j::Int, next_in_location::Vector{Vector{Int}}, coverage::BitArray{2})
    nlocations = length(next_in_location)
    earliest_time = Inf
    earliest_location = 0
    for i in 1:nlocations
        if !isempty(next_in_location[i]) && coverage[j,i]
            if next_in_location[i][1] < earliest_time
                earliest_location = i
                earliest_time = next_in_location[i][1]
            end
        end
    end
    @assert earliest_location != 0
    earliest_location
end

function form_queue(emergency_calls::DataFrame)
    pq = PriorityQueue{(Int,Symbol,Int,Int),Int}()
    for i in 1:nrow(emergency_calls)
        enqueue!(pq, (i, emergency_calls[i, :arrival_seconds],
                         emergency_calls[i, :neighborhood]),
                 emergency_calls[i, :arrival_seconds])
    end
    pq
end

function simulate_events!(problem::DispatchProblem,
                          model::DispatchModel,
                          model_name::ASCIIString="model",
                          verbose::Bool=false)
    
    amb_allocation = copy(problem.available) # for checking of invariance properties
    events = form_queue(problem.emergency_calls)
    ncalls = nrow(problem.emergency_calls)
    dispatch_col = symbol("$(model_name)_dispatch")
    delay_col = symbol("$(model_name)_delay")
    problem.emergency_calls[dispatch_col] = 0
    problem.emergency_calls[delay_name] = 0

    while !isempty(events)
        (id, t, value) = dequeue!(events)
        if id <= problem.nlocations # event == :call
            i = ambulance_for(model, value, problem, verbose)
            problem.emergency_calls[id, dispatch_col] = i
            verbose && println("time $t: incoming emergency call from region $value")
            verbose && println("time $t: requesting for ambulance from location $i")

            delay = request_for!(problem, i, t, verbose)
            problem.emergency_calls[id, delay_col] = delay
            verbose && (delay  > 0) && println("time $t:  send from $(problem.amb_queue[i])")
            verbose && (delay == 0) && println("time $t:  ambulance dispatched")

            t_end = t + delay + 3600
            push!(problem.amb_queue[i], t_end)
            enqueue!(events, (id+ncalls, t_end, i), t_end)
        else # event == :done
            verbose && println("time $t: 1 ambulance returned to location $value")
            returned_to!(problem, value, t, verbose)
            update_ambulances(model, value, 1)
        end
        # Sanity check
        @assert amb_allocation[i] == problem.available + length(problem.amb_queue) - problem.wait_queue[i]
    end

    for i in eachindex(problem.amb_queue)
        @assert isempty(problem.amb_queue[i])
    end
    @assert sum(problem.wait_queue) == 0
    @assert all(problem.available .== amb_allocation)
end


