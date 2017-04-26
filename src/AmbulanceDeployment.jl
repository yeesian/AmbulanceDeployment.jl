module AmbulanceDeployment

    import JuMP
    import Base.Collections: PriorityQueue, enqueue!, dequeue!
    import Colors: LCHab, Colorant
    import Compose: Polygon, UnitBox, context, compose, linewidth
    import Compose: stroke, fill, mm, circle
    import DataFrames: DataFrame, isna, nrow
    import Distributions: Poisson, LogNormal, quantile, sample
    import Gadfly: lab_gradient
    import GeoInterface: coordinates
    import GeoConverters: composeform
    import Gurobi: GurobiSolver
    import Winston: FramedPlot, Curve, Legend, setattr, add
    import Distributions
    import StatsBase

    export DeploymentProblem,
           RobustDeployment,
           StochasticDeployment,
           MALPDeployment,
           MEXCLPDeployment,

           DispatchProblem,
           ClosestDispatch,
           LPDispatchGreedy,
           LPDispatchBacklog,
           LPDispatchRandom,

           NoRedeployModel,

           solve,
           evaluate,
           deployment,
           convergence_plot,
           compose_neighborhoods,
           compose_locations,
           compose_chloropleth,
           performance,
           test_performance,
           plot_timings,
           simulate_events!,
           initialize!

    type DeploymentProblem{ IM <: AbstractMatrix{Int},
                            BM <: AbstractMatrix{Bool}}
        nambulances::Int
        nlocations::Int
        nregions::Int
        demand::IM      # nperiod x nregion
        train::Vector{Int}
        test::Vector{Int}
        coverage::BM    # nregion x nlocation
        adjacency::BM   # nregion x nregion
    end

    abstract DeploymentModel

    function naive_solution(p::DeploymentProblem)
        # evenly distribute the ambulances over all the locations
        x = zeros(Int, p.nlocations)
        for i in 0:p.nambulances-1
            x[1+ (i % p.nlocations)] += 1
        end
        x
    end

    type DispatchProblem
        emergency_calls::DataFrame
        hospitals::DataFrame
        stations::DataFrame
        coverage::Matrix{Bool} # (nbhd x stns)
        turnaround::LogNormal
        wait_queue::Vector{Vector{Int}} # length nbhd
        available::Vector{Int}
        deployment::Vector{Int}
        
        DispatchProblem{BM}(emergency_data::DataFrame,
                            hospitals::DataFrame,
                            stations::DataFrame,
                            coverage::BM,
                            turnaround::LogNormal = LogNormal(3.65, 0.3)) =
            new(emergency_data, hospitals, stations, coverage, turnaround)
    end

    function initialize!(problem::DispatchProblem,
                         deployment::Vector{Int})
        problem.wait_queue = [Int[] for i in 1:size(problem.coverage,1)]
        problem.available = copy(deployment)
        problem.deployment = deepcopy(deployment)
        problem
    end

    abstract DispatchModel

    function ambulance_for(model::DispatchModel,
                           id::Int,
                           problem::DispatchProblem;
                           verbose::Bool=false)
        i = available_for(model, id, problem, verbose=verbose)
        verbose && println("dispatching $i")
        if i == 0
            j = problem.emergency_calls[id, :neighborhood]
            return waiting_for(j, problem)
        else
            update_ambulances!(model, i, -1)
            return i
        end
    end

    function waiting_for(j::Int, problem::DispatchProblem)
        earliest_time = Inf
        earliest_location = 0
        for i in 1:length(problem.amb_queue)
            if problem.coverage[j,i]
                @assert problem.available[i] == 0
                if problem.deployment[i] > 0
                    @assert !isempty(problem.amb_queue[i])
                    if problem.amb_queue[i][1] < earliest_time
                        earliest_location = i
                        earliest_time = problem.amb_queue[i][1]
                    end
                end
            end
        end
        earliest_location
    end

    function returned_to!(problem::DispatchProblem, location::Int, t::Int)
        @assert problem.deployment[location] > 0
        @assert problem.available[location] >= 0
        problem.available[location] += 1
    end

    type Params
        α::Float64 # Probabilistic Guarantee
        ε::Float64 # Convergence
        δ::Float64 # Solver Tolerance

        nperiods::Int # for StochasticDeployment

        maxiter::Int # for RobustDeployment
    end
    Params() = Params(0.01, 0.5, 1e-6, 500, 50)
    params = Params()

    include("robust.jl")
    include("stochastic.jl")
    include("mexclp.jl")
    include("malp.jl")
    include("evaluate.jl")
    include("plot.jl")
    #include("dispatch.jl")

    type ClosestDispatch <: DispatchModel
        drivetime::DataFrame
        candidates::Vector{Vector{Int}}
    end

    function ambulance_for(model::ClosestDispatch,
                           id::Int,
                           problem::DispatchProblem; verbose=false)
        verbose && print("*")
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


    include("redeploy.jl")
    include("simulate.jl")
end