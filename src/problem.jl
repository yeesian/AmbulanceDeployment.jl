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
    deployment::Vector{Int}

    wait_queue::Vector{Vector{Int}} # length nbhd
    available::Vector{Int} # length stns
    
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
    # @assert problem.deployment[location] > 0
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
