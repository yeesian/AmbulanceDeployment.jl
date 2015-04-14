module AmbulanceDeployment

    using JuMP, Distributions, Gurobi

    export RobustDeployment,
           StochasticDeployment,
           DeploymentProblem

    type DeploymentProblem{ IM <: AbstractMatrix{Int},
                            BM <: AbstractMatrix{Bool}}
        nambulances::Int
        nlocations::Int
        nregions::Int
        demand::IM      # nperiod x nregion
        coverage::BM    # nregion x nlocation
        adjacency::BM   # nregion x nregion
    end

    abstract DeploymentModel
    model{T <: DeploymentModel}(m::T) = m.model
    deployment{T <: DeploymentModel}(m::T) = m.x

    include("robust.jl")
    include("stochastic.jl")
    include("evaluate.jl")
end