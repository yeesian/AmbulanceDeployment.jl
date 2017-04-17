using AmbulanceDeployment, DataFrames, Winston, JLD

hourly_calls = DataFrames.readtable("data/processed/2-weekday_calls.csv")
# weekend_hourly_calls = DataFrames.readtable("data/processed/2-weekend_calls.csv")
adjacent_nbhd = DataFrames.readtable("data/processed/2-adjacent_nbhd.csv")
coverage = JLD.load("data/processed/3-coverage.jld", "stn_coverage")
incidents = DataFrames.readtable("data/processed/3-incidents_drivetime.csv");

regions = Int[parse(Int,string(x)[2:end]) for x in names(hourly_calls[5:end])]
locations = collect(1:size(coverage,2))
adjacent = convert(Array, adjacent_nbhd[2:end])[regions,regions] .> 0.5
demand = convert(Array,hourly_calls[:,5:end]);

incidents = incidents[~DataFrames.isna(incidents[:stn1_min]), :] # drop 44 calls that were "unreachable" (because all stations are reachable from each other)
incidents[:interarrival_seconds] = [0; incidents[:arrival_seconds][2:end] - incidents[:arrival_seconds][1:end-1]]

incidents[:isweekday] = true
incidents[:isweekday][incidents[:dow] .== "Sun"] = false
incidents[:isweekday][incidents[:dow] .== "Sat"] = false;
incidents = incidents[incidents[:isweekday],:];

regions2index = Dict{Int,Int}(regions[i]=>i for i in 1:length(regions))
incidents[:neighborhood] = [regions2index[x] for x in incidents[:neighborhood]];

calls = incidents[:,[[:hour,:dow,:month,:year,:neighborhood,:interarrival_seconds];
                     Symbol[Symbol("stn$(i)_min") for i in locations]]]

# We focus on emergency calls during the "peak period" (8AM - 8PM),
# with the emergency calls from the first 3 month as our training set,
# and the subsequent emergency calls from the remaining months as our test set
peak_period = (hourly_calls[:hour] .>= 8) .* (hourly_calls[:hour] .<= 20)
indices = 1:DataFrames.nrow(hourly_calls);

train_filter = (hourly_calls[:year] .== 2012) .* (hourly_calls[:month] .<= 3)
test_filter  = !train_filter;

train_indices = indices[train_filter]
test_indices = indices[test_filter];

# Same as for the hourly calls; but this is for individual emergency calls
inc_peak_period = (calls[:hour] .>= 8) .* (calls[:hour] .<= 20)
inc_indices = 1:DataFrames.nrow(calls);

inc_train_filter = (calls[:year] .== 2012) .* (calls[:month] .<= 3)
inc_test_filter  = !inc_train_filter

inc_train_indices = inc_indices[inc_train_filter]
inc_test_indices = inc_indices[inc_test_filter];

# we distinguish between peak and offpeak hours

train_peak = indices[peak_period .* train_filter]
train_offpeak = indices[!peak_period .* train_filter]

test_peak = indices[peak_period .* test_filter]
test_offpeak = indices[!peak_period .* test_filter]

test_inc_peak = inc_peak_period .* inc_test_filter
test_inc_offpeak = !inc_peak_period .* inc_test_filter;

p = DeploymentProblem(30, length(locations), length(regions), demand, train_indices,
                      test_indices, coverage[regions,:], Array{Bool,2}(adjacent));

solverstats = JLD.load("data/processed/4-solve-stats.jld")
amb_deployment = solverstats["amb_deployment"]

peakstats = JLD.load("data/processed/4-peak-stats.jld")
peak_deployment = peakstats["peak_deployment"]

offpeakstats = JLD.load("data/processed/4-offpeak-stats.jld")
offpeak_deployment = offpeakstats["offpeak_deployment"]

const model_names = (:Stochastic, :Robust01, :Robust005, :Robust001, :Robust0001, :Robust00001, :MEXCLP, :MALP)

using Distributions

const turnard = Distributions.LogNormal(3.65, 0.3)
const volatile_turnard = Distributions.LogNormal(3.57, 0.5)
const steady_turnard = Distributions.LogNormal(3.69, 0.1)


const stn_names = [Symbol("stn$(i)_min") for i in 1:size(p.coverage,2)];

# overall simulation
@time begin
    p = DeploymentProblem(30, length(locations), length(regions), demand, indices[train_filter],
                                  indices[test_filter], coverage[regions,:], Array{Bool,2}(adjacent))
    for turnaround in (turnard,)
        problem = DispatchProblem(calls[inc_test_indices,:], p.coverage)
        problem.emergency_calls[:arrival_seconds] = 0
        arrival_sec = cumsum(problem.emergency_calls[:interarrival_seconds])
        problem.emergency_calls[:arrival_seconds] = arrival_sec;

        for name in model_names[1:2]
            print("$name: ")
            for namb in 45:5:50 # 10:5:20 #50
                print("$namb ")
                x = amb_deployment[name][namb]
                p.nambulances = namb

                initialize!(problem, x)
                dispatch = ClosestDispatch(p, problem.emergency_calls[:, stn_names], x)

                srand(1234) # reset seed
                simulate_events!(problem, dispatch, turnaround, "$(name)_$(namb)", false)
            end
            println()
        end
        #DataFrames.writetable("data/processed/4-overall_simulation.csv", problem.emergency_calls)
    end
end
