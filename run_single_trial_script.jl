include("startup.jl")

global a, dx, dt, Nsteps, LeftClickTimes, RightClickTimes, a_trace, c_trace, total_rate

ps = matread("Params/params.mat")
params = [ps["lam"], ps["sigma2_a"], ps["sigma2_s"], ps["sigma2_i"],  
          ps["B"], ps["phi"], ps["tau_phi"], ps["bias"], ps["lapse"]];

test_data   = matread("Data/test_data.mat")
RCs             = test_data["RightClickTimes"]
LCs             = test_data["LeftClickTimes"]
total_rate      = test_data["total_rate"]
maxT            = test_data["maxT"]
dx              = test_data["dx"]
dt              = test_data["dt"]
side_choice     = test_data["side_choice"]
Nsteps          = convert(Int, test_data["Nsteps"])

LL, LLgrad, LLhess = llikey(params, side_choice, maxT=maxT, dx=dx, dt=dt, 
                            RightPulseTimes=RCs, LeftPulseTimes=LCs)


my_c_trace = zeros(size(c_trace))
for i=1:length(c_trace) my_c_trace[i] = convert(Float64, c_trace[i]) end

matwrite("Results/results.mat", Dict("params" => ps, "LL" => LL, "LLgrad" => LLgrad, "LLhess" => LLhess,
                                     "a_trace" => a_trace, "c_trace" => my_c_trace))

matwrite("Results/results_single.mat", Dict("maxT" => maxT, "rat_choice" => side_choice, "RClickTimes" => RCs,
                                             "LClickTimes" => LCs, "params" => params, "LL" => LL, "LLgrad" => LLgrad))
