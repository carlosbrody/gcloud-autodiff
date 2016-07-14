include("startup.jl")

ratdata   = matread("Data/B069_10_trials.mat")
trial_ids = matread("Data/trial_ids.mat"); trial_ids = trial_ids["trial_ids"]

global total_rate
total_rate      = 40

ps = matread("Params/params.mat")
mv("Params/params.mat", "Old/params.mat", remove_destination=true)

params = [ps["lam"], ps["sigma2_a"], ps["sigma2_s"], ps["sigma2_i"],  
          ps["B"], ps["phi"], ps["tau_phi"], ps["bias"], ps["lapse"]];

RClickTimes, LClickTimes, maxT, rat_choice = trialdata(ratdata, 1)
# LL, LLgrad, LLhess = llikey(params, rat_choice, maxT=maxT, 
#                            RightPulseTimes=RClickTimes, LeftPulseTimes=LClickTimes)

test_data   = matread("Data/test_data.mat")
# RCs             = test_data["RightClickTimes"]
# LCs             = test_data["LeftClickTimes"]
# total_rate      = test_data["total_rate"]
# maxT            = test_data["maxT"]
# dx              = test_data["dx"]
# dt              = test_data["dt"]
# side_choice     = test_data["side_choice"]
# Nsteps          = convert(Int, test_data["Nsteps"])

# dx = 0.25
# dt = 0.02

LL, LLgrad, LLhess = multiLikey(ratdata, params, collect(1))


matwrite("Results/results_sen.mat", Dict("LL" => LL, "LLgrad" => LLgrad))
 
                                   #      "maxT" => maxT, "rat_choice" => rat_choice, "RClickTimes" => RClickTimes,
                                   #      "LClickTimes" => LClickTimes, "params" => params, "dx" => dx, "dt" => dt))


