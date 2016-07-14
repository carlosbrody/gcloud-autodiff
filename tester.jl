@everywhere include("startup.jl")

ratdata   = matread("Data/chrono_rawdata.mat")
trial_ids = matread("Data/trial_ids.mat")

ps = [];
params = [];

RClickTimes, LClickTimes, maxT, rat_choice = trialdata(ratdata, 1)   

ps = matread("Params/params.mat")
params = [ps["sigma_a"], ps["sigma_s"], ps["sigma_i"], ps["lam"], 
          ps["B"], ps["bias"], ps["phi"], ps["tau_phi"], ps["lapse"]];

LL = @spawn multiLikey(ratdata, params, collect(1))
# LL = multiLikey(ratdata, params, collect(1))

