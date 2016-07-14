@everywhere include("startup.jl")

display("---- reading rat data ----")
ratdata   = matread("Data/chrono_rawdata.mat")
display("    data read")
ratdata   = matread("Data/chrono_B069_rawdata.mat")
# ratdata   = matread("Data/B069_10_trials.mat")
trial_ids = matread("Data/trial_ids.mat"); trial_ids = trial_ids["trial_ids"]

global total_rate, dx
total_rate      = 40
dx              = 0.25

ps = [];
params = [];

d = readdir("Params")
if "params.mat" in d
    ps = matread("Params/params.mat")
    params = [ps["lam"], ps["sigma2_a"], ps["sigma2_s"], ps["sigma2_i"],  
              ps["B"], ps["phi"], ps["tau_phi"], ps["bias"], ps["lapse"]];
end


RClickTimes, LClickTimes, maxT, rat_choice = trialdata(ratdata, 1)   

testing = false; if testing 
    ps = matread("Params/params.mat")
    # mv("Params/params.mat", "Old/params.mat")

    params = [ps["lam"], ps["sigma2_a"], ps["sigma2_s"], ps["sigma2_i"],  
              ps["B"], ps["phi"], ps["tau_phi"], ps["bias"], ps["lapse"]];

    LL = run_trials(ratdata, params, trial_ids, parallel=true)
    fetch(LL[1])

    matwrite("Results/results.mat", Dict("params" => ps, "trial_ids" => trial_ids, "LL" => LL[1][1], 
                                             "LLgrad" => LL[1][2], "LLhess" => LL[1][3]))
    return
end



while true
    while length(readdir("Params"))==0    
        sleep(1)
        display("sleeping")
    end

    d = readdir("Params")

    if "params.mat" in d
        ps = matread("Params/params.mat")
        mv("Params/params.mat", "Old/params.mat", remove_destination=true)

        params = [ps["lam"], ps["sigma2_a"], ps["sigma2_s"], ps["sigma2_i"],  
          ps["B"], ps["phi"], ps["tau_phi"], ps["bias"], ps["lapse"]];

        LL = run_trials(ratdata, params, trial_ids, parallel=true)
        fetch(LL)

        matwrite("Results/results.mat", Dict("params" => ps, "trial_ids" => trial_ids, "LL" => LL[1], 
                                             "LLgrad" => LL[2], "LLhess" => LL[3])) 

        # RClickTimes, LClickTimes, maxT, rat_choice = trialdata(ratdata, 1)
        # LLz, LLgradz, LLhessz = llikey(params, rat_choice, maxT=maxT, 
        #                    RightPulseTimes=RClickTimes, LeftPulseTimes=LClickTimes)

        # matwrite("Results/results_sen.mat", Dict("LLz" => LLz, "LLgradz" => LLgradz, 
        #                                     "maxT" => maxT, "rat_choice" => rat_choice, "RClickTimes" => RClickTimes,
        #                                     "LClickTimes" => LClickTimes, "params" => params))

        display("emerged from an LL round")
    elseif "stop" in d
        mv("Params/stop", "Old/stop", remove_destination=true)
        quit()

    elseif "reread_ratdata" in d
        mv("Params/reread_ratdata", "Old/reread_ratdata", remove_destination=true)

        ratdata   = matread("Data/chrono_rawdata.mat")

        display("emerged from a round of reading ratdata")

    elseif "reread_trial_ids" in d
        mv("Params/reread_trial_ids", "Old/reread_trial_ids")

        trial_ids = matread("Data/trial_ids.mat"); trial_ids = trial_ids["trial_ids"]

        display("emerged from a round of reading trial_ids")

    elseif "reread_ratdata_and_trial_ids" in d
        mv("Params/reread_ratdata_and_trial_ids", "Old/reread_ratdata_and_trial_ids")

        ratdata   = matread("Data/chrono_rawdata.mat")
        trial_ids = matread("Data/trial_ids.mat"); trial_ids = trial_ids["trial_ids"]

        display("emerged from a round of reading ratdata and trial_ids")
    end

    sleep(0.5)
end


 
