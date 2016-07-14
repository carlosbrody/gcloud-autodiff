function trialdata(ratdata, trial)
    if ratdata["rawdata"]["pokedR"][trial] > 0
        rat_choice = "R"
    else
        rat_choice = "L"
    end;
    
    return ratdata["rawdata"]["rightbups"][trial], ratdata["rawdata"]["leftbups"][trial], 
    ratdata["rawdata"]["T"][trial], rat_choice
end

# RightClickTimes, LeftClickTimes, maxT, rat_choice = trialdata(ratdata, 1)

"""
LL, LLgrad, LLhesian = multiLikey(ratdata, params, set_of_trials)

Returns the summed log likelihood, LL gradient, and hessian for the trials of ratdata indicated by set_of_trials

"""
rat_choice = []

function multiLikey(ratdata, params, set_of_trials)

    LL        = zeros(1, length(set_of_trials))
    LLgrad    = zeros(length(params), length(set_of_trials))
    LLhess    = zeros(length(params), length(params), length(set_of_trials))
    
    for i in 1:length(set_of_trials)
        RClickTimes, LClickTimes, maxT, rat_choice = trialdata(ratdata, set_of_trials[i])   

        # return rat_choice
        # matwrite("trash.mat", Dict("rat_choice" => rat_choice))

        global dx
        dx = 0.25
        LLi, LLgradi, LLhessi = llikey(params, rat_choice, maxT=maxT, RightPulseTimes=RClickTimes,
                                          LeftPulseTimes = LClickTimes)

        LL[i]         = LLi;
        LLgrad[:,i]   = LLgradi;
        LLhess[:,:,i] = LLhessi;

        # println(i)
    end

    if length(set_of_trials)==1
        LLgrad = squeeze(LLgrad,2)
        LLhess = squeeze(LLhess,3)
    end
        
    return LL, LLgrad, LLhess
end

