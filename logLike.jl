"""
function logLike(params::Vector)

Needs to read globals 
    Nsteps number of timesteps to simulate
    RightClickTimes   vector with elements indicating times of right clicks
    LeftClickTimes    vector with elements indicating times of left clicks


Mutates globals
    a      (column vector representing distribution of values of accumulator a)

    a_trace (length(bin_centers)-by-Nsteps+1), a trace of the distribution of a as 
            a function of time
    c_trace (row vector Nsteps+1 long, effective value of c as 
            a function of time after adaptation)

Takes params
    sigma_a = params[1]; sigma_s = params[2]; sigma_i = params[3]; 
    lambda = params[4]; B = params[5]; bias = params[6]; 
    phi = params[7]; tau_phi = params[8]; lapse = params[9]; 
    rat_choice = params[10]

For rat_choice, +1 means went right, -1 means went left.

Returns the log of the probability that the agent did what the rat did. 
"""
function logLike(params::Vector)
    rat_choice = params[end]

    # display(rat_choice)
    
    return rat_choice*logProbRight(params[1:end-1]) + (1-rat_choice)*log(1 - exp(logProbRight(params[1:end-1])))
end
