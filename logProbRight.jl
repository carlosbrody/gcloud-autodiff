"""
function aL, aR = make_adapted_cat_clicks(phi, tau_phi)

Given the times of the clicks, computes their effective weights after adaptation. This function implements
across-stream adaptation (i.e., a pulse on one side affects pulses on the other side also) as used in the main text 
of Brunton et al. 2013.  Returns vectors that are the same size as RichtClickTimes and LeftClickTimes. Each entry
indicates the effective weight of the corresppinding pulse

Needs to read globals 
    RightClickTimes   vector with elements indicating times of right clicks
    LeftClickTimes    vector with elements indicating times of left clicks
    epsilon           smallest usable number, typically something like 1e-10

"""
function make_adapted_cat_clicks(phi, tau_phi; cross_side_suppression=0)

    global LeftClickTimes, RightClickTimes, epsilon

    if !isdefined(:epsilon) epsilon = 1e-10 end

    if abs(phi - 1) < epsilon
        aL = ones(size(LeftClickTimes))*phi*0.0 + 1   # weird typecasting error where ForwardDiff needs the *phi*0.0
        aR = ones(size(RightClickTimes))*phi*0.0 + 1

        return aL, aR
    else
        aL = ones(size(LeftClickTimes))*phi*0.0 + phi
        aR = ones(size(RightClickTimes))*phi*0.0 + phi
    end


    if length(LeftClickTimes) > 1
        lefts  = [LeftClickTimes[:]  -ones(length(LeftClickTimes), 1)*phi*0.0 - 1]
    else
        lefts  = [LeftClickTimes     -ones(length(LeftClickTimes), 1)*phi*0.0 - 1]
        
    end
    if length(RightClickTimes) > 1
        rights = [RightClickTimes[:]  ones(length(RightClickTimes), 1)*phi*0.0 + 1]
    else
        rights = [RightClickTimes     ones(length(RightClickTimes), 1)*phi*0.0 + 1]
    end


    allbups = sortrows([lefts; rights])'   # one pulse in each column, first row is time second row has side bup was on
    ici = diff(allbups[1,:]')'
    adapted = ones(1, size(allbups,2))*phi*0.0 + 1

    for i in 2:size(allbups,2)
        if ici[i-1] <= cross_side_suppression
            adapted[i-1] = 0;
            adapted[i] = 0;
        else
            last = tau_phi * log(1 - adapted[i-1]*phi);
            adapted[i] = 1 - exp((-ici[i-1] + last)/tau_phi)
        end
    end

    adapted = real(adapted);

    aL = adapted[allbups[2,:].==-1]';
    aR = adapted[allbups[2,:].==+1]';

    return aL, aR
end


""" 
function net_input, tot_input, nclicks = make_click_inputs35(t, leftbups, rightbups, clicks_L, clicks_R, NL, NR)

Compute the net input (right-left) and total input (right+left) of
clicks for all time steps.  


OUTPUTS
--------
net_input: array of length N with net input of clicks at each time step
    (right, weighted clicks - left, weighted clicks at each time step)

tot_input: array of length N with total input of clicks at each time step
    (right, weighted clicks + left, weighted clicks at each time step)

INPUTS
---------
t: time axis
leftbups:  times of left clicks
rightbups: times of right clicks
clicks_L:  adapted weights of left clicks
clicks_R:  adapted weights of right clicks
phi:       adaptation phi, one of the parameters to be differentiated w.r.t., here purely for weird hack purposes.
NL: 
NR:

"""
function make_click_inputs35(t, leftbups, rightbups, clicks_L, clicks_R, phi; NL=[], NR=[])

    if isempty(NL)   NL = ones(size(leftbups))*phi*0.0 + 1;  end
    if isempty(NR)   NR = ones(size(rightbups))*phi*0.0 + 1; end

    here_L = qfind(t, leftbups);
    here_R = qfind(t, rightbups);
    nclicks = zeros(length(t),1)*phi*0.0;
    for i = 1:length(here_L),
        nclicks[here_L[i]] = nclicks[here_L[i]] + NL[i];
    end;
    for i = 1:length(here_R),
        nclicks[here_R[i]] = nclicks[here_R[i]] + NR[i];
    end;


    net_input = zeros(length(t),1)*phi*0.0;
    tot_input = zeros(length(t),1)*phi*0.0;
    for i = 1:length(leftbups),
        net_input[here_L[i]] = net_input[here_L[i]] - clicks_L[i];
        tot_input[here_L[i]] = tot_input[here_L[i]] + clicks_L[i];
    end;
    for i = 1:length(rightbups),
        net_input[here_R[i]] = net_input[here_R[i]] + clicks_R[i];
        tot_input[here_R[i]] = tot_input[here_R[i]] + clicks_R[i];
    end;

    return net_input, tot_input, nclicks
end



"""
function lR = logProbRight(params::Vector)

Needs to read globals 
    Nsteps number of timesteps to simulate
    RightClickTimes   vector with elements indicating times of right clicks
    LeftClickTimes    vector with elements indicating times of left clicks
    dx                accumualtor bin width
    dt                timestep
    total_rate        Left + Right Poisson Clicks rate, sued to scale sigma2_s so that sensory
                      and accumulator sigmas are on simialr scales

Mutates globals
    a      (column vector representing distribution of values of accumulator a)

    a_trace (length(bin_centers)-by-Nsteps+1), a trace of the distribution of a as 
            a function of time
    c_trace (row vector Nsteps+1 long, effective value of c as 
            a function of time after adaptation)

Takes params
    sigma2_a = params[1]; sigma2_s = params[2]; sigma2_i = params[3]; 
    lambda = params[4]; B = params[5]; bias = params[6]; 
    phi = params[7]; tau_phi = params[8]; lapse = params[9]

Returns the log of the probability that the agent chose Right. 
"""
function logProbRight(params::Vector)
    lambda = params[1]; sigma2_a = params[2]; sigma2_s = params[3]; sigma2_i = params[4]; 
    B = params[5]; phi = params[6]; tau_phi = params[7]; bias = params[8]; lapse = params[9]
    global a, dx, dt, Nsteps, total_rate, LeftClickTimes, RightClickTimes, a_trace, c_trace
    
    RightClicks = zeros(1, Nsteps); if isempty(RightClickTimes) RightClickTimes = zeros(0) end;
    LeftClicks  = zeros(1, Nsteps); if isempty(LeftClickTimes ) LeftClickTimes  = zeros(0) end;
    for i in ceil((LeftClickTimes+epsilon)/dt)  LeftClicks[Int(i)]  = 1 end
    for i in ceil((RightClickTimes+epsilon)/dt) RightClicks[Int(i)] = 1 end
    
    my_B = convert(Float64, B) # my_B won't be differentiated; ForwardDiff can't do ceil()
    my_bias = convert(Float64, bias)  # my_bias won't be differentiated' FD can't do floor()
    binN = Int(ceil(my_B/dx))  
    binBias = Int(floor(my_bias/dx)) + binN+1  
    bin_centers = make_bins(B, dx, binN) 

    a_trace = zeros(length(bin_centers), Nsteps); c_trace = zeros(1, Nsteps)
    a = zeros(length(bin_centers),1)*sigma2_a*0.0; # That weirdo inexact error thing
    a[binN+1] = 1-lapse; a[1] = lapse/2; a[end] = lapse/2;
    # c_eff   = 1; c_trace[1] = c_eff;

    t = dt*(0:Nsteps-1)
    aL, aR = make_adapted_cat_clicks(phi, tau_phi)

    net_input, tot_input, nclicks = make_click_inputs35(t, LeftClickTimes, RightClickTimes, aL, aR, phi)
    c_trace = net_input

    Fi = Fmatrix(collect([sigma2_i/dt; lambda; 0.0]), bin_centers); a = Fi*a;
    a_trace[:,1] = a;

    F0 = Fmatrix(collect([sigma2_a; lambda; 0.0]), bin_centers)
    for i in 1:Nsteps-1 
        if (abs(net_input[i]) < epsilon)
            a = F0*a
        else
            net_sigma2 = sigma2_a + (sigma2_s*tot_input[i]/total_rate)/dt
            # net_sigma2 = sigma2_a + (sigma2_s*abs(net_input[i])/total_rate)/dt
            F = Fmatrix(collect([net_sigma2; lambda; net_input[i]/dt]), bin_centers)
            a = F*a
        end
        a_trace[:,i+1] = convert(Array{Float64}, a)
    end;

    pright = sum(a[binBias+1:end]) + 
        a[binBias]*0.5*(dx/2 - (bias - bin_centers[binBias]))/(dx/2)

    return log(pright)
end



