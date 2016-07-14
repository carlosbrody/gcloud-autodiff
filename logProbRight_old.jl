"""
function logProbRight(params::Vector)

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
    sigma2_a = params[1]; sigma2_s = params[2]; sigma2_i = params[3]; 
    lambda = params[4]; B = params[5]; bias = params[6]; 
    phi = params[7]; tau_phi = params[8]; lapse = params[9]
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
    c_eff   = 1; c_trace[1] = c_eff;
    
    Fi = Fmatrix(collect([sigma2_i/dt; lambda; 0.0]), bin_centers); a = Fi*a;
    a_trace[:,1] = a;

    F0 = Fmatrix(collect([sigma2_a; lambda; 0.0]), bin_centers)
    for i in 1:Nsteps-1 
        c_eff = 1 + (c_eff - 1)*exp(-dt/tau_phi)
        c_trace[i+1]   = convert(Float64, c_eff)
        if (RightClicks[i]==0) & (LeftClicks[i]==0)
            a = F0*a
        elseif (RightClicks[i]==1) & (LeftClicks[i]==1)
            c_eff = 0
            a = F0*a
        else
            net_sigma2 = sigma2_a + (sigma2_s*c_eff/total_rate)/dt
            F = Fmatrix(collect([net_sigma2; lambda; (RightClicks[i] - LeftClicks[i])*c_eff/dt]), bin_centers)
            a = F*a
            c_eff = c_eff*phi
        end
        a_trace[:,i+1] = convert(Array{Float64}, a)
    end;
    pright = sum(a[binBias+1:end]) + 
        a[binBias]*0.5*(dx/2 - (bias - bin_centers[binBias]))/(dx/2)
    return log(pright)
end



