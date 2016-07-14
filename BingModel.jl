
import ForwardDiff
import PyPlot
using PyPlot
import Base.convert

convert(::Type{Float64}, x::ForwardDiff.GradientNumber) = Float64(x.value)
function convert(::Array{Float64}, x::Array{ForwardDiff.GradientNumber}) 
    y = zeros(size(x)); 
    for i in 1:prod(size(x)) 
        y[i] = convert(Float64, x[i]) 
    end
    return y
end


convert(::Type{Float64}, x::ForwardDiff.HessianNumber) = Float64(x.gradnum.value)
function convert(::Array{Float64}, x::Array{ForwardDiff.HessianNumber}) 
    y = zeros(size(x)); 
    for i in 1:prod(size(x)) 
        y[i] = convert(Float64, x[i]) 
    end
    return y
end



"""
function bin_centers = make_bins(B, dx, binN)

Makes a series of points that will indicate bin centers. The first and
last points will indicate sticky bins. No "bin edges" are made-- the edge
between two bins is always implicity at the halfway point between their
corresponding centers. The center bin is always at x=0; bin spacing
(except for last and first bins) is always dx; and the position
of the first and last bins is chosen so that |B| lies exactly at the
midpoint between 1st (sticky) and 2nd (first real) bins, as well as
exactly at the midpoint between last but one (last real) and last
(sticky) bins.

Playing nice with ForwardDiff means that the *number* of bins must be predetermined.
So this function will not actually set the number of bins; what it'll do is determine their
locations. To accomplish this separation, the function uses as a third parameter binN,
which should be equal to the number of bins with bin centers > 0, as follows: 
   binN = ceil(B/dx)
and then the total number of bins will be 2*binN+1, with the center one always corresponding
to position zero. Use non-differentiable types for B and dx for this to work.
"""
function make_bins(B, dx, binN)
    bins = collect(1.0:binN)*B
    bins = dx*bins/B

    if bins[end] == B
        bins[end] = B + dx
    else
        bins[end] = 2*B - bins[end-1]
    end

    bins = [-bins[end:-1:1]; 0; bins]
    return bins
end


"""
function F = Fmatrix([sigma, lambda, c], bin_centers)

Uses globals
    dt
    epsilon       (=10.0^-10)

Returns a square Markov matrix of transition probabilities. 

sigma  should be in (accumulator units) per (second^(1/2))
lambda should be in s^-1
c      should be in accumulator units per second
bin_centers should be a vector of the centers of all the bins. Edges will be at midpoints
       between the centers, and the first and last bin will be sticky.

dx is not used inside Fmatrix, because bin_centers specifies all we need to know.
dt *is* used inside Fmatrix, to convert sigma, lambda, and c into timestep units
"""
function Fmatrix(params::Vector, bin_centers)
    sigma = params[1];
    lam   = params[2];
    c     = params[3];

    sigma = sqrt(dt)*sigma;
    sigma_sbin = convert(Float64, sigma)

    F = collect(1.0:length(bin_centers))*collect(1.0:length(bin_centers))';
    F = 0.0*sigma*F; # Multiplying by that sigma is needed, 
                     # for type casting reasons I do not understand...

    mus      = (bin_centers + c/lam)*exp(lam*dt) - c/lam
    sbinsize = 0.1*sigma_sbin
    swidth   = 4*sigma_sbin
    sbins    = collect(-swidth:sbinsize:swidth+epsilon)
    ps       = exp(-sbins.^2/(2*sigma^2)) / sqrt(2*sigma^2)
    ps       = ps/sum(ps);

    base_sbins = sbins;
    
    for j in 1:length(bin_centers)
        # sbins = [mus[j] - swidth : sbinsize : mus[j] + swidth+epsilon]
        sbins = collect(0:(length(base_sbins)-1))*sbinsize
        sbins = sbins + mus[j]-swidth
        for k in 1:length(sbins)
            if sbins[k] < (bin_centers[1] + bin_centers[2])/2
                F[1,j] = F[1,j] + ps[k]
            elseif (bin_centers[end]+bin_centers[end-1])/2 <= sbins[k]
                F[end,j] = F[end,j] + ps[k]
            else
                bot = find(bin_centers .<= sbins[k])[end]
                top = bot+1
                F[bot,j] = F[bot,j] + 
                    ps[k]*(bin_centers[top] - sbins[k])/(bin_centers[top] - bin_centers[bot])
                F[top,j] = F[top,j] + 
                    ps[k]*(sbins[k] - bin_centers[bot])/(bin_centers[top] - bin_centers[bot])
            end
        end
    end
    F[:,1] = 0; F[:,end] = 0; F[1,1] = 1; F[end,end] = 1;
    return F
end



"""
function logProbRight(params::Vector)

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
    phi = params[7]; tau_phi = params[8]; lapse = params[9]

Returns the log of the probability that the agent chose Right. 
"""

function logProbRight(params::Vector)
    sigma_a = params[1]; sigma_s = params[2]; sigma_i = params[3]; 
    lambda = params[4]; B = params[5]; bias = params[6]; 
    phi = params[7]; tau_phi = params[8]; lapse = params[9]
    global a, dx, dt, Nsteps, LeftClickTimes, RightClickTimes, a_trace, c_trace
    
    LeftClicks  = zeros(1, Nsteps); if isempty(RightClickTimes) RightClickTimes = zeros(0) end;
    RightClicks = zeros(1, Nsteps); if isempty(LeftClickTimes ) LeftClickTimes  = zeros(0) end;
    for i in ceil((LeftClickTimes+epsilon)/dt)  LeftClicks[Int(i)]  = 1 end
    for i in ceil((RightClickTimes+epsilon)/dt) RightClicks[Int(i)] = 1 end
    
    my_B = convert(Float64, B) # my_B won't be differentiated; ForwardDiff can't do ceil()
    my_bias = convert(Float64, bias)  # my_bias won't be differentiated' FD can't do floor()
    binN = Int(ceil(my_B/dx))  
    binBias = Int(floor(my_bias/dx)) + binN+1  
    bin_centers = make_bins(B, dx, binN) 

    a_trace = zeros(length(bin_centers), Nsteps+1); c_trace = zeros(1, Nsteps+1)
    a = zeros(length(bin_centers),1)*sigma_a*0.0; # That weirdo inexact error thing
    a[binN+1] = 1-2*lapse; a[1] = lapse; a[end] = lapse;
    c_eff   = 1;
    
    Fi = Fmatrix(collect([sigma_i/sqrt(dt); lambda; 0.0]), bin_centers); a = Fi*a;
    a_trace[:,1] = a;

    F0 = Fmatrix(collect([sigma_a; lambda; 0.0]), bin_centers)
    for i in 1:Nsteps 
        c_eff = 1 + (c_eff - 1)*exp(-dt/tau_phi)
        if (RightClicks[i]==0) & (LeftClicks[i]==0)
            a = F0*a
        elseif (RightClicks[i]==1) & (LeftClicks[i]==1)
            c_eff = 0
            a = F0*a
        else
            net_sigma = sqrt(sigma_a^2 + (sigma_s*c_eff)^2/dt)
            F = Fmatrix(collect([net_sigma; lambda; (RightClicks[i] - LeftClicks[i])*c_eff/dt]), bin_centers)
            a = F*a
            c_eff = c_eff*phi
        end
        c_trace[i+1]   = convert(Float64, c_eff)
        a_trace[:,i+1] = convert(Array{Float64}, a)
    end;
    pright = sum(a[binBias+1:end]) + 
        a[binBias]*0.5*(dx/2 - (bias - bin_centers[binBias]))/(dx/2)
    return log(pright)
end



function logLike(params::Vector)
    if rat_choice == "R"
        return logProbRight(params)
    elseif rat_choice == "L"
        return log(1 - exp(logProbRight(params)))
    else
        error("Rat did what?? It was neither R nor L")
    end
end


epsilon = 10.0^(-10); dx = 0.125; dt = 0.02; Nsteps = Int(ceil(1.0/dt))

""" 
function (LL, LLgrad, LLhessian, bin_centers, bin_times, a_trace) = 
    llikey(params, rat_choice, maxT=1, RightPulseTimes=[], LeftPulseTimes=[], dx=0.25, dt=0.02)

Computes the log likelihood according to Bing's model, and returns log likelihood, gradient, and hessian

params is a vector whose elements, in order, are
    sigma_a    square root of accumulator variance per unit time sqrt(click units^2 per second)
    sigma_s    standard deviation introduced with each click (will get scaled by click adaptation)
    sigma_i    square root of initial accumulator variance sqrt(click units^2)
    lambda     1/accumulator time constant (sec^-1). Positive means unstable, neg means stable
    B          sticky bound height (click units)
    bias       where the decision boundary lies (click units)
    phi        click adaptation/facilitation multiplication parameter
    tau_phi    time constant for recovery from click adaptation (sec)
    lapse      2*lapse fraction of trials are decided randomly

rat_choice     should be either "R" or "L"


RETURNS:


"""
function llikey(params::Vector, subject_choice; maxT=1, RightPulseTimes=[], LeftPulseTimes=[], dx=0.25, dt=0.02)

    global RightClickTimes, LeftClickTimes, Nsteps

    RightClickTimes = filter(x -> x<=maxT, collect(RightPulseTimes))
    LeftClickTimes  = filter(x -> x<=maxT, collect(LeftPulseTimes))    
    Nsteps = Int(ceil(maxT/dt))
    binN = Int(ceil(params[5]/dx)); bin_centers = make_bins(params[5], dx, binN); 
    bin_times = dt*collect(0:Nsteps);
    
    # LLhessian, allresults = ForwardDiff.hessian(logLike, params, ForwardDiff.AllResults)

    # LL     = ForwardDiff.value(allresults)
    # LLgrad = ForwardDiff.gradient(allresults)

    rat_choice = subject_choice
    LL = logLike(params)
    
    return LL, zeros(size(params)), zeros(length(params), length(params)) # LLgrad, LLhessian
end



sigma_a = 1; sigma_s = 0.1; sigma_i = 0.2; 
sigma_a_sbin = sigma_a  # remember we need this copy for Fmatrix
lam = -0.0005; B = 4.1; bias = 0.1; 
phi = 0.3; tau_phi = 0.1; lapse = 0.05;
params = [sigma_a, sigma_s, sigma_i, lam, B, bias, phi, tau_phi, lapse]   
rat_choice = "R"

# LL, LLgrad, LLhess = llikey(params, "R", RightPulseTimes=[0.2 0.4]);

LL, LLgrad, LLhess = llikey(params, "R", RightPulseTimes=[0.2 0.4]);

display(LL)


