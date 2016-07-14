""" 
function (LL, LLgrad, LLhessian, bin_centers, bin_times, a_trace) = 
llikey(params, subject_choice, maxT=1, RightPulseTimes=[], LeftPulseTimes=[], dx=0.25, dt=0.02)

Computes the log likelihood according to Bing's model, and returns log likelihood, gradient, and hessian. 
Mutates the globals bin_centers, bin_times, a_trace, RightClickTimes, LeftClickTimes, and Nsteps

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

subject_choice     should be either "R" or "L"


RETURNS:


"""
function llikey(params::Vector, subject_choice; maxT=1, RightPulseTimes=[], LeftPulseTimes=[], dx=0.25, dt=0.02)

    global RightClickTimes, LeftClickTimes, Nsteps

    B = params[5];

    RightClickTimes = filter(x -> x<=maxT, collect(RightPulseTimes))
    LeftClickTimes  = filter(x -> x<=maxT, collect(LeftPulseTimes))    
    Nsteps = Int(ceil(maxT/dt))
    binN = Int(ceil(B/dx)); bin_centers = make_bins(B, dx, binN); 
    bin_times = dt*collect(0:Nsteps);
 
    if subject_choice == "R"
        rat_choice = 1;
    elseif subject_choice == "L"
        rat_choice = 0;
    else
        error("Rat did what?? It was neither R nor L")
    end

    # introduce rat_choice as if it were another parameter, but we'll ignore the derivative with respect to it:
    params = [params ; rat_choice]

    result =  GradientResult(params)
    
    ForwardDiff.gradient!(result, logLike, params)

    LL        = ForwardDiff.value(result); 
    LLgrad    = ForwardDiff.gradient(result); 
    LLhessian = 0; # ForwardDiff.hessian(result);

    # Now delete the derivative w.r.t. rat choice from the gradient and hessian results
    LLgrad = LLgrad[1:length(LLgrad)-1]
    # LLhessian = LLhessian[1:length(LLgrad), 1:length(LLgrad)]

    return LL, LLgrad, LLhessian

    # rat_choice = subject_choice
    # LL = logLike(params)
    
    # return LL, zeros(size(params)), zeros(length(params), length(params)) # LLgrad, LLhessian
end

