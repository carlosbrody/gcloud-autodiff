"""
function F = Fmatrix([sigma, lambda, c], bin_centers)

Uses globals
    dt
    epsilon       (=10.0^-10)

Returns a square Markov matrix of transition probabilities. 

sigma2  should be in (accumulator units^2) per second
lambda should be in s^-1
c      should be in accumulator units per second
bin_centers should be a vector of the centers of all the bins. Edges will be at midpoints
       between the centers, and the first and last bin will be sticky.

dx is not used inside Fmatrix, because bin_centers specifies all we need to know.
dt *is* used inside Fmatrix, to convert sigma, lambda, and c into timestep units
"""
function Fmatrix(params::Vector, bin_centers)
    sigma2 = params[1];
    lam    = params[2];
    c      = params[3];

    sigma2 = dt*sigma2;
    sigma2_sbin = convert(Float64, sigma2)

    F = collect(1.0:length(bin_centers))*collect(1.0:length(bin_centers))';
    F = 0.0*sigma2*F; # Multiplying by that sigma is needed, 
                     # for type casting reasons I do not understand...

    mus      = (bin_centers + c/lam)*exp(lam*dt) - c/lam
    sbinsize = 0.1*sqrt(sigma2_sbin)
    swidth   = 5*sqrt(sigma2_sbin)
    sbins    = collect(-swidth:sbinsize:swidth+epsilon)
    ps       = exp(-sbins.^2/(2*sigma2)) / sqrt(2*sigma2)
    ps_o     = ps/sum(ps);

    base_sbins = sbins;
    sbins_o = collect(0:(length(base_sbins)-1))*sbinsize
    
    N = length(bin_centers)
    B_bottom = (bin_centers[1]   + bin_centers[2])/2
    B_top    = (bin_centers[end] + bin_centers[end-1])/2
    for j in 2:N-1
        sbins = sbins_o + mus[j]-swidth
        ps    = ps_o
        dd1   = bin_centers[2] - bin_centers[1]
        ddN   = bin_centers[N] - bin_centers[N-1]
        dd    = dx
        
        for k=1:length(sbins)
            if B_top <= sbins[k]
                F[end,j] = F[end,j] + ps[k]
            elseif sbins[k] < B_bottom
                F[1,j]   = F[1,j]   + ps[k]
            elseif (bin_centers[1] < sbins[k]) & (sbins[k] < bin_centers[2])
                F[1,j]   = F[1,j] + ps[k]*(bin_centers[2] - sbins[k])/dd1
                F[2,j]   = F[2,j] + ps[k]*(sbins[k] - bin_centers[1])/dd1
            elseif (bin_centers[N-1] < sbins[k]) & (sbins[k] < bin_centers[N])
                F[N-1,j]   = F[N-1,j] + ps[k]*(bin_centers[N] - sbins[k])/ddN
                F[N,  j]   = F[N,  j] + ps[k]*(sbins[k] - bin_centers[N-1])/ddN
            else
                bot = floor(Int, (sbins[k] - bin_centers[2])/dx) + 2
                top =  ceil(Int, (sbins[k] - bin_centers[2])/dx) + 2
                if bot == top
                    F[bot, j] = F[bot,j] + ps[k]
                else
                    F[bot,j] = F[bot,j] + ps[k]*(bin_centers[top] - sbins[k])/dd
                    F[top,j] = F[top,j] + ps[k]*(sbins[k] - bin_centers[bot])/dd
                end
            end
        end        
    end
    
    F[:,1] = 0; F[:,end] = 0; F[1,1] = 1; F[end,end] = 1;
    return F
end

