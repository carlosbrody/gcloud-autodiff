import ForwardDiff
# import PyPlot
# using PyPlot
import Base.convert
# using Optim

include("BM.jl")
include("Fmatrix.jl")
include("logProbRight.jl")
include("logLike.jl")
include("llikey.jl")
include("multiLikey.jl")
include("run_trials.jl")
include("data.jl")


epsilon = 10.0^(-10); dx = 0.125; dt = 0.02; Nsteps = Int(ceil(1.0/dt))

# sigma_a = 1; sigma_s = 0.1; sigma_i = 0.002; 
# sigma_a_sbin = sigma_a  # remember we need this copy for Fmatrix
# lambda = -0.0005; B = 4.1; bias = 0.1; 
# phi = 0.3; tau_phi = 0.1; lapse = 0.05;
# params = [sigma_a, sigma_s, sigma_i, lambda, B, bias, phi, tau_phi, lapse];   
rat_choice = "R"

a = 0
a_trace = 0

# LL, LLgrad, LLhess = llikey(params, "R", RightPulseTimes=[0.2 0.4])
# imshow(log(abs(LLhess)), interpolation="none")

# display(LL)
# display(LLgrad)
# display(LLhess)

using MAT
# if !isdefined(:ratdata)
#    ratdata = matread("Data/chrono_B069_rawdata.mat")
# end

sigma_a = 1; sigma_s = 0.1; sigma_i = 0.2; 
# sigma_a_sbin = sigma_a  # remember we need this copy for Fmatrix
lam = -0.0005; B = 4.1; bias = 0.1; 
phi = 0.3; tau_phi = 0.1; lapse = 0.05;
params = [sigma_a, sigma_s, sigma_i, lam, B, bias, phi, tau_phi, lapse];
   

params = Dict("sigma_a" => sigma_a, "sigma_s" => sigma_s, "sigma_i" => sigma_i, 
              "lam" => lam, "B" => B, "bias" => bias, "phi" => phi, "tau_phi" => tau_phi, 
              "lapse" => lapse);
