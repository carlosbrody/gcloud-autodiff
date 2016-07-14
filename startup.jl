import ForwardDiff
import Base.convert

rat_choice = "nothing";

include("qfind.jl")
include("BM.jl")
include("Fmatrix.jl")
include("logProbRight.jl")
include("logLike.jl")
include("llikey.jl")
include("multiLikey.jl")
include("run_trials.jl")
include("data.jl")


ps = [];
params =[];
LL = [];

epsilon = 10.0^(-10); dx = 0.125; dt = 0.02; 
global total_rate
total_rate = 40

using MAT
using ForwardDiff
