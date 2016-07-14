"""
LL = run_trials(ratdata, params, trial_ids=[])

Given ratdata (as read from a chronodata matfile) and a vector params (as defined in multiLikey), and 
a vector of trial numbers, spreads its calculation  over the local cores, and returns an array, 
same size as trial_ids, where each entry is LL, LLgrad, LLhess for the corresponding 
trial.

If trial_ids is not passed or is passed as an empty array, computes on all trials of ratdata.

"""
function run_trials(ratdata, params, trial_ids=[]; parallel=true)

    ntrials = length(trial_ids)
    if ntrials == 0
        trial_ids = collect(1:length(ratdata["rawdata"]["T"]))
        ntrials = length(trial_ids)
    end
    ncpus = nworkers()

    trial_vecs = []
    allocated_trial = 0;
    chunklen = Int(floor(ntrials/ncpus))
    remainder = rem(ntrials, ncpus)

    for i in 1:(ncpus-1)
      	my_n_trials = chunklen
	if remainder > 0
	    my_n_trials += 1
	    remainder    -= 1
	end

        push!(trial_vecs, trial_ids[(allocated_trial+1):(allocated_trial+my_n_trials)])
	allocated_trial += my_n_trials
    end

    push!(trial_vecs, trial_ids[allocated_trial+1:ntrials])

    display(trial_vecs)
    LLs = []
    for i in 1:length(trial_vecs)
        if ~parallel
            push!(LLs, multiLikey(ratdata, params, trial_vecs[i]));
        else
            @printf "spawning process from trial %d to trial %d\n" trial_vecs[i][1] trial_vecs[i][end]
            push!(LLs, @spawn multiLikey(ratdata, params, trial_vecs[i]));
        end
    end

    LL     = zeros(1, ntrials)
    LLgrad = zeros(length(params), ntrials)
    LLhess = zeros(length(params), length(params), ntrials)

    tidx = 1
    for i in 1:length(trial_vecs)
        LLi, LLgradi, LLhessi = fetch(LLs[i])
        n = length(LLi)
        LL[          tidx:tidx+n-1] = LLi
        LLgrad[:,    tidx:tidx+n-1] = LLgradi
        LLhess[:, :, tidx:tidx+n-1] = LLhessi
        tidx = tidx + n
    end

    # LL     = LLs[1][1];
    # LLgrad = LLs[1][2];
    # LLhess = LLs[1][3];

    # for i in collect(2:length(trial_vecs))
    # LL     += LLs[i][1];
    # LLgrad += LLs[i][2];
    # LLhess += LLs[i][3];
    # end

    return LL, LLgrad, LLhess
end


# display(LL[1])

