

ntrials =120

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

        push!(trial_vecs, collect((allocated_trial+1):(allocated_trial+my_n_trials)))
	allocated_trial += my_n_trials
end

push!(trial_vecs, collect(allocated_trial+1:ntrials))

display(trial_vecs)
LL = []


for i in 1:length(trial_vecs)
    push!(LL, @spawn multiLikey(ratdata, params, trial_vecs[i]))
end

for i in 1:length(trial_vecs)
    LL[i] = fetch(LL[i])
end

display(LL[1])

