"""
Given a ratdata variable, as read from a Matfile, and a vector of trial_ids, returns a variable just like
the ratdata, but only for the indicated trials.
"""
function slice_data(ratdata, trial_ids)

    slice = copy(ratdata)
    for gloop = ["rawdata", "avgdata"]
        slice[gloop] = copy(ratdata[gloop])
        for k in keys(slice[gloop])
            slice[gloop][k] = ratdata[gloop][k][trial_ids]
        end
    end
    
    slice["total_trials"] = length(trial_ids)
    return slice
end



    
