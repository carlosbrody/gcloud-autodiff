"""
vms_create_group(n_machines; ncores=2, groupname="julia-group")

"""
function vms_create_group(n_machines; ncores=2, groupname="julia-group")	 
    for i in 1:n_machines	
        instance_name = string(groupname, "-", i)
	machine_type  = string("n1-standard-", ncores)
	image         = "julia-and-packages-installed"

    	@async(run(`gcloud compute instances create $instance_name --machine-type $machine_type --image $image`))
    end    

    vl = vms_list(groupname=groupname, status="RUNNING")
    while length(vl) < n_machines
    	  sleep(1)
	  vl = vms_list(groupname=groupname, status="RUNNING")
    end
end

function vms_delete_group(groupname; wait_until_done=false)	 
    vl = vms_list(groupname=groupname)
    for i in 1:length(vl)	
        instance_name = vl[i]

    	@async(run(`gcloud compute instances delete --quiet $instance_name`))
    end    

    if wait_until_done
        vl = vms_list(groupname=groupname)
        while length(vl) > 0
    	  sleep(1)
	  vl = vms_list(groupname=groupname)
        end
    end
end


function vms_list(;groupname="julia-group", status="")
    a =[]
    try	 
        if isempty(status)
    	    a = split(readall(pipeline(`gcloud compute instances list`, `grep $groupname`, `awk '{print $1}'`)))
    	else
	    a = split(readall(pipeline(`gcloud compute instances list`, `grep $groupname`, `grep $status`, `awk '{print $1}'`)))
        end
    catch m
        a = []
    end

    return a
end


function vms_stop(vmlist)
    stopped_vms = []
    for i in 1:length(vmlist)
    	myvm =vmlist[i]
	try
            @async(run(`gcloud compute instances stop $myvm`))
	    push!(stopped_vms, myvm)
        catch
        end
    end
    return stopped_vms
end


function vms_start(vmlist)
    started_vms = []
    for i in 1:length(vmlist)
    	myvm =vmlist[i]
	try
            @async(run(`gcloud compute instances start $myvm`))
	    push!(started_vms, myvm)
        catch
        end
    end
    return started_vms
end
