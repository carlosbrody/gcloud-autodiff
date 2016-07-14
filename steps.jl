# 1. Make sure VMs exist; if they don't, create them.
# 2. Make sure VMs are started; if they are not, start them.
# 3. Make sure each VM has the right data -- send ratdata for each machine, including trial id list even if machine will run all the trials it is given
# On each machine:
#     3.1  Start a Julia sentinel that initializes by loading all packages and all the ratdata
#     3.2  Make a Params directory. If empty, machine is working or has finished. Sentinel checks this directory periodically.
#     3.3  When sentinel picks up params, it empties the Params directory, runs, and then writes results (including trial id list) in Results, then goes back to sentineling
# 4. Send parameters to all machines; check Results directory in all machines, when results appear read and remove them.
# 5. Combine results, set next params, then send again.
# When done stop all VMs
