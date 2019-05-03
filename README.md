# Readme

## Usage

This repository comes with a number of tools.

*N.B. You should run these tools from the repository root directory...*

- `./walk` - this is the data generating workhorse. It performs 1D and (constrained) 2D random walks according to the parameters passed as command line arguments (use `./walk -h` for more information).
    - Data is reported periodically according to the options given and printed on STDOUT
    - In addition, the data can be output to a .csv file passed with the -q option
        - this .csv file lists the estimate of the _current_ rather than the mfpt; the mfpt is simply the reciprocal of the current
        - this is because the current is more suitable for aggregation
        - the third column gives the 'weight' of the result based on the number of step iterations used to generate the result; this should be used when aggregating the data
    - If a .dat file is passed with -p then the current walker distribution will be periodically checkpointed there; the tool will automatically resume from this .dat file (unless it is empty)
    - The tool will make use of all available cores (via OpenMP), and has been observed to do well on NUMA architectures (when -n,-m and -s are high enough). For example, on a 10-node NUMA system consisting of 360 CPUs, it was found that for reasonable parameters a speedup of ~360x was achievable
    - If you are interested instead in collecting the positions of the walkers in order to observe the steady state distribution (i.e. full MCMC mode), pass -r and the tool will stream this data to STDOUT
    - if sent a SIGINT or SIGTERM signal, the tool will aim to wrap up computation within the next few seconds and will emit a partial result
        - to stop it immediately, send it SIGKILL
        - this is useful for HPC architectures where job runners may pass SIGTERM shortly before SIGKILL
        - thus this tool should hopefully never waste any CPU time
        - `job.sh` and `batch.py` should handle these signals properly too, passing them on to `walk`

- `./job.sh` - this tool automatically sets up appropriate checkpointing, logging and csv output, inferring appropriate filenames from -b, -d and -w
    - slightly nicer interface to `./walk`
    - assumes a 2D walk
    - to run multiple instances of one job with the same -bdw parameters, use the -S option to append a distinct suffix to the filenames
        - for example, one could first burn in the distribution for a given -bdw parameter set, then copy the .dat file for each suffix
        - then multiple computers can 

- `./batch.py` - this tool allows for automatically running a large number of jobs.
    - job descriptions are .job files in JSON format in the `./queue/` directory. (Possibly nested) subdirectories can be used to create 'jobsets', and the batch tools can be restricted to these
        - job descriptions will crucially define a target which describes how long to run `./walk` for, either giving the total number of desired outputs or giving the target precision (possibly with bounds on the number of outputs)
        - example description:
            ```json
            {
                "options": ["-b", "0.01", "-d", "5", ...], // -2v implicit
                "target": {
                    "skip": 2, // number of 'burn-in' outputs (optional)

                    // either... total count
                    "count": 30 // number of real outputs to generate

                    // or...     precision
                    "chunk": 10 // number of outputs to generate at a time
                    "min": 20 // minimum number of outputs (optional)
                    "max": 150 // maximum number of outputs (optional)
                    "prec": 0.01 // target standard error
                },
                "requirements": {
                    "cpu": 4 // don't run on computers with fewer than 4 cpus (optional)
                }
            }
            ```
    - `./batch.py enqueue -h` - use `enq` to automatically generate the supplied jobset. The jobs are generated from `./queue/<jobset>/template.py`
        - `template.py` files will be dynamically imported by batch.py; they should define at least one function, `generate()`, which should yield a sequence of tuples `(name,job)` where `name` is the name of the job (leading to a job file `./queue/<jobset>/<name>.job`) and `job` is a dictionary which will be output as JSON into this file
            - a template might, for example, generate a series of similar job descriptors where the distance option increases logarithmically
            - it might even adjust the target precision for more realistic goals for higher distances, and increase the `-x` and `-n` parameters appropriately
        - `enq` can be rerun as needed to update job files (note that it will not delete job files that are no longer named)
            - if a job file is currently being processed by `./batch.py run` then `enq` will refuse to overwrite it; rerun `enq` later to rectify this
    - `./batch.py run -h` - `run` is used to run all jobs, or a given set of jobsets. The tool will iterate over each job and execute it via `./job.sh`
        - the number of iterations to execute is calculated and appended to the options provided to `./job.sh`
        - for dynamic targets (i.e. precision-based), `./job.sh` will be executed repeatedly with a given chunk size until the precision is reached (or the output bounds are reached)
        - whilst executing a job, the job file is locked (via fcntl)
            - this means that `run` can be called multiply, both on the same computer and across different computers (when over nfs or some other file system supporting remote atomic fcntl locks)
            - therefore multiple computers can coöperate in parallel on jobsets, though each job will be executed by at most one computer
        - if the job is interrupted for whatever reason, or if the target outputs later increase, then `run` will rerun each affected job when called again
    - `./batch.py status -h` - `stat` is used to provide status on all jobs or some jobsets
        - if used with `-n N` it will give give a status update every `N` seconds
        - if used with `-e [email]` it will send these status updates to the supplied email address instead of displaying them to STDOUT
        - status updates are derived from .status files generated by `run`, and therefore `stat` can display the status even from different computers
            - it will also check the fcntl lock statuses to check whether the relevant job is actually actively running, or whether `./batch.py run` has exited or died

- `./rawfine.py` - this tool can be used to approximately reconstruct .csv output from .log output
    - this is useful because early versions did not support csv output
    - given a set of input files, it will create a .csv~ file with the inferred output

- `./refine.py` - this tool takes a .csv file(s) and reports the combined data (with lower error) for each file by aggregating all the outputs
    - pass `-s S` to skip the first `S` datums (e.g. for distribution burn-in purposes)

- `./distribution.py` takes piped output from `./walk -r` and builds 'histogram' data of the distribution
    - the data is printed (in csv format) when the stream ends or on SIGINT
        - or after `-n` results
    - is provided with a filename, output is stored there
        - in addition, said csv file will be updated on subsequent runs
        - that is, the counts will be combined
    - arguments can be used to bin the results, skip some initial data, and give progress updates

### Tips
- The MFPT is computed as the reciprocal of the _steady state_ current for the distribution of walkers (where walkers reaching the absorbing boundary are teleported back to the starting position)
    - Therefore, `./walk` should first be run long enough for the initial distribution to converge on the steady state distribution
    - This is 'burn in', and it may not be obvious when this is achieved
        - In practice, the standard error for the first few results is higher than subsequent results
        - The subsequent results tend to be fairly consistent, particularly in terms of standard error
        - Some minimal effort is applied to obtain a (poor) initial approximation of this distribution, which hopefully substantially reduces the burn-in time
            - The main tactic is to spread initial walkers enough to reduce any correlation between these walkers, as one of the primary purposes of burn in is to decorrelate the positions of the walker ensemble
            - As such, skipping one or two initial results is often sufficient


## Background

We are seeking the mean first passage time, the MFPT, for different quadrant systems with varying bias, constriction width, and starting position. We have determined the MFPT analytically for continuous phase space, but this is more challenging for the discrete spaces of interest. The purpose of this program, then, is to compute numerical answers for the discrete case in order to compare with the continous analytical result in order to justify its use as an approximation.

An additional continous/discrete consideration is that of time. In fact, as the time between each transition is governed by iid Poisson/exponential statistics, the continuous and discrete timestep cases are equivalent, and the moments of the FPT can be trivially converted between these cases. Therefore, we can compute the moments for the discrete timestep case (i.e. 'computational time'), which is much more efficient, and then recover the desired continuous time results. It is more efficient because computing in continuous time requires discretising the time with a sufficiently small timestep such that more than one Poisson event per timestep is statistically negligible. Such a timestep size will be very small, much smaller than the unit timestep for the discrete case, and therefore necessitate far greater computational times (and the approximation may itself lead to numerical errors). Therefore the discrete case is far superior.

There are perhaps a couple 'obvious' approaches to determining the discrete MFPT:

- We could attempt to compute it directly and exactly via a Fokker-Planck type approach. In this method, we compute the full probability distribution at each time, thereby tracking every walker path instance efficiently. The MFPT follows by placing an absorbing boundary at the constriction, and tracking the probability lost over this boundary at each time. The MFPT is then
    $$\sum_{t=0}^\infty t \Delta p_t$$
Unfortuantely, doing so correctly requires tracking an infinite state space. One can truncate this state space, perhaps by placing an appropriate reflective boundary, but this will obviously effect the results somewhat. Placing the boundary sufficiently far from the constriction alleviates this to some degree, but the memory usage then becomes impractical for even reasonable problems. Furthermore, optimising for time efficiency (which is also poor) leads to convoluted code that is prone to bugs, and initial attempts were inconclusive.

- The next option is to appeal to the problem's probabilistic origins and employ a Monte-Carlo approach. The 'obvious' Monte-Carlo approach here is to simulate a large number of walkers and simply time how long it takes for them to reach the constriction. If the MFPT is finite, then every walker is guaranteed to terminate. Unfortunately, we are running at low bias and the geometry of the system tends to 'repel' walkers, such that it is likely that some walkers will get 'lost' for a prolonged period of time, and therefore the program runtime tends to diverge in practice. We could set a timeout, but this will give a systematic underestimate to the MFPT, and so is strongly undesirable.

Instead, we use a subtly different MCMC approach...

- An equivalent way to determine the MFPT from a set of walkers is to compute the rate/current of walkers crossing the constriction, where we maintain a fixed population of walkers by teleporting absorbed walkers back to the starting point (cite: Hill's algorithm). In principle this is practically the same as our original MC approach, except that it is robust to 'lost' walkers. Rather, the steady state dynamics of this walker population leads to the correct constriction-current, with the MFPT its reciprocal. The estimated MFPT can be improved by observing this (stochastic) current over a longer period of time. The primary issue is in reaching the steady state, and so there is an initial burn-in period wherein we eliminate the original correlation between the walkers by letting them diffuse for a sufficiently long period of time. Thus we see the resemblance to a Markov-Chain Monte-Carlo approach -- we have established a Markov Chain whose steady state distribution gives the MFPT current; therefore, by simulating this Markov Chain for a sufficiently long period of time, we can approximate this MFPT to arbitrary accuracy (after a burn-in time wherein the seed distribution is 'forgotten').

### Algorithm

The core algorithm is fairly simple - simulate the walkers, count the number of constriction incidences, and divide by time and number of walkers to get the current. To get multiple independent estimates, we can observe over distinct windows, and we can use these to get an error estimate.

The primary issue is to introduce some heuristic to establish when we have achieved burn-in. Initially, all the walkers are correlated and start from the beginning. There is also a minimum time before they can possibly reach the constriction. Therefore, over a short initial time period, the MFPT is infinite; after some threshold, the current will start to rise, and so the MFPT will start to fall. Emprically however, we also find that during the burn-in period the MFPT tends to be a significant _underestimate_, and this may be attributed to the fact that the walkers are initially correlated and so evolve somewhat in tandem.

One option to establish burn in is thus:
- Pick a sample ensemble size, say 100, and a sample window size, say 1000.
- To record a sample ensemble, we record each of the 100 samples in turn; to record a sample, we evolve all the walkers for 1000 steps and record the current within the window; we then determine the mean and std-err for the ensemble of samples, and print it.
- We record sample ensembles ad infinitum, printing each one.
- Observing the printout, we should see some initial inconsistent data that eventually yields to a more consistent stream with valid data.

There are some issues, however
- We may want some more guarantees about burn-in; one approach is to start in a less correlated seed state. One simple choice of seed state is to determine the equilibrium distribution assuming the constriction is reflecting (easily done via detailed balance) and to distribute the walkers so, at least initially. Thus the walkers have no intrinsic correlation, and the burn in now serves instead to colour the distribution with the starting and ending positions.
- For a quadrant walk, to where should we teleport absorbed walkers? If there's a single starting position this is obvious, but we will probably average over a transverse line. The best approach is probably to teleport to a random position along this line, maintaining proper mixing.
- How to pick good sample windows and ensemble sizes? We could do this by trial and error, and pay attention to the standard error, but perhaps we can also infer these results? I.e., we could present an adaptive algorithm which tweaks these values over time to achieve a desired standard error. In addition, we probably want to increase the sample window size over time as our confidence of burn-in increases.
