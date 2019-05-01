# Todo

## `./batch.py`
- on sigint/sigterm, first try to pass onto current walk instance (then die)
- after multiple, sigkill it?

## `./walk`
- histogram/streaming mode for collecting many samples from ss distribution

# Done

## `./walk`
- safer checkpointing:
    - first, copy current checkpoint file to `./path/to/checkpoint.dat~`
    - then replace current checkpoint file...
- `-i` option: iteration limit
    - cease after this many ensemble measurements
    - `-i 0` places no limit
- `-q` option: clean output file
    - record fields: mean, error, #measurements
    - from these, we can reconstruct variance
    - can also improve estimates (always take mean in reciprocal...)
    - APPEND TO FILE
- progress meter
    - if iteration count exceeds some threshold...
    - place within outermost loop of ensemble_walk
    - will need a count for each thread
        - actually, should only report progress in master thread
        - otherwise won't be very useful - will tend to get a burst of
          progress ...... followed by nothing for a while
    - in master thread, print .s evenly up to some # of dots
    - when a thread completes its share, print a , (?)
- persistent data should end with a comment `# eof`
    - else assume error and die...
- hpc checkpointing
    - listen for sigterm and checkpoint early

## `./job.sh`
- pass through most options
- just extract certain values for filename purposes...
- and maybe automatically insert `-2v` ?
- pass sigint,sigterm signals onto ./walk properly

## `./refine.py`
- produce a refined measurement from the `-q` clean output
- `-s` skip option: skip first few measurements (burnin)
- computes a refined mean, error ~and variance~

## `./batch.py`
- build job batching utility
- jobs are stored in ./queue folder
    - can have subdirs for jobsets
    - each file describes a distinct job
    - job descriptor:
        - options to walk
        - target output count (-i adjusted to reach this goal)
        - cpu minimum (check against proc-cpuset)
            - alternatively, iteration time upper bound
        - target precision
            - supply a chunk size (-i) and a maximum output count
            - will keep running until either exceed output count or reach prec
- batch tool
    - verbs...
    - run
        - first, run `os.nice(40)`
        - optionally restrict to a jobset
        - place a fcntl lock on a job
            - fcntl locks should be nfs compatible
        - make a file ./queue/jobset/job-id.status
            - record hostname
            - periodically record psutil cpu_percent
            - job status info
    - status
        - loop over all jobs
            - if lock file held, processing
            - else maybe waiting, paused or done
                - no .status file => waiting
                - else confer with .status file
            - print relevant status file info
                - maybe .status should just be a single line of info to print...
    - status email
        - periodically send status update as email to designated email address
        - simply use sendmail
    - enqueue
        - generate jobs by some formula
        - don't overwrite existinc jobs
- job descriptor
    - options (job.sh)
        - string or array
    - target
        - _either_
            - count
        - _or_
            - chunk
            - limit
            - prec
    - requirements
        - cpu


# Var
## `./walk`
- More frequent checkpointing for long jobs...
- Progress monitor for long jobs...
- Explicit walker count?
    - If too few walkers, reseed (either randomly or from other walkers)
    - If too many, ignore excess (so that don't lose good burned data...)
- Target precision??
- Iteration count

Problems with more frequent checkpointing:
- it doesn't make a difference to current checkpointing system:
    - we assume the checkpoint is a random burned in state
    - therefore, recording this state data doesn't actually achieve anything..
    - furthermore, the parallelism makes recording a mid-ensemble state difficult

Alternative:
- Pick a desired sample window
- Then record as many ensembles as desired
- We can post-process to collate these measurements and obtain a more precise measurement...
- All we need to do is to discard the first few due to burn in
- We must also ensure sufficient precision