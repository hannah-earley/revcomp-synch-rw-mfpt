# Todo
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

### Plan
- progress meter
    - if iteration count exceeds some threshold...
    - place within outermost loop of ensemble_walk
    - will need a count for each thread
        - actually, should only report progress in master thread
        - otherwise won't be very useful - will tend to get a burst of
          progress ...... followed by nothing for a while
    - in master thread, print .s evenly up to some # of dots
    - when a thread completes its share, print a , (?)

# Done
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

## `./job.sh`
- pass through most options
- just extract certain values for filename purposes...
- and maybe automatically insert `-2v` ?

## `./refine.sh`
- produce a refined measurement from the `-q` clean output
- `-s` skip option: skip first few measurements (burnin)
- computes a refined mean, error ~and variance~