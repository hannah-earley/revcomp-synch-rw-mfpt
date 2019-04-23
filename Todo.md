# Todo
## `./walk`
- More frequent checkpointing for long jobs...
- Progress monitor for long jobs...
- Explicit walker count?
    - If too few walkers, reseed (either randomly or from other walkers)
    - If too many, ignore excess (so that don't lose good burned data...)
- Target precision??
- Iteration count

## `./job.sh`
- pass through most options
- just extract certain values for filename purposes...
- and maybe automatically insert `-2v` ?