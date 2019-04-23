# Experiment Plan

## Width=1 constant, bias x distance
- bias: 0.1 0.01 0.001
- dist: 1 3 10 32 100 316 1000 3162 10000 31623 100000

## Plateau investigation
- bias: 0.01
    - extra data points in ~[100,1000]
    - 178 427 562 759 1778

## Advice
- don't necessarily maximise sample window and measurement ensemble
- rather, collect results more frequently to improve checkpointing and get more fine-grained data
- can still obtain just as high quality results by combining the results

## Computers
- beehive.maths (144 core)
- subliminal.maths (64 core)
- mea.damtp (24 core)  -- personal computer???
- maxwell.damtp (12 core)
- carme.damtp (8 core) -- me
