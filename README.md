# Readme

We are seeking the mean first passage time, the MFPT, for different quadrant systems with varying bias, constriction width, and starting position. We have determined the MFPT analytically for continuous phase space, but this is more challenging for the discrete spaces of interest. The purpose of this program, then, is to compute numerical answers for the discrete case in order to compare with the continous analytical result in order to justify its use as an approximation.

An additional continous/discrete consideration is that of time. In fact, as the time between each transition is governed by iid Poisson/exponential statistics, the continuous and discrete timestep cases are equivalent, and the moments of the FPT can be trivially converted between these cases. Therefore, we can compute the moments for the discrete timestep case (i.e. 'computational time'), which is much more efficient, and then recover the desired continuous time results. It is more efficient because computing in continuous time requires discretising the time with a sufficiently small timestep such that more than one Poisson event per timestep is statistically negligible. Such a timestep size will be very small, much smaller than the unit timestep for the discrete case, and therefore necessitate far greater computational times (and the approximation may itself lead to numerical errors). Therefore the discrete case is far superior.

There are perhaps a couple 'obvious' approaches to determining the discrete MFPT:

- We could attempt to compute it directly and exactly via a Fokker-Planck type approach. In this method, we compute the full probability distribution at each time, thereby tracking every walker path instance efficiently. The MFPT follows by placing an absorbing boundary at the constriction, and tracking the probability lost over this boundary at each time. The MFPT is then
    $$\sum_{t=0}^\infty t \Delta p_t$$
Unfortuantely, doing so correctly requires tracking an infinite state space. One can truncate this state space, perhaps by placing an appropriate reflective boundary, but this will obviously effect the results somewhat. Placing the boundary sufficiently far from the constriction alleviates this to some degree, but the memory usage then becomes impractical for even reasonable problems. Furthermore, optimising for time efficiency (which is also poor) leads to convoluted code that is prone to bugs, and initial attempts were inconclusive.

- The next option is to appeal to the problem's probabilistic origins and employ a Monte-Carlo approach. The 'obvious' Monte-Carlo approach here is to simulate a large number of walkers and simply time how long it takes for them to reach the constriction. If the MFPT is finite, then every walker is guaranteed to terminate. Unfortunately, we are running at low bias and the geometry of the system tends to 'repel' walkers, such that it is likely that some walkers will get 'lost' for a prolonged period of time, and therefore the program runtime tends to diverge in practice. We could set a timeout, but this will give a systematic underestimate to the MFPT, and so is strongly undesirable.

Instead, we use a subtly different MCMC approach...

- An equivalent way to determine the MFPT from a set of walkers is to compute the rate/current of walkers crossing the constriction, where we maintain a fixed population of walkers by teleporting absorbed walkers back to the starting point (cite: Hill's algorithm). In principle this is practically the same as our original MC approach, except that it is robust to 'lost' walkers. Rather, the steady state dynamics of this walker population leads to the correct constriction-current, with the MFPT its reciprocal. The estimated MFPT can be improved by observing this (stochastic) current over a longer period of time. The primary issue is in reaching the steady state, and so there is an initial burn-in period wherein we eliminate the original correlation between the walkers by letting them diffuse for a sufficiently long period of time. Thus we see the resemblance to a Markov-Chain Monte-Carlo approach -- we have established a Markov Chain whose steady state distribution gives the MFPT current; therefore, by simulating this Markov Chain for a sufficiently long period of time, we can approximate this MFPT to arbitrary accuracy (after a burn-in time wherein the seed distribution is 'forgotten').

## Algorithm

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

## Log

- Stats for MFPT were giving systematic overestimate... this appears to be because we are taking average and variance of reciprocal of random variable (current); instead, we should compute stats for the rv itself (current mean, stdev, sterr) and then derive mfpt stats from these. This is because the distr for mfpt is a bit skewed, whereas the distr for current is more symmetric.

## Usage

- Run the tool with persistence (`-p filename`); that way, the current distribution is always stored on disk, and if the program is interrupted or crashes then we can resume from this save distribution. This means that we can burn the distribution in, and then always start from a burned in distribution, and we can also take more measurements whenever needed.
    - If we want to increase sample window or measurement ensemble, we can easily do so by adjusting the other parameters to ./walk
    - If we want to increase the number of walkers, we could manually edit the data file to replicate the distribution as needed, though this may require more burnin...