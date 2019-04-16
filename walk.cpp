#include <iostream>
#include <random>
#include <vector>
#include "pcg/pcg_random.hpp"

// template <
//     typename Walker,
//     std::size_t N
// > struct MFPT {
//     std::array<Walker, N> ensemble;
// };

template <typename S>
void mfpt1d(double bias, S init, size_t n) {
    S term = 0;

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    std::bernoulli_distribution fwd(0.5 * (1 + bias));

    std::vector<S> ensemble(n, init);
    size_t burnin = 0, burnt = 0;

    while (burnt++, burnin < n) {
        for (S& w : ensemble) {
            w += fwd(rng) ? 1 : -1;
            if (w == term) {
                w = init;
                // std::cout << "!";
                burnin++;
            }
            // std::cout << w << ",";
        }
        // std::cout << "\n";
    }
    std::cout << "BURNIN THRESHOLD\n";
    burnt*=100000;
    while (burnt--) {
        for (S& w : ensemble) {
            w += fwd(rng) ? 1 : -1;
            if (w == term) {
                w = init;
                // std::cout << "!";
            }
            // std::cout << w << ",";
        }
        // std::cout << "\n";
    }
    std::cout << "BURNIN COMPLETE\n";

    size_t t = 0, r = 0, t_max = 10000;
    while (++t) { // while (t++ < t_max) {
        for (S& w : ensemble) {
            w += fwd(rng) ? 1 : -1;
            if (w == term) {
                w = init;
                r++;
            }
        }

        if (t % 10000 == 0) {
            double j = (double)r / (double)n / (double)t;
            std::cout << t << " " << r << " : " << j << " " << (1/j) << "\n";
        }
    }

    /* TODO:

        - appropriate burn in time is perhaps hard to predict, and possibly far exceeds the time for most walkers to have completed a cycle

        - also, we would like to ascribe some sort of stdev or error to our results; to do so, we want to have a notion of samples of the mfpt estimate

        - we can achieve both of these goals by windowing over the current (j) calculation;

            - have an ensemble of walker ensembles
            - run the ensemble^2 for the window time, each walker ensemble then gives an mfpt estimate
            - combine the ensemble^2 mfpt estimate to get an error bound/stdev
            - by only computing j within the window, we will continue to burn in over time
            - can estimate burn in by comparing current estimate to previous estimate(s)
            - may also want to do some sort of adaptive window size...

    */
}

int main(int argc, char *argv[]) {
    std::cout << "Hello, world!\n";

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    std::cout << rng() << "\n";

    mfpt1d<int>(0.001, -5, 1000);
}