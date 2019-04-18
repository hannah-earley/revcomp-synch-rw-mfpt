#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include "pcg/pcg_random.hpp"

struct Stats {
    size_t n;
    double mean;
    double variance;
    double error;

    Stats(size_t n, double m, double v, double e)
        : n(n), mean(m), variance(v), error(e)
    {};

    Stats(std::vector<double> xs)
        : n(xs.size()), mean(0), variance(0), error(0)
    {
        for (double x : xs) {
            mean += x;
            variance += x*x;
        }

        mean /= n;
        variance -= n * mean * mean;
        variance /= n - 1;
        error = variance / n;

        variance = std::sqrt(variance);
        error = std::sqrt(error);
    }

    Stats pow(double exp) {
        double m, v, e;
        m = std::pow(mean, exp);
        e = std::abs(exp) * m * error / mean;
        v = e * std::sqrt((double)n);
        return Stats(n, m, v, e);
    }

    friend std::ostream& operator<< (std::ostream &os, const Stats &stats);
};

std::ostream& operator<<(std::ostream& os, const Stats &stats) {
    return os << "mean:" << stats.mean
              << " var:" << stats.variance
              << " err:" << stats.error
              << std::endl;
}

template <typename S>
void mfpt1d(double bias, S init, size_t n) {
    S term = 0;

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    std::bernoulli_distribution fwd(0.5 * (1 + bias));
    std::vector<S> ensemble(n, init);

    size_t sample_window = 10000;
    size_t ensemble_count = 10000;

    std::vector<double> js(ensemble_count, 0);
    while (1) {
        for (double& j : js) {
            j = 0;
            size_t r = 0;
            for (size_t t = 0; t < sample_window; t++) {
                for (S& w : ensemble) {
                    w += fwd(rng) ? 1 : -1;
                    if (w == term) {
                        w = init;
                        r++;
                    }
                }
            }
            j = (double)r / (double)n / (double)sample_window;
        }
        std::cout << Stats(js).pow(-1);
    }
}

int main(int argc, char *argv[]) {
    mfpt1d<int>(0.001, -5, 1000);
}