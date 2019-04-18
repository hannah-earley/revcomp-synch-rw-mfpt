#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include "pcg/pcg_random.hpp"

void stats(std::vector<double> xs) {
    double mean_ = 0, var_ = 0, err_;
    size_t n = 0;

    for (double x : xs) {
        mean_ += x;
        var_ += x*x;
        n++;
    }

    mean_ /= n;
    var_ = var_ - n * mean_ * mean_;
    var_ /= n - 1;
    err_ = var_ / n;

    var_ = std::sqrt(var_);
    err_ = std::sqrt(err_);

    std::cout << "mean:" << mean_;
    std::cout << " var:" << var_;
    std::cout << " err:" << err_ << "\n";
}

void stats_recip(std::vector<double> xs) {
    double mean_ = 0, var_ = 0, err_;
    size_t n = 0;

    for (double x : xs) {
        mean_ += x;
        var_ += x*x;
        n++;
    }

    mean_ /= n;
    var_ = var_ - n * mean_ * mean_;
    var_ /= n - 1;
    err_ = var_ / n;

    var_ = std::sqrt(var_);
    err_ = std::sqrt(err_);

    mean_ = 1 / mean_;
    err_ *= mean_ * mean_;
    var_ = err_ * std::sqrt((double)n);

    std::cout << "mean:" << mean_;
    std::cout << " var:" << var_;
    std::cout << " err:" << err_ << "\n";
}

void drift_stats(std::vector<double> xs) {
    size_t n = xs.size();
    std::vector<double> ys(n-1, 0);
    for (size_t i = 0; i < n-1; i++)
        ys[i] = xs[i+1] - xs[i];
    std::cout << "drift.. ";
    stats(ys);
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
        stats_recip(js);
    }
}

int main(int argc, char *argv[]) {
    mfpt1d<int>(0.001, -5, 1000);
}