#include <iostream>
#include <random>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include "pcg/pcg_random.hpp"
#include "walk.h"

struct Stats {
    size_t n;
    double mean;
    double error;

    Stats(size_t n, double m, double e)
        : n(n), mean(m), error(e)
    {}

    Stats(std::vector<double> xs)
        : n(xs.size()), mean(0), error(0)
    {
        for (double x : xs) {
            mean += x;
            error += x*x;
        }

        mean /= n;
        error /= n;
        error -= mean * mean;
        error /= n - 1;
        error = std::sqrt(error);
    }

    Stats pow(double exp) const {
        double m, v, e;
        m = std::pow(mean, exp);
        e = std::abs(exp) * m * error / mean;
        return Stats(n, m, e);
    }

    double variance() const {
        return error * std::sqrt((double)n);
    }

    friend std::ostream& operator<< (std::ostream &os, const Stats &stats);
};

std::ostream& operator<<(std::ostream& os, const Stats &stats) {
    return os << "mean:" << stats.mean
              << " (±" << stats.error << ")"
              << " var:" << stats.variance()
              << std::endl;
}

template <typename S>
void mfpt1d(double bias, S init, size_t n) {
    S term = 0;
    double p = 0.5 * (1 + bias);
    double q = 0.5 * (1 - bias);

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    std::bernoulli_distribution fwd(p);
    std::vector<S> ensemble(n, init);

    std::geometric_distribution<S> seed_distr(bias/p);
    for (S& w : ensemble)
        w = -seed_distr(rng);

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
                    if (w >= term) {
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

struct Coord2D {
    int x, y;
    Coord2D(int x, int y) : x(x), y(y) {}
};

void mfpt2d_seed(double bias, int width, pcg32& rng, std::vector<Coord2D>& ensemble) {
    std::uniform_real_distribution<double> unit_int(0,1);
    for (Coord2D& w : ensemble) {
        double x = unit_int(rng);
        int r = (width-1) + 0.5 * (1 + (1 + LambertW0((x - 1) / M_E)) / bias);
        std::uniform_int_distribution<int> col_dist(0,r);
        int c = col_dist(rng);
        w = Coord2D(-r, c);
        std::cout << r << "," << c << "\n";
    }
}

void mfpt2d(double bias, int init, int width, size_t n) {
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    std::vector<Coord2D> ensemble(n, Coord2D(init,0));
    mfpt2d_seed(bias, width, rng, ensemble);
}

int main(int argc, char *argv[]) {
    mfpt2d(0.01, 0, 5, 25);
    // mfpt1d<int>(0.001, -1, 1000);
}