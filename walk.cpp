#include <iostream>
#include <random>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <cstdlib>
#include "pcg/pcg_random.hpp"
#include "walk.h"

////

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

////

template <typename S>
void mfpt1d(double bias, S init, size_t n) {
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    S term = 0;
    double p = 0.5 * (1 + bias);
    double q = 0.5 * (1 - bias);
    std::bernoulli_distribution fwd(p);

    std::vector<S> ensemble(n, init);
    std::geometric_distribution<S> seed_distr(bias/p);
    for (S& w : ensemble)
        w = -seed_distr(rng);

    size_t sample_window = 10000;
    size_t ensemble_count = 10000;
    std::vector<double> js(ensemble_count, 0);
    while (1) {
        // auto c_start = std::clock();
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
        // auto c_end = std::clock();
        // auto t = (1.0 * (c_end - c_start)) / CLOCKS_PER_SEC;
        // std::cout << (1e9 * t / (n * sample_window * ensemble_count)) << " ns\n";
        std::cout << Stats(js).pow(-1);
    }
}

////

struct QuadWalkStepDistr32 {
    struct Step {
        bool right;
        bool up;

        Step(bool right, bool up) : right(right), up(up) {};
        friend std::ostream& operator<< (std::ostream &os, const Step &step);
    };

    uint32_t pt;

    QuadWalkStepDistr32(double p) {
        const uint64_t u32max = 4294967295;
        pt = 4294967295 * p;
        uint64_t pt_ = p * u32max;
        if (pt_ > u32max) pt = u32max;
        else if (pt_ < 0) pt = 0;
        else pt = pt_;
    }

    inline Step operator() (pcg32& g) {
        uint32_t r = g();
        return Step(pt >= r, (r & 1) == 1);
    }
};
typedef QuadWalkStepDistr32::Step QStep;
std::ostream& operator<<(std::ostream& os, const QStep &step) {
    return os << (step.right ? "r" : "l")
              << (step.up    ? "u" : "d");
}

struct Coord2D {
    int x, y;
    Coord2D(int x, int y) : x(x), y(y) {}
    friend std::ostream& operator<< (std::ostream &os, const Coord2D &coord);

    inline void move_ql(QStep dw) {
        // move within a pre-constriction quadrant...
        if (false) // below condition should never happen...
            if (y < 0 || y > -x) {
                std::cout << "!" << *this << "!";
                if (y < 0) y = 0;
                if (y > -x) y = -x;
            }

        if (dw.right) {
            if (y == 0 && dw.up)
                return;
            if (y == -x && !dw.up)
                return;
        }

        x += dw.right ?  1 : -1;
        y += dw.up    ? -1 :  0;
        y += dw.right ?  0 :  1;
    }
};
std::ostream& operator<< (std::ostream &os, const Coord2D &coord) {
    return os << "(" << coord.x << ", " << coord.y << ")";
}

void quad_step_test() {
    QStep ld(false,false), lu(false,true), rd(true,false), ru(true,true);
    std::vector<QStep> dirs = {lu, ld, ru, rd};
    std::vector<Coord2D> tests = {{-4,0}, {-4,1}, {-4,2}, {-4,3}, {-4,4},
                                {-5,0}, {-6,6}, {-7,4}};

    for (Coord2D& c : tests) {
        std::cout << c << " : ";
        for (QStep& dir : dirs) {
            Coord2D cc = Coord2D(c);
            cc.move_ql(dir);
            std::cout << cc << " ";
        }
        std::cout << std::endl;
    }
}

////

void mfpt2d_seed(double bias, int width, pcg32& rng, std::vector<Coord2D>& ensemble) {
    std::uniform_real_distribution<double> unit_int(0,1);
    for (Coord2D& w : ensemble) {
        double x = unit_int(rng);
        int r = (width-1) + 0.5 * (1 + (1 + LambertW0((x - 1) / M_E)) / bias);
        std::uniform_int_distribution<int> col_dist(0,r);
        int c = col_dist(rng);
        w = Coord2D(-r, c);
        // std::cout << w << "\n";
    }
}

void mfpt2d(double bias, int init, uint width, size_t n) {
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);

    init = 1 - init - width;
    int term = 1 - width;
    std::uniform_int_distribution<int> init_dist(0,-init);

    double p = 0.5 * (1 + bias);
    double q = 0.5 * (1 - bias);
    QuadWalkStepDistr32 step(p);

    std::vector<Coord2D> ensemble(n, Coord2D(init,0));
    mfpt2d_seed(bias, width, rng, ensemble);

    size_t sample_window = 1000;
    size_t ensemble_count = 1000;
    std::vector<double> js(ensemble_count, 0);
    while (1) {
        auto c_start = std::clock();
        for (double& j : js) {
            j = 0;
            size_t r = 0;
            for (size_t t = 0; t < sample_window; t++) {
                for (Coord2D& w : ensemble) {
                    w.move_ql(step(rng));
                    if (w.x >= term) {
                        w.x = init;
                        w.y = init_dist(rng);
                        r++;
                    }
                }
            }
            j = (double)r / (double)n / (double)sample_window;
        }
        auto c_end = std::clock();
        auto t = (1.0 * (c_end - c_start)) / CLOCKS_PER_SEC;
        std::cout << (1e9 * t / (n * sample_window * ensemble_count)) << " ns\n";
        std::cout << Stats(js).pow(-1);
    }
}

////

int main(int argc, char *argv[]) {
    // quad_step_test();
    mfpt2d(0.5, 80, 1, 1000);
    // mfpt1d<int>(0.001, -1, 1000);
}