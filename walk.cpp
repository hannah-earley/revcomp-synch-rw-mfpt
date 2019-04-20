#include <iostream>
#include <random>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <unistd.h>
#include <omp.h>
#include "pcg/pcg_random.hpp"
#include "walk.h"

bool VERBOSE_MODE = false;
#define VERBOSE if (VERBOSE_MODE)

//// statistical analysis

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

//// parallel walking routine

pcg32 fresh_rng() {
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    return pcg32(seed_source);
}

template <typename S, typename F>
Stats ensemble_walk(std::vector<S>& walkers,
                    size_t ensemble_count,
                    size_t sample_window,
                    std::vector<pcg32>& rngs,
                    F step_fn)
{
    size_t n = walkers.size();
    std::vector<size_t> rs(ensemble_count, 0);
    std::vector<double> js(ensemble_count, 0);

    #pragma omp parallel
    {
        // get a fresh rng for each thread
        #pragma omp master
        {
            size_t thread_count = omp_get_num_threads();
            while (rngs.size() < thread_count)
                rngs.push_back(fresh_rng());
        }
        #pragma omp barrier
        pcg32& rng = rngs[omp_get_thread_num()];

        // do the actual ensemble of walks
        #pragma omp for
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < ensemble_count; j++) {
                size_t r = 0;
                for (size_t t = 0; t < sample_window; t++)
                    step_fn(walkers[i], rng, r);
                #pragma omp atomic
                rs[j] += r;
            }
        } 
    }

    // statistical summary
    for (size_t i = 0; i < ensemble_count; i++)
        js[i] = (double)rs[i] / (double)n / (double) sample_window;
    return Stats(js).pow(-1);
}

//// 1d walk

template <typename S>
void mfpt1d(double bias, S init, size_t n) {
    S term = 0;
    double p = 0.5 * (1 + bias);
    std::bernoulli_distribution fwd(p);

    std::vector<S> ensemble(n, init);
    if (bias > 0) {
        pcg32 rng = fresh_rng();
        std::geometric_distribution<S> seed_distr(bias/p);
        for (S& w : ensemble)
            w = -seed_distr(rng);
    }

    size_t sample_window = 10000;
    size_t ensemble_count = 10000;
    std::vector<pcg32> rngs;

    while (1) {
        SimTimer bench;
        std::cout << ensemble_walk(ensemble, ensemble_count, sample_window, rngs,
        [&](S& w, pcg32& rng, size_t& r) {
            w += fwd(rng) ? 1 : -1;
            if (w >= term) {
                w = init;
                r++;
            }
        });
        VERBOSE bench.report(n * sample_window * ensemble_count);
    }
}

//// 2d walk helpers

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

//// 2d walk

void mfpt2d_seed(double bias, int width, std::vector<Coord2D>& ensemble) {
    pcg32 rng = fresh_rng();
    std::uniform_real_distribution<double> unit_int(0,1);
    for (Coord2D& w : ensemble) {
        double x = unit_int(rng);
        int r = (width-1) + 0.5 * (1 + (1 + LambertW0((x - 1) / M_E)) / bias);
        std::uniform_int_distribution<int> col_dist(0,r);
        int c = col_dist(rng);
        w = Coord2D(-r, c);
        VERBOSE std::cout << w << ", ";
    }
    VERBOSE std::cout << std::endl;
}

void mfpt2d(double bias, uint init, uint width, size_t n) {
    init = 1 - init - width;
    int term = 1 - width;
    std::uniform_int_distribution<int> init_dist(0,-init);

    double p = 0.5 * (1 + bias);
    QuadWalkStepDistr32 step(p);

    std::vector<Coord2D> ensemble(n, Coord2D(init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    size_t sample_window = 1000;
    size_t ensemble_count = 1000;
    std::vector<pcg32> rngs;

    while (1) {
        SimTimer bench;
        std::cout << ensemble_walk(ensemble, ensemble_count, sample_window, rngs,
        [&](Coord2D& w, pcg32& rng, size_t& r) {
            w.move_ql(step(rng));
            if (w.x >= term) {
                w.x = init;
                w.y = init_dist(rng);
                r++;
            }
        });
        VERBOSE bench.report(n * sample_window * ensemble_count);
    }
}

//// cli interface

enum Simulation { WALK_1D, WALK_2D, UNITEST };
std::ostream& operator<< (std::ostream &os, const Simulation &sim) {
    switch (sim) {
        case WALK_1D: return os << "MFPT - 1D Walk";
        case WALK_2D: return os << "MFPT - 2D Walk (Constrained/Quadrant)";
        case UNITEST: return os << "Unit Tests";
    }
}

void help(int argc, char *argv[]) {
    char progn_def[] = "./walk";
    char *progn = progn_def;
    if (argc > 0) progn = argv[0];

    fprintf(stderr, "Usage: %s [-1|-2|-t] [-v] [-b bias] [-d distance] [-w width] [-n count]\n", progn);
    fprintf(stderr, "    -1           Compute 1D walk MFPT\n");
    fprintf(stderr, "    -2           Compute 2D walk MFPT\n");
    fprintf(stderr, "    -t           Perform unit tests\n");
    fprintf(stderr, "    -v           Verbose/debug mode\n");
    fprintf(stderr, "    -b bias      Biased walk, bias \\in [-1,1]\n");
    fprintf(stderr, "    -d distance  Starting point, [nat]\n");
    fprintf(stderr, "    -w width     Constriction width (2D only), [nat]\n");
    fprintf(stderr, "    -n count     Number of walkers in ensemble, [nat]\n");
}

int unit_tests() {
    std::cout << "# Quadrant Step Test\n";
    quad_step_test();
    std::cout << "\n";

    return 0;
}

int main(int argc, char *argv[]) {
    double bias = 0;
    uint width = 1;
    uint dist = 1;
    size_t n = 1000;
    Simulation sim = UNITEST;

    int c;
    while ((c = getopt(argc, argv, "12tvb:w:n:d:")) != -1) switch(c) {
        case '1':
            sim = WALK_1D;
            break;
        case '2':
            sim = WALK_2D;
            break;
        case 't':
            sim = UNITEST;
        case 'v':
            VERBOSE_MODE = true;
            break;
        case 'b':
            bias = std::stod(optarg);
            if (bias < -1 || bias > 1) {
                std::cerr << "Bias outside range..\n\n";
                goto help;
            }
            break;
        case 'w':
            width = stou(optarg);
            break;
        case 'n':
            n = std::stoul(optarg);
            break;
        case 'd':
            dist = stou(optarg);
            if (dist > 2147483647) { 
                std::cerr << "Distance outside range..\n\n";
                goto help;
            }
            break;
        case 'h':
        case '?':
        default:
            goto help;
    }

    std::cout << "Running simulation: " << sim << std::endl;
    std::cout << "  Bias:               " << bias << std::endl;
    std::cout << "  Distance:           " << dist << std::endl;
  if (sim == WALK_2D)
    std::cout << "  Constriction Width: " << width << std::endl;
    std::cout << "  Ensemble Count:     " << n << std::endl;
  VERBOSE
    std::cout << "Verbose output..." << std::endl;
    std::cout << std::endl;

    switch (sim) {
        case WALK_1D:
            mfpt1d<int>(bias, -(int)dist, n);
            break;

        case WALK_2D:
            mfpt2d(bias, dist, width, n);
            break;

        case UNITEST:
            unit_tests();
            break;
    }

    return 0;
help:
    help(argc, argv);
    return 1;
}