#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <unistd.h>
#include <omp.h>
#include <exception>
#include <cctype>
#include <cstring>
#include <string>
#include <signal.h>
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

pcg32 &omp_thread_rng(std::vector<pcg32>& rngs) {
    // ensures each thread has a distinct rng
    #pragma omp master
    {
        size_t thread_count = omp_get_num_threads();
        while (rngs.size() < thread_count)
            rngs.push_back(fresh_rng());
    }
    #pragma omp barrier
    return rngs[omp_get_thread_num()];
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
        pcg32& rng = omp_thread_rng(rngs);
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

    for (size_t i = 0; i < ensemble_count; i++)
        js[i] = (double)rs[i] / (double)n / (double) sample_window;
    return Stats(js).pow(-1);
}

//// walk miscellanea

struct WalkConfig {
    size_t n;
    std::string pv_filename;
    size_t sample_window;
    size_t ensemble_count;
    std::string cmd;

    double iterations() {
        return (double)n * (double)sample_window * (double)ensemble_count;
    }
};

//// persistence

template<typename S>
struct PersistentVector {
    WalkConfig wc;

    PersistentVector(WalkConfig wc) : wc(wc) {}

    void maybe_load(std::vector<S>& vec) {
        if (wc.pv_filename.empty())
            return;
        std::fstream fs;
        std::string line;
        fs.open(wc.pv_filename, std::ios_base::in);
        if (!fs.is_open() || fs.eof())
            return;

        vec.clear();
        while (std::getline(fs, line)) {
            if (line[0] == '#')
                continue;
            vec.push_back(read1(line));
        }
    }

    void store(const std::vector<S>& vec) {
        // prevent interrupts during storage
        struct sigaction new_act, old_act;
        new_act.sa_handler = SIG_IGN;
        sigemptyset(&new_act.sa_mask);
        new_act.sa_flags = 0;
        sigaction(SIGINT, &new_act, &old_act);
        // interrupts...

        if (wc.pv_filename.empty())
            return;
        std::fstream fs;
        fs.open(wc.pv_filename, std::ios_base::out|std::ios_base::trunc);

        fs << "# " << wc.cmd << std::endl;
        for (const S& x : vec)
            write1(fs, x);

        // uninstall interrupt handler
        sigaction(SIGINT, &old_act, NULL);
        // interrupts...
    }

    static S read1(std::string);
    static std::ostream& write1(std::ostream&, const S&);

};

template<>
int PersistentVector<int>::read1(std::string s) {
    return std::stoi(s);
}
template<>
std::ostream& PersistentVector<int>::write1(std::ostream& os, const int& n) {
    return os << n << std::endl;
}

//// 1d walk

template<typename S>
void mfpt1d_seed(double bias, std::vector<S>& ensemble) {
    pcg32 rng = fresh_rng();
    std::geometric_distribution<S> seed_distr(2*bias/(1+bias));
    for (S& w : ensemble)
        w = -seed_distr(rng);
}

template <typename S>
void mfpt1d(double bias, S init, WalkConfig wc) {
    S term = 0;
    double p = 0.5 * (1 + bias);
    std::bernoulli_distribution fwd(p);

    std::vector<S> ensemble(wc.n, init);
    if (bias > 0) mfpt1d_seed(bias, ensemble);
    std::vector<pcg32> rngs;

    PersistentVector<S> pv(wc);
    pv.maybe_load(ensemble);

    while (1) {
        SimTimer bench;
        std::cout << ensemble_walk(ensemble, wc.ensemble_count, wc.sample_window, rngs,
        [&](S& w, pcg32& rng, size_t& r) {
            w += fwd(rng) ? 1 : -1;
            if (w >= term) {
                w = init;
                r++;
            }
        });
        VERBOSE bench.report(wc.iterations());
        pv.store(ensemble);
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
                std::cerr << "!" << *this << "!";
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
template<>
Coord2D PersistentVector<Coord2D>::read1(std::string s) {
    auto n = s.find(',');
    if (n == std::string::npos)
        throw std::invalid_argument("Corrupt persistent state... (invalid Coord2D)");

    int x = std::stoi(s);
    int y = std::stoi(s.substr(n));
    return Coord2D(x, y);
}
template<>
std::ostream& PersistentVector<Coord2D>::write1(std::ostream& os, const Coord2D& coord) {
    return os << coord.x << "," << coord.y << std::endl;
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
        // use an approximate inverse cmf method to obtain the row...
        double x = unit_int(rng);
        int r = (width-1) + 0.5 * (1 + (1 + LambertW0((x - 1) / M_E)) / bias);
        // the columns are then equidistributed...
        std::uniform_int_distribution<int> col_dist(0,r);
        int c = col_dist(rng);
        w = Coord2D(-r, c);
    }
}

void mfpt2d(double bias, uint init, uint width, WalkConfig wc) {
    init = 1 - init - width;
    int term = 1 - width;
    std::uniform_int_distribution<int> init_dist(0,-init);

    double p = 0.5 * (1 + bias);
    QuadWalkStepDistr32 step(p);

    std::vector<Coord2D> ensemble(wc.n, Coord2D(init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);
    std::vector<pcg32> rngs;

    PersistentVector<Coord2D> pv(wc);
    pv.maybe_load(ensemble);

    while (1) {
        SimTimer bench;
        std::cout << ensemble_walk(ensemble, wc.ensemble_count, wc.sample_window, rngs,
        [&](Coord2D& w, pcg32& rng, size_t& r) {
            w.move_ql(step(rng));
            if (w.x >= term) {
                w.x = init;
                w.y = init_dist(rng);
                r++;
            }
        });
        VERBOSE bench.report(wc.iterations());
        pv.store(ensemble);
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

    fprintf(stderr, "Usage: %s\n", progn);
    fprintf(stderr, "         [-1|-2|-t] [-v] [-b bias] [-d distance]\n");
    fprintf(stderr, "         [-w width] [-n count] [-p filename]\n");
    fprintf(stderr, "         [-m count] [-s window] [-x count]\n");
    fprintf(stderr, "    -1           Compute 1D walk MFPT\n");
    fprintf(stderr, "    -2           Compute 2D walk MFPT\n");
    fprintf(stderr, "    -t           Perform unit tests\n");
    fprintf(stderr, "    -v           Verbose/debug mode\n");
    fprintf(stderr, "    -b bias      Biased walk, bias \\in [-1,1]\n");
    fprintf(stderr, "    -d distance  Starting point, [nat]\n");
    fprintf(stderr, "    -w width     Constriction width (2D only), [nat]\n");
    fprintf(stderr, "    -n count     Number of walkers in ensemble, [nat]\n");
    fprintf(stderr, "    -p filename  Persistence file\n");
    fprintf(stderr, "    -m count     Number of measurements, [nat]\n");
    fprintf(stderr, "    -s window    Sample time window, [nat]\n");
    fprintf(stderr, "    -x count     Sets both m and s, [nat]\n");
}

std::string args2cmd(int argc, char *argv[]) {
    // roughly reconstruct the original command line
    // good for a simple diagnostic header for persistence files
    std::string cmd;
    for (int i = 0; i < argc; i++) {
        char *arg = argv[i];

        bool quote = strchr(argv[i], ' ') != NULL;
        if (quote) cmd.push_back('"');

        for (char *arg = argv[i]; *arg != '\0'; arg++) {
            if (*arg == '\\' || *arg == '"')
                cmd.push_back('\\');
            cmd.push_back(std::isprint(*arg) ? *arg : '?');
        }

        if (quote) cmd.push_back('"');
        if (i < argc - 1)
            cmd.push_back(' ');
    }
    return cmd;
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
    Simulation sim = UNITEST;
    WalkConfig wc = {
        .n = 1000,
        .pv_filename = "",
        .sample_window = 1000,
        .ensemble_count = 1000,
        .cmd = args2cmd(argc, argv)
    };

    int c;
    while ((c = getopt(argc, argv, "12tvb:w:n:d:p:m:s:x:")) != -1) switch(c) {
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
            wc.n = std::stoul(optarg);
            break;
        case 'd':
            dist = stou(optarg);
            if (dist > 2147483647) { 
                std::cerr << "Distance outside range..\n\n";
                goto help;
            }
            break;
        case 'p':
            wc.pv_filename = optarg;
            break;
        case 'm':
            wc.ensemble_count = std::stoul(optarg);
            break;
        case 's':
            wc.sample_window = std::stoul(optarg);
            break;
        case 'x':
            wc.ensemble_count = std::stoul(optarg);
            wc.sample_window = std::stoul(optarg);
            break;
        case 'h':
        case '?':
        default:
            goto help;
    }

    VERBOSE {
        std::cerr << "Running simulation: " << sim << std::endl;
        std::cerr << "  Bias:               " << bias << std::endl;
        std::cerr << "  Distance:           " << dist << std::endl;
    if (sim == WALK_2D)
        std::cerr << "  Constriction Width: " << width << std::endl;
        std::cerr << "  Walker Count:       " << wc.n << std::endl;
        std::cerr << "  Measurement Count:  " << wc.ensemble_count << std::endl;
        std::cerr << "  Sample Window:      " << wc.sample_window << std::endl;
    if (!wc.pv_filename.empty())
        std::cerr << "  Persisting to:      " << wc.pv_filename << std::endl;
        std::cerr << std::endl;
    }

    switch (sim) {
        case WALK_1D:
            mfpt1d<int>(bias, -(int)dist, wc);
            break;

        case WALK_2D:
            mfpt2d(bias, dist, width, wc);
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