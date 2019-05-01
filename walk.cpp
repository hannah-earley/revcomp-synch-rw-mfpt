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
#include <limits>
#include "pcg/pcg_random.hpp"
#include "walk.h"

bool VERBOSE_MODE = false;
#define VERBOSE if (VERBOSE_MODE)

bool DYING = false;
size_t DEATH_CHECK = 100000000;

//// walk config

struct WalkConfig {
    size_t n;
    std::string pv_filename;
    std::string out_filename;
    size_t sample_window;
    size_t ensemble_count;
    size_t loop_count;
    std::string cmd;
    bool streaming;

    double iterations() {
        return (double)n * (double)sample_window * (double)ensemble_count;
    }
};

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
        double m, e;
        m = std::pow(mean, exp);
        e = std::abs(exp) * m * error / mean;
        return Stats(n, m, e);
    }

    double variance() const {
        return error * std::sqrt((double)n);
    }

    friend std::ostream& operator<< (std::ostream &os, const Stats &stats);

    template<typename T>
    void store(std::string filename, T weight) {
        if (filename.empty())
            return;

        std::fstream fs;
        fs.open(filename, std::ios_base::out|std::ios_base::app);
        if (fs.tellp() <= 0) // write header
            fs << "#mean,stderr,weight" << std::endl;

        fs.precision(16);
        fs << mean << "," << error << "," << weight << std::endl;
    }

    void store(std::string filename) {
        store(filename, n);
    }
};

std::ostream& operator<<(std::ostream& os, const Stats &stats) {
    return os << "mean:" << stats.mean
              << " (±" << stats.error << ")"
              << " var:" << stats.variance()
              << std::endl;
}

//// persistence

const char PV_EOF[] = "# eof";

template<typename S>
struct PersistentVector {
    std::string filename_bak;
    WalkConfig wc;

    PersistentVector(WalkConfig wc) : wc(wc) {
        filename_bak = wc.pv_filename + "~";
    }

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
            if (line[0] == '#') {
                if (line == PV_EOF)
                    return;
                continue;
            }
            vec.push_back(read1(line));
        }

        std::cerr << "Corrupt persistent data file!" << std::endl;
        std::exit(1);
    }

    void store(const std::vector<S>& vec) {
        Uninterruptable unint;

        if (wc.pv_filename.empty())
            return;
        backup();

        std::fstream fs;
        fs.open(wc.pv_filename, std::ios_base::out|std::ios_base::trunc);

        fs << "# " << wc.cmd << std::endl;
        for (const S& x : vec)
            write1(fs, x);
        fs << PV_EOF << std::endl;
    }

    void backup() {
        std::fstream orig, bak;
        std::string line;

        orig.open(wc.pv_filename, std::ios_base::in);
        if (!orig.is_open() || orig.eof())
            return;

        bak.open(filename_bak, std::ios_base::out|std::ios_base::trunc);
        while (std::getline(orig, line))
            bak << line << std::endl;
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

#define POLLSIGS if ((sig = unterm.poll()) > 0) {\
    std::cerr << "...got interrupt(" << sig << ")...";\
    DYING = true;\
}

#define if_POLLSIGS POLLSIGS if (sig > 0)

template <typename S, typename E, typename F>
Stats ensemble_walk(std::vector<S>& walkers,
                    WalkConfig wc,
                    std::vector<pcg32>& rngs,
                    size_t& its,
                    F step_fn,
                    E step_env)
{
    // catch sigint, sigterm...
    Unterminable unterm;
    int sig;

    its = 0;
    size_t n = walkers.size();
    std::vector<size_t> rs(wc.ensemble_count, 0);
    std::vector<double> js(wc.ensemble_count, 0);
    std::vector<size_t> ns(wc.ensemble_count, 0);
    size_t threads_done = 0;

    #pragma omp parallel
    {
        pcg32& rng = omp_thread_rng(rngs);
        std::vector<size_t> rs_loc(wc.ensemble_count, 0);
        std::vector<size_t> ns_loc(wc.ensemble_count, 0);
        size_t its_loc = 0,
               ensemble_count = wc.ensemble_count,
               sample_window = wc.sample_window,
               death_check = DEATH_CHECK;

        size_t num_threads = omp_get_num_threads();
        size_t threads_done_loc;
        int thread = omp_get_thread_num();

        // progress
        size_t progress_its = 0;
        size_t progress_10pc = 0;
        size_t progress_chunk = (n / num_threads) + 1;

        #pragma omp for nowait
        for (size_t i = 0; i < n; i++) {
            if (DYING) continue;

            // progress (only from master)
            VERBOSE if (thread == 0) {
                if (progress_its == 0) {
                    std::cerr << "0%.." << std::flush;
                }
                progress_its++;
                size_t n10pc = (progress_its * 10) / progress_chunk;
                while (progress_10pc < n10pc && progress_10pc < 9) {
                    progress_10pc++;
                    std::cerr << progress_10pc << "0%.." << std::flush;
                }
            }

            S walker = walkers[i];
            for (size_t j = 0; j < ensemble_count; j++) {
                size_t r = 0;
                for (size_t t = 0; t < sample_window; t++, its_loc++)
                    step_fn(walker, rng, r, step_env);
                rs_loc[j] += r;
                ns_loc[j]++;

                // death
                if (its_loc % death_check == 0) {
                    if (thread == 0) POLLSIGS;
                    if (DYING) break;
                }
            }
            walkers[i] = walker;
        }

        #pragma omp atomic
        threads_done++;

        #pragma omp master
        {
            if (!DYING) {
                #pragma omp critical
                threads_done_loc = threads_done;

                while (threads_done_loc < num_threads) {
                    if_POLLSIGS break;
                    sleep(1);
                    #pragma omp critical
                    threads_done_loc = threads_done;
                }
            }
        }

        for (size_t j = 0; j < ensemble_count; j++) {
            #pragma omp atomic
            rs[j] += rs_loc[j];
            #pragma omp atomic
            ns[j] += ns_loc[j];
        }

        #pragma omp atomic
        its += its_loc;

        // progress (only from master)
        VERBOSE if (thread == 0 && !DYING) {
            std::cerr << "100%" << std::flush;
        }
    }

    // death
    POLLSIGS;

    // progress
    VERBOSE std::cerr << std::endl;

    for (size_t i = 0; i < wc.ensemble_count; i++)
        js[i] = (double)rs[i] / (double)ns[i] / (double) wc.sample_window;

    Stats j(js);
    j.store(wc.out_filename, its);
    std::cerr << "its:  " << its << std::endl;
    return j.pow(-1);
}

template <typename S, typename E, typename F>
void streaming_walk(std::vector<S>& walkers,
                    WalkConfig wc,
                    std::vector<pcg32>& rngs,
                    F step_fn,
                    E step_env)
{
    Unterminable unterm;
    int sig;
    size_t n = walkers.size();
    size_t threads_done = 0;

    #pragma omp parallel
    {
        pcg32& rng = omp_thread_rng(rngs);
        size_t its_loc = 0,
               ensemble_count = wc.ensemble_count,
               sample_window = wc.sample_window,
               death_check = DEATH_CHECK;

        size_t num_threads = omp_get_num_threads();
        size_t threads_done_loc;
        int thread = omp_get_thread_num();

        std::vector<S> history;
        history.reserve(ensemble_count);

        #pragma omp for nowait
        for (size_t i = 0; i < n; i++) {
            if (DYING) continue;

            S walker = walkers[i];

            for (size_t j = 0; j < ensemble_count; j++) {
                size_t r = 0;
                for (size_t t = 0; t < sample_window; t++, its_loc++)
                    step_fn(walker, rng, r, step_env);
                history.push_back(walker);

                // death
                if (its_loc % death_check == 0) {
                    if (thread == 0) POLLSIGS;
                    if (DYING) break;
                }
            }

            #pragma omp critical
            {
                for (const S& x : history)
                    PersistentVector<S>::write1(std::cout, x);
            }

            history.clear();
            walkers[i] = walker;
        }

        #pragma omp atomic
        threads_done++;

        #pragma omp master
        {
            if (!DYING) {
                #pragma omp critical
                threads_done_loc = threads_done;

                while (threads_done_loc < num_threads) {
                    if_POLLSIGS break;
                    sleep(1);
                    #pragma omp critical
                    threads_done_loc = threads_done;
                }
            }
        }
    }

    POLLSIGS;
}

template <typename S, typename E, typename F>
void walk_loop(std::vector<S>& ensemble,
               WalkConfig wc,
               F step_fn,
               E step_env) {

    size_t its;
    std::vector<pcg32> rngs;
    PersistentVector<S> pv(wc);
    pv.maybe_load(ensemble);

    for (size_t loopi = 0; !DYING && loopi < wc.loop_count; loopi++) {
        if (wc.streaming) {
            streaming_walk(ensemble, wc, rngs, step_fn, step_env);
        } else {
            SimTimer bench;
            std::cout << ensemble_walk(ensemble, wc, rngs, its, step_fn, step_env);
            VERBOSE bench.report(its);
        }
        
        pv.store(ensemble);
    }
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
struct Params1D {
    S init, term;
    double bias, p;
    std::bernoulli_distribution fwd;

    Params1D(S init, double bias) :
        init(init), term(0),
        bias(bias), p(0.5 * (1+bias)),
        fwd(p)
    {}
};

template <typename S>
void mfpt1d(double bias, S init, WalkConfig wc) {
    Params1D<S> params(init, bias);
    std::vector<S> ensemble(wc.n, params.init);
    if (bias > 0) mfpt1d_seed(bias, ensemble);

    walk_loop(ensemble, wc, [](S& w, pcg32& rng, size_t& r, Params1D<S>& params) {
        w += params.fwd(rng) ? 1 : -1;
        if (w >= params.term) {
            w = params.init;
            r++;
        }
    }, params);
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

    inline Step operator() (pcg32& g) const {
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

    int x = std::stoi(s.substr(0,n));
    int y = std::stoi(s.substr(n+1));
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

struct Params2D {
    const uint init;
    const int term;
    std::uniform_int_distribution<int> init_dist;
    double bias, p;
    QuadWalkStepDistr32 step;

    Params2D(uint init, int term, double bias) :
        init(init), term(term),
        init_dist(0,-init),
        bias(bias), p(0.5 * (1+bias)),
        step(p)
    {}
};

void mfpt2d(double bias, uint init_, uint width, WalkConfig wc) {
    Params2D params(1 - init_ - width, 1 - width, bias);
    std::vector<Coord2D> ensemble(wc.n, Coord2D(params.init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    walk_loop(ensemble, wc, [](Coord2D& w, pcg32& rng, size_t& r, Params2D& params) {
        w.move_ql(params.step(rng));
        if (w.x >= params.term) {
            w.x = params.init;
            w.y = params.init_dist(rng);
            r++;
        }
    }, params);
}

//// cli interface

enum Simulation { WALK_1D, WALK_2D, UNITEST, TESTBED };
std::ostream& operator<< (std::ostream &os, const Simulation &sim) {
    switch (sim) {
        case WALK_1D: return os << "MFPT - 1D Walk";
        case WALK_2D: return os << "MFPT - 2D Walk (Constrained/Quadrant)";
        case UNITEST: return os << "Unit Tests";
        case TESTBED: return os << "Test Bed...";
    }
    return os << "???";
}

void help(int argc, char *argv[]) {
    char progn_def[] = "./walk";
    char *progn = progn_def;
    if (argc > 0) progn = argv[0];

    fprintf(stderr, "Usage: %s [options]\n", progn);
    fprintf(stderr, "    -1           Compute 1D walk MFPT\n");
    fprintf(stderr, "    -2           Compute 2D walk MFPT\n");
    fprintf(stderr, "    -t           Perform unit tests\n");
    fprintf(stderr, "    -h           Print this help message\n");
    fprintf(stderr, "    -v           Verbose/debug mode\n");
    fprintf(stderr, "    -b bias      Biased walk, bias \\in [-1,1]\n");
    fprintf(stderr, "    -d distance  Starting point, [nat]\n");
    fprintf(stderr, "    -w width     Constriction width (2D only), [nat]\n");
    fprintf(stderr, "    -n count     Number of walkers in ensemble, [nat]\n");
    fprintf(stderr, "    -p filename  Persistence file\n");
    fprintf(stderr, "    -q filename  Clean output file\n");
    fprintf(stderr, "    -m count     Number of measurements, [nat]\n");
    fprintf(stderr, "    -s window    Sample time window, [nat]\n");
    fprintf(stderr, "    -x count     Sets both m and s, [nat]\n");
    fprintf(stderr, "    -i count     Number of outputs (0=unlimited), [nat]\n");
    fprintf(stderr, "    -r           Streaming mode (e.g. for histograms)\n");
}

std::string args2cmd(int argc, char *argv[]) {
    // roughly reconstruct the original command line
    // good for a simple diagnostic header for persistence files
    std::string cmd;
    for (int i = 0; i < argc; i++) {
        char *arg = argv[i];

        bool quote = strchr(arg, ' ') != NULL;
        if (quote) cmd.push_back('"');

        for (; *arg != '\0'; arg++) {
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

#include "testbed.cpp"

int main(int argc, char *argv[]) {
    std::cout.precision(12);
    std::cerr.precision(12);

    double bias = 0;
    uint width = 1;
    uint dist = 1;
    Simulation sim = UNITEST;
    WalkConfig wc = {
        .n = 1000,
        .pv_filename = "",
        .out_filename = "",
        .sample_window = 1000,
        .ensemble_count = 1000,
        .loop_count = 0,
        .cmd = args2cmd(argc, argv),
        .streaming = false
    };

    int c;
    while ((c = getopt(argc, argv, "129tvh?b:w:n:d:p:q:m:s:x:i:r")) != -1) switch(c) {
        case '1':
            sim = WALK_1D;
            break;
        case '2':
            sim = WALK_2D;
            break;
        case 't':
            sim = UNITEST;
        case '9':
            sim = TESTBED;
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
        case 'q':
            wc.out_filename = optarg;
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
        case 'i':
            wc.loop_count = std::stoul(optarg);
            break;
        case 'r':
            wc.streaming = true;
            break;
        case 'h':
        case '?':
        default:
            goto help;
    }

    VERBOSE {
        std::cerr << "Running simulation: " << sim << std::endl;
    if (wc.streaming) {
        std::cerr << "Streaming mode active" << std::endl; }
        std::cerr << "  Bias:               " << bias << std::endl;
        std::cerr << "  Distance:           " << dist << std::endl;
    if (sim == WALK_2D) {
        std::cerr << "  Constriction Width: " << width << std::endl; }
        std::cerr << "  Walker Count:       " << wc.n << std::endl;
        std::cerr << "  Measurement Count:  " << wc.ensemble_count << std::endl;
        std::cerr << "  Sample Window:      " << wc.sample_window << std::endl;
    if (wc.loop_count > 0) {
        std::cerr << "  Output Count:       " << wc.loop_count << std::endl; }
    if (!wc.pv_filename.empty()) {
        std::cerr << "  Persisting to:      " << wc.pv_filename << std::endl; }
        std::cerr << std::endl;
    }

    if (wc.loop_count == 0)
        wc.loop_count = std::numeric_limits<size_t>::max();

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

        case TESTBED:
            testbed();
            break;
    }

    return 0;
help:
    help(argc, argv);
    return 1;
}
