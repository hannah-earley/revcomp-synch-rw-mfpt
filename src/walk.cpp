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
#include "walk.hpp"

bool VERBOSE_MODE = false;
#define VERBOSE if (VERBOSE_MODE)

bool DYING = false;
#define DEATH_POLL if ((its_loc & (~(~0ULL << 24))) == 0)

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
    static void printn(FILE* file, const std::vector<S>&);
    static void printn(const std::vector<S>& vec) {
        printn(stdout, vec);
    }

};

template<>
int PersistentVector<int>::read1(std::string s) {
    return std::stoi(s);
}
template<>
std::ostream& PersistentVector<int>::write1(std::ostream& os, const int& n) {
    return os << n << std::endl;
}
template<>
void PersistentVector<int>::printn(FILE* file, const std::vector<int>& ns) {
    for (const int n : ns)
        fprintf(file, "%d\n", n);
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
                    const WalkConfig wc,
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
        size_t its_loc = 0;
        const size_t ensemble_count = wc.ensemble_count,
                     sample_window = wc.sample_window;

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
                for (size_t t = 0; t < sample_window; t++)
                    step_fn(walker, rng, r, step_env);
                rs_loc[j] += r;
                ns_loc[j]++;
                its_loc += sample_window;

                // death
                DEATH_POLL {
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
                    const WalkConfig wc,
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
        size_t its_loc = 0;
        const size_t ensemble_count = wc.ensemble_count,
                     sample_window = wc.sample_window;

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
                for (size_t t = 0; t < sample_window; t++)
                    step_fn(walker, rng, r, step_env);
                history.push_back(walker);
                its_loc += sample_window;

                // death
                DEATH_POLL {
                    if (thread == 0) POLLSIGS;
                    if (DYING) break;
                }
            }

            #pragma omp critical
            PersistentVector<S>::printn(history);

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

// convert real in [0,1] to int in [0,max]
uint32_t prob2int(double p) {
    const uint64_t u32max = 4294967295;
    uint32_t pt = 4294967295 * p;
    uint64_t pt_ = p * u32max;
    if (pt_ > u32max) pt = u32max;
    else if (pt_ < 0) pt = 0;
    else pt = pt_;
    return pt;
}

//// 1d walk

template<typename S>
void mfpt1d_seed(double bias, std::vector<S>& ensemble) {
    const double BIAS_THRESH = 0.0001;
    if (bias < BIAS_THRESH)
        bias = BIAS_THRESH;

    pcg32 rng = fresh_rng();
    std::geometric_distribution<S> seed_distr(2*bias/(1+bias));

    // 1D walks from an even distance conserve parity,
    // therefore the parity ratio (even:odd) established here
    // may give rise to separate populations which can affect
    // results (most prominent in histograms where it leads
    // to jaggedness).
    //
    // To prevent this, we guarantee equal parity here. This
    // is achieved by forcing even entries to be even and odd
    // odd

    int parity = 0;
    for (S& w : ensemble) {
        w = -seed_distr(rng);
        if ((parity - w) % 2 != 0)
            w--;
        parity = 1 - parity;
    }
}

template <typename S>
struct Params1D {
    S init, term;
    double bias, p;
    uint32_t pt;

    Params1D(S init, double bias) :
        init(init), term(0),
        bias(bias), p(0.5 * (1+bias)),
        pt(prob2int(p))
    {}
};

template <typename S>
void mfpt1d(double bias, S init, WalkConfig wc) {
    Params1D<S> params(init, bias);
    std::vector<S> ensemble(wc.n, params.init);
    mfpt1d_seed(bias, ensemble);

    walk_loop(ensemble, wc, [](S& w, pcg32& rng, size_t& r, Params1D<S>& params) {
        bool right = params.pt >= rng();
        w += 2 * right - 1;
        if (w >= params.term) {
            w = params.init;
            r++;
        }
    }, params);
}

//// 2d walk helpers

struct Coord2D {
    int x, y;
    Coord2D(int x, int y) : x(x), y(y) {}
    friend std::ostream& operator<< (std::ostream &os, const Coord2D &coord);

    inline void move_ql(bool right, bool down) {
        // move within a pre-constriction quadrant...

        /*  // below condition should never happen...
            if (y < 0 || y > -x) {
                std::cerr << "!" << *this << "!";
                if (y < 0) y = 0;
                if (y > -x) y = -x;
            }
        */

        int a = right ? y : -1;
        int b = x * down;
        if (a == -b) return;

        x += 2 * right - 1;
        y += down - right;

    }

    inline void move_rng(uint32_t pt, pcg32& g) {
        uint32_t r = g();
        bool right = pt >= r;
        bool down = (r & 1);
        move_ql(right, down);
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
template<>
void PersistentVector<Coord2D>::printn(FILE* file, const std::vector<Coord2D>& coords) {
    for (const Coord2D &coord : coords)
        fprintf(file, "%d,%d\n", coord.x, coord.y);
}

struct QStep {
    bool right:1, up:1;
    QStep(bool right, bool up) : right(right), up(up) {};
    friend std::ostream& operator<< (std::ostream &os, const QStep &step) {
        return os << (step.right ? "r" : "l") << (step.up ? "u" : "d");
    }
};
void quad_step_test() {
    QStep ld(false,false), lu(false,true), rd(true,false), ru(true,true);
    std::vector<QStep> dirs = {lu, ld, ru, rd};
    std::vector<Coord2D> tests = {{-4,0}, {-4,1}, {-4,2}, {-4,3}, {-4,4},
                                {-5,0}, {-6,6}, {-7,4}};

    for (Coord2D& c : tests) {
        std::cout << c << " : ";
        for (QStep& dir : dirs) {
            Coord2D cc = Coord2D(c);
            cc.move_ql(dir.right, !dir.up);
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
    const int init;
    const int term;
    const unsigned int col;
    std::uniform_int_distribution<int> init_dist;
    double bias, p;
    uint32_t pt;

    Params2D(int init, int term, double bias, unsigned int col = 0) :
        init(init), term(term), col(col),
        init_dist(0,-init),
        bias(bias), p(0.5 * (1+bias)),
        pt(prob2int(p))
    {}

    Params2D(int init, int term, double bias, bool cusp, unsigned int col = 0) :
        init(init), term(term), col(col),
        init_dist(0,-2*init),
        bias(bias), p(0.5 * (1+bias)),
        pt(prob2int(p))
    {}
};

void mfpt2d(double bias, unsigned int init_, unsigned int width, WalkConfig wc) {
    Params2D params(1 - init_ - width, 1 - width, bias);
    std::vector<Coord2D> ensemble(wc.n, Coord2D(params.init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    walk_loop(ensemble, wc, [](Coord2D& w, pcg32& rng, size_t& r, Params2D& params) {
        w.move_rng(params.pt, rng);
        if (w.x >= params.term) {
            w.x = params.init;
            w.y = params.init_dist(rng);
            r++;
        }
    }, params);
}

void mfpt2d_cusp(double bias, unsigned int init_, unsigned int width, WalkConfig wc) {
    Params2D params(1 - init_ - width, 1 - width, bias, true);
    std::vector<Coord2D> ensemble(wc.n, Coord2D(params.init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    walk_loop(ensemble, wc, [](Coord2D& w, pcg32& rng, size_t& r, Params2D& params) {
        w.move_rng(params.pt, rng);
        if (w.x >= params.term) {
            w.x = params.init;
            w.y = params.init_dist(rng);
            if (w.y < -w.x) {
                w.x -= w.y;
            } else {
                w.y += w.x;
                w.x -= w.y;
                w.y = -w.x-w.y;
            }
            r++;
        }
    }, params);
}

void mfpt2d_single(double bias, unsigned int init_, unsigned int width, unsigned int col, WalkConfig wc) {
    Params2D params(1 - init_ - width, 1 - width, bias, col);
    std::vector<Coord2D> ensemble(wc.n, Coord2D(params.init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    walk_loop(ensemble, wc, [](Coord2D& w, pcg32& rng, size_t& r, Params2D& params) {
        w.move_rng(params.pt, rng);
        if (w.x >= params.term) {
            w.x = params.init;
            w.y = params.col;
            r++;
        }
    }, params);
}

void mfpt2d_cusp_single(double bias, unsigned int init_, unsigned int width, unsigned int col, WalkConfig wc) {
    Coord2D w(1 - init_ - width, col);
    if (w.y < -w.x) {
        w.x -= w.y;
    } else {
        w.y += w.x;
        w.y = -w.x-w.y;
        w.x -= w.y;
        w.y = -w.x-w.y;
    }

    Params2D params(w.x, 1 - width, bias, true, w.y);
    std::vector<Coord2D> ensemble(wc.n, Coord2D(params.init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    walk_loop(ensemble, wc, [](Coord2D& w, pcg32& rng, size_t& r, Params2D& params) {
        w.move_rng(params.pt, rng);
        if (w.x >= params.term) {
            w.x = params.init;
            w.y = params.col;
            r++;
        }
    }, params);
}

void mfpt2d_refl(double bias, unsigned int init_, unsigned int width, WalkConfig wc) {
    Params2D params(2 - init_ - width, 1 - width, bias);
    std::vector<Coord2D> ensemble(wc.n, Coord2D(params.init,0));
    if (bias > 0) mfpt2d_seed(bias, width, ensemble);

    walk_loop(ensemble, wc, [](Coord2D& w, pcg32& rng, size_t& r, Params2D& params) {
        w.move_rng(params.pt, rng);
        if (w.x >= params.term) {
            w.x = params.init - 1;
            w.y = params.init_dist(rng) + (rng() % 2);
            r++;
        }
    }, params);
}

//// 2d gessel walk

struct Coord2DGessel {
    int x, y;
    Coord2DGessel(int x, int y) : x(x), y(y) {}
    friend std::ostream& operator<< (std::ostream &os, const Coord2DGessel &coord);

    inline void move_ql(bool up, bool right) {
        right = x == 0 ? 1 : right;
        y -= up - right;
        x += 2 * right - 1;
    }

    inline void move_rng(uint32_t pt, pcg32& g) {
        uint32_t r = g();
        bool up = pt >= r;
        bool right = (r & 1);
        move_ql(up, right);
    }
};
std::ostream& operator<< (std::ostream &os, const Coord2DGessel &coord) {
    return os << "(" << coord.x << ", " << coord.y << ")";
}
template<>
Coord2DGessel PersistentVector<Coord2DGessel>::read1(std::string s) {
    auto n = s.find(',');
    if (n == std::string::npos)
        throw std::invalid_argument("Corrupt persistent state... (invalid Coord2DGessel)");

    int x = std::stoi(s.substr(0,n));
    int y = std::stoi(s.substr(n+1));
    return Coord2DGessel(x, y);
}
template<>
std::ostream& PersistentVector<Coord2DGessel>::write1(std::ostream& os, const Coord2DGessel& coord) {
    return os << coord.x << "," << coord.y << std::endl;
}
template<>
void PersistentVector<Coord2DGessel>::printn(FILE* file, const std::vector<Coord2DGessel>& coords) {
    for (const Coord2DGessel &coord : coords)
        fprintf(file, "%d,%d\n", coord.x, coord.y);
}

void mfpt2d_gessel_seed(double bias, int init, int col, std::vector<Coord2DGessel>& ensemble) {
    for (Coord2DGessel& w : ensemble) {
        w = Coord2DGessel(col, init);
    }
}

struct Params2DGessel {
    const unsigned int init, col;
    double bias, p;
    uint32_t pt;

    Params2DGessel(unsigned int init, unsigned int col, double bias) :
        init(init), col(col),
        bias(bias), p(0.5 * (1+bias)),
        pt(prob2int(p))
    {}
};

void mfpt2d_gessel(double bias, unsigned int init, unsigned int col, WalkConfig wc) {
    Params2DGessel params(init, col, bias);
    std::vector<Coord2DGessel> ensemble(wc.n, Coord2DGessel(col, init));
    // if (bias > 0) mfpt2d_gessel_seed(bias, init, col, ensemble);

    walk_loop(ensemble, wc, [](Coord2DGessel& w, pcg32& rng, size_t& r, Params2DGessel& params) {
        w.move_rng(params.pt, rng);
        if (w.y <= 0) {
            w.y = params.init;
            w.x = params.col;
            r++;
        }
    }, params);
}

//// cli interface

enum Simulation { WALK_1D, WALK_2D, WALK_2G, UNITEST, TESTBED };
enum DistrShape { DIST_FLAT, DIST_CUSP, DIST_REFL };

std::ostream& operator<< (std::ostream &os, const Simulation &sim) {
    switch (sim) {
        case WALK_1D: return os << "MFPT - 1D Walk";
        case WALK_2D: return os << "MFPT - 2D Walk (Constrained/Quadrant)";
        case WALK_2G: return os << "MFPT - 2D Walk (Gessel)";
        case UNITEST: return os << "Unit Tests";
        case TESTBED: return os << "Test Bed...";
    }
    return os << "???";
}

std::ostream& operator<< (std::ostream &os, const DistrShape &shape) {
    switch (shape) {
        case DIST_FLAT: return os << "Flat";
        case DIST_CUSP: return os << "Cusp";
        case DIST_REFL: return os << "Reflective";
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
    fprintf(stderr, "    -g           Compute 2D gessel walk MFPT\n");
    fprintf(stderr, "    -t           Perform unit tests\n");
    fprintf(stderr, "    -h           Print this help message\n");
    fprintf(stderr, "    -v           Verbose/debug mode\n");
    fprintf(stderr, "    -b bias      Biased walk, bias \\in [-1,1]\n");
    fprintf(stderr, "    -d distance  Starting point, [nat]\n");
    fprintf(stderr, "    -c column    Single point distribution with given column, [nat]\n");
    fprintf(stderr, "    -C           Use a cusp-shape distribution\n");
    fprintf(stderr, "    -R           Use a 'reflective' (1/2,1,1,...,1,1,1/2) distribution\n");
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

    DistrShape shape = DIST_FLAT;
    double bias = 0;
    unsigned int width = 1;
    unsigned int dist = 1;
    int column = -1;
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
    while ((c = getopt(argc, argv, "12g9tvh?b:w:n:d:c:CRp:q:m:s:x:i:r")) != -1) switch(c) {
        case '1':
            sim = WALK_1D;
            break;
        case '2':
            sim = WALK_2D;
            break;
        case 'g':
            sim = WALK_2G;
            break;
        case 't':
            sim = UNITEST;
            VERBOSE_MODE = true;
            break;
        case '9':
            sim = TESTBED;
            VERBOSE_MODE = true;
            break;
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
        case 'c':
            column = stou(optarg);
            break;
        case 'C':
            shape = DIST_CUSP;
            break;
        case 'R':
            shape = DIST_REFL;
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
        std::cerr << "  Constriction Width: " << width << std::endl;
        std::cerr << "  Distribution Shape: " << shape << std::endl; }
    if (column >= 0) {
        std::cerr << "  Initial Column:     " << column << std::endl; }
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
            switch (shape) {
                case DIST_FLAT:
                    if (column < 0)      mfpt2d(bias, dist, width, wc);
                    else          mfpt2d_single(bias, dist, width, column, wc);
                    break;
                case DIST_CUSP:
                    if (column < 0) mfpt2d_cusp(bias, dist, width, wc);
                    else     mfpt2d_cusp_single(bias, dist, width, column, wc);
                    break;
                case DIST_REFL:
                    mfpt2d_refl(bias, dist, width, wc);
                    break;
            }
            break;

        case WALK_2G:
            mfpt2d_gessel(bias, dist, column, wc);
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
