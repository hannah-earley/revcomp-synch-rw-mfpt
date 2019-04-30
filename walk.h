#include <chrono>
#include <signal.h>
#include <vector>

double LambertW0(double x, const double acc=1e-15);
unsigned stou(std::string const& str, size_t* idx=0, int base=10);
std::string format_si(double quantity, std::string unit);
std::string format_time(double seconds);
std::string format_timepoint(std::chrono::time_point<std::chrono::system_clock> tp,
                             const char* fmt);

struct SimTimer {
    std::clock_t start_cpu;
    std::chrono::high_resolution_clock::time_point start_wall;
    std::chrono::time_point<std::chrono::system_clock> start_real;

    SimTimer() {
        start_cpu = std::clock();
        start_wall = std::chrono::high_resolution_clock::now();
        start_real = std::chrono::system_clock::now();
    }

    void report(double its=1) {
        std::clock_t end_cpu = std::clock();
        std::chrono::high_resolution_clock::time_point end_wall
            = std::chrono::high_resolution_clock::now();
        std::chrono::time_point<std::chrono::system_clock> end_real
            = std::chrono::system_clock::now();

        double t_cpu = (1.0 * (end_cpu - start_cpu)) / CLOCKS_PER_SEC;
        std::chrono::duration<double> diff_wall = end_wall - start_wall;
        double t_wall = diff_wall.count();

        const char fmt[] = "%H:%M:%S %d/%m/%y";
        std::cerr << "CPU:  " << format_time(t_cpu/its) << " (" << format_time(t_cpu) << ")\n";
        std::cerr << "Wall: " << format_time(t_wall/its) << " (" << format_time(t_wall) << ")\n";
        std::cerr << "Real: " << format_timepoint(start_real, fmt) << " - "
                              << format_timepoint(end_real, fmt) << "\n";
    }
};

// temporarily install an ignoring signal handler for ^C interrupt
class Uninterruptable {
    struct sigaction old_act;

public:
    Uninterruptable() {
        struct sigaction new_act;
        new_act.sa_handler = SIG_IGN;
        sigemptyset(&new_act.sa_mask);
        new_act.sa_flags = 0;
        sigaction(SIGINT, &new_act, &old_act);
    }

    ~Uninterruptable() {
        sigaction(SIGINT, &old_act, NULL);
    }
};

class Unterminable {
    sigset_t oldmask, newmask;
    std::vector<int> signals;

public:
    Unterminable(std::vector<int> signals) : signals(signals) {
        sigemptyset(&newmask);
        for (int signal : signals)
            sigaddset(&newmask, signal);
        sigprocmask(SIG_BLOCK, &newmask, &oldmask);
    }

    Unterminable() : Unterminable({SIGINT, SIGTERM}) {}

    int poll() {
        sigset_t sigpend;
        sigpending(&sigpend);
        for (int signal : signals) {
            if (sigismember(&sigpend, signal)) {
                int sigret;
                if (sigwait(&newmask, &sigret) == 0)
                    return sigret;
                break;
            }
        }
        return false;
    }

    ~Unterminable() {
        sigprocmask(SIG_SETMASK, &oldmask, NULL);
    }
};