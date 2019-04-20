#include <chrono>
#include <signal.h>

double LambertW0(double x, const double acc=1e-15);
unsigned stou(std::string const& str, size_t* idx=0, int base=10);
std::string format_si(double quantity, std::string unit);
std::string format_time(double seconds);

struct SimTimer {
    std::clock_t start_cpu;
    std::chrono::high_resolution_clock::time_point start_wall;

    SimTimer() {
        start_cpu = std::clock();
        start_wall = std::chrono::high_resolution_clock::now();
    }

    void report(double its=1) {
        std::clock_t end_cpu = std::clock();
        std::chrono::high_resolution_clock::time_point end_wall
            = std::chrono::high_resolution_clock::now();

        double t_cpu = (1.0 * (end_cpu - start_cpu)) / CLOCKS_PER_SEC;
        std::chrono::duration<double> diff_wall = end_wall - start_wall;
        double t_wall = diff_wall.count();

        std::cerr << "CPU:  " << format_time(t_cpu/its) << " (" << format_time(t_cpu) << ")\n";
        std::cerr << "Wall: " << format_time(t_wall/its) << " (" << format_time(t_wall) << ")\n";
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