#include <chrono>

double LambertW0(double x, const double acc=1e-15);
unsigned stou(std::string const& str, size_t* idx=0, int base=10);

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

        std::cerr << "CPU:  " << (1e9 * t_cpu  / its) << " ns\n";
        std::cerr << "Wall: " << (1e9 * t_wall / its) << " ns\n";
    }
};