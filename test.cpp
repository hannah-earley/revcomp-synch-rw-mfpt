#include <iostream>
#include <random>
#include <string>
#include <cstdint>
#include <cstdio>
#include <chrono>
#include "pcg/pcg_random.hpp"

int main(int argc, char *argv[]) {
    std::cout << "Hello, world!\n";

    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32 rng(seed_source);


    size_t n = 100;
    if (argc >= 2)
        n = std::stol(argv[1]);
    std::cout << n << " random nums\n";

    // low bit base64...
    // char b64[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    // for (int i=0; i<n; i++) {
    //     uint32_t n = 0;
    //     for (int i=0; i<6; i++) {
    //         n <<= 1;
    //         n += rng() & 1;
    //     }

    //     std::cout << b64[n];
    // }

    // raw throughput...
    // for (int i=0; i<n; i++) {
    //     uint32_t r = rng();
    //     char a,b,c,d;
    //     a = r;
    //     b = r >> 8;
    //     c = r >> 16;
    //     d = r >> 24;
    //     printf("%c%c%c%c", a, b, c, d);
    //     // std::cout << a << b << c << d;
    // }


    auto c_start = std::clock();
    uint32_t r = 0;
    for (int i=0; i<n; i++) {
        r += rng();
    }
    auto c_end = std::clock();

    auto t = (1.0 * (c_end - c_start)) / CLOCKS_PER_SEC;
    t += (r * 1.0) - (r * 1.0);
    // std::cout << r << "\n";
    std::cout << (1e9 * t / n) << " ns\n";

    // std::cout << "\n";
}