#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include "walk.h"

double LambertW0_neg_(double x) {
    // Estimate LambertW (0 branch) using two approaches (5 terms)
    // for non-positive inputs (within the range [-1/e, 0]).
    // Maximum error is ~0.005.
    // Use LambertW0_neg for better accuracy (Newton-Raphson)
    // We also give a (poor) approximation for positive x,
    // to act as a seed for the Newton-Raphson method.

    const double x0 = -0.24522772368742732;
    double r;

    if (x < x0) {
        // https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf - (4.22-4.24)
        double p2 = 1 + M_E * x;
        if (p2 < 0) p2 = 0; // clamp to range
        double p = std::sqrt(p2);
        double pp = p;

        r = -1;
        r += pp;
        pp *= p;
        r += (-1/3) * pp;
        pp *= p;
        r += (11/72) * pp;
        pp *= p;
        r += (-43/540) * pp;
        pp *= p;
        r += (769/17280) * pp;
    } else if (x < 0.5) {
        // series expansion about x=0
        double xx = x;
        r = xx;
        xx *= x;
        r += -xx;
        xx *= x;
        r += (3/2)*xx;
        xx *= x;
        r += (-8/3)*xx;
        xx *= x;
        r += (125/24)*xx;

        if (r > 0) r = 0;
    } else if (x < 10) {
        r = 0.7 * std::log(1 + x);
    } else {
        double l = std::log(x);
        double ll = std::log(l);
        double lll = std::log(1 - ll / l);
        r = l - ll - lll;
    }

    return r;
}

double LambertW0(double x, const double acc) {
    size_t n = 0;
    double w0, w1 = LambertW0_neg_(x);
    // std::cout << "w(" << x << ")~" << w1 << "\n";
    do {
        w0 = w1;
        double ew = std::exp(w0);
        w1 = w0 - (w0*ew - x) / (ew * (1 + w0));
        if (std::isnan(w1)) {
            if (-0.37 < x && x < 0) {
                w1 = -1;
                break;
            }
        }
    } while(n++ < 10 && std::abs((w1-w0)/w1) > acc);
    // std::cout << "w(" << x << ")=" << w1 << ", " << n << " iterations\n";
    return w1;
}

// int main(int argc, char *argv[]) {
//     double x = std::stod(argv[1]);
//     LambertW0(x);
// }