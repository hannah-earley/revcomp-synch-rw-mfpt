#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include "walk.hpp"

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
    return w1;
}

unsigned stou(std::string const& str, size_t* idx, int base) {
    unsigned long result = std::stoul(str, idx, base);
    if (result > std::numeric_limits<unsigned>::max()) {
        throw std::out_of_range("stou");
    }
    return result;
}

std::string format_si(double quantity, std::string unit) {
    const char* prefixes_g[] = {"", "k", "M", "G", "T", "P", "E", "Z", "Y", NULL};
    const char* prefixes_l[] = {"", "m", "u", "n", "p", "f", "a", "z", "y", NULL};

    const char **prefixes = prefixes_g;
    if (quantity == 0) {
    } else if (quantity >= 1) {
        while (quantity >= 1000 && *(prefixes + 1)) {
            quantity /= 1000;
            prefixes++;
        }
    } else {
        prefixes = prefixes_l;
        while (quantity < 1 && *(prefixes + 1)) {
            quantity *= 1000;
            prefixes++;
        }
    }

    char buf[50];
    std::sprintf(buf, "%g", quantity);

    std::string out = "";
    out += buf;
    out += " ";
    out += *prefixes;
    out += unit;
    return out;
}

std::string format_time(double seconds) {
    std::string out = "";
    unsigned long long secs = seconds;
    double subseconds = seconds - secs;

    if (secs > 365 * 86400) {
        unsigned long long years = secs / (365 * 86400);
        secs -= years * (365 * 86400);
        out += std::to_string(years);
        out += " y ";
    }

    if (secs > 86400) {
        unsigned long long days = secs / 86400;
        secs -= days * 86400;
        out += std::to_string(days);
        out += " d ";
    }

    if (secs > 3600) {
        unsigned long long hours = secs / 3600;
        secs -= hours * 3600;
        out += std::to_string(hours);
        out += " h ";
    }

    if (secs > 60) {
        unsigned long long mins = secs / 60;
        secs -= mins * 60;
        out += std::to_string(mins);
        out += " m ";
    }

    out += format_si(subseconds + secs, "s");
    return out;
}

std::string format_timepoint(std::chrono::time_point<std::chrono::system_clock> tp,
                             const char* fmt)
{
    std::stringstream ss;
    std::time_t time = std::chrono::system_clock::to_time_t(tp);
    ss << std::put_time(std::localtime(&time), fmt);
    return ss.str();
}