#pragma once
#include <random>

namespace qcd {
    class RNG {
    public:
        RNG(unsigned seed = 12345) : gen(seed), dist(0.0, 1.0) {}
        double uniform() { return dist(gen); }
        double uniform(double a, double b) { return a + (b-a)*uniform(); }
    private:
        std::mt19937 gen;
        std::uniform_real_distribution<double> dist;
    };
}