#pragma once
#include <random>

namespace bwmu {

    // Simple MT19937 wrapper for reproducible sampling
    class RNG {
    public:
        explicit RNG(unsigned seed = 12345);
        double uniform();                     // [0, 1)
        double uniform(double a, double b);   // [a, b)
        double gaussian();                    // standard normal
        double gaussian(double mu, double sigma);
    private:
        std::mt19937 gen_;
        std::uniform_real_distribution<double> udist_;
        std::normal_distribution<double>       ndist_;
    };

} // namespace bwmu
