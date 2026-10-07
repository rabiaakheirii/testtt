#include "utils/random.hpp"

namespace bwmu {

    RNG::RNG(unsigned seed)
        : gen_(seed), udist_(0.0, 1.0), ndist_(0.0, 1.0) {}

    double RNG::uniform() {
        return udist_(gen_);
    }

    double RNG::uniform(double a, double b) {
        return a + (b - a) * udist_(gen_);
    }

    double RNG::gaussian() {
        return ndist_(gen_);
    }

    double RNG::gaussian(double mu, double sigma) {
        return mu + sigma * ndist_(gen_);
    }

} // namespace bwmu
