#ifndef BLASTWAVE_CONSTANTS_H
#define BLASTWAVE_CONSTANTS_H

#include <cmath>
#include <string>

namespace bw {

constexpr int kMaxParticles = 6;
constexpr int kMaxCentralities = 20;

// Particle indices: 0=pi+, 1=pi-, 2=K+, 3=K-, 4=p, 5=anti-p
constexpr double kMasses[kMaxParticles] = {
    0.13957061, 0.13957061,
    0.493667,   0.493667,
    0.938272,   0.938272
};

constexpr const char* kParticleNames[kMaxParticles] = {
    "pi+", "pi-", "K+", "K-", "p", "pbar"
};

constexpr const char* kParticleTitles[kMaxParticles] = {
    "#pi^{+}", "#pi^{#minus}", "K^{+}", "K^{#minus}", "p", "#bar{p}"
};

// Species group for shared fit ranges: 0=pi, 1=K, 2=p
inline int SpeciesGroup(int particle)
{
    return particle / 2;
}

inline double MtMinusM(int particle, double pT)
{
    const double m = kMasses[particle];
    return std::sqrt(pT * pT + m * m) - m;
}

}  // namespace bw

#endif
