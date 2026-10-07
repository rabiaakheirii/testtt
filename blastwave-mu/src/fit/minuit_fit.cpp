#include "fit/minuit_fit.hpp"
#include "fit/tbeta_penalty.hpp"

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace bwmu {

    namespace {

        std::vector<int> build_norm_groups(const std::vector<Species>& species) {
            std::vector<int> g(species.size(), -1);
            std::vector<double> seen;
            for (size_t i = 0; i < species.size(); ++i) {
                int found = -1;
                for (size_t k = 0; k < seen.size(); ++k)
                    if (std::abs(seen[k] - species[i].mass) < 1e-6)
                        { found = (int)k; break; }
                if (found < 0) { found = (int)seen.size(); seen.push_back(species[i].mass); }
                g[i] = found;
            }
            return g;
        }

        class Chi2FCN : public ROOT::Minuit2::FCNBase {
        public:
            Chi2FCN(const std::vector<Species>& sp,
                    const std::vector<Spectrum>& dt,
                    const BWParams& bwp,
                    Statistics st,
                    double lam)
                : species_(sp), data_(dt), bw_(bwp), stats_(st), lambda_(lam)
            {
                groups_ = build_norm_groups(sp);
            }

            double Up() const override { return 1.0; }

            double operator()(const std::vector<double>& x) const override {
                FitParams fp;
                fp.T = x[0]; fp.beta_s = x[1];
                fp.mu_B = x[2]; fp.mu_Q = x[3]; fp.mu_S = x[4];

                auto r = compute_chi2(species_, data_, bw_, fp, stats_, groups_);
                if (r.ndof <= 0) return 1e12;

                const double pen = tbeta_penalty(fp.T, fp.beta_s, lambda_);
                return r.chi2 + pen;
            }

            Chi2Result compute_at(const std::vector<double>& x) const {
                FitParams fp;
                fp.T = x[0]; fp.beta_s = x[1];
                fp.mu_B = x[2]; fp.mu_Q = x[3]; fp.mu_S = x[4];
                return compute_chi2(species_, data_, bw_, fp, stats_, groups_);
            }

        private:
            const std::vector<Species>& species_;
            const std::vector<Spectrum>& data_;
            const BWParams& bw_;
            Statistics stats_;
            double lambda_;
            std::vector<int> groups_;
        };

    } // anonymous namespace

    JointFitResult joint_fit_minuit(
        const std::vector<Species>& species,
        const std::vector<Spectrum>& data,
        const BWParams& bw,
        const FitBounds& bounds,
        Statistics stats,
        double lambda_tbeta,
        bool verbose)
    {
        JointFitResult result;
        result.converged = false;
        result.evaluations = 0;

        Chi2FCN fcn(species, data, bw, stats, lambda_tbeta);

        ROOT::Minuit2::MnUserParameters upar;
        // (name, start, step, min, max)
        upar.Add("T",      0.120, 0.005, bounds.T_min,   bounds.T_max);
        upar.Add("beta_s", 0.600, 0.020, 0.55, 0.70);   // restricted to physical range
        upar.Add("mu_B",   0.400,  0.020, 0.300, 0.700);   // NICA
        upar.Add("mu_Q",   0.010,  0.003, 0.000, 0.030);   // NICA positive
        upar.Add("mu_S",   0.020,  0.003, 0.000, 0.050);   // NICA

        if (verbose)
            std::cout << "[minuit] starting, lambda_tbeta = " << lambda_tbeta << "\n";

        ROOT::Minuit2::MnMigrad migrad(fcn, upar);
        auto min = migrad(0, 1e-4);

        auto st = min.UserState();

        if (verbose) {
            std::cout << "[minuit] chi2_tot = " << min.Fval()
                      << "  valid = " << (min.IsValid() ? "yes" : "no") << "\n";
            std::cout << "[minuit] T      = " << st.Value(0)*1000 << " MeV\n";
            std::cout << "[minuit] beta_s = " << st.Value(1) << "\n";
            std::cout << "[minuit] mu_B   = " << st.Value(2)*1000 << " MeV\n";
            std::cout << "[minuit] mu_Q   = " << st.Value(3)*1000 << " MeV\n";
            std::cout << "[minuit] mu_S   = " << st.Value(4)*1000 << " MeV\n";
        }

        std::vector<double> x = { st.Value(0), st.Value(1), st.Value(2),
                                  st.Value(3), st.Value(4) };
        auto chi = fcn.compute_at(x);

        result.params.T      = x[0];
        result.params.beta_s = x[1];
        result.params.mu_B   = x[2];
        result.params.mu_Q   = x[3];
        result.params.mu_S   = x[4];
        result.norms        = chi.norms;
        result.chi2         = chi.chi2;   // pure data chi2
        result.ndof         = chi.ndof;
        result.chi2_ndof    = chi.chi2_ndof;
        result.converged    = min.IsValid();

        return result;
    }

} // namespace bwmu
