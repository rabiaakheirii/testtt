#include "percolation/percolation.hpp"
#include "utils/random.hpp"
#include <cmath>
#include <fstream>
#include <iostream>

namespace qcd {

    PercolationResult run_percolation(int N_strings,
                                      double R_fm,
                                      double r_string_fm,
                                      double overlap_prob,
                                      unsigned seed)
    {
        PercolationResult res;
        RNG rng(seed);

        res.strings.reserve(N_strings);
        for (int i = 0; i < N_strings; ++i) {
            const double u   = rng.uniform();
            const double rad = R_fm * std::sqrt(u);
            const double phi = rng.uniform(0.0, 2.0 * M_PI);
            String s;
            s.x = rad * std::cos(phi);
            s.y = rad * std::sin(phi);
            s.r = r_string_fm;
            if (rng.uniform() < overlap_prob)
                res.strings.push_back(s);
        }
        res.n_strings = static_cast<int>(res.strings.size());

        UnionFind uf(res.n_strings);
        for (int i = 0; i < res.n_strings; ++i) {
            for (int j = i + 1; j < res.n_strings; ++j) {
                const double dx = res.strings[i].x - res.strings[j].x;
                const double dy = res.strings[i].y - res.strings[j].y;
                const double rr = res.strings[i].r + res.strings[j].r;
                if (dx*dx + dy*dy <= rr*rr) uf.unite(i, j);
            }
        }
        res.n_clusters           = uf.cluster_count();
        res.largest_cluster_size = uf.largest_cluster_size();
        res.percolation_fraction = (res.n_strings > 0)
            ? static_cast<double>(res.largest_cluster_size) / res.n_strings : 0.0;

        res.grid_nx = 100;
        res.grid_ny = 100;
        res.grid_xmin = -R_fm; res.grid_xmax = R_fm;
        res.grid_ymin = -R_fm; res.grid_ymax = R_fm;
        res.density.assign(res.grid_nx * res.grid_ny, 0.0);

        const double dxg = (res.grid_xmax - res.grid_xmin) / (res.grid_nx - 1);
        const double dyg = (res.grid_ymax - res.grid_ymin) / (res.grid_ny - 1);

        double sum = 0.0;
        for (int iy = 0; iy < res.grid_ny; ++iy) {
            for (int ix = 0; ix < res.grid_nx; ++ix) {
                const double x = res.grid_xmin + ix * dxg;
                const double y = res.grid_ymin + iy * dyg;
                const int c = coverage(res.strings, x, y);
                res.density[iy * res.grid_nx + ix] = static_cast<double>(c);
                sum += c;
            }
        }
        res.energy_density_avg = sum / (res.grid_nx * res.grid_ny);

        std::cout << "[percolation] N=" << res.n_strings
                  << " clusters=" << res.n_clusters
                  << " largest=" << res.largest_cluster_size
                  << " P=" << res.percolation_fraction
                  << " <n>=" << res.energy_density_avg << "\n";

        return res;
    }

    void write_density_map(const PercolationResult& res,
                           const std::string& filename)
    {
        std::ofstream out(filename);
        out << "# x[fm]  y[fm]  n_strings\n";
        const double dxg = (res.grid_xmax - res.grid_xmin) / (res.grid_nx - 1);
        const double dyg = (res.grid_ymax - res.grid_ymin) / (res.grid_ny - 1);
        for (int iy = 0; iy < res.grid_ny; ++iy) {
            const double y = res.grid_ymin + iy * dyg;
            for (int ix = 0; ix < res.grid_nx; ++ix) {
                const double x = res.grid_xmin + ix * dxg;
                out << x << "  " << y << "  "
                    << res.density[iy * res.grid_nx + ix] << "\n";
            }
            out << "\n";
        }
        std::cout << "[percolation] wrote " << filename << "\n";
    }
}
