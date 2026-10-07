#pragma once
#include "percolation/cluster.hpp"
#include <string>
#include <vector>

namespace qcd {

    struct PercolationResult {
        std::vector<String> strings;
        int    n_strings;
        int    largest_cluster_size;
        int    n_clusters;
        double percolation_fraction;   // largest / total
        double energy_density_avg;     // <n(x,y)>
        int    grid_nx, grid_ny;
        double grid_xmin, grid_xmax;
        double grid_ymin, grid_ymax;
        std::vector<double> density;   // grid_nx * grid_ny
    };

    PercolationResult run_percolation(int N_strings,
                                      double R_fm,
                                      double r_string_fm,
                                      double overlap_prob,
                                      unsigned seed);

    void write_density_map(const PercolationResult& res,
                           const std::string& filename);
}
