#include "percolation/cluster.hpp"
#include <algorithm>
#include <cmath>

namespace qcd {

    UnionFind::UnionFind(int n) : parent(n), rank_(n, 0), n_(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int UnionFind::find(int i) const {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }

    void UnionFind::unite(int i, int j) {
        int ri = find(i);
        int rj = find(j);
        if (ri == rj) return;
        if (rank_[ri] < rank_[rj]) std::swap(ri, rj);
        parent[rj] = ri;
        if (rank_[ri] == rank_[rj]) ++rank_[ri];
    }

    int UnionFind::cluster_count() const {
        int c = 0;
        for (int i = 0; i < n_; ++i)
            if (parent[i] == i) ++c;
        return c;
    }

    int UnionFind::largest_cluster_size() const {
        std::vector<int> sz(n_, 0);
        for (int i = 0; i < n_; ++i) sz[find(i)]++;
        return *std::max_element(sz.begin(), sz.end());
    }

    int coverage(const std::vector<String>& strings, double x, double y) {
        int n = 0;
        for (const auto& s : strings) {
            const double dx = x - s.x;
            const double dy = y - s.y;
            if (dx*dx + dy*dy <= s.r * s.r) ++n;
        }
        return n;
    }

} // namespace qcd
