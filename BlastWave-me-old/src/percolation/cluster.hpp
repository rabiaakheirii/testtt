#pragma once
#include <vector>

namespace qcd {

    struct String {
        double x, y;
        double r;
    };

    class UnionFind {
    public:
        explicit UnionFind(int n);
        int  find(int i) const;
        void unite(int i, int j);
        int  cluster_count() const;
        int  largest_cluster_size() const;
    private:
        mutable std::vector<int> parent;
        std::vector<int> rank_;
        int n_;
    };

    int coverage(const std::vector<String>& strings, double x, double y);

} // namespace qcd
