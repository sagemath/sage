struct PathCandidate {
    double cost;
    int path_idx;
    int dev_idx;
    bool is_simple;
    PathCandidate() : cost(0), path_idx(0), dev_idx(0), is_simple(false) {}
    PathCandidate(double c, int p, int d, bool s)
        : cost(c), path_idx(p), dev_idx(d), is_simple(s) {}
    // compare for "min"-heap
    bool operator<(const PathCandidate& other) const {
        return cost > other.cost;
    }
};