#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using Point = std::array<double, 3>;
using Points = std::vector<Point>;

namespace {

using Clock = std::chrono::steady_clock;

struct Rotation {
    double a[3][3];
};

struct Candidate {
    Points q;
    double estimate = 0.0;
};

Clock::time_point start_time;
std::mt19937_64 rng;

double elapsed() {
    return std::chrono::duration<double>(Clock::now() - start_time).count();
}

uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

Point rotate(const Rotation& r, const Point& p) {
    return {
        r.a[0][0] * p[0] + r.a[0][1] * p[1] + r.a[0][2] * p[2],
        r.a[1][0] * p[0] + r.a[1][1] * p[1] + r.a[1][2] * p[2],
        r.a[2][0] * p[0] + r.a[2][1] * p[1] + r.a[2][2] * p[2]
    };
}

Rotation identity_rotation() {
    Rotation r{};
    for (int i = 0; i < 3; ++i) r.a[i][i] = 1.0;
    return r;
}

Rotation random_rotation() {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double u1 = unit(rng);
    const double u2 = unit(rng);
    const double u3 = unit(rng);
    const double pi = std::acos(-1.0);
    const double x = std::sqrt(1.0 - u1) * std::sin(2.0 * pi * u2);
    const double y = std::sqrt(1.0 - u1) * std::cos(2.0 * pi * u2);
    const double z = std::sqrt(u1) * std::sin(2.0 * pi * u3);
    const double w = std::sqrt(u1) * std::cos(2.0 * pi * u3);

    Rotation r{};
    r.a[0][0] = 1.0 - 2.0 * (y * y + z * z);
    r.a[0][1] = 2.0 * (x * y - z * w);
    r.a[0][2] = 2.0 * (x * z + y * w);
    r.a[1][0] = 2.0 * (x * y + z * w);
    r.a[1][1] = 1.0 - 2.0 * (x * x + z * z);
    r.a[1][2] = 2.0 * (y * z - x * w);
    r.a[2][0] = 2.0 * (x * z - y * w);
    r.a[2][1] = 2.0 * (y * z + x * w);
    r.a[2][2] = 1.0 - 2.0 * (x * x + y * y);
    return r;
}

void tighten_box(Points& q) {
    for (int axis = 0; axis < 3; ++axis) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (const Point& p : q) {
            lo = std::min(lo, p[axis]);
            hi = std::max(hi, p[axis]);
        }
        const double span = hi - lo;
        if (span < 1e-14) {
            for (Point& p : q) p[axis] = 0.5;
        } else {
            for (Point& p : q) {
                p[axis] = std::clamp((p[axis] - lo) / span, 0.0, 1.0);
            }
        }
    }
}

double exact_min_distance(const Points& q) {
    double best2 = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < q.size(); ++i) {
        for (size_t j = i + 1; j < q.size(); ++j) {
            const double dx = q[i][0] - q[j][0];
            const double dy = q[i][1] - q[j][1];
            const double dz = q[i][2] - q[j][2];
            best2 = std::min(best2, dx * dx + dy * dy + dz * dz);
        }
    }
    return std::sqrt(best2);
}

void retain_candidate(std::vector<Candidate>& pool, Points q, double estimate, int limit) {
    tighten_box(q);
    Candidate c{std::move(q), estimate};
    pool.push_back(std::move(c));
    std::sort(pool.begin(), pool.end(), [](const Candidate& x, const Candidate& y) {
        return x.estimate > y.estimate;
    });
    if (static_cast<int>(pool.size()) > limit) pool.resize(limit);
}

Points simple_cubic(int n) {
    int s = 1;
    while (1LL * s * s * s < n) ++s;

    struct Item {
        uint64_t key;
        Point p;
    };
    std::vector<Item> all;
    all.reserve(1LL * s * s * s);
    for (int x = 0; x < s; ++x) {
        for (int y = 0; y < s; ++y) {
            for (int z = 0; z < s; ++z) {
                const uint64_t id = (static_cast<uint64_t>(x) * s + y) * s + z;
                const uint64_t key = splitmix64(id ^ (static_cast<uint64_t>(s) << 48));
                all.push_back({key, Point{double(x), double(y), double(z)}});
            }
        }
    }
    std::sort(all.begin(), all.end(), [](const Item& a, const Item& b) {
        return a.key < b.key;
    });
    Points q;
    q.reserve(n);
    for (int i = 0; i < n; ++i) q.push_back(all[i].p);
    tighten_box(q);
    return q;
}

void add_corner_candidate(int n, std::vector<Candidate>& pool, int limit) {
    if (n > 8) return;
    double best = -1.0;
    Points best_q;
    for (int mask = 0; mask < 256; ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) != n) continue;
        Points q;
        for (int v = 0; v < 8; ++v) {
            if (mask & (1 << v)) {
                q.push_back(Point{double((v >> 2) & 1), double((v >> 1) & 1), double(v & 1)});
            }
        }
        tighten_box(q);
        const double d = exact_min_distance(q);
        if (d > best) {
            best = d;
            best_q = std::move(q);
        }
    }
    if (!best_q.empty()) retain_candidate(pool, std::move(best_q), best, limit);
}

Points generate_sc_patch(int k) {
    Points raw;
    raw.reserve(1LL * (2 * k + 1) * (2 * k + 1) * (2 * k + 1));
    for (int i = -k; i <= k; ++i)
        for (int j = -k; j <= k; ++j)
            for (int z = -k; z <= k; ++z)
                raw.push_back(Point{double(i), double(j), double(z)});
    return raw;
}

Points generate_fcc_patch(int k) {
    Points raw;
    for (int i = -k; i <= k; ++i) {
        for (int j = -k; j <= k; ++j) {
            for (int z = -k; z <= k; ++z) {
                if (((i + j + z) & 1) == 0) raw.push_back(Point{double(i), double(j), double(z)});
            }
        }
    }
    return raw;
}

Points rectangular_fcc(int n) {
    const int limit = 2 * static_cast<int>(std::ceil(std::cbrt(2.0 * n))) + 8;
    int best_a = 1, best_b = 1, best_c = 1;
    double best_distance = -1.0;
    for (int a = 1; a <= limit; ++a) {
        for (int b = a; b <= limit; ++b) {
            for (int c = b; c <= limit; ++c) {
                const long long count = (1LL * (a + 1) * (b + 1) * (c + 1) + 1) / 2;
                if (count < n) continue;
                const double da = 1.0 / a;
                const double db = 1.0 / b;
                const double dc = 1.0 / c;
                const double distance = std::min({
                    2.0 * da, 2.0 * db, 2.0 * dc,
                    std::hypot(da, db), std::hypot(da, dc), std::hypot(db, dc)
                });
                if (distance > best_distance + 1e-15) {
                    best_distance = distance;
                    best_a = a;
                    best_b = b;
                    best_c = c;
                }
            }
        }
    }

    struct Item {
        double shell;
        double radial;
        Point p;
    };
    std::vector<Item> items;
    for (int x = 0; x <= best_a; ++x) {
        for (int y = 0; y <= best_b; ++y) {
            for (int z = 0; z <= best_c; ++z) {
                if ((x + y + z) & 1) continue;
                const Point p{double(x) / best_a, double(y) / best_b, double(z) / best_c};
                const double dx = std::abs(p[0] - 0.5);
                const double dy = std::abs(p[1] - 0.5);
                const double dz = std::abs(p[2] - 0.5);
                items.push_back({std::max({dx, dy, dz}), dx * dx + dy * dy + dz * dz, p});
            }
        }
    }
    auto less = [](const Item& left, const Item& right) {
        if (left.shell != right.shell) return left.shell < right.shell;
        if (left.radial != right.radial) return left.radial < right.radial;
        return left.p < right.p;
    };
    if (static_cast<int>(items.size()) > n) {
        std::nth_element(items.begin(), items.begin() + n, items.end(), less);
        items.resize(n);
    }
    Points q;
    q.reserve(n);
    for (const Item& item : items) q.push_back(item.p);
    tighten_box(q);
    return q;
}

Points tetrahedral_cluster(int n) {
    int side = 5;
    Points root;
    for (int family_index = 1;; ++family_index) {
        side = (family_index & 1) ? 6 * ((family_index + 1) / 2) - 1
                                  : 6 * (family_index / 2) + 2;
        const int residue = (side + 1) % 6;
        Points candidate;
        for (int x = 0; x <= side; ++x) {
            for (int y = 0; y <= side; ++y) {
                for (int z = 0; z <= side; ++z) {
                    bool in_corner = false;
                    const int corners[4][3] = {
                        {0, 0, 0}, {0, side, side},
                        {side, 0, side}, {side, side, 0}
                    };
                    for (const auto& corner : corners) {
                        const int dx = std::abs(x - corner[0]);
                        const int dy = std::abs(y - corner[1]);
                        const int dz = std::abs(z - corner[2]);
                        const int taxicab = dx + dy + dz;
                        if (dx % 3 == 0 && dy % 3 == 0 && dz % 3 == 0
                            && taxicab <= side - 5 && taxicab % 6 == residue) {
                            in_corner = true;
                            break;
                        }
                    }

                    const int sum = x + y + z;
                    const bool in_center = x % 3 == 1 && y % 3 == 1 && z % 3 == 1
                        && sum >= side + 1 && sum % 6 == residue
                        && x + y - z <= side - 1
                        && x - y + z <= side - 1
                        && -x + y + z <= side - 1;
                    if (in_corner || in_center) {
                        candidate.push_back(Point{double(x), double(y), double(z)});
                    }
                }
            }
        }
        if (static_cast<int>(candidate.size()) >= n) {
            root = std::move(candidate);
            break;
        }
    }

    if (static_cast<int>(root.size()) > n) {
        std::sort(root.begin(), root.end(), [side](const Point& a, const Point& b) {
            auto shell = [side](const Point& p) {
                return std::max({std::abs(2.0 * p[0] - side),
                                 std::abs(2.0 * p[1] - side),
                                 std::abs(2.0 * p[2] - side)});
            };
            const double shell_a = shell(a);
            const double shell_b = shell(b);
            if (shell_a != shell_b) return shell_a > shell_b;
            const double radius_a = (2.0 * a[0] - side) * (2.0 * a[0] - side)
                                  + (2.0 * a[1] - side) * (2.0 * a[1] - side)
                                  + (2.0 * a[2] - side) * (2.0 * a[2] - side);
            const double radius_b = (2.0 * b[0] - side) * (2.0 * b[0] - side)
                                  + (2.0 * b[1] - side) * (2.0 * b[1] - side)
                                  + (2.0 * b[2] - side) * (2.0 * b[2] - side);
            if (radius_a != radius_b) return radius_a > radius_b;
            return a < b;
        });
        root.resize(n);
    }
    tighten_box(root);
    return root;
}

Points diagonal_fcc_cluster(int n) {
    Points root;
    for (int t = 1; root.empty(); ++t) {
        const double side = t + 1.0 / 3.0;
        const double cut = 1.5 * side;
        for (int parity = 0; parity < 2; ++parity) {
            Points candidate;
            for (int x = 0; x <= t; ++x) {
                for (int y = 0; y <= t; ++y) {
                    for (int z = 0; z <= t; ++z) {
                        const int sum = x + y + z;
                        if ((sum & 1) == parity && sum < cut) {
                            candidate.push_back(Point{double(x), double(y), double(z)});
                        }
                        if ((sum & 1) != parity && sum + 1.0 > cut) {
                            candidate.push_back(Point{x + 1.0 / 3.0,
                                                      y + 1.0 / 3.0,
                                                      z + 1.0 / 3.0});
                        }
                    }
                }
            }
            if (static_cast<int>(candidate.size()) >= n
                && (root.empty() || candidate.size() < root.size())) {
                root = std::move(candidate);
            }
        }
    }

    if (static_cast<int>(root.size()) > n) {
        double side = 0.0;
        for (const Point& p : root)
            for (double x : p) side = std::max(side, x);
        std::sort(root.begin(), root.end(), [side](const Point& a, const Point& b) {
            auto shell = [side](const Point& p) {
                return std::max({std::abs(2.0 * p[0] - side),
                                 std::abs(2.0 * p[1] - side),
                                 std::abs(2.0 * p[2] - side)});
            };
            const double shell_a = shell(a);
            const double shell_b = shell(b);
            if (shell_a != shell_b) return shell_a > shell_b;
            return a < b;
        });
        root.resize(n);
    }
    tighten_box(root);
    return root;
}

Points staggered_rows(int n) {
    const int root = static_cast<int>(std::ceil(std::cbrt(double(n))));
    const int limit = 2 * root + 2;
    int best_x = 2, best_y = 2, best_z = 2;
    double best_step = 0.5;
    double best_distance = -1.0;
    for (int nx = 2; nx <= limit; ++nx) {
        const int intervals = nx - 1;
        for (int ny = 2; ny <= limit; ++ny) {
            const double dy = 1.0 / (ny - 1);
            for (int nz = 2; nz <= limit; ++nz) {
                if (1LL * nx * ny * nz < n) continue;
                const double dz = 1.0 / (nz - 1);
                const double transverse = std::min(dy, dz);
                const double fixed = std::min({2.0 * dy, 2.0 * dz,
                                               std::hypot(dy, dz)});
                auto distance = [&](double step) {
                    const double phase = 1.0 - intervals * step;
                    const double offset = std::min(phase, step - phase);
                    return std::min({step, std::hypot(offset, transverse), fixed});
                };
                double lo = 1.0 / nx;
                double hi = 1.0 / intervals;
                for (int iteration = 0; iteration < 70; ++iteration) {
                    const double left = (2.0 * lo + hi) / 3.0;
                    const double right = (lo + 2.0 * hi) / 3.0;
                    if (distance(left) < distance(right)) lo = left;
                    else hi = right;
                }
                const double step = 0.5 * (lo + hi);
                const double d = distance(step);
                if (d > best_distance + 1e-14) {
                    best_distance = d;
                    best_step = step;
                    best_x = nx;
                    best_y = ny;
                    best_z = nz;
                }
            }
        }
    }

    const double phase = 1.0 - (best_x - 1) * best_step;
    Points points;
    points.reserve(1LL * best_x * best_y * best_z);
    for (int y = 0; y < best_y; ++y) {
        for (int z = 0; z < best_z; ++z) {
            const double shift = ((y + z) & 1) ? phase : 0.0;
            for (int x = 0; x < best_x; ++x) {
                points.push_back(Point{shift + x * best_step,
                                       double(y) / (best_y - 1),
                                       double(z) / (best_z - 1)});
            }
        }
    }
    if (static_cast<int>(points.size()) > n) {
        std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
            auto shell = [](const Point& p) {
                return std::max({std::abs(p[0] - 0.5),
                                 std::abs(p[1] - 0.5),
                                 std::abs(p[2] - 0.5)});
            };
            const double shell_a = shell(a);
            const double shell_b = shell(b);
            if (shell_a != shell_b) return shell_a > shell_b;
            return a < b;
        });
        points.resize(n);
    }
    tighten_box(points);
    return points;
}

Points generate_bct_patch(int k, double aspect) {
    Points raw;
    raw.reserve(2LL * (2 * k + 1) * (2 * k + 1) * (2 * k + 1));
    for (int i = -k; i <= k; ++i) {
        for (int j = -k; j <= k; ++j) {
            for (int z = -k; z <= k; ++z) {
                raw.push_back(Point{double(i), double(j), aspect * z});
                raw.push_back(Point{i + 0.5, j + 0.5, aspect * (z + 0.5)});
            }
        }
    }
    return raw;
}

Points generate_barlow_patch(int k, bool abc) {
    constexpr double sqrt3 = 1.7320508075688772935;
    const double layer = std::sqrt(2.0 / 3.0);
    const Point offsets[3] = {
        Point{0.0, 0.0, 0.0},
        Point{0.5, sqrt3 / 6.0, 0.0},
        Point{0.0, sqrt3 / 3.0, 0.0}
    };
    Points raw;
    raw.reserve(1LL * (2 * k + 1) * (2 * k + 1) * (2 * k + 1));
    for (int z = -k; z <= k; ++z) {
        const int phase = abc ? ((z % 3) + 3) % 3 : ((z % 2) + 2) % 2;
        for (int j = -k; j <= k; ++j) {
            for (int i = -k; i <= k; ++i) {
                raw.push_back(Point{
                    i + 0.5 * (j & 1) + offsets[phase][0],
                    0.5 * sqrt3 * j + offsets[phase][1],
                    layer * z
                });
            }
        }
    }
    return raw;
}

Points lower_fcc_repair_pattern(int n) {
    Points base;
    int intervals = 0;
    for (int m = 1;; ++m) {
        Points lattice;
        for (int x = 0; x <= m; ++x)
            for (int y = 0; y <= m; ++y)
                for (int z = 0; z <= m; ++z)
                    if (((x + y + z) & 1) == 0)
                        lattice.push_back(Point{double(x) / m, double(y) / m, double(z) / m});
        if (static_cast<int>(lattice.size()) > n) break;
        base = std::move(lattice);
        intervals = m;
    }
    if (static_cast<int>(base.size()) == n || intervals == 0) return base;

    const int grid = 2 * intervals;
    Points candidates;
    for (int x = 0; x <= grid; ++x)
        for (int y = 0; y <= grid; ++y)
            for (int z = 0; z <= grid; ++z)
                candidates.push_back(Point{double(x) / grid, double(y) / grid, double(z) / grid});
    std::vector<double> clearance2(candidates.size(), std::numeric_limits<double>::infinity());
    for (int x = 0; x <= grid; ++x) {
        for (int y = 0; y <= grid; ++y) {
            for (int z = 0; z <= grid; ++z) {
                double nearest = std::numeric_limits<double>::infinity();
                const int x0 = x / 2;
                const int y0 = y / 2;
                const int z0 = z / 2;
                for (int bx = std::max(0, x0 - 1); bx <= std::min(intervals, x0 + 2); ++bx) {
                    for (int by = std::max(0, y0 - 1); by <= std::min(intervals, y0 + 2); ++by) {
                        for (int bz = std::max(0, z0 - 1); bz <= std::min(intervals, z0 + 2); ++bz) {
                            if ((bx + by + bz) & 1) continue;
                            const double dx = double(x - 2 * bx) / grid;
                            const double dy = double(y - 2 * by) / grid;
                            const double dz = double(z - 2 * bz) / grid;
                            nearest = std::min(nearest, dx * dx + dy * dy + dz * dz);
                        }
                    }
                }
                clearance2[(static_cast<size_t>(x) * (grid + 1) + y) * (grid + 1) + z]
                    = nearest;
            }
        }
    }
    while (static_cast<int>(base.size()) < n) {
        const size_t index = std::max_element(clearance2.begin(), clearance2.end())
                           - clearance2.begin();
        const Point point = candidates[index];
        base.push_back(point);
        clearance2[index] = -1.0;
        for (size_t c = 0; c < candidates.size(); ++c) {
            if (clearance2[c] < 0.0) continue;
            const double dx = candidates[c][0] - point[0];
            const double dy = candidates[c][1] - point[1];
            const double dz = candidates[c][2] - point[2];
            clearance2[c] = std::min(clearance2[c], dx * dx + dy * dy + dz * dz);
        }
    }
    return base;
}

void crop_candidate(const Points& raw, int n, double lattice_distance, int sample,
                    std::vector<Candidate>& pool, int limit) {
    struct Item {
        double key;
        double tie;
        int index;
        Point p;
    };

    const bool aligned = sample < 8;
    const Rotation rotation = aligned ? identity_rotation() : random_rotation();
    std::uniform_real_distribution<double> phase_dist(-0.45, 0.45);
    std::uniform_real_distribution<double> shape_dist(0.90, 1.10);
    Point center{0.0, 0.0, 0.0};
    Point shape{1.0, 1.0, 1.0};
    if (sample != 0) {
        for (int a = 0; a < 3; ++a) {
            center[a] = phase_dist(rng);
            if (!aligned) shape[a] = shape_dist(rng);
        }
    }

    std::vector<Item> items;
    items.reserve(raw.size());
    for (int index = 0; index < static_cast<int>(raw.size()); ++index) {
        const Point p = rotate(rotation, raw[index]);
        const double dx = std::abs(p[0] - center[0]) / shape[0];
        const double dy = std::abs(p[1] - center[1]) / shape[1];
        const double dz = std::abs(p[2] - center[2]) / shape[2];
        const double key = std::max({dx, dy, dz});
        const double tie = dx * dx + dy * dy + dz * dz;
        items.push_back({key, tie, index, p});
    }

    auto less = [](const Item& a, const Item& b) {
        if (a.key != b.key) return a.key < b.key;
        if (a.tie != b.tie) return a.tie < b.tie;
        return a.index < b.index;
    };
    if (static_cast<int>(items.size()) > n) {
        std::nth_element(items.begin(), items.begin() + n, items.end(), less);
        items.resize(n);
    }

    Point lo{1e100, 1e100, 1e100};
    Point hi{-1e100, -1e100, -1e100};
    Points selected;
    selected.reserve(n);
    for (const Item& item : items) {
        selected.push_back(item.p);
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], item.p[a]);
            hi[a] = std::max(hi[a], item.p[a]);
        }
    }
    const double max_span = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
    const double estimate = lattice_distance / std::max(max_span, 1e-12);
    retain_candidate(pool, std::move(selected), estimate, limit);
}

struct HashGrid {
    int g = 1;
    std::vector<int> head;
    std::vector<int> next;

    int cell(double x) const {
        return std::clamp(static_cast<int>(x * g), 0, g - 1);
    }

    int index(int x, int y, int z) const {
        return (x * g + y) * g + z;
    }

    void build(const Points& q, double distance) {
        g = std::max(1, static_cast<int>(std::floor(1.0 / std::max(distance, 1e-9))));
        g = std::min(g, 256);
        head.assign(static_cast<size_t>(g) * g * g, -1);
        next.assign(q.size(), -1);
        for (int i = 0; i < static_cast<int>(q.size()); ++i) {
            const int id = index(cell(q[i][0]), cell(q[i][1]), cell(q[i][2]));
            next[i] = head[id];
            head[id] = i;
        }
    }
};

void jitter(Points& q, double amplitude);

bool has_violation(const Points& q, double distance) {
    HashGrid grid;
    grid.build(q, distance);
    const double distance2 = distance * distance;
    for (int i = 0; i < static_cast<int>(q.size()); ++i) {
        const int cx = grid.cell(q[i][0]);
        const int cy = grid.cell(q[i][1]);
        const int cz = grid.cell(q[i][2]);
        for (int dx = -1; dx <= 1; ++dx) {
            const int x = cx + dx;
            if (x < 0 || x >= grid.g) continue;
            for (int dy = -1; dy <= 1; ++dy) {
                const int y = cy + dy;
                if (y < 0 || y >= grid.g) continue;
                for (int dz = -1; dz <= 1; ++dz) {
                    const int z = cz + dz;
                    if (z < 0 || z >= grid.g) continue;
                    for (int j = grid.head[grid.index(x, y, z)]; j != -1; j = grid.next[j]) {
                        if (j <= i) continue;
                        const double ax = q[i][0] - q[j][0];
                        const double ay = q[i][1] - q[j][1];
                        const double az = q[i][2] - q[j][2];
                        if (ax * ax + ay * ay + az * az < distance2) return true;
                    }
                }
            }
        }
    }
    return false;
}

struct ForceState {
    double max_overlap = 0.0;
};

ForceState overlap_forces(const Points& q, double target, Points& force) {
    const int n = static_cast<int>(q.size());
    force.assign(n, Point{0.0, 0.0, 0.0});
    HashGrid grid;
    grid.build(q, target);
    const double target2 = target * target;
    ForceState state;
    for (int i = 0; i < n; ++i) {
        const int cx = grid.cell(q[i][0]);
        const int cy = grid.cell(q[i][1]);
        const int cz = grid.cell(q[i][2]);
        for (int dx = -1; dx <= 1; ++dx) {
            const int x = cx + dx;
            if (x < 0 || x >= grid.g) continue;
            for (int dy = -1; dy <= 1; ++dy) {
                const int y = cy + dy;
                if (y < 0 || y >= grid.g) continue;
                for (int dz = -1; dz <= 1; ++dz) {
                    const int z = cz + dz;
                    if (z < 0 || z >= grid.g) continue;
                    for (int j = grid.head[grid.index(x, y, z)]; j != -1; j = grid.next[j]) {
                        if (j <= i) continue;
                        double ax = q[i][0] - q[j][0];
                        double ay = q[i][1] - q[j][1];
                        double az = q[i][2] - q[j][2];
                        double d2 = ax * ax + ay * ay + az * az;
                        if (d2 >= target2) continue;
                        if (d2 < 1e-28) {
                            const uint64_t h = splitmix64((static_cast<uint64_t>(i) << 32) ^ j);
                            ax = ((h & 1023) + 1) / 1024.0 - 0.5;
                            ay = (((h >> 10) & 1023) + 1) / 1024.0 - 0.5;
                            az = (((h >> 20) & 1023) + 1) / 1024.0 - 0.5;
                            d2 = ax * ax + ay * ay + az * az;
                        }
                        const double d = std::sqrt(d2);
                        const double overlap = target - d;
                        const double magnitude = overlap / d;
                        force[i][0] += magnitude * ax;
                        force[i][1] += magnitude * ay;
                        force[i][2] += magnitude * az;
                        force[j][0] -= magnitude * ax;
                        force[j][1] -= magnitude * ay;
                        force[j][2] -= magnitude * az;
                        state.max_overlap = std::max(state.max_overlap, overlap);
                    }
                }
            }
        }
    }
    return state;
}

void fire_minimize(Points& q, double target, int iterations, double deadline) {
    const int n = static_cast<int>(q.size());
    Points velocity(n, Point{0.0, 0.0, 0.0});
    Points force;
    double dt = 0.08;
    double alpha = 0.15;
    int positive_steps = 0;
    ForceState state = overlap_forces(q, target, force);
    for (int iteration = 0; iteration < iterations && elapsed() < deadline; ++iteration) {
        for (int i = 0; i < n; ++i)
            for (int a = 0; a < 3; ++a) velocity[i][a] += dt * force[i][a];

        double power = 0.0;
        double velocity_norm2 = 0.0;
        double force_norm2 = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int a = 0; a < 3; ++a) {
                power += velocity[i][a] * force[i][a];
                velocity_norm2 += velocity[i][a] * velocity[i][a];
                force_norm2 += force[i][a] * force[i][a];
            }
        }
        if (power > 0.0) {
            if (++positive_steps > 5) {
                dt = std::min(0.8, dt * 1.1);
                alpha *= 0.99;
            }
        } else {
            positive_steps = 0;
            dt *= 0.5;
            alpha = 0.15;
            std::fill(velocity.begin(), velocity.end(), Point{0.0, 0.0, 0.0});
        }
        if (velocity_norm2 > 0.0 && force_norm2 > 0.0) {
            const double scale = std::sqrt(velocity_norm2 / force_norm2);
            for (int i = 0; i < n; ++i)
                for (int a = 0; a < 3; ++a)
                    velocity[i][a] = (1.0 - alpha) * velocity[i][a]
                                   + alpha * scale * force[i][a];
        }
        for (int i = 0; i < n; ++i) {
            for (int a = 0; a < 3; ++a) {
                q[i][a] += dt * velocity[i][a];
                if (q[i][a] < 0.0) {
                    q[i][a] = 0.0;
                    if (velocity[i][a] < 0.0) velocity[i][a] = 0.0;
                } else if (q[i][a] > 1.0) {
                    q[i][a] = 1.0;
                    if (velocity[i][a] > 0.0) velocity[i][a] = 0.0;
                }
            }
        }
        state = overlap_forces(q, target, force);
        if (state.max_overlap < 1e-9 * target) break;
    }
}

void repair_fcc_vacancies(Points start, Points& best, double& best_distance,
                          double deadline) {
    double local_best_distance = exact_min_distance(start);
    Points local_best = start;
    rng.seed(0x9e3779b97f4a7c15ULL);
    jitter(start, 0.05 * local_best_distance);
    tighten_box(start);
    double target = exact_min_distance(start) * 1.001;
    int failures = 0;
    int batches = 0;
    const int batch_limit = start.size() <= 128 ? 500 : 452;
    while (batches < batch_limit && elapsed() < deadline) {
        fire_minimize(start, target, 1000, deadline);
        if (elapsed() >= deadline) break;
        tighten_box(start);
        const double distance = exact_min_distance(start);
        if (distance > local_best_distance) {
            local_best_distance = distance;
            local_best = start;
        }
        if (distance >= target * (1.0 - 2e-4)) {
            target *= 1.0009;
            failures = 0;
        } else if (++failures >= 4) {
            target = std::max(distance * 1.001, target * 0.9995);
            failures = 0;
        }
        ++batches;
    }
    if (local_best_distance > best_distance) {
        best_distance = local_best_distance;
        best = std::move(local_best);
    }
    if (std::getenv("SPHERE_DEBUG")) {
        std::cerr << "repair batches=" << batches << " distance=" << local_best_distance
                  << " time=" << elapsed() << '\n';
    }
}

double projection_sweep(Points& q, double distance, std::vector<int>& order,
                        std::vector<int>& rank) {
    HashGrid grid;
    grid.build(q, distance);
    std::shuffle(order.begin(), order.end(), rng);
    for (int pos = 0; pos < static_cast<int>(order.size()); ++pos) rank[order[pos]] = pos;

    const double distance2 = distance * distance;
    double max_overlap = 0.0;
    std::uniform_real_distribution<double> random_component(-1.0, 1.0);
    for (int pos = 0; pos < static_cast<int>(order.size()); ++pos) {
        const int i = order[pos];
        const int cx = grid.cell(q[i][0]);
        const int cy = grid.cell(q[i][1]);
        const int cz = grid.cell(q[i][2]);
        for (int dx = -1; dx <= 1; ++dx) {
            const int x = cx + dx;
            if (x < 0 || x >= grid.g) continue;
            for (int dy = -1; dy <= 1; ++dy) {
                const int y = cy + dy;
                if (y < 0 || y >= grid.g) continue;
                for (int dz = -1; dz <= 1; ++dz) {
                    const int z = cz + dz;
                    if (z < 0 || z >= grid.g) continue;
                    for (int j = grid.head[grid.index(x, y, z)]; j != -1; j = grid.next[j]) {
                        if (rank[j] <= pos) continue;
                        double ax = q[i][0] - q[j][0];
                        double ay = q[i][1] - q[j][1];
                        double az = q[i][2] - q[j][2];
                        double d2 = ax * ax + ay * ay + az * az;
                        if (d2 >= distance2) continue;
                        if (d2 < 1e-24) {
                            ax = random_component(rng);
                            ay = random_component(rng);
                            az = random_component(rng);
                            d2 = ax * ax + ay * ay + az * az;
                        }
                        const double d = std::sqrt(d2);
                        const double overlap = distance - d;
                        max_overlap = std::max(max_overlap, overlap);
                        const double move = 0.44 * overlap / d;
                        for (int a = 0; a < 3; ++a) {
                            const double delta = move * (a == 0 ? ax : (a == 1 ? ay : az));
                            q[i][a] = std::clamp(q[i][a] + delta, 0.0, 1.0);
                            q[j][a] = std::clamp(q[j][a] - delta, 0.0, 1.0);
                        }
                    }
                }
            }
        }
    }
    return max_overlap;
}

void jitter(Points& q, double amplitude) {
    std::uniform_real_distribution<double> delta(-amplitude, amplitude);
    for (Point& p : q) {
        for (double& x : p) x = std::clamp(x + delta(rng), 0.0, 1.0);
    }
}

Point random_in_ball(double radius) {
    std::uniform_real_distribution<double> coord(-1.0, 1.0);
    Point v{};
    double norm2;
    do {
        norm2 = 0.0;
        for (double& x : v) {
            x = coord(rng);
            norm2 += x * x;
        }
    } while (norm2 > 1.0 || norm2 < 1e-20);
    for (double& x : v) x *= radius;
    return v;
}

int random_walk(Points& q, double distance, double magnitude, int attempts,
                double deadline) {
    const double threshold2 = distance * distance * (1.0 - 2e-14);
    int accepted = 0;
    for (int move = 0; move < attempts; ++move) {
        if ((move & 255) == 0 && elapsed() >= deadline) break;
        const int i = static_cast<int>(rng() % q.size());
        const Point delta = random_in_ball(magnitude);
        Point proposal{
            std::clamp(q[i][0] + delta[0], 0.0, 1.0),
            std::clamp(q[i][1] + delta[1], 0.0, 1.0),
            std::clamp(q[i][2] + delta[2], 0.0, 1.0)
        };
        bool valid = true;
        for (int j = 0; j < static_cast<int>(q.size()); ++j) {
            if (j == i) continue;
            const double dx = proposal[0] - q[j][0];
            const double dy = proposal[1] - q[j][1];
            const double dz = proposal[2] - q[j][2];
            if (dx * dx + dy * dy + dz * dz < threshold2) {
                valid = false;
                break;
            }
        }
        if (valid) {
            q[i] = proposal;
            ++accepted;
        }
    }
    return accepted;
}

double stochastic_billiard(Points& q, double magnitude, double minimum_magnitude,
                           int attempts, double deadline) {
    double distance = exact_min_distance(q);
    int stagnant = 0;
    while (magnitude > minimum_magnitude && elapsed() < deadline) {
        const int accepted = random_walk(q, distance, magnitude, attempts, deadline);
        if (elapsed() >= deadline) break;
        const double next_distance = exact_min_distance(q);
        if (next_distance > distance * (1.0 + 2e-12)) {
            distance = next_distance;
            magnitude = std::min(0.25, magnitude * 1.8);
            stagnant = 0;
        } else {
            magnitude *= accepted == 0 ? 0.35 : 0.55;
            ++stagnant;
        }
        if (stagnant > 12) break;
    }
    return exact_min_distance(q);
}

void billiard_relax(Points& best, double& best_distance, double deadline) {
    const int n = static_cast<int>(best.size());
    const int attempt_factor = n <= 128 ? 7 : 3;
    double perturbation = std::min(0.16, best_distance * 0.35);
    int rounds = 0;
    int improvements = 0;
    while (elapsed() < deadline) {
        Points trial = best;
        for (Point& p : trial) {
            const Point delta = random_in_ball(perturbation);
            for (int a = 0; a < 3; ++a) p[a] = std::clamp(p[a] + delta[a], 0.0, 1.0);
        }
        tighten_box(trial);

        const double remaining = deadline - elapsed();
        if (remaining <= 0.0) break;
        const double inner_deadline = elapsed()
                                    + std::min(remaining, std::max(0.012, remaining * 0.40));
        const int attempts = std::max(300, n * attempt_factor);
        const double d = stochastic_billiard(trial, perturbation,
                                             std::max(perturbation / 400.0, 1e-6),
                                             attempts, inner_deadline);
        ++rounds;
        if (d > best_distance * (1.0 + 2e-10)) {
            best_distance = d;
            best = std::move(trial);
            perturbation = std::min(0.18, perturbation * 1.35);
            ++improvements;
        } else {
            perturbation *= 0.68;
            if (perturbation < best_distance * 0.002) {
                perturbation = best_distance * (0.12 + 0.03 * (rounds % 4));
            }
        }
    }
    if (std::getenv("SPHERE_DEBUG")) {
        std::cerr << "billiard rounds=" << rounds << " improvements=" << improvements
                  << " distance=" << best_distance << " time=" << elapsed() << '\n';
    }
}

void relax(Points& best, double& best_distance, double deadline) {
    const int n = static_cast<int>(best.size());
    Points current = best;
    std::vector<int> order(n), rank(n);
    std::iota(order.begin(), order.end(), 0);

    jitter(current, 0.012 * best_distance);
    for (int sweep = 0; sweep < 20 && elapsed() < deadline; ++sweep) {
        projection_sweep(current, best_distance, order, rank);
        if ((sweep & 7) == 7) tighten_box(current);
    }
    if (!has_violation(current, best_distance * (1.0 - 2e-7))) current = best;

    double accepted = best_distance;
    double step = n <= 32 ? 0.010 : (n <= 256 ? 0.006 : (n <= 1024 ? 0.004 : 0.003));
    int failures = 0;
    int successes = 0;
    int attempts = 0;
    while (elapsed() < deadline) {
        ++attempts;
        const double target = accepted * (1.0 + step);
        Points trial = current;
        if (n <= 64 && attempts % 4 == 0) {
            std::uniform_real_distribution<double> unit(0.0, 1.0);
            for (Point& p : trial) for (double& x : p) x = unit(rng);
        } else if (failures > 0) {
            const double shake = (n <= 64 ? 0.10 : (n <= 256 ? 0.045 : 0.018)) * accepted;
            jitter(trial, shake);
        }
        const int sweeps = n <= 64 ? 220 : (n <= 512 ? 70 : 28);
        for (int sweep = 0; sweep < sweeps && elapsed() < deadline; ++sweep) {
            const double overlap = projection_sweep(trial, target, order, rank);
            if ((sweep & 7) == 7) tighten_box(trial);
            if (overlap < target * 2e-7) break;
        }
        if (elapsed() >= deadline) break;
        tighten_box(trial);

        if (!has_violation(trial, target * (1.0 - 3e-6))) {
            current = std::move(trial);
            accepted = target * (1.0 - 3e-6);
            best = current;
            best_distance = accepted;
            failures = 0;
            ++successes;
            step = std::min(step * 1.08, n <= 64 ? 0.015 : 0.007);
        } else {
            ++failures;
            step *= 0.55;
            if (failures >= 3 || step < 2e-5) {
                current = best;
                jitter(current, (0.010 + 0.004 * (failures & 3)) * best_distance);
                for (int sweep = 0; sweep < 12 && elapsed() < deadline; ++sweep) {
                    projection_sweep(current, best_distance, order, rank);
                }
                tighten_box(current);
                if (has_violation(current, best_distance * (1.0 - 2e-6))) current = best;
                failures = 0;
                step = std::max(step, 0.00012);
            }
        }

        if ((successes & 7) == 7 && n <= 768 && elapsed() + 0.015 < deadline) {
            const double exact = exact_min_distance(best);
            if (exact > best_distance) best_distance = exact;
            accepted = best_distance;
        }
    }
    if (std::getenv("SPHERE_DEBUG")) {
        std::cerr << "relax attempts=" << attempts << " successes=" << successes
                  << " distance=" << best_distance << " time=" << elapsed() << '\n';
    }
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    start_time = Clock::now();

    int n;
    if (!(std::cin >> n) || n < 2) return 0;
    rng.seed(15);

    constexpr int pool_limit = 21;
    std::vector<Candidate> pool;
    Points fcc_repair_start;
    Points cubic = simple_cubic(n);
    const int s = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(n)) - 1e-12));
    retain_candidate(pool, cubic, 1.0 / std::max(1, s - 1), pool_limit);
    add_corner_candidate(n, pool, pool_limit);
    Points rectangular = rectangular_fcc(n);
    const double rectangular_distance = exact_min_distance(rectangular);
    retain_candidate(pool, std::move(rectangular), rectangular_distance, pool_limit);
    Points tetrahedral = tetrahedral_cluster(n);
    const double tetrahedral_distance = exact_min_distance(tetrahedral);
    retain_candidate(pool, std::move(tetrahedral), tetrahedral_distance, pool_limit);
    Points diagonal = diagonal_fcc_cluster(n);
    const double diagonal_distance = exact_min_distance(diagonal);
    retain_candidate(pool, std::move(diagonal), diagonal_distance, pool_limit);
    Points rows = staggered_rows(n);
    const double rows_distance = exact_min_distance(rows);
    retain_candidate(pool, std::move(rows), rows_distance, pool_limit);
    if (n <= 768) {
        int lower_side = 2;
        auto lower_count = [](int side) { return (side * side * side + 1) / 2; };
        while (lower_count(lower_side + 1) <= n) ++lower_side;
        const int base_count = lower_count(lower_side);
        const double excess_ratio = double(n - base_count) / base_count;
        if (n >= 40 && n > base_count && excess_ratio <= 0.035) {
            fcc_repair_start = lower_fcc_repair_pattern(n);
        }
    }

    const int k = static_cast<int>(std::ceil(1.10 * std::cbrt(static_cast<double>(n)))) + 7;
    const int samples = n <= 128 ? 48 : 24;

    Points raw = generate_sc_patch(k);
    for (int sample = 0; sample < std::max(8, samples / 2); ++sample) {
        crop_candidate(raw, n, 1.0, sample, pool, pool_limit);
    }

    raw = generate_fcc_patch(k);
    for (int sample = 0; sample < samples; ++sample) {
        crop_candidate(raw, n, std::sqrt(2.0), sample, pool, pool_limit);
    }

    raw = generate_barlow_patch(k, false);
    for (int sample = 0; sample < samples; ++sample) {
        crop_candidate(raw, n, 1.0, sample, pool, pool_limit);
    }

    raw = generate_barlow_patch(k, true);
    for (int sample = 0; sample < std::max(8, samples / 2); ++sample) {
        crop_candidate(raw, n, 1.0, sample, pool, pool_limit);
    }

    const double aspects[] = {0.78, 0.90, 1.00, 1.16, 1.4142135623730951, 1.58};
    for (double aspect : aspects) {
        raw = generate_bct_patch(k, aspect);
        const double d0 = std::min({1.0, aspect, std::sqrt(0.5 + 0.25 * aspect * aspect)});
        for (int sample = 0; sample < 12; ++sample) {
            crop_candidate(raw, n, d0, sample, pool, pool_limit);
        }
        if (elapsed() > 0.30) break;
    }

    if (std::getenv("SPHERE_DEBUG")) {
        std::cerr << "constructed pool=" << pool.size() << " time=" << elapsed() << '\n';
    }

    Points best = cubic;
    double best_distance = exact_min_distance(best);
    for (Candidate& candidate : pool) {
        if (elapsed() > 0.42 && n > 512) break;
        const double d = exact_min_distance(candidate.q);
        candidate.estimate = d;
        if (std::getenv("SPHERE_DEBUG")) {
            std::cerr << "candidate estimate=" << candidate.estimate << " exact=" << d << '\n';
        }
        if (d > best_distance) {
            best_distance = d;
            best = candidate.q;
        }
    }

    if (std::getenv("SPHERE_DEBUG")) {
        std::cerr << "selected distance=" << best_distance << " time=" << elapsed() << '\n';
    }

    std::sort(pool.begin(), pool.end(), [](const Candidate& a, const Candidate& b) {
        return a.estimate > b.estimate;
    });
    if (n <= 64) {
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        Points random_start(n);
        for (Point& p : random_start) for (double& x : p) x = unit(rng);
        tighten_box(random_start);
        pool.push_back(Candidate{std::move(random_start), 0.0});
        pool.back().estimate = exact_min_distance(pool.back().q);
        std::sort(pool.begin(), pool.end(), [](const Candidate& a, const Candidate& b) {
            return a.estimate > b.estimate;
        });
    }

    constexpr double relaxation_deadline = 0.84;
    if (best_distance >= std::sqrt(3.0) - 1e-12) {
        // The cube diagonal is an absolute upper bound for two centers.
    } else if (!fcc_repair_start.empty()) {
        repair_fcc_vacancies(std::move(fcc_repair_start), best, best_distance,
                             relaxation_deadline);
    } else if (n <= 512) {
        if (n <= 32 && pool.size() > 3) {
            Points alternate = pool[3].q;
            double alternate_distance = exact_min_distance(alternate);
            billiard_relax(alternate, alternate_distance, 0.18);
            if (alternate_distance > best_distance) {
                best_distance = alternate_distance;
                best = std::move(alternate);
            }
        }
        billiard_relax(best, best_distance, relaxation_deadline);
    } else {
        relax(best, best_distance, relaxation_deadline);
    }

    tighten_box(best);
    best_distance = exact_min_distance(best);
    const double radius = best_distance / (2.0 * (1.0 + best_distance));
    const double scale = 1.0 / (1.0 + best_distance);

    std::cout << std::setprecision(17);
    for (const Point& q : best) {
        std::cout << radius + scale * q[0] << ' '
                  << radius + scale * q[1] << ' '
                  << radius + scale * q[2] << '\n';
    }
    return 0;
}
