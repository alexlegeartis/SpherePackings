#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

constexpr int DIM = 24;
constexpr int WIDTH = 4096;
constexpr int FRAME_NORM = 8;
using Row = std::array<int, DIM>;
// Two bit planes hold the low and high bits of 24 coordinates modulo four.
using Word = std::pair<std::uint32_t, std::uint32_t>;

Word add(Word a, Word b) {
    return {a.first ^ b.first, a.second ^ b.second ^ (a.first & b.first)};
}

Word word(const Row& row) {
    Word w{0, 0};
    for (int i = 0; i < DIM; ++i) {
        int digit = static_cast<unsigned>(row[i]) & 3U;
        w.first |= std::uint32_t(digit & 1) << i;
        w.second |= std::uint32_t(digit >> 1) << i;
    }
    return w;
}

struct Geometry {
    std::array<Row, DIM> frame{}, forward{};
    std::vector<Row> image;
    std::array<int, WIDTH> common_lo{}, common_hi{}, sum_lo{}, sum_hi{};

    explicit Geometry(const char* filename) {
        std::ifstream f(filename);
        for (auto* matrix : {&frame, &forward})
            for (auto& row : *matrix)
                for (int& x : row)
                    if (!(f >> x) || x < -1000 || x > 1000)
                        throw std::runtime_error("invalid frame matrix");
        std::array<Row, 2> masks{};
        for (auto& row : masks)
            for (int& x : row)
                if (!(f >> x) || x < 0 || x >= WIDTH)
                    throw std::runtime_error("invalid projection mask");
        for (int x = 0; x < WIDTH; ++x)
            for (int i = 0; i < 12; ++i)
                if ((x >> i) & 1) {
                    common_lo[x] ^= masks[0][i];
                    common_hi[x] ^= masks[0][i + 12];
                    sum_lo[x] ^= masks[1][i];
                    sum_hi[x] ^= masks[1][i + 12];
                }
        int count;
        if (!(f >> count) || count < 1 || count > 65536)
            throw std::runtime_error("invalid residue count");
        image.resize(count);
        for (auto& row : image)
            for (int& x : row)
                if (!(f >> x) || x < 0 || x > 3)
                    throw std::runtime_error("invalid frame residue");
    }
};

struct Decoder {
    const Geometry& g;
    const Row& anchor;
    int bound;
    int minimum = 100;
    std::vector<std::uint8_t> cost;
    bool save_costs;
    std::uint64_t nodes = 0, leaves = 0;
    std::array<std::vector<std::pair<int, Word>>, DIM> options;
    std::array<int, DIM+1> suffix{};

    Decoder(const Geometry& g, const Row& anchor, int bound, bool save_costs)
        : g(g), anchor(anchor), bound(bound), cost(save_costs ? WIDTH * WIDTH : 0, 100),
          save_costs(save_costs) {}

    void visit(int pos, int space, Word state) {
        ++nodes;
        if (pos == DIM) {
            // The tracked coordinates are y - anchor. The high bits give
            // (y - anchor) / 2 modulo two when the low bits vanish.
            if (state.first) return;
            const int squared = bound*FRAME_NORM - space;
            if (squared % FRAME_NORM) throw std::runtime_error("nonintegral Leech norm");
            const int norm = squared / FRAME_NORM;
            ++leaves;
            minimum = std::min(minimum, norm);
            if (save_costs) {
                int low = state.second & (WIDTH-1), high = state.second >> 12;
                int n = g.common_lo[low] ^ g.common_hi[high];
                int m = g.sum_lo[low] ^ g.sum_hi[high];
                auto& value = cost[n * WIDTH + m];
                value = std::min(value, std::uint8_t(norm));
            }
            return;
        }
        for (auto [value, increment] : options[pos]) {
            if (value + suffix[pos+1] > space) break;
            visit(pos+1, space-value, add(state, increment));
        }
    }

    void walk() {
        Row projected{};
        for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j) projected[i] += anchor[j] * g.forward[j][i];

        std::array<std::array<Word, 4>, DIM> increments{};
        for (int i = 0; i < DIM; ++i)
            for (int k = 1; k < 4; ++k)
                increments[i][k] = add(increments[i][k-1], word(g.frame[i]));

        for (const auto& image : g.image) {
            Row base{};
            suffix[DIM] = 0;
            for (int i = DIM-1; i >= 0; --i) {
                const int r = static_cast<unsigned>(projected[i] + 2*image[i]) & 7U;
                base[i] = r > 4 ? r-8 : r;
                suffix[i] = suffix[i+1] + base[i]*base[i];
            }
            if (suffix[0] > FRAME_NORM*bound) continue;
            Row coordinates{};
            // Convert the closest representative back to lattice coordinates.
            // Moving a frame coordinate by 8 adds the corresponding frame row.
            for (int i = 0; i < DIM; ++i) {
                for (int j = 0; j < DIM; ++j)
                    coordinates[i] += (base[j]-projected[j])*g.frame[j][i];
                if (coordinates[i] % FRAME_NORM)
                    throw std::runtime_error("frame residue outside the Leech half-coset image");
                coordinates[i] /= FRAME_NORM;
            }
            for (int i = 0; i < DIM; ++i) {
                options[i].clear();
                // The largest parser budget is 60, with 8 times that in frame coordinates.
                // Therefore |base + 8k| <= sqrt(480) < 22 and |k| <= 3.
                for (int k = -3; k <= 3; ++k) {
                    int v = base[i]+8*k;
                    if (v*v <= bound*FRAME_NORM)
                        options[i].emplace_back(v*v, increments[i][static_cast<unsigned>(k)&3U]);
                }
                std::sort(options[i].begin(), options[i].end());
            }
            visit(0, FRAME_NORM*bound, word(coordinates));
        }
    }
};

void decode(const Geometry& g, const std::array<Row, 3>& anchors, int budget) {
    std::array<int, 3> minimum{};
    int sum = 0;
    for (int i = 0; i < 3; ++i) {
        Decoder local(g, anchors[i], 8, false);
        local.walk();
        minimum[i] = std::min(9, local.minimum);
        sum += minimum[i];
    }
    if (sum > budget) {
        std::cout << "EXCLUDED " << budget+1 << " 0 0:0 0:0 0:0" << std::endl;
        return;
    }
    std::array<Decoder, 3> blocks{
        Decoder(g, anchors[0], budget - sum + minimum[0], true),
        Decoder(g, anchors[1], budget - sum + minimum[1], true),
        Decoder(g, anchors[2], budget - sum + minimum[2], true)
    };
    for (auto& block : blocks) block.walk();
    unsigned long long pairs = 0;
    int best = budget + 1;
    for (int common = 0; common < WIDTH; ++common) {
        std::array<std::vector<std::pair<int, int>>, 2> terms;
        int start = common * WIDTH;
        for (int i = 0; i < 2; ++i) {
            for (int m = 0; m < WIDTH; ++m) {
                int norm = blocks[i].cost[start + m];
                if (norm <= blocks[i].bound) terms[i].emplace_back(norm, m);
            }
            std::sort(terms[i].begin(), terms[i].end());
        }
        for (auto [a_cost, a] : terms[0])
            for (auto [b_cost, b] : terms[1]) {
                if (a_cost + b_cost + minimum[2] >= best) break;
                ++pairs;
                // All three first labels agree and the second labels sum to zero.
                best = std::min(best,
                    a_cost + b_cost + blocks[2].cost[start + (a ^ b)]);
            }
    }
    std::cout << (best <= budget ? "FOUND" : "EXCLUDED") << ' ' << best
              << ' ' << pairs;
    for (auto& block : blocks) std::cout << ' ' << block.nodes << ':' << block.leaves;
    std::cout << std::endl;
}
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const Geometry geometry(argv[1]);
    std::ifstream f(argv[2]);
    int budget;
    while (f >> budget) {
        if (budget < 0 || budget > 60) return 2;
        std::array<Row, 3> anchors{};
        for (auto& a : anchors)
            for (auto& v : a) {
                f >> v;
                if (!f || v < 0 || v > 3) return 2;
            }
        decode(geometry, anchors, budget);
    }
    return 0;
}
