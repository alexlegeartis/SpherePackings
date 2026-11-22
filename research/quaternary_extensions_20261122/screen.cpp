#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

using Word = std::uint32_t;
using Syndromes = std::uint64_t;
constexpr Syndromes all_syndromes = std::numeric_limits<Syndromes>::max();

struct Residue {
    Word low_bits;
    unsigned high_syndrome;
};

// For a mask of possible syndromes, shift every syndrome by XOR with value.
Syndromes shift_syndromes(Syndromes mask, unsigned value) {
    constexpr std::array<Syndromes, 6> low_masks{
        0x5555555555555555ULL, 0x3333333333333333ULL,
        0x0f0f0f0f0f0f0f0fULL, 0x00ff00ff00ff00ffULL,
        0x0000ffff0000ffffULL, 0x00000000ffffffffULL,
    };
    while (value) {
        const unsigned bit = std::countr_zero(value);
        value &= value - 1;
        const auto low = low_masks[bit];
        const unsigned distance = 1U << bit;
        mask = ((mask & low) << distance) | ((mask & ~low) >> distance);
    }
    return mask;
}

struct Frame {
    Word active;
    Word double_weight;
    std::array<unsigned, 32> columns{};
    std::array<unsigned, 1U << 16> lower_syndrome{};
    std::array<unsigned, 1U << 16> upper_syndrome{};

    Frame(Word active, Word double_weight) : active(active), double_weight(double_weight) {
        for (unsigned j = 0; j < columns.size(); ++j) {
            if ((active >> j) & 1U) {
                columns[j] = (j == 2 && double_weight) ? 2 : 1 | (j << 1);
            }
        }
        for (Word bits = 1; bits < lower_syndrome.size(); ++bits) {
            const unsigned j = std::countr_zero(bits);
            const Word rest = bits & (bits - 1);
            lower_syndrome[bits] = lower_syndrome[rest] ^ columns[j];
            upper_syndrome[bits] = upper_syndrome[rest] ^ columns[j + 16];
        }
    }

    unsigned syndrome(Word bits) const {
        return lower_syndrome[bits & 65535] ^ upper_syndrome[bits >> 16];
    }

    int weight(Word bits) const {
        return std::popcount(bits) + std::popcount(bits & double_weight);
    }

    // Squared costs are in units of 1/4 of the pre-scaling squared length.
    // A set bit in offset adds 1/2 to the corresponding pre-scaling coordinate.
    Syndromes allowed(Word odd, Word offset, int minimum) const {
        const int base = 4 * weight(odd & ~offset) + weight(offset);
        if (base >= minimum) {
            return all_syndromes;
        }
        // An odd unshifted coordinate can change its high bit for free.
        Syndromes free = 1;
        for (Word bits = odd & ~offset; bits; bits &= bits - 1) {
            free |= shift_syndromes(free, columns[std::countr_zero(bits)]);
            if (free == all_syndromes) {
                return 0;
            }
        }
        // Each remaining high-bit change costs 8 or 16, times the axis weight.
        const int budget = (minimum - 1 - base) / 8;
        std::array<Syndromes, 8> reachable{};
        std::fill_n(reachable.begin(), budget + 1, free);
        for (Word bits = active & (~odd | offset); bits; bits &= bits - 1) {
            const Word bit = bits & -bits;
            const int cost = ((bit & offset) ? 1 : 2) * weight(bit);
            const unsigned column = columns[std::countr_zero(bits)];
            for (int i = budget; i >= cost; --i) {
                reachable[i] |= shift_syndromes(reachable[i - cost], column);
            }
            if (reachable[budget] == all_syndromes) {
                return 0;
            }
        }
        // The cheaper high bit is 1 exactly at the odd shifted positions.
        return shift_syndromes(~reachable[budget], syndrome(odd & offset));
    }
};

Syndromes translate_syndromes(const Frame& frame, const std::vector<Residue>& residues,
                              Word low, Word offset, int minimum) {
    Syndromes choices = all_syndromes;
    for (const auto& residue : residues) {
        const auto allowed = frame.allowed(low ^ residue.low_bits, offset, minimum);
        const auto carry_syndrome = frame.syndrome(low & residue.low_bits);
        choices &= shift_syndromes(allowed, residue.high_syndrome ^ carry_syndrome);
        if (!choices) {
            break;
        }
    }
    return choices;
}

void print_pairs(const Frame& frame, const std::unordered_map<Word, Syndromes>& groups) {
    std::vector<std::pair<Word, Syndromes>> ordered(groups.begin(), groups.end());
    std::sort(ordered.begin(), ordered.end());
    const auto even = groups.find(0);
    std::uint64_t edges = 0, index_four_pairs = 0, cyclic = 0;
    std::array<unsigned, 4> first{}, first_index_four{};
    std::array<unsigned, 2> first_cyclic{};
    for (std::size_t i = 0; i < ordered.size(); ++i) {
        const auto [x, left] = ordered[i];
        if (even != groups.end() && ((even->second >> frame.syndrome(x)) & 1)) {
            if (!cyclic) {
                first_cyclic = {x, static_cast<unsigned>(std::countr_zero(left))};
            }
            cyclic += std::popcount(left);
        }
        for (std::size_t j = i; j < ordered.size(); ++j) {
            const auto [y, right] = ordered[j];
            const auto difference = groups.find(x ^ y);
            if (difference == groups.end()) {
                continue;
            }
            const unsigned borrow = frame.syndrome(y & ~x);
            for (auto bits = left; bits; bits &= bits - 1) {
                const unsigned h = std::countr_zero(bits);
                auto neighbors = right & shift_syndromes(difference->second, h ^ borrow);
                if (i == j) {
                    neighbors &= all_syndromes << h;
                }
                if (!neighbors) {
                    continue;
                }
                const std::array<unsigned, 4> pair{x, h, y,
                    static_cast<unsigned>(std::countr_zero(neighbors))};
                const auto count = std::popcount(neighbors);
                if (!edges) {
                    first = pair;
                }
                edges += count;
                if (!frame.syndrome(x) && !frame.syndrome(y)) {
                    if (!index_four_pairs) {
                        first_index_four = pair;
                    }
                    index_four_pairs += count;
                }
            }
        }
    }
    std::cout << edges << ' ' << index_four_pairs << ' ' << cyclic;
    for (auto witness : {first, first_index_four}) {
        for (auto value : witness) {
            std::cout << ' ' << value;
        }
    }
    for (auto value : first_cyclic) {
        std::cout << ' ' << value;
    }
    std::cout << '\n';
}

int main(int argc, char** argv) {
    Word active, double_weight, offset;
    int minimum, generators;
    if (!(std::cin >> active >> double_weight >> offset >> minimum >> generators)
        || (offset & ~active) || (double_weight != 0 && double_weight != 4)
        || minimum < 1 || minimum > 64 || generators < 1 || generators > 10) {
        throw std::runtime_error("invalid frame or search parameters");
    }
    const bool scan_pairs = argc == 2 && std::string_view(argv[1]) == "--pairs";
    const bool stack = argc == 2 && std::string_view(argv[1]) == "--stack";
    if ((argc != 1 && !scan_pairs && !stack) || (scan_pairs && offset)) {
        throw std::runtime_error("invalid search mode");
    }
    int second_minimum = 0;
    if (stack && (!(std::cin >> second_minimum) || second_minimum < 1 || second_minimum > 64)) {
        throw std::runtime_error("invalid second-layer minimum");
    }
    const Frame frame(active, double_weight);
    std::vector<Residue> residues{{0, 0}};
    for (int i = 0; i < generators; ++i) {
        Word low, high;
        if (!(std::cin >> low >> high) || ((low | high) & ~active) || frame.syndrome(low)) {
            throw std::runtime_error("invalid lifted generator");
        }
        const unsigned size = residues.size();
        for (unsigned j = 0; j < size; ++j) {
            const auto [bits, syndrome] = residues[j];
            residues.push_back({bits ^ low, syndrome ^ frame.syndrome(high ^ (bits & low))});
        }
    }

    // Put one binary residue in each coset by clearing a pivot set.
    std::array<Word, 32> basis{};
    Word pivots = 0;
    for (const auto& residue : residues) {
        Word bits = residue.low_bits;
        while (bits) {
            const unsigned j = std::countr_zero(bits);
            if (!basis[j]) {
                basis[j] = bits;
                pivots |= 1U << j;
                break;
            }
            bits ^= basis[j];
        }
    }
    if (std::popcount(pivots) != generators) {
        throw std::runtime_error("dependent binary residues");
    }

    Word free_positions = active & ~pivots;
    const unsigned count = 1U << std::popcount(free_positions);
    std::uint64_t accepted = 0, index_two = 0;
    std::unordered_map<Word, Syndromes> groups;
    std::vector<Residue> stack_candidates;
    constexpr std::size_t stack_candidate_limit = 65536;
    Word first = 0, first_index_two = 0;
    unsigned first_syndrome = 0, first_index_two_syndrome = 0;
    const auto doubled_residue = std::find_if(residues.begin(), residues.end(),
        [offset](const Residue& residue) { return residue.low_bits == offset; });
    const auto stack_doubles = stack ? translate_syndromes(frame, residues, offset, 0, second_minimum)
                                     : all_syndromes;
    // Compare against each residue. High-bit syndrome represents the remaining
    // cosets modulo the binary torsion code. Adding a seed word is sufficient
    // for distance checks because the seed is an additive group.
    Word low = 0;
    do {
        const auto choices = ((stack_doubles >> frame.syndrome(low)) & 1)
            ? translate_syndromes(frame, residues, low, offset, minimum) : 0;
        if (choices) {
            if (scan_pairs) {
                groups.emplace(low, choices);
            }
            if (stack && !offset) {
                for (auto bits = choices; bits && stack_candidates.size() < stack_candidate_limit;
                     bits &= bits - 1) {
                    stack_candidates.push_back({low, static_cast<unsigned>(std::countr_zero(bits))});
                }
            }
            if (!accepted) {
                first = low;
                first_syndrome = std::countr_zero(choices);
            }
            const auto size = std::popcount(choices);
            accepted += size;
            // Twice the translate has low bits `offset` and high bits `low`.
            if (doubled_residue != residues.end()
                && frame.syndrome(low) == doubled_residue->high_syndrome) {
                if (!index_two) {
                    first_index_two = low;
                    first_index_two_syndrome = std::countr_zero(choices);
                }
                index_two += size;
            }
        }
        low = (low - free_positions) & free_positions;
    } while (low);

    std::cout << count << ' ' << accepted << ' ' << index_two << ' '
              << first << ' ' << first_syndrome << ' '
              << first_index_two << ' ' << first_index_two_syndrome << '\n';
    if (scan_pairs) {
        print_pairs(frame, groups);
    }
    for (const auto& candidate : stack_candidates) {
        std::cout << candidate.low_bits << ' ' << candidate.high_syndrome << '\n';
    }
}
