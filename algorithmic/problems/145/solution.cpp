// Number Loop Construction - problem 145
// Constructive solver: outputs a precomputed unique-solution puzzle
// (maximized 1-clues) for the given input mode.
//
// Mode 0 (clues in {0,1,2,3}): 35 ones, 1 solution  -> score 81.25
// Mode 1 (clues in {1,2,3}):   31 ones, 1 solution  -> score 77.68
//
// The grids were found by loop-based local search (unbump/pinch/2-opt/slide
// moves) with an embedded port of the checker's loop counter to verify
// uniqueness. Each grid is a fixed 12x12 assignment matching the required
// template (blank positions = ' ', clue positions = digit).

#include <cstdio>
#include <cstring>

static const char* GRID0[12] = {
    "1   1   111 ",
    "11 31  3   1",
    "1 1 3  3   1",
    "1 2 1  3133 ",
    "1 1 2  3    ",
    "1   3  3    ",
    "            ",
    "1  3   11101",
    "2 0      1  ",
    "21   1 0 1  ",
    "1 3  3 1 1  ",
    "1  1 111 0  "
};

static const char* GRID1[12] = {
    "1   3   111 ",
    "11 21  3   1",
    "1 1 1  3   1",
    "1 1 3  3133 ",
    "1 1 3  3    ",
    "2   3  3    ",
    "            ",
    "3  1   31311",
    "1 1      2  ",
    "13   1 3 1  ",
    "1 3  3 1 2  ",
    "3  1 311 3  "
};

int main() {
    int w = 0;
    if (scanf("%d", &w) != 1) w = 0;
    const char* const* g = (w == 1) ? GRID1 : GRID0;
    for (int i = 0; i < 12; ++i) {
        // Each row is exactly 12 chars; preserve trailing spaces.
        fwrite(g[i], 1, 12, stdout);
        fputc('\n', stdout);
    }
    return 0;
}
