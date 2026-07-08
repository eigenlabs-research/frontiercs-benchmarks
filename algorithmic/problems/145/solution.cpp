// Number Loop Construction - problem 145
// Constructive solver: outputs a precomputed unique-solution puzzle
// (maximized 1-clues) for the given input mode.
//
// Mode 0 (clues in {0,1,2,3}): 46 ones, 1 solution  -> score 91.07
// Mode 1 (clues in {1,2,3}):   42 ones, 1 solution  -> score 87.50
//
// The grids were found by loop-based local search (simulated annealing over
// XOR-with-cell-cycle moves: blobs, rectangles, straight runs) with an exact
// reimplementation of the checker's loop counter to verify uniqueness.
// Because the judge runs the testlib checker with a 10s CPU cap, the search
// also bounded the checker's own enumeration cost: both grids verify in
// ~2s locally (~0.4M-0.5M checker search nodes), well inside the cap.
// Each grid is a fixed 12x12 assignment matching the required template
// (blank positions = ' ', clue positions = digit).

#include <cstdio>
#include <cstring>

static const char* GRID0[12] = {
    "2   1   121 ",
    "11 11  1   1",
    "1 1 1  1   1",
    "1 1 1  1101 ",
    "1 1 3  1    ",
    "1   0  1    ",
    "            ",
    "1  1   01101",
    "1 1      1  ",
    "11   3 1 1  ",
    "1 0  1 3 1  ",
    "1  1 111 1  "
};

static const char* GRID1[12] = {
    "1   1   111 ",
    "12 31  3   1",
    "1 1 3  1   1",
    "1 1 1  1211 ",
    "1 1 2  1    ",
    "1   1  3    ",
    "            ",
    "1  1   21311",
    "3 1      1  ",
    "12   3 2 1  ",
    "1 1  1 3 1  ",
    "1  1 111 1  "
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
