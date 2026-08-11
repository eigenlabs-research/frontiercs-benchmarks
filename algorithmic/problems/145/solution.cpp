#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

// Number Loop Construction (problem 145).
//
// Construction is loop-first (Lean NumberLoop / Sable notes):
//   - Embed a single closed loop as a set of grid-edge indices (0..311),
//     matching the official checker edge order (156 horizontal, then 156 vertical).
//   - Induce each template clue as the number of selected edges around that cell.
//   - Digits are therefore consistent with a witness loop by construction.
//
// The embedded loops are the best locally verified UNIQUE constructions found
// so far under face-flip / rectangle-flip search:
//   mode 0: 46 ones (zeros allowed), unique
//   mode 1: 42 ones (digits in 1..3), unique
// Local score estimate: average of (50 + 50 * ones / 56) / 100 ≈ 0.8929.
//
// Neighborhood search finds 43–48-one loops nearby, but every improving
// neighbor examined so far has k>=2 solutions. Further gains need a move set
// or certification method that escapes this uniqueness plateau.

static const char *T[12] = {
    "?   ?   ??? ", "?? ??  ?   ?", "? ? ?  ?   ?", "? ? ?  ???? ",
    "? ? ?  ?    ", "?   ?  ?    ", "            ", "?  ?   ?????",
    "? ?      ?  ", "??   ? ? ?  ", "? ?  ? ? ?  ", "?  ? ??? ?  "};

// Edge lists recovered from the unique seed puzzles.
static const int LOOP0[] = {
    0,   1,   2,   3,   9,   10,  11,  14,  15,  17,  30,  31,  32,  41,  42,  52,
    54,  56,  58,  63,  65,  69,  74,  80,  84,  85,  86,  95,  97,  98,  113, 124,
    127, 136, 138, 145, 146, 147, 148, 152, 156, 160, 165, 168, 169, 171, 174, 175,
    178, 181, 182, 184, 187, 194, 195, 197, 202, 207, 208, 210, 212, 213, 214, 216,
    217, 218, 219, 220, 221, 223, 224, 229, 232, 233, 234, 243, 245, 246, 250, 256,
    261, 269, 274, 278, 279, 282, 287, 290, 292, 293, 294, 295, 300, 304, 307, 308};
static const int LOOP0_N = 96;

static const int LOOP1[] = {
    1,   6,   7,   8,   9,   10,  11,  15,  16,  17,  19,  27,  29,  30,  40,  53,
    56,  57,  58,  63,  64,  67,  68,  69,  70,  74,  77,  79,  81,  82,  85,  88,
    90,  91,  93,  95,  96,  97,  100, 101, 103, 104, 105, 108, 109, 110, 113, 115,
    116, 121, 122, 124, 127, 136, 138, 145, 146, 147, 148, 152, 157, 158, 162, 168,
    170, 171, 172, 176, 177, 181, 183, 184, 186, 187, 190, 194, 196, 197, 203, 207,
    209, 210, 213, 214, 219, 220, 222, 223, 224, 227, 228, 233, 235, 239, 242, 243,
    245, 246, 249, 251, 253, 257, 260, 267, 276, 278, 279, 282, 287, 290, 292, 293,
    294, 295, 300, 304, 307, 308};
static const int LOOP1_N = 118;

static int Hidx[13][12], Vidx[12][13];

static void init_index() {
    int id = 0;
    for (int r = 0; r <= 12; r++)
        for (int c = 0; c < 12; c++) Hidx[r][c] = id++;
    for (int r = 0; r < 12; r++)
        for (int c = 0; c <= 12; c++) Vidx[r][c] = id++;
}

static void induce_and_print(const int *loop, int n) {
    bool sel[312];
    memset(sel, 0, sizeof(sel));
    for (int i = 0; i < n; i++) sel[loop[i]] = true;

    for (int r = 0; r < 12; r++) {
        string row(12, ' ');
        for (int c = 0; c < 12; c++) {
            if (T[r][c] != '?') continue;
            int d = (sel[Hidx[r][c]] ? 1 : 0) + (sel[Hidx[r + 1][c]] ? 1 : 0) +
                    (sel[Vidx[r][c]] ? 1 : 0) + (sel[Vidx[r][c + 1]] ? 1 : 0);
            row[c] = char('0' + d);
        }
        cout << row << '\n';
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int mode = 0;
    if (!(cin >> mode)) mode = 0;
    init_index();
    if (mode == 1)
        induce_and_print(LOOP1, LOOP1_N);
    else
        induce_and_print(LOOP0, LOOP0_N);
    return 0;
}
