#include "common.hpp"

using namespace std;
using lint = long long;
using pint = std::pair<int, int>;
using plint = std::pair<lint, lint>;

struct fast_ios {
    fast_ios() {
        std::cin.tie(nullptr), std::ios::sync_with_stdio(false), std::cout << std::fixed << std::setprecision(20);
    };
} fast_ios_;

int main(int argc, char *argv[]) {
    int X = 0;
    if (argc >= 2) { X = std::stoi(argv[1]); }

    int x;
    std::cin >> x;
    jdump("stdin", x);

    dump_onlinejudge("solution");

    jdump("score", -(X - 6) * (X - 6));
}
