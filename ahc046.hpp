#pragma once

#include <array>
#include <bitset>
#include <vector>

/* ★ 問題は N = 20,  M = 40 固定 ★ */
constexpr int N = 20;
constexpr int M = 40;

/* 座標 → 1 次元 id への簡易射影 */
using Point = int;
inline Point f(int x, int y) { return x * N + y; }
inline bool isin(int x, int y) { return 0 <= x and x < N and 0 <= y and y < N; }

const std::array<int, 4> dx{1, -1, 0, 0};
const std::array<int, 4> dy{0, 0, 1, -1};

using AlterPlan = std::vector<Point>;
using AlterPlans = std::array<AlterPlan, M>;

using BanStates = std::bitset<N * N>;
