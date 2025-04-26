#pragma once
#include "ahc046.hpp"

#include <queue>

/* 方向テーブル */
const int DX[4] = {-1, 1, 0, 0};
const int DY[4] = {0, 0, -1, 1};
const char DIRC[4] = {'U', 'D', 'L', 'R'};

/* -------------------------------------------------------------------------- */
/* 1 区間 (start → goal,  AlterPlan = plan) を最短で組み立てる BFS            */
/* -------------------------------------------------------------------------- */
namespace detail {

struct Prev {
    int prev;
    char act, dir;
};

std::vector<std::pair<char, char>>
bfs_segment(Point start, Point goal, const AlterPlan &plan, const std::array<std::array<uint8_t, N>, N> &board0) {
    using namespace std;
    const int sx = start / N, sy = start % N;
    const int gx = goal / N, gy = goal % N;
    const int K = plan.size();

    /* ---- 盤面スナップショット (0 … K 回目の Alter 終了時) ---- */
    vector<array<array<uint8_t, N>, N>> snap(K + 1);
    snap[0] = board0;
    for (int i = 0; i < K; ++i) {
        snap[i + 1] = snap[i];
        int x = plan[i] / N, y = plan[i] % N;
        snap[i + 1][x][y] ^= 1; // toggle
    }

    const int SZ = (K + 1) * N * N; // 全状態数
    auto idx = [&](int x, int y, int s) { return s * N * N + x * N + y; };

    vector<char> seen(SZ, 0);
    vector<Prev> prv(SZ);
    queue<tuple<int, int, int>> que;

    int sidx = idx(sx, sy, 0);
    seen[sidx] = 1;
    que.emplace(sx, sy, 0);

    int gidx = -1;

    while (!que.empty()) {
        auto [x, y, s] = que.front();
        que.pop();

        if (x == gx && y == gy && s == K) { // 目的地＋Alter 完了
            gidx = idx(x, y, s);
            break;
        }

        /* ---- Move ---- */
        for (int d = 0; d < 4; ++d) {
            int nx = x + DX[d], ny = y + DY[d];
            if (!isin(nx, ny) || snap[s][nx][ny]) continue;
            int ni = idx(nx, ny, s);
            if (!seen[ni]) {
                seen[ni] = 1;
                prv[ni] = {idx(x, y, s), 'M', DIRC[d]};
                que.emplace(nx, ny, s);
            }
        }

        /* ---- Slide ---- */
        for (int d = 0; d < 4; ++d) {
            int nx = x, ny = y;
            while (true) {
                int tx = nx + DX[d], ty = ny + DY[d];
                if (!isin(tx, ty) || snap[s][tx][ty]) break;
                nx = tx;
                ny = ty;
            }
            if (nx == x && ny == y) continue; // 動けなかった
            int ni = idx(nx, ny, s);
            if (!seen[ni]) {
                seen[ni] = 1;
                prv[ni] = {idx(x, y, s), 'S', DIRC[d]};
                que.emplace(nx, ny, s);
            }
        }

        /* ---- Alter (順番厳守) ---- */
        if (s < K) {
            int ax = plan[s] / N, ay = plan[s] % N;
            for (int d = 0; d < 4; ++d) {
                if (x + DX[d] == ax && y + DY[d] == ay) {
                    int ni = idx(x, y, s + 1);
                    if (!seen[ni]) {
                        seen[ni] = 1;
                        prv[ni] = {idx(x, y, s), 'A', DIRC[d]};
                        que.emplace(x, y, s + 1);
                    }
                    break; // 隣接は高々 1 方向
                }
            }
        }
    }

    /* ---- 経路復元 ---- */
    vector<pair<char, char>> ret;
    if (gidx == -1) return ret; // 到達不可能 (設計ミス)

    for (int cur = gidx; cur != sidx;) {
        auto p = prv[cur];
        ret.emplace_back(p.act, p.dir);
        cur = p.prev;
    }
    reverse(ret.begin(), ret.end());
    return ret;
}

} // namespace detail

/* -------------------------------------------------------------------------- */
/* 要求関数：全区間まとめて文字列 (動作列) を生成                              */
/* -------------------------------------------------------------------------- */
std::vector<std::pair<char, char>> RetrieveAll(const std::array<Point, M> &points, const AlterPlans &plans) {
    using namespace std;
    vector<pair<char, char>> answer;

    /* 現在の盤面 (true = block) ―― 区間をまたいで持ち越し */
    array<array<uint8_t, N>, N> board{};
    for (auto &row : board) row.fill(0);

    for (int i = 1; i < M; ++i) {
        /* 区間 i-1 → i を最短で組み立てる */
        auto seg = detail::bfs_segment(points[i - 1], points[i], plans[i], board);
        answer.insert(answer.end(), seg.begin(), seg.end());

        /* 盤面を実際に更新 (トグル適用) */
        for (Point p : plans[i]) {
            int x = p / N, y = p % N;
            board[x][y] ^= 1;
        }
    }
    return answer; // vector<pair<action, dir>>
}
