#include "common.hpp"
#include <array>

using namespace std;
using lint = long long;
using pint = std::pair<int, int>;
using plint = std::pair<lint, lint>;

struct fast_ios {
    fast_ios() {
        std::cin.tie(nullptr), std::ios::sync_with_stdio(false), std::cout << std::fixed << std::setprecision(20);
    };
} fast_ios_;

#include "ahc046.hpp"
#include "llm.hpp"
#include "rnd.hpp"

constexpr int inf = 1e7;
int EvaluateSteps(Point s, Point t, const AlterPlan &plan, BanStates ban) {
    std::array<int, N * N> dp, dp_next;

    if (ban.test(s)) return inf;
    if (ban.test(t)) return inf;

    dp.fill(inf);
    dp.at(s) = 0;

    auto is_ok = [&](int x, int y) {
        return isin(x, y) and !ban.test(f(x, y));
    };

    std::vector<Point> st{s};
    int left = 0;

    auto func = [&](Point p, bool neighbor) {
        const int px = p / N, py = p % N;

        int rem = 1;

        if (neighbor) {
            rem = 0;
            REP(dir, 4) {
                const int nx = px + dx[dir], ny = py + dy[dir];
                if (is_ok(nx, ny)) ++rem;
            }
        }

        while (left < (int)st.size()) {
            const Point now = st.at(left++);
            const int x = now / N, y = now % N;

            if (!neighbor and now == p) break;
            if (neighbor and abs(x - px) + abs(y - py) == 1) {
                --rem;
                if (rem <= 0) break;
            }

            REP(dir, 4) {
                int nx = x + dx[dir], ny = y + dy[dir];
                if (!is_ok(nx, ny)) continue;

                if (const int np = f(nx, ny); chmin(dp.at(np), dp.at(now) + 1)) {
                    st.push_back(np);
                }

                while (is_ok(nx + dx[dir], ny + dy[dir])) nx += dx[dir], ny += dy[dir];
                if (const int np = f(nx, ny); chmin(dp.at(np), dp.at(now) + 1)) {
                    st.push_back(np);
                }
            }
        }

        dp_next.fill(inf);
        if (neighbor) {
            const int x = p / N, y = p % N;
            st.clear();
            left = 0;

            REP(dir, 4) {
                const int nx = x + dx[dir];
                const int ny = y + dy[dir];
                const int np = f(nx, ny);
                if (is_ok(nx, ny) and chmin(dp_next.at(np), dp.at(np) + 1)) st.push_back(np);
            }
        } else {
            if (is_ok(px, py) and chmin(dp_next.at(p), dp.at(p))) {
                st.push_back(p);
            }
        }

        swap(dp, dp_next);
    };

    for (Point p : plan) {
        func(p, true);
        ban.set(p);
    }

    func(t, false);
    return dp.at(t);
}

std::array<int, N * N> ForwardDP(Point s, const BanStates &ban) {

    auto is_ok = [&](int x, int y) { return isin(x, y) and !ban.test(f(x, y)); };

    std::array<int, N * N> ret;
    ret.fill(inf);
    
    std::array<Point, N * N> st;
    int l = 0, r = 0;
    if (is_ok(s / N, s % N)) {
        ret.at(s) = 0;
        st[r++] = s;
    }

    ret.at(s) = 0;
    st[0] = s;
    while (l < r) {
        const Point now = st.at(l++);
        const int x = now / N, y = now % N;
        REP(dir, 4) {
            int nx = x + dx[dir], ny = y + dy[dir];
            if (!is_ok(nx, ny)) continue;
            if (chmin(ret.at(f(nx, ny)), ret.at(now) + 1)) { st[r++] = f(nx, ny); }

            while (is_ok(nx + dx[dir], ny + dy[dir])) {
                nx += dx[dir];
                ny += dy[dir];
            }

            if (chmin(ret.at(f(nx, ny)), ret.at(now) + 1)) { st[r++] = f(nx, ny); }
        }
    }

    return ret;
}

std::array<int, N * N> BackwardDP(Point t, const BanStates &ban) {
    auto is_ok = [&](int x, int y) { return isin(x, y) and !ban.test(f(x, y)); };

    std::array<int, N * N> ret;
    ret.fill(inf);
    std::array<Point, N * N> st;
    int l = 0, r = 0;
    if (is_ok(t / N, t % N)) {
        st[r++] = t;
        ret.at(t) = 0;
    }

    while (l < r) {
        const Point now = st.at(l++);
        const int x = now / N, y = now % N;
        REP(dir, 4) {
            int nx = x + dx[dir], ny = y + dy[dir];
            if (!is_ok(nx, ny)) continue;
            if (chmin(ret.at(f(nx, ny)), ret.at(now) + 1)) st[r++] = f(nx, ny);

            if (!is_ok(x + dx[dir ^ 1], y + dy[dir ^ 1])) {
                while (is_ok(nx, ny)) {
                    if (chmin(ret.at(f(nx, ny)), ret.at(now) + 1)) st[r++] = f(nx, ny);
                    nx += dx[dir];
                    ny += dy[dir];
                }
            }
        }
    }

    return ret;
}

std::array<int, N * N> GetReductionIf(Point s, Point t, const AlterPlan &plan, BanStates ban) {
    std::array<int, N * N> ret;
    ret.fill(0);

    Point now = s;

    auto eval = [&](Point nxt) {
        const auto forward = ForwardDP(now, ban);
        const auto backward = BackwardDP(nxt, ban);
        if (forward.at(nxt) != backward.at(now)) {
            dbg(make_tuple(s, t, plan, ban, forward.at(nxt), backward.at(now)));
        }
        assert(forward.at(nxt) == backward.at(now));

        const int base_cost = forward.at(nxt);


        REP(i, N * N) {
            if (ban.test(i)) continue;
            const int x = i / N, y = i % N;
            int best = 0;
            REP(dir, 4) {
                const int x1 = x + dx[dir], y1 = y + dy[dir];
                if (!isin(x1, y1)) continue;
                if (ban.test(f(x1, y1))) continue;

                int x2 = x1 + dx[dir], y2 = y1 + dy[dir];
                while (isin(x2, y2) and !ban.test(f(x2, y2))) {
                    int cost = forward.at(f(x2, y2)) + 1 + backward.at(f(x1, y1));
                    chmax(best, base_cost - cost);
                    x2 += dx[dir];
                    y2 += dy[dir];
                }
            }
            ret.at(i) += best;
        }
    };

    for (auto nxt : plan) {
        eval(nxt);
        ban.set(nxt);
    }
    eval(t);

    REP(i, N * N) {
        if (ban.test(i)) ret.at(i) = -inf;
    }
    ret.at(t) = -inf;

    return ret;
}

int FastEvaluateAll(std::array<Point, M> points, AlterPlans &plans) {

    int ret = 0;
    BanStates ban;
    FOR(i, 1, M) {
        ret += EvaluateSteps(points.at(i - 1), points.at(i), plans.at(i), ban);
        for (auto p : plans.at(i)) {
            ban.set(p);
        }
    }

    return ret;
}

void tspimprove(std::array<Point, M> points, AlterPlans &state) {
    BanStates ban;
    FOR(tick, 1, M) {
        auto &vec = state.at(tick);
        int eval = EvaluateSteps(points.at(tick - 1), points.at(tick), vec, ban);
        if (vec.size() > 1) {
            FOR(i, 1, vec.size()) {
                swap(vec.at(i - 1), vec.at(i));
                if (chmin(eval, EvaluateSteps(points.at(tick - 1), points.at(tick), vec, ban))) {
                } else {
                    swap(vec.at(i - 1), vec.at(i));
                }
            }
        }
        for (auto i : state.at(tick)) ban.set(i);
    }

    int eval = FastEvaluateAll(points, state);

    FOR(tick, 1, M - 1) {
        if (state.at(tick).empty()) continue;
        Point v = state.at(tick).back();
        state.at(tick).pop_back();
        state.at(tick + 1).insert(state.at(tick + 1).begin(), v);
        if (chmin(eval, FastEvaluateAll(points, state))) {

        } else {
            state.at(tick).push_back(v);
            state.at(tick + 1).erase(state.at(tick + 1).begin());
        }
    }
}

#include <chrono>

class timer_ {
    std::chrono::system_clock::time_point start_;

public:
    timer_() : start_(now()) {}

    static std::chrono::system_clock::time_point now() { return std::chrono::system_clock::now(); }

    int spent_ms() const {
        auto diff = now() - start_;
        return std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    }
} timer;


int main(int argc, char *argv[]) {
    {
        int N_, M_;
        cin >> N_ >> M_;
        assert(N_ == N and M_ == M);
    }
    array<int, M> xs, ys;
    REP(i, M) cin >> xs.at(i) >> ys.at(i);
    dbg(xs);
    dbg(ys);

    array<Point, M> points;
    REP(i, M) points.at(i) = f(xs.at(i), ys.at(i));
    dbg(points);

    auto greedy_insert = [&](AlterPlans state) -> pair<int, AlterPlans> {
        vector<tuple<int, int, int, Point>> cands;

        BanStates ban;
        FOR(t, 1, M) for (int i : state.at(t)) ban.set(i);

        std::array<int, N * N> gri_sum;
        gri_sum.fill(0);

        IFOR(t, 1, M) {

            IREP(d, state.at(t).size() + 1) {
                const Point s = d ? state.at(t).at(d - 1) : points.at(t - 1);
                const Point g = d < (int)state.at(t).size() ? state.at(t).at(d) : points.at(t);

                const auto backward = BackwardDP(g, ban);
                if (d) ban.reset(state.at(t).at(d - 1));
                const auto forward = ForwardDP(s, ban);
                const int base_cost = backward.at(s);
                if (base_cost >= inf / 2) continue;
                REP(i, N * N) {
                    if (gri_sum.at(i) < 0) continue;
                    int max_reward = -1;

                    REP(dir, 4) {
                        const int x = i / N + dx[dir];
                        const int y = i % N + dy[dir];
                        if (!isin(x, y) or ban.test(f(x, y))) continue;
                        const int cost_inc = forward.at(f(x, y)) + 1 + backward.at(f(x, y)) - base_cost;
                        chmax(max_reward, gri_sum.at(i) - cost_inc);
                    }

                    if (max_reward > 0) {
                        cands.emplace_back(-max_reward, t, d, i);
                    }
                }
            }

            const auto gri = GetReductionIf(points.at(t - 1), points.at(t), state.at(t), ban);
            REP(i, N * N) gri_sum.at(i) += gri.at(i);
        }

        sort(cands.begin(), cands.end());
        while (cands.size() > 100) cands.pop_back();
        // dbg(cands);
        int best_e = inf;
        int best_t = -1, best_d = -1, best_idx = -1;
        for (auto [_, t, d, idx] : cands) {
            state.at(t).insert(state.at(t).begin() + d, idx);
            const int e = FastEvaluateAll(points, state);
            state.at(t).erase(state.at(t).begin() + d);
            if (chmin(best_e, e)) {
                best_t = t;
                best_d = d;
                best_idx = idx;
                // dbg(make_tuple(best_e, best_t, best_d, best_idx));
            }
        }

        if (best_e < inf) {
            state.at(best_t).insert(state.at(best_t).begin() + best_d, best_idx);
            return {best_e, state};
        } else {
            return {inf, state};
        }
    };

    auto greedy_delete = [&](AlterPlans state) -> pair<int, AlterPlans> {
        int best_e = inf;
        int best_t = -1, best_d = -1;
        FOR(t, 1, M) {
            REP(d, state.at(t).size()) {
                const int v = state.at(t).at(d);
                state.at(t).erase(state.at(t).begin() + d);
                if (chmin(best_e, FastEvaluateAll(points, state))) {
                    best_t = t;
                    best_d = d;
                }
                state.at(t).insert(state.at(t).begin() + d, v);
            }
        }

        if (best_e < inf) {
            state.at(best_t).erase(state.at(best_t).begin() + best_d);
            return {best_e, state};
        } else {
            return {inf, state};
        }
    };

    auto try_move = [&](AlterPlans state) -> pair<int, AlterPlans> {
        int best_e = inf;
        AlterPlans best_state;
        FOR(t, 1, M) {
            auto &v = state.at(t);
            REP(d, v.size()) {
                const Point s = v.at(d);
                REP(dir, 4) {
                    const int x = s / N + dx[dir];
                    const int y = s % N + dy[dir];
                    if (!isin(x, y)) continue;
                    v.at(d) = f(x, y);
                    if (chmin(best_e, FastEvaluateAll(points, state))) { best_state = state; }
                    v.at(d) = s;
                }
            }
        }

        if (best_e < inf) {
            return {best_e, best_state};
        } else {
            return {inf, state};
        }
    };

    AlterPlans state;
    int opt = FastEvaluateAll(points, state);

    while (true) {
        auto [e, new_state] = greedy_insert(state);
        if (e < opt) {
            opt = e;
            state = new_state;
            dbg(make_tuple(opt, e));
        } else {
            break;
        }
    }

    while (true) {
        auto [e, new_state] = greedy_delete(state);
        if (chmin(opt, e)) {
            opt = e;
            state = new_state;
            dbg(make_tuple(-2, opt, e));
        } else {
            auto [e, newnew_state] = greedy_insert(new_state);
            if (chmin(opt, e)) {
                opt = e;
                state = newnew_state;
                dbg(make_tuple(-3, opt, e));
            } else {
                break;
            }
        }
    }

    FOR(tick, 1, M) {
        for (int d = 0; d < (int)state.at(tick).size(); d++) {
            auto nxt_state = state;
            nxt_state.at(tick).erase(nxt_state.at(tick).begin() + d);
            auto [e, new_state] = greedy_insert(nxt_state);
            if (chmin(opt, e)) {
                opt = e;
                state = new_state;
                dbg(make_tuple(-5, tick, d, opt));
            }
        }
    }

    REP(_, 0) {
        FOR(tick, 1, M) {
            if (state.at(tick).size()) {
                auto tmp = state;
                tmp.at(tick).clear();
                auto eval = FastEvaluateAll(points, state);
                while (true) {
                    auto [e, new_state] = greedy_insert(tmp);
                    tspimprove(points, new_state);
                    if (chmin(eval, e)) {
                        tmp = new_state;
                    } else {
                        break;
                    }
                }

                if (chmin(opt, eval)) {
                    dbg(make_tuple(_, tick, opt));
                    state = tmp;
                }
            }
        }
    }

    REP(_, 0) {
        const int choice = rand_int() % 2;
        if (choice == 0) {
            auto [e, new_state] = greedy_insert(state);
            tspimprove(points, new_state);
            if (e < opt) {
                opt = e;
                state = new_state;
                dbg(make_tuple(choice, opt, e));
            }
        } else if (choice == 1) {
            auto [e, new_state] = greedy_delete(state);
            tspimprove(points, new_state);
            if (e < opt) {
                opt = e;
                state = new_state;
                dbg(make_tuple(choice, opt, e));
            }
        } else {
            auto [e, new_state] = try_move(state);
            if (e < opt) {
                opt = e;
                state = new_state;
                dbg(make_tuple(choice, opt, e));
            }
        }
    }

    {
        std::string retstr;
        int next_target = 1;
        int x = xs.at(0), y = ys.at(0);
        vector is_blocked(N, vector<int>(N));

        int cost = 0;

        dbg(state);

        for (auto [a, d] : RetrieveAll(points, state)) {
            ++cost;
            // dbg(make_tuple(a, d, x, y));
            retstr += a;
            retstr += ' ';
            retstr += d;
            retstr += '\n';

            int dx = 0, dy = 0;
            if (d == 'U') {
                dx = -1;
            } else if (d == 'D') {
                dx = 1;
            } else if (d == 'L') {
                dy = -1;
            } else if (d == 'R') {
                dy = 1;
            } else {
                dbg(d);
                assert(false);
            }

            if (a == 'M') {
                const int nx = x + dx, ny = y + dy;
                assert(isin(nx, ny));
                assert(!is_blocked.at(nx).at(ny));
                x = nx;
                y = ny;
            } else if (a == 'A') {
                const int nx = x + dx, ny = y + dy;
                assert(isin(nx, ny));
                is_blocked.at(nx).at(ny) ^= 1;
            } else if (a == 'S') {
                while (isin(x + dx, y + dy) and !is_blocked.at(x + dx).at(y + dy)) {
                    x += dx;
                    y += dy;
                }
            } else {
                assert(false);
            }

            if (next_target < M and xs.at(next_target) == x and ys.at(next_target) == y) { ++next_target; }
        }

        if (next_target < M) {
            dbg(make_tuple(next_target));
        }

        assert(next_target == M);
        dump_onlinejudge(retstr);

        const int score = 1640 - cost;
        jdump("cost", cost);
        jdump("score", score);
        jdump("score150", score * 150);
    }
    // if (argc >= 2) { X = std::stoi(argv[1]); }
}
