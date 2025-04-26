とりあえず、途中状態として 
```
using Point = int;
Point f(int x, int y) { return x * N + y; }
bool isin(int x, int y) { return 0 <= x and x < N and 0 <= y and y < N; }

using AlterPlan = std::vector<Point>;
using AlterPlans = std::array<AlterPlan, M>;
```
のようなクラスを導入しました。 ここで i = 1, ..., M - 1 について、 AlterPlans[i] に含まれる Point は、i - 1 番目の目的地から i 番目の目的地まで移動するまでの間にこの順にブロックが置かれる / 撤去される場所です。
AlterPlans と 入力で与えられた目的地たちの情報をもとに、本問題で出力すべき文字列を出力する関数 std::vector<std::pair<char, char>> RetrieveAll(const std::array<Point, M> &points, const AlterPlans &plans) を C++ で書いてください。

なお、実装方針は以下です。
- これは、単純な BFS を繰り返し適用することで求められるはずです。ただし、ブロックは複数の方向から設置することが可能なので、それらを全て考慮した上で最短のケースを構築してください。単に N*N 頂点のグラフを考えるのではなく、 ブロックの設置個数を K として、 N * N * (M + K) 頂点くらいのグラフを考えて拡張 Dijkstra 法のようなことをするとよいはずです。
- 既に設置したブロックの上を通るような経路を出力しないよう気をつけてください。
