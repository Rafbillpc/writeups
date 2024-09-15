#include "header.hpp"

const f64 TL = 1.996;
timer TIMER;

const i64 N = 1000;
const i64 NUM_NODES = N + N-1;

struct pt {
  i32 x, y;

  auto operator<=>(pt const& o) const = default;
};

i64 combine_cost(pt const& a, pt const& b) {
  return abs(a.x-b.x) + abs(a.y-b.y);
}

pt combine(pt const& a, pt const& b) {
  return pt { min(a.x,b.x), min(a.y, b.y) };
}

pt input[N];

void read(){
  { i64 n; cin>>n; runtime_assert(n == N); }
  FOR(i, N) {
    cin>>input[i].x>>input[i].y;
  }
}

struct node {
  bool          is_leaf;
  pt            point;
  array<i32, 2> children;
  i32           parent;
  i64           cost;

  FORCE_INLINE i64 other(i64 x) const{
    FOR(j, 2) if(children[j] != x) return children[j];
    return -1;
  }
};

struct state {
  node nodes[NUM_NODES];
  i64  cost;

  void greedy_init() {
    cost = 0;

    i32 next_node = 0;
    vector<i32> I;
    FOR(i, N) {
      nodes[next_node] = node {
        .is_leaf = true,
        .point = input[i],
        .children = {-1,-1},
        .parent = -1,
        .cost = 0,
      };
      I.eb(next_node);
      next_node += 1;
    }

    while(I.size() >= 2) {
      i32 bi = -1, bj = -1; i64 bc = 1e18;
      FOR(i0, I.size()) FOR(j0, i0) {
        i32 i = I[i0], j = I[j0];
        auto c = combine_cost(nodes[i].point, nodes[j].point);
        if(c < bc) {
          bi = i;
          bj = j;
          bc = c;
        }
      }

      nodes[bi].parent = next_node;
      nodes[bj].parent = next_node;
      cost += combine_cost(nodes[bi].point, nodes[bj].point);
      
      nodes[next_node] = node {
        .is_leaf = false,
        .point = combine(nodes[bi].point, nodes[bj].point),
        .children = {bi, bj},
        .parent = -1,
        .cost = combine_cost(nodes[bi].point, nodes[bj].point),
      };
      I.eb(next_node);
      I.erase(find(all(I), bi));
      I.erase(find(all(I), bj));
      next_node += 1;
    }

    runtime_assert(next_node == NUM_NODES);
    debug(TIMER.elapsed());
  }

  vector<array<pt, 2>> reconstruct() const {
    i32 root = 0;
    FOR(i, NUM_NODES) if(nodes[i].parent == -1) { root = i; break; }
    vector<array<pt, 2>> solution;
    auto dfs = [&](auto self, i32 i, pt par) -> void {
      solution.pb({ par, nodes[i].point });
      if(!nodes[i].is_leaf) {
        FOR(j, 2) self(self, nodes[i].children[j], nodes[i].point);
      }
    };
    dfs(dfs, root, {0,0});

    return solution;
  }

  FORCE_INLINE bool is_ancestor_of(i32 i, i32 j) const {
    while(i != j) {
      if(j == -1) return false;
      if(nodes[j].point.x < nodes[i].point.x ||
         nodes[j].point.y < nodes[j].point.y) return false;
      j = nodes[j].parent;
    }
    return true;
  }
  
  FORCE_INLINE void cut(i32 i) {
    i32 p = nodes[i].parent;
    if(p == -1) return;
    nodes[i].parent = -1;
    FOR(j, 2) if(nodes[p].children[j] == i) {
      nodes[p].children[j] = -1;
      return;
    }
  }

  FORCE_INLINE void link(i32 i, i32 p) {
    if(p != -1) {
      nodes[i].parent = p;
      FOR(j, 2) if(nodes[p].children[j] == -1) {
        nodes[p].children[j] = i;
        break;
      }
    }
  }

  void update(i32 p) {
    while(1) {
      if(p == -1) return;
      i32 x = nodes[p].children[0], y = nodes[p].children[1];
      if(x == -1 || y == -1) return;
      
      i64 new_cost = combine_cost(nodes[x].point, nodes[y].point);
      cost += new_cost - nodes[p].cost;
      nodes[p].cost = new_cost;
      
      auto new_point = combine(nodes[x].point, nodes[y].point);
      if(new_point == nodes[p].point) return;
      nodes[p].point = new_point;
      
      p = nodes[p].parent;
    }
  }

  // List of transitions:
  // - Rotation
  // - Swap two nodes (BAD)
  // - Insert node as neighbor of other node

  template<class ACC>
  void transition_swap(ACC&& acc) {
    i32 x1 = rng.random32(NUM_NODES);
    if(nodes[x1].parent == -1) return;
    i32 x2 = rng.random32(NUM_NODES);
    if(nodes[x2].parent == -1) return;

    if(is_ancestor_of(x1, x2) || is_ancestor_of(x2, x1)) return;

    i32 p1 = nodes[x1].parent;
    i32 p2 = nodes[x2].parent;

    cut(x1);
    cut(x2);
    link(x1, p2);
    link(x2, p1);
    update(p1);
    update(p2);

    if(!acc()) {
      cut(x1);
      cut(x2);
      link(x1, p1);
      link(x2, p2);
      update(p1);
      update(p2);
    }
  }
  
  template<class ACC>
  void transition_rotate(ACC&& acc) {
    i32 x1 = rng.random32(NUM_NODES);
    if(nodes[x1].parent == -1) return;
    i32 p1 = nodes[x1].parent;
    if(nodes[p1].parent == -1) return;
    i32 q1 = nodes[p1].parent;

    auto r1 = nodes[q1].parent; // can be (-1);

    i32 x2 = nodes[p1].other(x1);
    // i32 p2 = nodes[q1].other(p1);

    // DO
    cut(x1);
    cut(x2);
    cut(p1);
    cut(q1);
    link(x1, p1);
    link(x2, q1);
    link(p1, r1);
    link(q1, p1);
    update(q1);
    update(p1);

    if(!acc()) {
      // UNDO
      cut(x1);
      cut(x2);
      cut(p1);
      cut(q1);
      link(x1, p1);
      link(x2, p1);
      link(p1, q1);
      link(q1, r1);
      update(p1);
      update(q1);
    }
  }

  template<class ACC>
  void transition_move_node(ACC&& acc) {
    i32 x1 = rng.random32(NUM_NODES);
    if(nodes[x1].parent == -1) return;
    i32 x2 = rng.random32(NUM_NODES);
    if(nodes[x2].parent == -1) return;
    
    i32 p1 = nodes[x1].parent;
    i32 p2 = nodes[x2].parent;

    i32 y1 = nodes[p1].other(x1);
    i32 y2 = nodes[p2].other(x2);

    if(x1 == y2) return;

    i64 delta = 0;
    delta -= combine_cost(nodes[x1].point, nodes[y1].point);
    delta -= combine_cost(nodes[x2].point, nodes[y2].point);
    delta += combine_cost(nodes[x1].point, nodes[y2].point);
    delta += combine_cost(nodes[x2].point, nodes[y1].point);
    if(delta > 2.5e8) return;

    if(is_ancestor_of(x1, x2) || is_ancestor_of(x2, x1)) {
      return;
    }
    
    i32 q1 = nodes[p1].parent; // possibly (-1)
    
    // DO
    cut(x1);
    cut(y1);
    cut(x2);
    cut(p1);
    link(x1, p1);
    link(x2, p1);
    link(p1, p2);
    link(y1, q1);
    update(p1);
    update(p2);
    update(q1);

    if(!acc()) {
      // UNDO
      cut(y1);
      cut(p1);
      cut(x2);
      cut(x1);
      link(x1, p1);
      link(y1, p1);
      link(p1, q1);
      link(x2, p2);
      update(p1);
      update(p2);
      update(q1);
    }
  }

};

void print_solution(vector<array<pt, 2>> const& sol) {
  cout << sol.size() << '\n';
  for(auto p : sol) {
    FOR(i, 2) cout << p[i].x << ' ' << p[i].y << ' ';
    cout << '\n';
  }
  cout << flush;
}

struct union_find {
  vector<i64> A, S;

  union_find(i64 n = 0) : A(n), S(n,0) {
    iota(all(A), 0);
  }

  i64 find(i64 a) {
    return A[a] == a ? a : A[a] = find(A[a]);
  }

  void unite(i64 a, i64 b) {
    a = find(a);
    b = find(b);
    if(a == b) {
      return;
    }
    A[a] = b;
    S[b] += S[a];
  }
};

//
//  Doesn't yield any improvement.
//
// void combine_solutions(state const& S, state const& T){
//   map<pt, i32> M_vertex;
//   auto get_vertex = [&](pt p) {
//     if(M_vertex.count(p)) return M_vertex[p];
//     i32 ix = M_vertex.size();
//     M_vertex[p] = ix;
//     return ix;
//   };
//   FOR(i, NUM_NODES) get_vertex(S.nodes[i].point);
//   FOR(i, NUM_NODES) get_vertex(T.nodes[i].point);

//   map<tuple<i32, i32>, i64> S_edges;
//   auto add_edge = [&](pt p, pt q, i64 sign)  {
//     i64 d = abs(p.x-q.x)+abs(p.y-q.y);
//     S_edges[mt(get_vertex(p),get_vertex(q))] += sign * d;
//   };
  
//   FOR(i, NUM_NODES) {
//     if(S.nodes[i].parent != -1) {
//       add_edge(S.nodes[S.nodes[i].parent].point, S.nodes[i].point, 1);
//     }
//   }
//   FOR(i, NUM_NODES) {
//     if(T.nodes[i].parent != -1) {
//       add_edge(T.nodes[T.nodes[i].parent].point, T.nodes[i].point, -1);
//     }
//   }

//   union_find uf(M_vertex.size());
  
//   for(auto p : S_edges) {
//     if(p.second != 0) {
//       i32 x = get<0>(p.first), y = get<1>(p.first);
//       uf.unite(x,y);
//       uf.S[uf.find(x)] += p.second;
//     }
//   }

//   FOR(i, M_vertex.size()) {
//     if(uf.find(i) == i && uf.S[i] > 0) debug(i, uf.S[i]);
//   }
// }

void solve(){
  state S; S.greedy_init();
  debug(S.cost);

  state best_state = S;

  f64 temp0 = 5e8;
  f64 temp1 = 5e2;
  f64 temp = temp0;
  f64 done = 0.0;
  
  i64 niter = 0;
  while(1) {
    niter += 1;
    if(niter % (1<<13) == 0) {
      done = TIMER.elapsed() / TL;
      temp = temp0 * pow(temp1/temp0, done);
      if(done > 1.0) break;
    }

    auto old_cost = S.cost;
    auto accept = [&]() ALWAYS_INLINE -> bool {
      auto delta = S.cost - old_cost;
      return delta <= 0 || delta <= temp * rng.randomDouble();
    };
    
    i32 ty = rng.random32(40);
    if(ty <= 34) {
      S.transition_move_node(accept);
    }else if(ty <= 39) {
      S.transition_rotate(accept);
    }

    if(done > 0.8 && S.cost < best_state.cost) best_state = S;
  }

  debug(niter);
  debug(S.cost);

  { auto sol = best_state.reconstruct();
    print_solution(sol);
  }
}

int main(){
  TIMER.reset();
  read();
  solve();

  return 0;
}
