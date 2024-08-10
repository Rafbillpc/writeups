// -*- compile-command: "g++ -std=c++20 -O3 -march=native -DLOCAL -o Arrows Arrows.cpp" -*-
#include "header.hpp"

const f64 TL_BEAM  = 0.9; 
const f64 TL_POST  = 0.975;
const f64 TL_POST2 = 1.0;

const f64 TL       = 9.75;

timer TIMER;

const i32 N = 30;
const i32 NN = N*N;

i32 n, nn;
i32 grid_dir[NN];
i32 grid_mult[NN];

i32 bitset_size;
i32 beam_state_size;

i32 nrepeat = 2;

struct llist {
  i32 data;
  i32 prev;
};

const i32 LLIST_POOL_SIZE = 1<<25;
llist* llist_pool;
i32 llist_pool_next = 0;

i32 make_llist(i32 data, i32 prev) {
  runtime_assert(llist_pool_next < LLIST_POOL_SIZE);
  llist_pool[llist_pool_next] = {
    .data = data,
    .prev = prev,
  };
  return llist_pool_next++;
}

const u32 dx[8] = {(u32)(-1),(u32)(-1),0,1,1,1,0,(u32)(-1)};
const u32 dy[8] = {0,1,1,1,0,(u32)(-1),(u32)(-1),(u32)(-1)};

const i32 BEAM_DATA_SIZE = 1 << 28;
u8 *beam_data0, *beam_data1;

struct beam_state {
  u32 u;
  i32 history;
  i32 score;
  i64 heuristic;
  u64 hash;
};

struct graph_array {
  i32 size;
  u32 data[N];

  graph_array() : size(0) { }

  void push(u32 x) {
    data[size] = x;
    size += 1;
  }

  template<class F>
  void iter(F&& f) const {
    FOR(i, size) f(data[i]);
  }

  template<class F>
  void iterr(F&& f) const {
    FORD(i, size-1, 0) f(data[i]);
  }
};

graph_array graph[NN];
bitset<NN> graph_adj[NN];

struct cgraph {
  i32 n;
  vector<i32> pos;
  vector<i32> data;

  void init(vector<vector<i32>> G) {
    i32 cur = 0;
    n = G.size();
    pos.resize(n+1, 0);
    i32 ne = 0;
    FOR(i, n) ne += G[i].size();
    data.resize(ne);
    FOR(i, n) {
      for(auto j : G[i]) {
        data[cur] = j;
        cur += 1;
      }
      pos[i+1] = cur;
    }
  }

  template<class F>
  void iter(i32 i, F&& f) const {
    FORU(x, pos[i], pos[i+1]-1) f(data[x]);
  }
};

cgraph cgraph, crgraph;

void read() {
  cin>>n; nn = n * n;
  cerr << "[DATA] n = " << n << endl;
  
  TIMER.reset();
  FOR(i, n) FOR(j, n) cin>>grid_dir[i*n+j]>>grid_mult[i*n+j];

  bitset_size = (nn+63)/64 * sizeof(u64);
  beam_state_size = sizeof(beam_state) + bitset_size;

  FOR(i, NN) {
    graph[i].size = 0;
    graph_adj[i] = 0;
  }
  
  vector<vector<i32>> G(nn), RG(nn);
  FOR(x0, n) FOR(y0, n) {
    i32 d = grid_dir[x0*n+y0];
    u32 x = x0 + dx[d];
    u32 y = y0 + dy[d];
    while(x < n && y < n) {
      auto a = x0*n+y0;
      auto b = x*n+y;
      graph[a].push(b);
      G[a].pb(b);
      RG[b].pb(a);
      graph_adj[a][b] = 1;
      x += dx[d]; y += dy[d];
    }
  }
  cgraph.init(G);
  crgraph.init(RG);

  i32 n0 = 0, n1 = 0;
  FOR(i, nn) {
    if(RG[i].size() == 1) {
      n0 += 1;
    }
    if(G[i].size() == 1) {
      n1 += 1;
    }
  }
  debug(n0, n1);
}

beam_state &get_beam_state(u8 *beam_data, i32 id){
  return *(beam_state*)(beam_data + id * beam_state_size);
}

u8 *get_beam_bitset(u8* beam_data, i32 id) {
  return beam_data + id * beam_state_size + sizeof(beam_state);
}

void fill_bitset(u8 *b) { memset(b,0,bitset_size); }
void copy_bitset(u8 *b, u8 const* c) { memcpy(b,c,bitset_size); }
FORCE_INLINE
bool get_bitset(u8 const* b, i32 id) { return ((u64*)b)[id / 64] & (1ull<<(id % 64)); }
FORCE_INLINE
void set_bitset(u8* b, i32 id) { ((u64*)b)[id / 64] |= (1ull<<(id % 64)); }
FORCE_INLINE
void unset_bitset(u8* b, i32 id) { ((u64*)b)[id / 64] &= ~(1ull<<(id % 64)); }

void get_beam_cutoff
(i32 nb, i32 width, i64 &cutoff, f64 &cutoff_prob)
{
  vector<i64> scores;
  FOR(ib, nb) {
    beam_state const& sb = get_beam_state(beam_data0, ib);
    scores.eb(sb.heuristic);
  }
  if(scores.size() <= width) {
    cutoff = 0;
    cutoff_prob = 1.0;
    return;
  }
  nth_element(begin(scores), begin(scores)+width, end(scores), greater<>());
  cutoff = scores[width];
  i32 above = 0, at = 0;
  FOR(ib, nb) {
    if(scores[ib] > cutoff) above += 1;
    if(scores[ib] == cutoff) at += 1;
  }
  cutoff_prob = (f64)(width-above) / at;
}

f64 heuristic_table[NN][N];

f64 heuristic_mult = 1.05;
const i32 heuristic_in_ncoeff = 4;
f64 heuristic_in_coeff[heuristic_in_ncoeff] =
  { 0.2, 0.4, 0.4, 0.4 };
const i32 heuristic_out_ncoeff = 4;
f64 heuristic_out_coeff[heuristic_out_ncoeff] =
  { 0.3, 0.3, 0.3, 0.3 };

void init_heuristic() {
  FOR(indeg, NN) FOR(outdeg, N) {
    f64 pin = 1.0, pout = 1.0;
    FOR(i, indeg) pin = pin * heuristic_in_coeff[min(i, heuristic_in_ncoeff-1)];
    FOR(i, outdeg) pout = pout * heuristic_out_coeff[min(i, heuristic_out_ncoeff-1)];
    heuristic_table[indeg][outdeg]
      = heuristic_mult * (1 - pin) * (1 - pout);
  }
}

void calc_heuristic(i32 step, beam_state& sa, u8 const* ta) {
  static i32 last_vis[NN];
  static i32 in_deg[NN];
  static i32 out_deg[NN];
  static i32 cur_date = 0;
  cur_date += 1;
  static u32 Q[NN], Q2[NN];

  f64 vis_count[6] = {0,0,0,0,0,0};
  
  i32 nq = 0, nq2 = 0;
  Q[nq++] = sa.u; last_vis[sa.u] = cur_date; 

  FOR(iq, nq) {
    auto p = Q[iq];
    out_deg[p] = 0;

    FOR(ito, graph[p].size) {
      u32 to = graph[p].data[ito];
      if(get_bitset(ta, to)) continue;
      if(last_vis[to] != cur_date) {
        in_deg[to] = 0;
        last_vis[to] = cur_date;
        Q[nq++] = to;
        Q2[nq2++] = to;
      }
      in_deg[to] += 1;
      out_deg[p] += 1;
    }
  }
  in_deg[sa.u] = 0;

  FOR(iq, nq2) {
    auto p = Q2[iq];

    vis_count[grid_mult[p]] +=
      heuristic_table[in_deg[p]][out_deg[p]];
  }

  f64 heuristic = sa.score;
  f64 fstep = step;
  FORU(i, 1, 5) {
    heuristic += i * (fstep - 0.5) * vis_count[i] * 1.0;
    heuristic += i * vis_count[i] * vis_count[i] * 0.5;
    fstep += vis_count[i];
  }

  sa.heuristic = 1000000000.0 * heuristic;
}

const i32 HASH_BITS = 16;
const i32 HASH_SIZE = 1<<HASH_BITS;
const i32 HASH_MASK = HASH_SIZE-1;

struct {
  u64 hash;
  i32 pos;
} hash_table[HASH_SIZE];

u64 hash_vis[NN];
u64 hash_pos[NN];

void init_hash(){
  FOR(i, NN) hash_vis[i] = rng.randomInt64();
  FOR(i, NN) hash_pos[i] = rng.randomInt64();
}

vector<i32> solve_beam(i32 irepeat, f64 tl){
  i32 na = 0;
  FOR(xy, nn) if(xy % nrepeat == irepeat) {
    beam_state& sa = get_beam_state(beam_data0, na);
    u8* ta = get_beam_bitset(beam_data0, na);
    na += 1;
    sa.u = xy;
    sa.score = grid_mult[xy];
    sa.history = (-1);
    sa.heuristic = sa.score;
    sa.hash = hash_vis[xy] ^ hash_pos[xy];
    fill_bitset(ta);
    set_bitset(ta, xy);
  }

 
  i32 width = 1.0 * 100 * pow(1.0 * N / n, 3.0);
  i32 best_score = 0;
  i32 best_history = (-1);
  i64 total_calls = 0;

  FORU(step, 2, nn) {
    f64 t0 = TIMER.elapsed();
    
    debug(step, width);
    
    i64 cutoff = 0;
    f64 cutoff_prob = 1.0;
    get_beam_cutoff(na, width, cutoff, cutoff_prob);

    i32 nb = 0;

    FOR(ia, na) {
      beam_state& sa = get_beam_state(beam_data0, ia);
      u8 const* ta = get_beam_bitset(beam_data0, ia);

      if(sa.heuristic > cutoff || (sa.heuristic == cutoff && rng.randomFloat() < cutoff_prob)) {
        sa.history = make_llist(sa.u, sa.history);

        i32 count = 0;
        graph[sa.u].iter([&](u32 to) ALWAYS_INLINE -> void {
          if(count == 8) return;
          if(get_bitset(ta, to)) return;

          count += 1;

          beam_state &sb = get_beam_state(beam_data1, nb);
          u8* tb = get_beam_bitset(beam_data1, nb);

          nb += 1;
          sb.u = to;
          sb.score = sa.score + step * grid_mult[to];
          sb.history = sa.history;
          sb.hash = sa.hash ^ hash_vis[to] ^ hash_pos[sa.u] ^ hash_pos[to];
          copy_bitset(tb, ta);
          set_bitset(tb, to);

          total_calls += 1;
          calc_heuristic(step, sb, tb);

          auto& he = hash_table[sb.hash & HASH_MASK];
          if(he.hash == sb.hash) {
            auto& sc = get_beam_state(beam_data1, he.pos);
            u8* tc = get_beam_bitset(beam_data1, he.pos);
            if(sb.heuristic > sc.heuristic) {
              sc = sb;
              copy_bitset(tc, tb);
              nb -= 1;
            }else{
              nb -= 1;
            }
          }else{
            he.hash = sb.hash;
            he.pos = nb;
          }
        });
      }
    }
 
    swap(beam_data0, beam_data1);
    na = nb;

    FOR(ia, na) {
      beam_state const& sa = get_beam_state(beam_data0, ia);
      if(sa.score > best_score) {
        best_score = sa.score;
        best_history = make_llist(sa.u, sa.history);
      }
    }

    f64 t1 = TIMER.elapsed();

    f64 cur = t1 - t0;
    f64 goal = 1.2 * (tl - t0) / (nn - step + 1); 
    width = ceil(width * sqrt(goal / cur));
    width = max(width, 5);
    width = min(width, 1'000'000);
  }

  debug(best_score);
  
  vector<i32> ans;
  for(auto h = best_history; h != -1; h = llist_pool[h].prev) {
    auto xy = llist_pool[h].data;
    ans.pb(xy);
  }
  reverse(all(ans));

  cerr << "calls: " << total_calls << endl;

  return ans;
}

vector<i32> postprocess(vector<i32> ans, f64 tl) {
  { vector<i32> seen(nn,0);
    for(auto p : ans) seen[p] = 1;
    FOR(i, nn) if(!seen[i]) ans.eb(i);
  }

  auto sim = [&](i32& sz) -> i32 {
    i32 u = ans[0];
    i32 score = grid_mult[u];
    sz = nn;
    FORU(i, 1, nn-1) {
      i32 v = ans[i];
      if(!graph_adj[u][v]) { sz = i; break; }
      score += (i+1) * grid_mult[v];
      u = v;
    }
    return score;
  };

  i32 sz;
  i32 score = sim(sz);
  debug(score);

  vector<i32> best_ans = ans;
  i32 best_score = score;
  
  i32 niter = 0;
  f64 t0 = TIMER.elapsed();
  f64 temp0 = sqrt(score);
  f64 temp1 = 1.0;
  f64 temp = temp0;
  
  while(1) {
    niter += 1;
    if(niter % 1024 == 0) {
      f64 done = (TIMER.elapsed() - t0) / (tl - t0);
      if(done > 1.0) break;
      temp = temp0 * pow(temp1 / temp0, done);
    }

    i32 ty = rng.random32(3);
    if(ty == 0) {
      i32 i = rng.random32(nn+1);
      i32 j = rng.random32(nn+1);
      i32 k = rng.random32(nn+1);
      if(i > j) swap(i,j);
      if(j > k) swap(j,k);
      if(i > j) swap(i,j);
      if(i == j || j == k) continue;
      if(i > 0 && !graph_adj[ans[i-1]][ans[j]]) continue;
      if(!graph_adj[ans[k-1]][ans[i]]) continue;
      rotate(begin(ans)+i, begin(ans)+j, begin(ans)+k);
      i32 new_score = sim(sz);
      if(new_score >= score || score-new_score <= temp * rng.randomDouble()) {
        score = new_score;
      }else{
        rotate(begin(ans)+i, begin(ans)+(i+(k-j)), begin(ans)+k);
      }
    }else if(ty == 1) {
      i32 i = rng.random32(nn+1);
      i32 j = rng.random32(nn);
      i32 k = j+1;
      if(i >= j) continue;
      if(i > 0 && !graph_adj[ans[i-1]][ans[j]]) continue;
      if(!graph_adj[ans[k-1]][ans[i]]) continue;
      rotate(begin(ans)+i, begin(ans)+j, begin(ans)+k);
      i32 new_score = sim(sz);
      if(new_score >= score || score-new_score <= temp * rng.randomDouble()) {
        score = new_score;
      }else{
        rotate(begin(ans)+i, begin(ans)+(i+(k-j)), begin(ans)+k);
      }
    }else{
      i32 i = rng.random32(nn+1);
      i32 j = i+1;
      i32 k = rng.random32(nn+1);
      if(j >= k) continue;
      if(i > 0 && !graph_adj[ans[i-1]][ans[j]]) continue;
      if(!graph_adj[ans[k-1]][ans[i]]) continue;
      rotate(begin(ans)+i, begin(ans)+j, begin(ans)+k);
      i32 new_score = sim(sz);
      if(new_score >= score || score-new_score <= temp * rng.randomDouble()) {
        score = new_score;
      }else{
        rotate(begin(ans)+i, begin(ans)+(i+(k-j)), begin(ans)+k);
      }
    }

    if(score > best_score) {
      best_ans = ans;
      best_score = score;
      debug(best_score);
    }
  }

  debug(niter);

  ans = best_ans;
  debug(best_score);
  i32 r = sim(sz);
  debug(r);
  ans.resize(sz);
  return ans;
}

vector<i32> postpostprocess(vector<i32> ans, f64 tl) {
  static i32 selected_in[NN];
  static i32 selected_out[NN];
  static i32 cur_date = 0;

  i32 cur_score = 0;
  
  while(TIMER.elapsed() < tl) {
    cur_date += 1;

    i32 i = rng.random32(nn);
    selected_in[i] = cur_date;
    i32 d = grid_dir[i];
    cgraph.iter(i, [&](i32 j) {
      if(grid_dir[j] != d) return;
      selected_in[j] = cur_date;
    });
    crgraph.iter(i, [&](i32 j) {
      if(grid_dir[j] != d) return;
      selected_in[j] = cur_date;
    });

    i32 score = 0;
    i32 steps = 0;
    i32 total = 0;
    i32 from = 0;

    struct range {
      i32 fr, to;
      i32 score, steps, total;
    };
    vector<range> ranges;
    
    FOR(j, ans.size()) {
      steps += 1;
      score += steps * grid_mult[ans[j]];
      total += grid_mult[ans[j]];
      if(selected_in[ans[j]] == cur_date) {
        ranges.pb(range {
            .fr = from,
            .to = j,
            .score = score,
            .steps = steps,
            .total = total,
          });
        
        if(j+1 < ans.size()) selected_out[ans[j+1]] = 1;
        steps = 0;
        score = 0;
        total = 0;
        from = j+1;
      }
    }
    ranges.pb(range {
        .fr = from,
        .to = (i32)ans.size()-1,
        .score = score,
        .steps = steps,
        .total = total,
      });
    i32 rsz = ranges.size();
    if(rsz <= 2) continue;
    
    vector<i32> I(rsz); iota(all(I), 0);
    vector<i32> BI = I;
    i32 c = 0;
    i32 best_score = 0;
    i32 iter = 0;
    do {
      iter += 1;
      if(iter % 1024 && TIMER.elapsed() > TL_POST2) {
        break;
      }
      { bool ok = 1;
        FOR(j, rsz-1) if(ranges[I[j+1]].fr < ans.size() &&
                         !graph_adj[ans[ranges[I[j]].to]][ans[ranges[I[j+1]].fr]]) {
          ok = 0;
          break;
        }
        if(!ok) continue;
      }
      i32 score = 0;
      i32 steps = 0;
      FOR(j, rsz) {
        score += ranges[I[j]].score;
        score += ranges[I[j]].total * steps;
        steps += ranges[I[j]].steps;
      }
      if(score >= best_score) {
        if(score > best_score) {
          best_score = score;
          c = 0;
        }
        c += 1;
        if(rng.random32(c) == 0) {
          BI = I;
        }
      }
    } while(next_permutation(begin(I)+1, begin(I)+rsz-1));

    if(best_score > cur_score) {
      cur_score = best_score;
      debug(cur_score);
    }
    
    vector<i32> nans;
    for(auto j : BI) {
      FORU(x, ranges[j].fr, ranges[j].to) nans.eb(ans[x]);
    }
    ans = nans;
  }
  return ans;
}

vector<i32> best_ans;
i32 best_score;

void try1(i32 i, f64 tl) {
  debug(tl);
  f64 t0 = TIMER.elapsed();
  auto ans = solve_beam(i, t0 + (tl - t0) * TL_BEAM);
  ans = postprocess(ans, t0 + (tl - t0) * TL_POST);
  ans = postpostprocess(ans, t0 + (tl - t0) * TL_POST2);
  
  i32 score = 0;
  i32 steps = 0;
  for(auto i : ans) {
    steps += 1;
    score += grid_mult[i] * steps;
  }
  if(score > best_score) {
    best_ans = ans;
    best_score = score;
  }
}

void solve() {
  FOR(i, nrepeat){
    try1(i, (i+1) * TL / nrepeat);
  }
  
  cerr << "[DATA] elapsed = " << TIMER.elapsed() << endl;
  
  cout << best_ans.size() << endl;
  for(auto p : best_ans) {
    cout << p/n << " " << p%n << endl;
  }
}

int main(){
  read();

  // ofstream os(argv[1]);
  // os << n << endl;
  // FOR(i, n) FOR(j, n) {
  //   os << grid_dir[i][j] << ' ' << grid_mult[i*n+j] << endl;
  // }
  // os.close();

  beam_data0 = new u8[BEAM_DATA_SIZE];
  beam_data1 = new u8[BEAM_DATA_SIZE];
  llist_pool = new llist[LLIST_POOL_SIZE];
  init_hash();
  init_heuristic();
  
  solve();

  return 0;
}
