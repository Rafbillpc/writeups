#include "header.hpp"
// #include "params.hpp"
const f64 temps0[4] = {0.7114597547328911, 0.72060661552706, 0.35759210358282834, 0.15};

timer TIMER;
#ifdef LOCAL
const f64 TIME_FACTOR = 0.5;
#else
const f64 TIME_FACTOR = 1.0;
#endif
const f64 TIME_LIMIT_CYCLES = TIME_FACTOR * 8.6;
const f64 TIME_LIMIT_END    = TIME_FACTOR * 9.7;

const i32 N = 100;
const i32 C = 6;
const i32 V = 10;

i32 n,nc,maxv;

struct pt {
  i32 x = 0, y = 0;
  pt() = default;
  pt(i32 x_, i32 y_) : x(x_), y(y_){ }

  bool operator==(pt const& o) const { return tie(x,y)==tie(o.x,o.y); }
};

struct star_t {
  i32 id;
  pt  pos;
  i32 v;
};

star_t stars[C][N];

FORCE_INLINE
i64 ccw(pt const& a, pt const& b, pt const& c) {
  return (c.y-a.y)*(b.x-a.x) - (b.y-a.y)*(c.x-a.x);
}

FORCE_INLINE
i64 ccw2(pt const& b, pt const& c) {
  return c.y*b.x - c.x*b.y;
}

FORCE_INLINE
bool intersects(pt const& a, pt const& b, pt const& c, pt const& d) {
  pt ab(b.x-a.x, b.y-a.y);
  pt ac(c.x-a.x, c.y-a.y);
  pt ad(d.x-a.x, d.y-a.y);
  return ccw2(ac,ad) * ccw(b,c,d) <= 0 && ccw2(ab,ac) * ccw2(ab,ad) <= 0;
}

void read() {
  cin>>n>>nc>>maxv;
  cerr << "[DATA] N=" << n << endl;
  cerr << "[DATA] C=" << nc << endl;
  cerr << "[DATA] V=" << maxv << endl;
  i32 nstar[C];
  FOR(i, nc) nstar[i] = 0;
  FOR(i, n*nc) {
    i32 r,c,col,v; cin>>r>>c>>col>>v;
    stars[col][nstar[col]++] = { .id = i, .pos = pt(r,c), .v = v };
  }
}

struct seg_t {
  i32 id;
  i32 col, a, b;
};

const i32 MAX_SEGS = 29'700;
i32 nsegs;
seg_t segs[MAX_SEGS];
i32 seg_id[C][N][N];

bitset<MAX_SEGS> seg_intersects[MAX_SEGS];

f64 logs[200001];
i32 nn2;

i32 seg_idf(i32 c, i32 a, i32 b) {
  if(a<b) swap(a,b);
  return c*nn2+a*(a-1)/2+b;
}

pt operator-(pt const& a, pt const& b) {
  return pt(a.x-b.x,a.y-b.y);
}

i64 area(pt p, pt q) {
  return abs(p.x*q.y-p.y*q.x);
}

i64 area(pt p, pt q, pt r) {
  return area(q-p,r-p);
}

void init(){
  nn2=n*(n-1)/2;
  FORU(i, 1, 200000) logs[i] = log(i+1);
  
  nsegs = 0;
  FOR(c1, nc) FOR(i1, n) FOR(j1, i1) {
    seg_id[c1][i1][j1] = nsegs;
    seg_id[c1][j1][i1] = nsegs;

    segs[nsegs] = { .id = nsegs, .col = c1, .a = i1, .b = j1 };
    nsegs += 1;
  }

  FOR(i, nsegs) seg_intersects[i][i] = 1;
  i32 id1 = 0;
  FOR(c1, nc) FOR(i1, n) FOR(j1, i1) {
    i32 id2 = 0;
    FOR(c2, c1) FOR(i2, n) FOR(j2, i2) {
      if(intersects(stars[c1][i1].pos, stars[c1][j1].pos,
                    stars[c2][i2].pos, stars[c2][j2].pos)) {
        seg_intersects[id1][id2] = true;
        seg_intersects[id2][id1] = true;
      }
      id2 += 1;
    }
    id1 += 1;
  }
  { i32 id1 = 0;
    FOR(c1, nc) {
      FOR(i1, n) FOR(j1, i1) {
        i32 id2 = c1*n*(n-1)/2;
        FOR(i2, i1) FOR(j2, i2){
          if(j1 != i2 && j1 != j2
             && intersects(stars[c1][i1].pos, stars[c1][j1].pos,
                           stars[c1][i2].pos, stars[c1][j2].pos)) {
            seg_intersects[id1][id2] = true;
            seg_intersects[id2][id1] = true;
          }
          id2 += 1;
        }
        id1 += 1;
      }
    }
  }
  cerr << "Intersects time: " << TIMER.elapsed() << endl;
}

struct segs_set {
  i32 nsegs;
  i32 segs[N*C];

  void reset(){
    nsegs = 0;
  }

  void push(i32 i) {
    FOR(j, nsegs) if(segs[j] == i) return;
    segs[nsegs++] = i;
  }
};

struct state_cycles {
  i32 ncycle[C];
  i32 cycle[C][N+1];
  bool in_cycle[C][N];

  i32 sum_cycle[C];

  i32 nsegs;
  i32 segs[N*C];
  i32 seg_where[MAX_SEGS];
  
  void reset() {
    nsegs = 0;
    FOR(i, nc) {
      ncycle[i] = 0;
      FOR(j, n) in_cycle[i][j] = 0;
      sum_cycle[i] = 0;
    }
    FOR(i, ::nsegs) seg_where[i] = -1;
  }

  void blocking(segs_set& set, i32 id, i32 max) {
    FOR(i, nsegs) if(seg_intersects[id][segs[i]]) {
      set.push(segs[i]);
      if(set.nsegs > max) return;
    }
  }
  
  bool can_add_segment(i32 id) {
    FOR(i, nsegs) if(seg_intersects[id][segs[i]]) return false;
    return true;
  }

  void add_seg(i32 id) {
    seg_where[id] = nsegs;
    segs[nsegs++] = id;
  }

  void rem_seg(i32 id) {
    i32 i = seg_where[id];
    swap(segs[i], segs[nsegs-1]);
    nsegs -= 1;
    seg_where[segs[i]] = i;
    seg_where[id] = -1;
  }

  bool has_seg(i32 id) const { return seg_where[id] != -1; }

  f64 value() const
  {
    f64 x = 0;
    FOR(i, nc) x += logs[2*ncycle[i]*sum_cycle[i]];
    return x;
  }
};

void greedy_cycles(state_cycles& best, f64 tl, i32 iter) {
  state_cycles S = best;

  auto try_add = [&](i32 c, i32 a, i32 b) -> bool {
    auto id = seg_idf(c,a,b);
    if(!S.can_add_segment(id)) return false;
    S.add_seg(id);
    return true;
  };
  
  if(S.ncycle[0] == 0) {
    FOR(i, nc) {
      vector<i32> I(n); iota(all(I),0);
      rng.shuffle(I);
      if(!try_add(i,I[0],I[1])) return;
      if(!try_add(i,I[1],I[2])) return;
      if(!try_add(i,I[2],I[0])) return;

      S.ncycle[i] = 3;
      FOR(j,3) {
        S.cycle[i][j] = I[j];
        S.in_cycle[i][I[j]] = 1;
        S.sum_cycle[i] += stars[i][I[j]].v;
      }
      S.cycle[i][3] = I[0];
    }
  }
  
  segs_set set;
  
  const i32 MAX_ITER =
    iter == 0 ? 200'000 :
    iter == 1 ? 500'000 :
    1'000'000;
  bool enable_trans = rng.random32(10) != 0;

  f64 done = 0;
  f64 temp0 = temps0[iter];
  f64 temp = temp0;
  f64 time0 = TIMER.elapsed();
  i64 niter = 0;
  bool use_time = tl > 0;
  
  while(1) {
    niter += 1;
    if(niter % 1024 == 0) {
      FOR(c, nc) {
        i32 x = rng.random32(S.ncycle[c]);
        rotate(S.cycle[c], S.cycle[c]+x, S.cycle[c]+S.ncycle[c]);
        S.cycle[c][S.ncycle[c]] = S.cycle[c][0];
      }
      done = use_time ? (TIMER.elapsed() - time0) / tl : (f64)niter / MAX_ITER;
      if(done > 1.0) break;
      temp = temp0*(1-done);
    }

    
    i32 ty = rng.random32(70);
    
    if(ty <= 19) {
      i32 c = rng.random32(nc);
      i32 i = rng.random32(n);
      if(S.in_cycle[c][i]) continue;
      i32 j = rng.random32(S.ncycle[c]);
      i32 a = S.cycle[c][j];
      i32 b = S.cycle[c][j+1];
      set.reset();
      S.rem_seg(seg_idf(c,a,b));
      S.blocking(set, seg_idf(c,i,a), 2);
      if(set.nsegs > 2) {
        S.add_seg(seg_idf(c,a,b));
        continue;
      }
      S.blocking(set, seg_idf(c,i,b), 2);
      if(set.nsegs != 0 && set.nsegs != 2) {
        S.add_seg(seg_idf(c,a,b));
        continue;
      }
      if(set.nsegs == 2) {
        auto id1 = set.segs[0], id2 = set.segs[1];
        auto const &s1 = segs[id1];
        auto const &s2 = segs[id2];
        if(s1.col == c || s1.col != s2.col || S.ncycle[s1.col] <= 5) {
          S.add_seg(seg_idf(c,a,b));
          continue;
        }
        i32 x,y;
        bool ok = 0;
        if(s1.a == s2.a || s1.a == s2.b) { x = s1.a; ok = 1; }
        if(s1.b == s2.a || s1.b == s2.b) { x = s1.b; ok = 1; }
        if(!ok) {
          if(S.has_seg(seg_idf(s1.col,s1.a,s2.a))) { x=s1.a; y=s2.a; ok = 1; }
          else if(S.has_seg(seg_idf(s1.col,s1.a,s2.b))) { x=s1.a; y=s2.b; ok = 1; }
          else if(S.has_seg(seg_idf(s1.col,s1.b,s2.a))) { x=s1.b; y=s2.a; ok = 1; }
          else if(S.has_seg(seg_idf(s1.col,s1.b,s2.b))) { x=s1.b; y=s2.b; ok = 1; }
          if(!ok
             || !S.can_add_segment(seg_id[s1.col][s1.a^s1.b^x][s2.a^s2.b^y])
             || seg_intersects[seg_id[s1.col][s1.a^s1.b^x][s2.a^s2.b^y]][seg_id[c][i][a]]
             || seg_intersects[seg_id[s1.col][s1.a^s1.b^x][s2.a^s2.b^y]][seg_id[c][i][b]])
            {
              S.add_seg(seg_id[c][a][b]);
              continue;
            }

          f64 delta =
            - logs[2*S.sum_cycle[c]*S.ncycle[c]]
            - logs[2*S.sum_cycle[s1.col]*S.ncycle[s1.col]]
            + logs[2*(S.sum_cycle[c] + stars[c][i].v)*(S.ncycle[c]+1)]
            + logs[2*(S.sum_cycle[s1.col] - stars[s1.col][x].v - stars[s1.col][y].v)
                   *(S.ncycle[s1.col]-2)];

          if(!enable_trans || (delta < 0 && -delta > temp * rng.randomDouble())) { // TODO SA
            S.add_seg(seg_id[c][a][b]);
            continue;
          }

          S.sum_cycle[c] += stars[c][i].v;
          S.sum_cycle[s1.col] -= stars[s1.col][x].v;
          S.sum_cycle[s1.col] -= stars[s1.col][y].v;
          
          S.rem_seg(id1);
          S.rem_seg(id2);
          S.rem_seg(seg_idf(s1.col,x,y));
          S.add_seg(seg_idf(s1.col,s1.a^s1.b^x,s2.a^s2.b^y));
          S.add_seg(seg_idf(c,i,a));
          S.add_seg(seg_idf(c,i,b));
          FORD(k, S.ncycle[c], j+1) S.cycle[c][k+1] = S.cycle[c][k];
          S.cycle[c][j+1] = i;
          S.ncycle[c] += 1;
          S.in_cycle[c][i] = 1;
        
          { i32 z = 0; while(S.cycle[s1.col][z] != x) z += 1;
            FORU(k, z, S.ncycle[s1.col]) S.cycle[s1.col][k] = S.cycle[s1.col][k+1];
            S.ncycle[s1.col] -= 1;
          }
          { i32 z = 0; while(S.cycle[s1.col][z] != y) z += 1;
            FORU(k, z, S.ncycle[s1.col]) S.cycle[s1.col][k] = S.cycle[s1.col][k+1];
            S.ncycle[s1.col] -= 1;
          }
          S.cycle[s1.col][S.ncycle[s1.col]] = S.cycle[s1.col][0];
          S.in_cycle[s1.col][x] = 0;
          S.in_cycle[s1.col][y] = 0;
        }else{
          if(!S.can_add_segment(seg_id[s1.col][s1.a^s1.b^x][s2.a^s2.b^x])
             || seg_intersects[seg_id[s1.col][s1.a^s1.b^x][s2.a^s2.b^x]][seg_id[c][i][a]]
             || seg_intersects[seg_id[s1.col][s1.a^s1.b^x][s2.a^s2.b^x]][seg_id[c][i][b]] )
            {
              S.add_seg(seg_id[c][a][b]);
              continue;
            }
          f64 delta =
            - logs[2*S.sum_cycle[c]*S.ncycle[c]]
            - logs[2*S.sum_cycle[s1.col]*S.ncycle[s1.col]]
            + logs[2*(S.sum_cycle[c] + stars[c][i].v)*(S.ncycle[c]+1)]
            + logs[2*(S.sum_cycle[s1.col] - stars[s1.col][x].v)*(S.ncycle[s1.col]-1)];
        
          if(!enable_trans || (delta < 0 && -delta > temp * rng.randomDouble())) { // TODO SA
            S.add_seg(seg_idf(c,a,b));
            continue;
          }
        
          S.sum_cycle[c] += stars[c][i].v;
          S.sum_cycle[s1.col] -= stars[s1.col][x].v;
        
          S.rem_seg(id1);
          S.rem_seg(id2);
          S.add_seg(seg_idf(s1.col,s1.a^s1.b^x,s2.a^s2.b^x));
          S.add_seg(seg_idf(c,i,a));
          S.add_seg(seg_idf(c,i,b));
          FORD(k, S.ncycle[c], j+1) S.cycle[c][k+1] = S.cycle[c][k];
          S.cycle[c][j+1] = i;
          S.ncycle[c] += 1;
          S.in_cycle[c][i] = 1;
        
          y = 0; while(S.cycle[s1.col][y] != x) y += 1;
          FORU(k, y, S.ncycle[s1.col]) S.cycle[s1.col][k] = S.cycle[s1.col][k+1];
          S.ncycle[s1.col] -= 1;
          S.cycle[s1.col][S.ncycle[s1.col]] = S.cycle[s1.col][0];
          S.in_cycle[s1.col][x] = 0;
        }
      }else{
        S.add_seg(seg_idf(c,i,a));
        S.add_seg(seg_idf(c,i,b));
        FORD(k, S.ncycle[c], j+1) S.cycle[c][k+1] = S.cycle[c][k];
        S.cycle[c][j+1] = i;
        S.ncycle[c] += 1;
        S.in_cycle[c][i] = 1;
        S.sum_cycle[c] += stars[c][i].v;
      }
    }else if(ty <= 69){
      i32 c = rng.random32(nc);
      i32 i = rng.random32(S.ncycle[c]);
      i32 j = rng.random32(S.ncycle[c]);
      i32 k = rng.random32(S.ncycle[c]);
      if(i > j) swap(i,j);
      if(j > k) swap(j,k);
      if(i > j) swap(i,j);
      if(i == j) continue;
      if(j == k) continue;

      i32 a1 = S.cycle[c][i];
      i32 a2 = S.cycle[c][i+1];
      i32 b1 = S.cycle[c][j];
      i32 b2 = S.cycle[c][j+1];
      i32 c1 = S.cycle[c][k];
      i32 c2 = S.cycle[c][k+1];
      
      S.rem_seg(seg_idf(c,a1,a2));
      S.rem_seg(seg_idf(c,b1,b2));
      S.rem_seg(seg_idf(c,c1,c2));

      if(!S.can_add_segment(seg_idf(c,a1,b2))) {
        S.add_seg(seg_idf(c,a1,a2));
        S.add_seg(seg_idf(c,b1,b2));
        S.add_seg(seg_idf(c,c1,c2));
        continue;
      }
      S.add_seg(seg_idf(c,a1,b2));

      if(!S.can_add_segment(seg_idf(c,b1,c2))) {
        S.rem_seg(seg_idf(c,a1,b2));
        S.add_seg(seg_idf(c,a1,a2));
        S.add_seg(seg_idf(c,b1,b2));
        S.add_seg(seg_idf(c,c1,c2));
        continue;
      }
      S.add_seg(seg_idf(c,b1,c2));

      if(!S.can_add_segment(seg_idf(c,c1,a2))) {
        S.rem_seg(seg_idf(c,a1,b2));
        S.rem_seg(seg_idf(c,b1,c2));
        S.add_seg(seg_idf(c,a1,a2));
        S.add_seg(seg_idf(c,b1,b2));
        S.add_seg(seg_idf(c,c1,c2));
        continue;
      }
      S.add_seg(seg_idf(c,c1,a2));

      rotate(S.cycle[c]+i+1, S.cycle[c]+j+1, S.cycle[c]+k+1);
    }else{
      i32 c = rng.random32(nc);
      if(S.ncycle[c] <= 3) continue;
      i32 i = rng.random32(S.ncycle[c]-1);
      i32 x = S.cycle[c][i];
      i32 y = S.cycle[c][i+1];
      i32 z = S.cycle[c][i+2];
      if(!S.can_add_segment(seg_idf(c,x,z))) continue;

      f64 delta =
        - logs[2*S.sum_cycle[c]*S.ncycle[c]]
        + logs[2*(S.sum_cycle[c] - stars[c][y].v)*(S.ncycle[c]-1)];
      
      if(!enable_trans || (delta < 0 && -delta > temp * rng.randomDouble())) { 
        continue;
      }

      S.sum_cycle[c] -= stars[c][y].v;
        
      S.rem_seg(seg_idf(c,x,y));
      S.rem_seg(seg_idf(c,y,z));
      S.add_seg(seg_idf(c,x,z));
      S.ncycle[c] -= 1;
      S.in_cycle[c][y] = 0;
      FORU(k, i+1, S.ncycle[c]) S.cycle[c][k] = S.cycle[c][k+1];
      S.cycle[c][S.ncycle[c]] = S.cycle[c][0];
    }

    if(S.value() > best.value()) {
      best = S;
    }
  }
  debug(niter);
}

struct state_climb {
  i32 nsegs;
  i32 segs[5*N*C+2];

  i32 seg_where[MAX_SEGS];

  i32 nksegs;
  i32 ksegs[5*N*C+2];
  
  i32 deg[N*C];
  i32 uf[N*C], uf_v[N*C], uf_sz[N*C], uf_cycle[N*C];

  i32 nlinks[C][N];
  i32 links[C][N][2];

  i32 npath;
  i32 path[N+1];
  
  void reset(){
    nsegs = 0;
  }

  i32 find(i32 i) {
    return uf[i]==i?i:uf[i]=find(uf[i]);
  }

  void compute_links() {
    FOR(i, nc) FOR(j, n) nlinks[i][j] = 0;
    FOR(i, nsegs) {
      auto const& s = ::segs[segs[i]];
      links[s.col][s.a][nlinks[s.col][s.a]++] = s.b;
      links[s.col][s.b][nlinks[s.col][s.b]++] = s.a;
    }
  }

  void compute_where() {
    FOR(i, nsegs) seg_where[segs[i]] = i;
  }
  bool can_add_segment(i32 id) {
    FOR(i, nsegs) if(seg_intersects[id][segs[i]]) return false;
    return true;
  }

  void add_seg(i32 id) {
    seg_where[id] = nsegs;
    segs[nsegs++] = id;
  }

  void rem_seg(i32 id) {
    i32 i = seg_where[id];
    swap(segs[i], segs[nsegs-1]);
    nsegs -= 1;
    seg_where[segs[i]] = i;
    seg_where[id] = -1;
  }

  void get_path(i32 c, i32 i){
    npath = 0;
    if(nlinks[c][i] == 0) return;
    path[npath++] = i;
    i32 dir = rng.random32(nlinks[c][i]);
    i32 last = i;
    i32 j = links[c][i][dir];
    path[npath++] = j;
    for(i32 k = j; nlinks[c][k] != 1; ) {
      auto tmp = k;
      k = links[c][k][0]^links[c][k][1] ^ last;
      last = tmp;
      path[npath++] = k;
      if(k == i) break;
    }
  }
  
  bool unite(i32 i, i32 j) {
    i=find(i);j=find(j);
    if(i == j) {
      uf_cycle[i] = 1;
      return false;
    }
    uf[i] = j;
    uf_v[j] += uf_v[i];
    uf_sz[j] += uf_sz[i];
    return true;
  }
  
  void eval_reset() {
    nksegs = 0;
    FOR(i, nc) FOR(j, n) {
      i32 id = i*n+j;
      uf[id] = id;
      uf_v[id] = stars[i][j].v;
      uf_sz[id] = 1;
      uf_cycle[id] = 0;
      deg[id] = 0;
    }
  }

  void eval_single(i32 seg_id) {
    FOR(j, nksegs) if(seg_intersects[seg_id][ksegs[j]]) return;
    auto const& s = ::segs[seg_id];
    auto id1 = s.col*n+s.a;
    auto id2 = s.col*n+s.b;
    if(deg[id1] == 2 || deg[id2] == 2) return;
    deg[id1] += 1;
    deg[id2] += 1;
    unite(id1,id2);
    ksegs[nksegs++] = seg_id;
  }

  void eval_list() {
    FOR(i, nsegs) eval_single(segs[i]);
  }
  
  f64 eval_end() {
    f64 score = 0;
    FOR(i, nc) {
      i32 sum = 0;
      FOR(j, n) {
        i32 id = i*n+j;
        if(uf[id] == id && deg[id] > 0) {
          sum += (uf_cycle[id]?2:1)*uf_v[id]*uf_sz[id];
        }
      }
      score += logs[sum];
    }
    return score;
  }
};

void solve_climb(state_cycles const& cycles, f64 tl) {
  f64 time0 = TIMER.elapsed();
  
  state_climb S;
  S.reset();
  f64 score = 0;

  auto best_state = S;
  f64 best_score = score;

  f64 temp0 = 0.001;
  f64 temp1 = 0.00001;
  f64 done = 0.0;
  f64 temp = temp0;
  
  auto accept = [&]() -> bool {
    auto new_score = S.eval_end();
    if(new_score > best_score) {
      best_score = new_score;
      best_state = S;
      best_state.nsegs = best_state.nksegs;
      FOR(i, best_state.nsegs) best_state.segs[i] = best_state.ksegs[i];
      debug(best_score + 1e-9, exp(best_score/nc), S.nksegs);
    }
    auto delta = score-new_score;
    if(delta <= 0 || delta <= temp * rng.randomDouble()) {
      score = new_score;
      S.nsegs = S.nksegs;
      FOR(i, S.nsegs) S.segs[i] = S.ksegs[i];
      return true;
    }else{
      return false;
    }
  };
  
  S.eval_reset(); 
  FOR(i, nc) FOR(j, cycles.ncycle[i]) {
    auto a = cycles.cycle[i][j];
    auto b = cycles.cycle[i][j+1];
    S.eval_single(seg_idf(i,a,b));
  }
  accept();
    
  i64 niter = 0;
  while(1) {
    niter += 1;
    if(niter%128 == 0) {
      done = (TIMER.elapsed()-time0) / tl;
      if(done > 1.0) break;
      temp = temp0 * pow(temp1/temp0, done);
    }

    i32 ty = rng.random32(31);
    
    if(ty <= 9) {
      i32 id = rng.random32(nsegs);
      rng.shuffle(S.segs, S.segs+S.nsegs);
      S.eval_reset();
      S.eval_single(id);
      S.eval_list();
      
      accept();
    } else if(ty <= 19) {
      if(S.nsegs == 0) continue;
      rng.shuffle(S.segs, S.segs+S.nsegs);
      i32 id1 = S.segs[rng.random32(S.nsegs)];
      auto const& s = segs[id1];
      i32 c = rng.random32(n);
      if(s.a == c || s.b == c) continue;
      S.eval_reset();
      S.eval_single(seg_idf(s.col,s.a,c));
      S.eval_single(seg_idf(s.col,s.b,c));
      S.eval_list();
      
      accept();
    }else if(ty <= 29) {
      S.compute_links();
      i32 c = rng.random32(nc);
      i32 a = rng.random32(n);
      S.get_path(c,a);
      if(S.npath <= 3) continue;
      i32 i = rng.random32(S.npath-1);
      i32 j = rng.random32(S.npath-1);
      i32 k = rng.random32(S.npath-1);
      if(i > j) swap(i,j);
      if(j > k) swap(j,k);
      if(i > j) swap(i,j);
      if(i == j) continue;
      if(j == k) continue;

      i32 a1 = S.path[i];
      i32 a2 = S.path[i+1];
      i32 b1 = S.path[j];
      i32 b2 = S.path[j+1];
      i32 c1 = S.path[k];
      i32 c2 = S.path[k+1];

      S.compute_where();

      S.rem_seg(seg_idf(c,a1,a2));
      S.rem_seg(seg_idf(c,b1,b2));
      S.rem_seg(seg_idf(c,c1,c2));

      if(!S.can_add_segment(seg_idf(c,a1,b2))) {
        S.add_seg(seg_idf(c,a1,a2));
        S.add_seg(seg_idf(c,b1,b2));
        S.add_seg(seg_idf(c,c1,c2));
        continue;
      }
      S.add_seg(seg_idf(c,a1,b2));

      if(!S.can_add_segment(seg_idf(c,b1,c2))) {
        S.rem_seg(seg_idf(c,a1,b2));
        S.add_seg(seg_idf(c,a1,a2));
        S.add_seg(seg_idf(c,b1,b2));
        S.add_seg(seg_idf(c,c1,c2));
        continue;
      }
      S.add_seg(seg_idf(c,b1,c2));

      if(!S.can_add_segment(seg_idf(c,c1,a2))) {
        S.rem_seg(seg_idf(c,a1,b2));
        S.rem_seg(seg_idf(c,b1,c2));
        S.add_seg(seg_idf(c,a1,a2));
        S.add_seg(seg_idf(c,b1,b2));
        S.add_seg(seg_idf(c,c1,c2));
        continue;
      }
      S.add_seg(seg_idf(c,c1,a2));

      S.eval_reset(); S.eval_list(); 
      accept();
    }else if(ty == 30) {
      S.eval_reset();
      S.eval_list();
      FOR(id, nsegs) S.eval_single(rng.random32(nsegs));
      accept();
    }

  }

  { S.eval_reset();
    S.eval_list();
    FOR(id, nsegs) S.eval_single(id);
    accept();
  }

  debug(niter);

  S = best_state;
  vector<vector<vector<i32>>> G(nc, vector<vector<i32>>(n));
  FOR(i, S.nsegs) {
    auto id = S.segs[i];
    auto const& s = segs[id];
    G[s.col][s.a].eb(s.b);
    G[s.col][s.b].eb(s.a);
  }
  vector<vector<i32>> E(nc, vector<i32>(n,0));
  vector<vector<i32>> ans;
  FOR(i, nc) FOR(j, n) if(!E[i][j] && G[i][j].size() == 1) {
    ans.eb();
    ans.back().eb(stars[i][j].id);
    E[i][j] = 1;
    i32 last = j;
    ans.back().eb(stars[i][G[i][j][0]].id);
    E[i][G[i][j][0]] = 1;
    for(i32 k = G[i][j][0]; G[i][k].size() == 2; ) {
      auto tmp = k;
      k = G[i][k][0]^G[i][k][1]^last;
      last = tmp;

      ans.back().eb(stars[i][k].id);
      E[i][k] = 1;
    }
  }
  FOR(i, nc) FOR(j, n) if(!E[i][j] && G[i][j].size() == 2) {
    ans.eb();
    ans.back().eb(stars[i][j].id);
    E[i][j] = 1;
    i32 last = j;
    ans.back().eb(stars[i][G[i][j][0]].id);
    E[i][G[i][j][0]] = 1;
    for(i32 k = G[i][j][0]; k != j; ) {
      auto tmp = k;
      k = G[i][k][0]^G[i][k][1]^last;
      last = tmp;

      ans.back().eb(stars[i][k].id);
      E[i][k] = 1;
    }
  }
  
  // cerr << ans.size() << endl;
  // for(auto p : ans) {
  //   for(auto x : p) cerr << x << ' ';
  //   cerr << endl;
  // }

  cout << ans.size() << endl;
  for(auto p : ans) {
    cout << p.size() << endl;
    for(auto x : p) cout << x << endl;
  }
}

void solve() {
  init();

  i32 ITERS = 3;
  
  vector<state_cycles> states;

  FOR(iter, ITERS) {
    f64 t0 = TIMER.elapsed();
    f64 tl = (TIME_LIMIT_CYCLES-t0)/(ITERS-iter);
    if(iter == 0) {
      while(TIMER.elapsed() < t0 + tl) {
        states.eb();
        states.back().reset();
        greedy_cycles(states.back(), 0, iter);
        if(states.back().ncycle[0] == 0) states.pop_back();
      }
    }else{
      for(auto& s : states) greedy_cycles(s, (f64)0.9 * tl / states.size(), iter);
    }
    sort(all(states), [&](auto const& a, auto const& b) { return a.value()>b.value(); });
    i32 sz = 1<<max(0,ITERS-1-iter);
    if((i32)states.size() > sz) states.resize(sz);
    for(auto const& s : states) cerr << s.value() << ' ';
    cerr << endl;
    while(states.size() >= 2 && states.back().value() < states[0].value() - 1.0) states.pop_back();
  }
  
  solve_climb(states[0], TIME_LIMIT_END-TIMER.elapsed());
}

int main(int, char**) {
  ios::sync_with_stdio(false); cin.tie(0);
  TIMER.reset();

  read();
  // ofstream os(argv[1]);
  // os << n << ' ' << nc << ' ' << maxv << endl;
  // FOR(i, nc) FOR(j, n) {
  //   os << stars[i][j].pos.x << ' ' << stars[i][j].pos.y << ' ' << i << ' ' << stars[i][j].v << endl;
  // }
  // cout << 0 << endl;
  // FOR(i, n*nc) {
  //   i32 r,c,col,v; cin>>r>>c>>col>>v;
  //   stars[col][nstar[col]++] = { .id = i, .pos = pt(r,c), .v = v };
  // }
  solve();
  cerr << "[DATA] time = " << TIMER.elapsed() << endl;

  return 0;
}
