// reduction_functions.stan — reduction catalog for rho_by_cat_qu (per-participant paired + group-mean).
// Included inside the model's `functions` block (alongside credit_model_functions.stan).
// Each reduction takes BOTH representations and a by_unit flag:
//   p1/p2 = PAIRED per-participant category probs (rows = units present at both questions, cols = K);
//   g1/g2 = the two GROUP MEANS of unit_summary (over all participants at cmp1 / cmp2).
// by_unit = 1 uses p1/p2 (g1/g2 ignored); by_unit = 0 uses g1/g2 (p1/p2 ignored). Stan has no null, so the
// caller passes both and the unused side is ignored. Reduction-specific constants (e.g. SCR cut/floor) are
// defined INTERNALLY here, NOT exposed as data.

// GMFR = 2^(Δ mean log2 level), E[y] = sum_k (k-1) p[k] over the 0-INDEXED log2-titre ladder (category k is
// 1-indexed in Stan but carries ladder value k-1: 1:5->0, 1:10->1, 1:40->3). by_unit: Δ = mean_u (E[y2_u]-
// E[y1_u]) from the paired matrices p1/p2 (g1/g2 unused); else Δ = g2 - g1 from the group means (p1/p2 unused).
real reduce_gmfr(matrix p1, matrix p2, real g1, real g2, int by_unit) {
  if (by_unit) {
    int M = rows(p1);
    int K = cols(p1);
    if (M == 0) {
      return 0;
    }
    vector[K] cats = linspaced_vector(K, 0, K - 1); // ladder values [0, 1, ..., K-1]
    return pow(2.0, mean((p2 - p1) * cats));
  }
  return pow(2.0, g2 - g1);
}

// SCR: mean_u P(seroconversion_u) = mean_u P(d >= scr_cut), d = (post - pre) - 1{pre == scr_floor}
//   (the '-1' censoring limb for baseline < 1:10). scr_cut = 2 (>= 4-fold = +2 cats); scr_floor = 1 (<1:10).
// G[a,b] = sum_u p1[u,a] p2[u,b] (cross-moment, K x K); the rate = sum_u P(d>=cut) = sum_{a,b: d>=cut} G[a,b],
// so we accumulate the qualifying G entries directly in the threshold loop (no mask matrix, no M x K matmul).
real reduce_scr(matrix p1, matrix p2) {
  int M = rows(p1);
  int K = cols(p1);
  int scr_cut = 2;
  int scr_floor = 1;
  if (M == 0) {
    return 0;
  }
  matrix[K, K] G = p1' * p2;
  real s = 0;
  for (a in 1 : K) {
    for (b in 1 : K) {
      int d = (b - a) - (a == scr_floor ? 1 : 0);
      if (d >= scr_cut) {
        s += G[a, b];
      }
    }
  }
  return s / M;
}

// Per-participant scalar summary of a (participants x K) prob matrix, mirroring the R get_endpoints pipeline
// (y is 0-INDEXED: category k carries value k-1). thr == 0 -> expected level E[y] = sum_k (k-1) p[k] (out-of-7
// means); thr >= 1 -> P(y >= thr) = sum_{k: k-1 >= thr} p[k] (categorical high-category aggregation).
vector unit_summary(matrix p, int thr) {
  int K = cols(p);
  if (thr == 0) {
    return p * linspaced_vector(K, 0, K - 1); // E[y], 0-indexed values
  }
  vector[K] mask = rep_vector(0., K);
  for (k in 1 : K) {
    if (k - 1 >= thr) {
      mask[k] = 1;
    }
  }
  return p * mask; // P(y >= thr)
}

// Relative % change = (s1 - s2)/s1 * 100 (baseline - endline over baseline; LOWER-is-better reads positive --
// flip in R for higher-is-better). by_unit: mean over participants of the per-unit ratio (s = unit_summary(.,
// thr) of the paired matrices; softmax probs > 0 so s1 > 0, no divide-by-zero); else ratio of the group means
// g1/g2 (already summarised, thr unused). p1/p2 unused when !by_unit, g1/g2 unused when by_unit.
real reduce_relpct(matrix p1, matrix p2, real g1, real g2, int by_unit, int thr) {
  if (by_unit) {
    int M = rows(p1);
    if (M == 0) {
      return 0;
    }
    vector[M] s1 = unit_summary(p1, thr);
    vector[M] s2 = unit_summary(p2, thr);
    return mean((s1 - s2) ./ s1) * 100;
  }
  return (g1 - g2) / g1 * 100;
}

// Difference = s1 - s2 (baseline - endline). by_unit: mean over participants (= group diff, LINEAR) of
// unit_summary(., thr) of the paired matrices; else g1 - g2 from the group means. Unused args as above.
real reduce_diff(matrix p1, matrix p2, real g1, real g2, int by_unit, int thr) {
  if (by_unit) {
    int M = rows(p1);
    if (M == 0) {
      return 0;
    }
    return mean(unit_summary(p1, thr) - unit_summary(p2, thr));
  }
  return g1 - g2;
}

// Unified reduction catalog. One red_fn id = ONE specific reduction. p1/p2 are the PAIRED per-participant
// prob matrices (units present at BOTH cmp questions); g1/g2 are the GROUP MEANS of unit_summary over all
// participants at cmp1 / cmp2 (no pairing). thr is the per-comparison level/threshold selector (0 = E[y],
// >=1 = P(y>=thr)); reductions that don't use it ignore it.
//   per-participant (paired):
//     1 GMFR   2^(mean_u (E[y2]-E[y1]))
//     2 SCR    mean_u P(d>=cut)
//     3 relpct mean_u (s1_u - s2_u)/s1_u * 100        (within-participant rel % change)
//     4 diff   mean_u (s1_u - s2_u)                   (within-participant difference)
//   group-mean (no pairing):
//     5 base   g1                                     (baseline group level)
//     6 end    g2                                     (endline group level)
//     7 relpct (g1 - g2)/g1 * 100                     (ratio-of-group-means)
//     8 diff   g1 - g2                                (difference-of-group-means)
//     9 GMFR   2^(g2 - g1)                            (group-mean fold rise)
real reduce_dispatch(int red_fn, matrix p1, matrix p2, real g1, real g2, int thr) {
  if (red_fn == 1) {
    return reduce_gmfr(p1, p2, g1, g2, 1); // per-participant
  }
  if (red_fn == 2) {
    return reduce_scr(p1, p2); // per-participant only (joint)
  }
  if (red_fn == 3) {
    return reduce_relpct(p1, p2, g1, g2, 1, thr); // per-participant
  }
  if (red_fn == 4) {
    return reduce_diff(p1, p2, g1, g2, 1, thr); // per-participant
  }
  if (red_fn == 5) {
    return g1; // baseline group level
  }
  if (red_fn == 6) {
    return g2; // endline group level
  }
  if (red_fn == 7) {
    return reduce_relpct(p1, p2, g1, g2, 0, thr); // group-mean (ratio-of-means)
  }
  if (red_fn == 8) {
    return reduce_diff(p1, p2, g1, g2, 0, thr); // group-mean (diff-of-means)
  }
  if (red_fn == 9) {
    return reduce_gmfr(p1, p2, g1, g2, 0); // group-mean
  }
  reject("reduce_dispatch: unknown red_fn = ", red_fn);
}

// Orchestrator: for each comparison jc, build BOTH the paired per-participant matrices (units present at both
// questions) and the two group means (over ALL participants at cmp1 / cmp2), then dispatch each reduction r.
// prob_obs[c][n, 1:K[c]] = category probs for observation n of category type c. Returns R x Qcmp.
matrix compute_rho(array[] matrix prob_obs, int R, int Qcmp, int U,
                   array[,] int question_cmp, array[] int cmp_cat_type,
                   array[] int red_fn, array[] int cmp_reduce_thr,
                   array[] int K, array[] int N, array[] int cat_obs_start,
                   array[] int unit_of_obs, array[] int question_of_obs) {
  matrix[R, Qcmp] out = rep_matrix(0., R, Qcmp);
  // decide ONCE (red_fn is the same for every jc) which representations are actually requested:
  // ids 1-4 (per-participant) need the paired matrices; ids 5-9 (levels / group-mean) need the group means
  int need_unit = 0;
  int need_group = 0;
  for (r in 1 : R) {
    if (red_fn[r] >= 1 && red_fn[r] <= 4) {
      need_unit = 1;
    }
    if (red_fn[r] >= 5 && red_fn[r] <= 9) {
      need_group = 1;
    }
  }
  for (jc in 1 : Qcmp) {
    int c = cmp_cat_type[jc];
    int Kc = K[c];
    int thr = cmp_reduce_thr[jc];
    int q1 = question_cmp[jc, 1];
    int q2 = question_cmp[jc, 2];
    int s = cat_obs_start[c];
    real g1 = 0;
    real g2 = 0;
    array[U] int row1 = rep_array(0, U);
    array[U] int row2 = rep_array(0, U);
    int M = 0;
    // group means (only if a group-mean reduction is requested): mean over ALL participants at q1 / q2
    if (need_group > 0) {
      int n1 = 0;
      int n2 = 0;
      for (n in 1 : N[c]) {
        int obs = s + n - 1;
        if (question_of_obs[obs] == q1) {
          n1 += 1;
        }
        if (question_of_obs[obs] == q2) {
          n2 += 1;
        }
      }
      if (n1 > 0) {
        matrix[n1, Kc] m1 = rep_matrix(0., n1, Kc);
        int m = 0;
        for (n in 1 : N[c]) {
          int obs = s + n - 1;
          if (question_of_obs[obs] == q1) {
            m += 1;
            m1[m,  : ] = prob_obs[c][n, 1 : Kc];
          }
        }
        g1 = mean(unit_summary(m1, thr));
      }
      if (n2 > 0) {
        matrix[n2, Kc] m2 = rep_matrix(0., n2, Kc);
        int m = 0;
        for (n in 1 : N[c]) {
          int obs = s + n - 1;
          if (question_of_obs[obs] == q2) {
            m += 1;
            m2[m,  : ] = prob_obs[c][n, 1 : Kc];
          }
        }
        g2 = mean(unit_summary(m2, thr));
      }
    }
    // paired per-participant matrices (only if a per-participant reduction is requested): units at BOTH q1, q2
    if (need_unit > 0) {
      for (n in 1 : N[c]) {
        int obs = s + n - 1;
        int u = unit_of_obs[obs];
        if (question_of_obs[obs] == q1) {
          row1[u] = n;
        }
        if (question_of_obs[obs] == q2) {
          row2[u] = n;
        }
      }
      for (u in 1 : U) {
        if (row1[u] > 0 && row2[u] > 0) {
          M += 1;
        }
      }
    }
    {
      matrix[M, Kc] p1 = rep_matrix(0., M, Kc);
      matrix[M, Kc] p2 = rep_matrix(0., M, Kc);
      int m = 0;
      for (u in 1 : U) {
        if (row1[u] > 0 && row2[u] > 0) {
          m += 1;
          p1[m,  : ] = prob_obs[c][row1[u], 1 : Kc];
          p2[m,  : ] = prob_obs[c][row2[u], 1 : Kc];
        }
      }
      for (r in 1 : R) {
        out[r, jc] = reduce_dispatch(red_fn[r], p1, p2, g1, g2, thr);
      }
    }
  }
  return out;
}
