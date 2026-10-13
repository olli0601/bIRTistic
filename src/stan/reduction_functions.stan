// reduction_functions.stan — per-participant PAIRED reduction functions for participant_reduced_rho.
// Included inside the model's `functions` block (alongside credit_model_functions.stan).
// Each reduction takes the SET of participants present at BOTH paired questions:
//   p1 = baseline category probs (rows = participants, cols = K categories), p2 = endline (same rows/cols).
// It averages the per-participant functional over participants AND applies the final transform (finalise).
// Reduction-specific constants (e.g. SCR cut/floor) are defined INTERNALLY here, NOT exposed as data.

// GMFR: 2^( mean_u (E[y2_u] - E[y1_u]) ), E[y] = sum_k (k-1) p[k] over the 0-INDEXED log2-titre ladder
// (category k is 1-indexed in Stan but carries ladder value k-1: 1:5->0, 1:10->1, 1:40->3).
real reduce_gmfr(matrix p1, matrix p2) {
  int M = rows(p1);
  int K = cols(p1);
  if (M == 0) {
    return 0;
  }
  vector[K] cats = linspaced_vector(K, 0, K - 1); // ladder values [0, 1, ..., K-1]
  return pow(2.0, mean((p2 - p1) * cats)); // E[y2]-E[y1] per unit, mean over units, then 2^(.)
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

// Relative % change, mean over participants: per unit (s1 - s2) / s1 * 100 (baseline - endline over baseline,
// so LOWER-is-better items read positive -- match reduce_diff's orientation; flip in R for higher-is-better).
// s = unit_summary(., thr). softmax probs are strictly positive so s1 > 0 (no divide-by-zero).
real reduce_relpct(matrix p1, matrix p2, int thr) {
  int M = rows(p1);
  if (M == 0) {
    return 0;
  }
  vector[M] s1 = unit_summary(p1, thr);
  vector[M] s2 = unit_summary(p2, thr);
  return mean((s1 - s2) ./ s1) * 100;
}

// Difference, mean over participants: per unit s1 - s2 (baseline - endline). LINEAR, so per-participant mean
// equals the group difference. s = unit_summary(., thr).
real reduce_diff(matrix p1, matrix p2, int thr) {
  int M = rows(p1);
  if (M == 0) {
    return 0;
  }
  return mean(unit_summary(p1, thr) - unit_summary(p2, thr));
}

// Dispatch the reduction by id (Stan has no function-pointer / function-name data). thr is the per-comparison
// threshold for the level/threshold family (ids 3/4); ids 1/2 (GMFR/SCR) ignore it.
real reduce_apply(int red_fn, matrix p1, matrix p2, int thr) {
  if (red_fn == 1) {
    return reduce_gmfr(p1, p2);
  }
  if (red_fn == 2) {
    return reduce_scr(p1, p2);
  }
  if (red_fn == 3) {
    return reduce_relpct(p1, p2, thr);
  }
  if (red_fn == 4) {
    return reduce_diff(p1, p2, thr);
  }
  reject("reduce_apply: unknown red_fn = ", red_fn);
}

// GROUP-MEAN family: reduce two pre-formed group means g1 (cmp1 group) and g2 (cmp2 group). This matches the
// legacy get_endpoints ratio-/diff-of-group-means. For a paired pre/post design this is NOT the right estimand
// (use the per-participant reduce_apply); it is the correct one for a disjoint between-group comparison.
real reduce_group(int red_fn, real g1, real g2) {
  if (red_fn == 1) {
    return pow(2.0, g2 - g1); // GMFR fold rise (g = group-mean log2 level)
  }
  if (red_fn == 2) {
    reject("reduce_group: SCR (red_fn 2) is a joint functional, not defined on group means");
  }
  if (red_fn == 3) {
    return (g1 - g2) / g1 * 100; // rel % change (baseline - endline)/baseline
  }
  if (red_fn == 4) {
    return g1 - g2; // diff (baseline - endline)
  }
  reject("reduce_group: unknown red_fn = ", red_fn);
}

// Orchestrator: for each comparison pair jc, find units present at BOTH question_cmp[jc,1] and [jc,2]
// (within cat type cmp_cat_type[jc]), stack their paired probs into p1/p2 (M x K[c]), then call each
// reduction once. prob_obs[c][n, 1:K[c]] = category probs for observation n of category type c.
matrix compute_reduced_rho(array[] matrix prob_obs, int R, int Qcmp, int U,
                           array[,] int question_cmp, array[] int cmp_cat_type,
                           array[] int red_fn, array[] int cmp_reduce_thr,
                           array[] int K, array[] int N,
                           array[] int cat_obs_start, array[] int unit_of_obs,
                           array[] int question_of_obs) {
  matrix[R, Qcmp] out = rep_matrix(0., R, Qcmp);
  for (jc in 1 : Qcmp) {
    int c = cmp_cat_type[jc];
    int Kc = K[c];
    int q1 = question_cmp[jc, 1];
    int q2 = question_cmp[jc, 2];
    int s = cat_obs_start[c];
    // locate each unit's obs row (1-indexed within cat c) at q1 and q2; 0 = absent
    array[U] int row1 = rep_array(0, U);
    array[U] int row2 = rep_array(0, U);
    int M = 0;
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
    if (M == 0) {
      continue; // out[,jc] stays 0
    }
    // stack the paired probability rows for the both-present units
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
      out[r, jc] = reduce_apply(red_fn[r], p1, p2, cmp_reduce_thr[jc]);
    }
  }
  return out;
}

// GROUP-MEAN levels: for each comparison jc, the mean over ALL participants at question_cmp[jc,1] (col 1) and
// at question_cmp[jc,2] (col 2) of the per-participant summary (E[y] or P(y>=thr)). These are the Baseline /
// Endline group levels. No pairing, so this also serves a disjoint between-group comparison.
matrix compute_group_means(array[] matrix prob_obs, int Qcmp,
                           array[,] int question_cmp, array[] int cmp_cat_type,
                           array[] int cmp_reduce_thr, array[] int K, array[] int N,
                           array[] int cat_obs_start, array[] int question_of_obs) {
  matrix[Qcmp, 2] gm = rep_matrix(0., Qcmp, 2);
  for (jc in 1 : Qcmp) {
    int c = cmp_cat_type[jc];
    int Kc = K[c];
    int thr = cmp_reduce_thr[jc];
    int s = cat_obs_start[c];
    int n1 = 0;
    int n2 = 0;
    for (n in 1 : N[c]) {
      int obs = s + n - 1;
      if (question_of_obs[obs] == question_cmp[jc, 1]) {
        n1 += 1;
      }
      if (question_of_obs[obs] == question_cmp[jc, 2]) {
        n2 += 1;
      }
    }
    if (n1 > 0) {
      // stack the cmp1 group's rows, then group mean = mean over participants of unit_summary
      matrix[n1, Kc] g1 = rep_matrix(0., n1, Kc);
      int m = 0;
      for (n in 1 : N[c]) {
        int obs = s + n - 1;
        if (question_of_obs[obs] == question_cmp[jc, 1]) {
          m += 1;
          g1[m,  : ] = prob_obs[c][n, 1 : Kc];
        }
      }
      gm[jc, 1] = mean(unit_summary(g1, thr));
    }
    if (n2 > 0) {
      matrix[n2, Kc] g2 = rep_matrix(0., n2, Kc);
      int m = 0;
      for (n in 1 : N[c]) {
        int obs = s + n - 1;
        if (question_of_obs[obs] == question_cmp[jc, 2]) {
          m += 1;
          g2[m,  : ] = prob_obs[c][n, 1 : Kc];
        }
      }
      gm[jc, 2] = mean(unit_summary(g2, thr));
    }
  }
  return gm;
}

// GROUP-MEAN reductions: apply reduce_group to the pre-formed group means (gmean[jc,1], gmean[jc,2]) for each
// reduction r and comparison jc. Mirrors the legacy ratio-/diff-of-group-means in get_endpoints.
matrix compute_group_reduced(matrix gmean, int R, int Qcmp, array[] int red_fn) {
  matrix[R, Qcmp] out = rep_matrix(0., R, Qcmp);
  for (jc in 1 : Qcmp) {
    for (r in 1 : R) {
      out[r, jc] = reduce_group(red_fn[r], gmean[jc, 1], gmean[jc, 2]);
    }
  }
  return out;
}
