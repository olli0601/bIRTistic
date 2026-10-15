// reduction_functions.stan — reduction catalog for rho_by_cat_qu (per-participant paired + group-mean).
// Included inside the model's `functions` block (alongside credit_model_functions.stan).
//
// Pooling convention: a comparison names a SET of baseline questions and a SET of endline questions. For each
// participant we average their category probabilities over the questions (in the set) that they answered,
// giving ONE pooled baseline prob vector and one pooled endline prob vector per participant. Reductions then
// act on those pooled per-participant probs:
//   p1/p2 = pooled per-participant prob matrices (rows = the participant set, cols = K) baseline / endline;
//   g1/g2 = the two GROUP MEANS (over participants) of the per-participant summary.
// by_unit = 1 uses p1/p2 (paired participants), by_unit = 0 uses g1/g2; the unused side is ignored.
// obs_of_unit[u, c, q] (transformed data) = local obs index (1..N[c]) of unit u at question q of category c,
// 0 if absent. Reduction constants (e.g. SCR cut/floor) are INTERNAL, not data.

// Per-participant scalar summary of a (participants x K) prob matrix, 0-INDEXED y (category k carries value
// k-1). thr == 0 -> E[y] = sum_k (k-1) p[k]; thr >= 1 -> P(y >= thr) = sum_{k: k-1 >= thr} p[k].
vector summarize_probs_py_participants_for_one_item(matrix p, int thr) {
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

// Pool probs per participant over a set of questions: for each unit in 
// `units`, average prob_obs[c][.,1:Kc] over the `questions` they answered 
//(via obs_of_unit); returns a (|units| x Kc) matrix. One question is just
// that question's probs (no averaging).
matrix average_probs_over_questions_for_each_participant_and_level(array[] matrix prob_obs, int c, int Kc,
                                     array[,,] int obs_of_unit, array[] int units, array[] int questions) {
  int nu = size(units);
  int nq = size(questions);
  matrix[nu, Kc] out = rep_matrix(0., nu, Kc);
  for (i in 1 : nu) {
    int u = units[i];
    row_vector[Kc] acc = rep_row_vector(0., Kc);
    int cnt = 0;
    for (j in 1 : nq) {
      int n = obs_of_unit[u, c, questions[j]];
      if (n > 0) {
        acc += prob_obs[c][n, 1 : Kc];
        cnt += 1;
      }
    }
    if (cnt > 0) {
      out[i,  : ] = acc / cnt;
    }
  }
  return out;
}

// Per-participant summary after pooling over the questions (average probs then summarise). For one question
// this is just summarize_probs_py_participants_for_one_item of that question's probs.
vector summarize_probs_py_participants(array[] matrix prob_obs, int c, int Kc, int thr,
                                       array[,,] int obs_of_unit, array[] int units, array[] int questions) {
  return summarize_probs_py_participants_for_one_item(
           average_probs_over_questions_for_each_participant_and_level(prob_obs, c, Kc, obs_of_unit, units, questions), thr);
}

// GMFR = 2^(Δ mean log2 level), E[y] on the 0-indexed log2 ladder. by_unit: Δ = mean_u (E[y2]-E[y1]) from the
// pooled paired matrices; else Δ = g2 - g1 from the group means.
real compute_rho_gmfr(matrix p1, matrix p2, real g1, real g2, int by_unit) {
  if (by_unit) {
    int M = rows(p1);
    int K = cols(p1);
    if (M == 0) {
      return not_a_number();
    }
    return pow(2.0, mean((p2 - p1) * linspaced_vector(K, 0, K - 1)));
  }
  return pow(2.0, g2 - g1);
}

// SCR = mean_u P(seroconversion) = mean_u P(d >= scr_cut), d = (post - pre) - 1{pre == scr_floor} (scr_cut = 2
// >= 4-fold; scr_floor = 1, baseline < 1:10), from the pooled paired matrices. G[a,b] = sum_u p1[u,a] p2[u,b]
// (K x K cross-moment); rate = sum_{a,b: d>=cut} G[a,b] / M. Joint-only (per-participant); no group form.
real compute_rho_scr(matrix p1, matrix p2) {
  int M = rows(p1);
  int K = cols(p1);
  int scr_cut = 2;
  int scr_floor = 1;
  if (M == 0) {
    return not_a_number();
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

// Relative % change = (s1 - s2)/s1 * 100 (baseline - endline over baseline; LOWER-is-better reads positive --
// flip in R for higher-is-better). by_unit: mean over participants of the per-unit ratio from the pooled paired
// matrices; else from the group means (g1 - g2)/g1.
real compute_rho_relpct(matrix p1, matrix p2, real g1, real g2, int by_unit, int thr) {
  if (by_unit) {
    int M = rows(p1);
    if (M == 0) {
      return not_a_number();
    }
    vector[M] s1 = summarize_probs_py_participants_for_one_item(p1, thr);
    vector[M] s2 = summarize_probs_py_participants_for_one_item(p2, thr);
    return mean((s1 - s2) ./ s1) * 100;
  }
  return (g1 - g2) / g1 * 100;
}

// Difference = s1 - s2 (baseline - endline). by_unit: mean over participants (= group diff, LINEAR) from the
// pooled paired matrices; else g1 - g2 from the group means.
real compute_rho_diff(matrix p1, matrix p2, real g1, real g2, int by_unit, int thr) {
  if (by_unit) {
    int M = rows(p1);
    if (M == 0) {
      return not_a_number();
    }
    return mean(summarize_probs_py_participants_for_one_item(p1, thr)
                - summarize_probs_py_participants_for_one_item(p2, thr));
  }
  return g1 - g2;
}

// Unified catalog. One red_fn id = ONE reduction. p1/p2 = pooled per-participant paired prob matrices,
// g1/g2 = group means (both already pooled over the comparison's question sets).
//   per-participant (paired): 1 GMFR, 2 SCR, 3 relpct, 4 diff
//   group-mean:               5 baseline level g1, 6 endline level g2, 7 relpct, 8 diff, 9 GMFR
real compute_rho_dispatcher(int red_fn, matrix p1, matrix p2, real g1, real g2, int thr) {
  if (red_fn == 1) {
    return compute_rho_gmfr(p1, p2, g1, g2, 1);
  }
  if (red_fn == 2) {
    return compute_rho_scr(p1, p2);
  }
  if (red_fn == 3) {
    return compute_rho_relpct(p1, p2, g1, g2, 1, thr);
  }
  if (red_fn == 4) {
    return compute_rho_diff(p1, p2, g1, g2, 1, thr);
  }
  if (red_fn == 5) {
    return g1;
  }
  if (red_fn == 6) {
    return g2;
  }
  if (red_fn == 7) {
    return compute_rho_relpct(p1, p2, g1, g2, 0, thr);
  }
  if (red_fn == 8) {
    return compute_rho_diff(p1, p2, g1, g2, 0, thr);
  }
  if (red_fn == 9) {
    return compute_rho_gmfr(p1, p2, g1, g2, 0);
  }
  reject("compute_rho_dispatcher: unknown red_fn = ", red_fn);
}

// Orchestrator. For each comparison jc, the baseline question set = non-zero question_cmp[jc,1,:], the endline
// set = non-zero question_cmp[jc,2,:] (a single item = one question per side; an item-group = several). Build,
// via obs_of_unit, the pooled per-participant baseline/endline prob matrices (paired participants = answered
// >=1 question on BOTH sides) and the group-mean levels g1/g2 (over participants with that side), then dispatch
// each reduction. Returns R x Qcmp.
// The per-comparison participant sets Qcmp_units_have_q1 / _q2 / _q12 (unit ids with >=1 answered baseline /
// endline / both question, padded with 0; counts n_q1 / n_q2 / n_q12) are precomputed in transformed data.
matrix compute_rho(array[] matrix prob_obs, int R, int Qcmp,
                   array[,,] int question_cmp, array[] int cmp_cat_type,
                   array[] int red_fn, array[,] int cmp_reduce_thr,
                   array[,,] int obs_of_unit, array[] int K,
                   array[,] int Qcmp_units_have_q1, array[,] int Qcmp_units_have_q2,
                   array[,] int Qcmp_units_have_q12,
                   array[] int n_q1, array[] int n_q2, array[] int n_q12) {
  matrix[R, Qcmp] out = rep_matrix(0., R, Qcmp);
  int Scmp = dims(question_cmp)[3];
  for (jc in 1 : Qcmp) {
    int c = cmp_cat_type[jc];
    int Kc = K[c];
    int nb = n_q1[jc];
    int ne = n_q2[jc];
    int np = n_q12[jc];
    // any per-participant reduction (ids 1-4) -> build the pooled paired matrices (thr-independent) once
    int anyu = 0;
    for (r in 1 : R) {
      if (red_fn[r] < 5) {
        anyu = 1;
      }
    }
    // baseline / endline question sets to pool (non-zero, contiguous)
    int nqb = 0;
    int nqe = 0;
    for (j in 1 : Scmp) {
      if (question_cmp[jc, 1, j] > 0) {
        nqb += 1;
      }
      if (question_cmp[jc, 2, j] > 0) {
        nqe += 1;
      }
    }
    array[nqb] int Qb = question_cmp[jc, 1, 1 : nqb];
    array[nqe] int Qe = question_cmp[jc, 2, 1 : nqe];
    // pooled per-participant paired matrices (participants with BOTH sides); probs only, so thr-independent
    matrix[np, Kc] p1 = rep_matrix(0., np, Kc);
    matrix[np, Kc] p2 = rep_matrix(0., np, Kc);
    if (anyu && np > 0) {
      p1 = average_probs_over_questions_for_each_participant_and_level(
             prob_obs, c, Kc, obs_of_unit, Qcmp_units_have_q12[jc, 1 : np], Qb);
      p2 = average_probs_over_questions_for_each_participant_and_level(
             prob_obs, c, Kc, obs_of_unit, Qcmp_units_have_q12[jc, 1 : np], Qe);
    }
    // each reduction uses its own threshold cmp_reduce_thr[r, jc]; group-mean reductions need g1/g2 at that thr
    for (r in 1 : R) {
      int thr = cmp_reduce_thr[r, jc];
      real g1 = 0;
      real g2 = 0;
      if (red_fn[r] >= 5) {
        if (nb > 0) {
          g1 = mean(summarize_probs_py_participants(prob_obs, c, Kc, thr, obs_of_unit,
                                                    Qcmp_units_have_q1[jc, 1 : nb], Qb));
        }
        if (ne > 0) {
          g2 = mean(summarize_probs_py_participants(prob_obs, c, Kc, thr, obs_of_unit,
                                                    Qcmp_units_have_q2[jc, 1 : ne], Qe));
        }
      }
      out[r, jc] = compute_rho_dispatcher(red_fn[r], p1, p2, g1, g2, thr);
    }
  }
  return out;
}
