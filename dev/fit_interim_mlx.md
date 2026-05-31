# MLX-Metal port of the batched-SMC interim sampler

Branch: `interim-mlx`. Experiment: port the batched-vmap SMC (Case A, fixed `x`)
to Apple Metal via MLX, to see whether the M4 Max 40-core GPU can beat the
CPU per-sample-parallel SMC baseline at `pps_z_total = 200`.

The CPU baseline (`fit_interim_SMC_PPS_of_posterior_xz_from_x`, 12 spawn
workers): **253 min for all 8 interims at S=200** — 3× slower than the HMC
refit (~10 min/interim → ~80 min total), which is the figure we wanted to close.

## What was built (on the branch)

- [`python/mlx_pcm.py`](../python/mlx_pcm.py) — partial-credit `loglik` + prior
  in pure MLX: multi-category eta + `take_along_axis` gather + `logsumexp` for
  the categorical lpmf (**no scatter** — that lives only in the ordered-prob
  generated quantities, which we keep on JAX/CPU for the final ratio). The
  prior matches the model exactly (`Normal(0, 1/sqrt(1-1/U))` marginal +
  `Normal(0,1)` + `Normal(0,3.5)` + `FoldedStudentT(3,0,1)` with hard-coded
  log-gamma constants).
- [`python/fit_interim_mlx.py`](../python/fit_interim_mlx.py) —
  `fit_interim_SMC_batched_GPU_PPS_of_posterior_xz_from_x`. All `S` future
  samples live together as one `(S, K, D)` particle tensor on the GPU, moved
  by a vmapped MALA kernel. **Per-sample β masking** advances each `β_s`
  independently (bisection in vectorised numpy on a pulled `llz` matrix), so
  easy samples finish in a few temperatures — the over-tempering that hurt the
  JAX-CPU batched prototype is gone. Final per-sample labels go through the
  existing JAX/CPU `_ratio_per_draw_from_params` (uses the ordered-prob
  generated quantity, which has scatter).
- [`scripts-py/Ukraine_interim_analysis_SMC_mlx.py`](../scripts-py/Ukraine_interim_analysis_SMC_mlx.py)
  — full all-interim analysis mirroring the other SMC scripts.

The `main` branch is untouched (no commits made; the MLX env install is
additive in `.pixi/`, not in any tracked manifest).

## Measurements

### Probe (synthetic, single MLX grad eval)

```
default device: Device(gpu, 0)
single: logpost finite, grad_shape=(152,)
vmap-grad over 25 600 particles: ~700 ms  (multi-cat PCM loglik + real prior, double-vmap)
[mlx-cpu] same: ~17 s   → GPU ≈ 46× MLX-CPU per eval
```

The Metal GPU + MLX autodiff handle our op pattern (gather, matmul, cumsum,
logsumexp) at full 200×128 budget in one fused kernel.

### Smoke test at `pps_z_total = 12` (real Ukraine data, all 8 interims)

```
interim  mins  n_temps   (per-sample β; no over-tempering)
1        6.85    39
2        4.83    28
3        3.47    19
4        2.37    12
5        1.59     7
6        1.43     6
7        1.44     6
8        0.97     3
                    total 22.95 min ; p_h1_xz valid (1920 rows in [0,1])
```

Comparable to CPU per-s SMC (18.4 min) and correctness validated.

### `pps_z_total = 200` (full HMC-comparable budget)

- First run: died with `Metal Command buffer execution failed
  (kIOGPUCommandBufferCallbackErrorPageFault)` — 20 lazy MALA steps × 25 600
  particles built too deep a Metal graph. Fix: `mx.eval(u)` inside the move
  loop to keep the graph shallow.
- Two subsequent runs hung in `STAT=U` after the laptop slept during travel —
  the Metal GPU context dies on suspend and the Python process is stuck
  on an unrecoverable GPU op. Solved by wrapping in `caffeinate -i` (idle
  sleep only — lid-close still triggers).
- The wrapped run reached **~100 min wall-clock with no output** before being
  killed for taking too long. Process was healthy (state R, CPU oscillating
  85% ↔ 7%, caffeinate attached) — orchestration overhead simply dominates at
  S=200 on this machine.

## Why the GPU per-eval win did not translate at S=200

- **Sync overhead per tempering step.** Every step pulls `llz` to numpy, runs
  the per-sample bisection + systematic resample, gathers particles, and
  pushes back to MLX. With shapes `(S=200, K=128, D=133)` that's `~14 MB`
  round-tripping every step.
- **MALA loop depth.** `n_move_steps = 20` × `(grad + grad + lp + lp)` over
  25 600 particles is heavy even on the GPU.
- **Per-sample β masking** kills over-tempering (good for total work) but
  preserves long tempering schedules for the worst samples — interim 1
  still needs ~40 temperatures, and the move dominates each one.
- **MLX-Metal maturity.** Stalls / page faults under deep lazy graphs;
  Metal context dies on macOS suspend. Production-fragile on a laptop that
  sleeps.

## Net

Per-eval the MLX-Metal kernel is decisively faster than CPU. **End-to-end at
S=200 it did not beat the 12-worker CPU per-s SMC** on this machine — the
gain is eaten by orchestration sync, the deep MALA loop, and the
suspend-kills-GPU-context issue. The CPU per-s SMC stays the production path.

## If revisiting

Highest-impact MLX optimizations, roughly in order:

1. **GPU-side per-s resample.** Replace the per-step numpy↔MLX round-trip
   with `mx.take_along_axis` over the particle axis so the particle tensor
   never leaves the GPU. Eliminates the 14 MB transfer + cumulative sync.
2. **Drop `n_move_steps` 20 → 6–8.** 20 MALA sweeps per temperature is
   overkill to decorrelate; the smoke test already shows the schedule is
   determined by `Δβ` adaptation, not move quality.
3. **Chunk `S` into ~3 groups** (e.g. `S=64` × 3 passes) so each MLX kernel
   stays inside the Metal scheduler's comfortable size, and the lazy graph
   never grows deep enough to page-fault.
4. **Keep the run short.** Per-interim ≤ ~30 min so the laptop doesn't
   suspend mid-run (or always wrap in `caffeinate -di` and leave the lid
   open).
5. **Defer the final ratio computation** off the critical path — batch the
   per-sample `_ratio_per_draw_from_params` calls.

If those land, the headline measurement to re-take is: **wall-clock minutes
per interim at S=200**, against the CPU 31.7 min/interim and the HMC
~10 min/interim.
