"""Static R0 PRCC at equal-width (+/-20 % and +/-50 %) parameter ranges.

The algebraic R0 needs no ODE solves, so this completes the 3x3 rank matrix
(responses {static R0, dynamic peak, dynamic cumulative} x ranges {published,
+/-20 %, +/-50 %}) at negligible cost.

Design pairing: the LHS construction, sample size (n=2000) and seed
(20260816) are identical to tasks A/D/E, so the uniform variates u are paired
sample-for-sample across every analysis in this folder.  R0 is computed with
the same K02 * K20 algebra used by sensitivity.samples_to_params (verified
equal to r0.compute_r0's spectral radius to 1e-6).

One deliberate deviation from a pure +/-50 % design: the upper bound of
beta_vh would be 1.05, which is not a valid probability; it is capped at 1.0
and flagged in the JSON (beta_vh_capped_at_1).  All other bounds stay at the
literal +/-50 % definition.
"""

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

import dv_common as C
import model as M
import sensitivity as S
from r0 import compute_r0, next_generation_matrix

SEED = 20260816
N_SAMPLES = 2000
HALF_WIDTHS = [0.20, 0.50]

BASE_VALUES = {
    'nv_ratio': float(M.NV_RATIO),
    'beta_hv': float(M.BETA_HV),
    'beta_vh': float(M.BETA_VH),
    'mu_v': float(M.MU_V),
    'gamma_h': float(M.GAMMA_H),
    'b': C.B_REF,
    'eip_days': C.EIP_REF_DAYS,
    'sigma_h': float(M.SIGMA_H),
}


def equal_width_ranges(half_width):
    rng = {k: (v * (1.0 - half_width), v * (1.0 + half_width))
           for k, v in BASE_VALUES.items()}
    capped = False
    if rng['beta_vh'][1] > 1.0:
        rng['beta_vh'] = (rng['beta_vh'][0], 1.0)
        capped = True
    return rng, capped


def r0_of_sample(row):
    """R0 with b and sigma_v as sampled constants (same algebra as
    sensitivity.samples_to_params; equals the NGM spectral radius)."""
    p = {
        'nv_ratio': float(row[0]),
        'beta_hv': float(row[1]),
        'beta_vh': float(row[2]),
        'mu_v': float(row[3]),
        'gamma_h': float(row[4]),
        'sigma_h': float(row[7]),
        'mu_h': M.MU_H,
    }
    b = float(row[5])
    sv = 1.0 / float(row[6])
    m = p['nv_ratio']
    K02 = p['beta_vh'] * b * sv / ((sv + p['mu_v']) * p['mu_v'])
    K20 = p['beta_hv'] * b * m * p['sigma_h'] / ((p['sigma_h'] + p['mu_h']) * (p['gamma_h'] + p['mu_h']))
    return float(np.sqrt(max(K02 * K20, 0.0)))


def r0_via_ngm(row):
    """Independent check: patch the temperature functions to the sampled
    constants and use r0.compute_r0 (spectral radius of the NGM).

    Note: r0.py binds biting_rate/eip_rate at import time, so both the model
    and the r0 module namespaces must be patched (and restored).
    """
    import r0 as r0_mod
    saved_m = (M.biting_rate, M.eip_rate)
    saved_r = (r0_mod.biting_rate, r0_mod.eip_rate)
    b, sv = float(row[5]), 1.0 / float(row[6])
    M.biting_rate = (lambda T, v=b: v)
    M.eip_rate = (lambda T, v=sv: v)
    r0_mod.biting_rate = M.biting_rate
    r0_mod.eip_rate = M.eip_rate
    try:
        p = M.baseline_params(
            nv_ratio=float(row[0]), beta_hv=float(row[1]),
            beta_vh=float(row[2]), mu_v=float(row[3]),
            gamma_h=float(row[4]), sigma_h=float(row[7]))
        return compute_r0(p, 28.0)
    finally:
        M.biting_rate, M.eip_rate = saved_m
        r0_mod.biting_rate, r0_mod.eip_rate = saved_r


def five(x):
    return {'min': float(np.min(x)), 'p05': float(np.percentile(x, 5)),
            'median': float(np.median(x)), 'p95': float(np.percentile(x, 95)),
            'max': float(np.max(x)), 'mean': float(np.mean(x))}


def main():
    t_start = time.perf_counter()
    C.reset_rates()

    # sanity: the K02*K20 algebra reproduces the published R0(28) values
    base = M.baseline_params()
    assert abs(compute_r0(base, 28.0) - 2.928116) < 1e-5, "R0(28) mismatch"

    blocks = {}
    for hw in HALF_WIDTHS:
        rng, capped = equal_width_ranges(hw)
        samples = C.sample_lhs_ranges(S.PARAM_NAMES, rng,
                                      n_samples=N_SAMPLES, seed=SEED)
        r0 = np.array([r0_of_sample(row) for row in samples])

        # independent verification on 100 random samples
        idx = np.random.default_rng(7).choice(N_SAMPLES, 100, replace=False)
        diff = np.array([r0_of_sample(samples[i]) - r0_via_ngm(samples[i])
                         for i in idx])
        verify = {
            'n_checked': int(idx.size),
            'max_abs_diff_vs_ngm_spectral_radius': float(np.max(np.abs(diff))),
            'passed': bool(np.max(np.abs(diff)) < 1e-9),
        }

        rows, dropped = C.prcc_table(samples, r0, S.PARAM_NAMES)
        blocks[f'half_width_{int(hw * 100)}'] = {
            'half_width': hw,
            'ranges': {k: {'lo': v[0], 'hi': v[1]} for k, v in rng.items()},
            'beta_vh_capped_at_1': capped,
            'prcc_table': rows,
            'n_dropped_nonfinite': dropped,
            'r0_distribution': five(r0),
            'verification_ngm_crosscheck': verify,
        }
        print(f"--- static PRCC, +/-{int(hw * 100)} % "
              f"(R0 {r0.min():.3f}-{r0.max():.3f}, capped={capped}, "
              f"ngm-check {verify['passed']}) ---")
        for r in rows:
            sig = '***' if r['p_value'] < 0.001 else (
                '**' if r['p_value'] < 0.01 else ('*' if r['p_value'] < 0.05 else 'n.s.'))
            print(f"{r['rank']:<3}{r['param']:<12}{r['prcc']:>9.3f}"
                  f"{C.fmt_p(r['p_value']):>12}  {sig}")

    wall = time.perf_counter() - t_start
    payload = {
        'task': 'F - static R0 PRCC at equal-width ranges (+/-20 % and +/-50 %)',
        'method': (
            'Algebraic R0 (K02*K20, identical to sensitivity.samples_to_params '
            'and verified against the NGM spectral radius), LHS n=2000, seed '
            '20260816, PRCC per Marino et al. (2008) with df = n - 9. No ODE '
            'solves.'
        ),
        'pairing': (
            'Same LHS construction and seed as tasks A/D/E, so the uniform '
            'variates are paired sample-for-sample across all analyses.'
        ),
        'outputs': blocks,
        'meta': C.make_meta(
            script='run_static_equalwidth.py', seed=SEED,
            n_samples=N_SAMPLES, wall_seconds=wall,
            per_run_seconds=0.0,
            extra={'n_ode_solves': 0, 'half_widths': HALF_WIDTHS,
                   'beta_vh_capped_at_1': True}),
    }
    C.write_json(os.path.join(HERE, 'static_equalwidth_prcc.json'), payload)


if __name__ == '__main__':
    main()
