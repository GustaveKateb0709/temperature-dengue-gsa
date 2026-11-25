"""Task G: dynamic PRCC at +/-50 % equal-width parameter ranges.

Same design as task D (run_equalwidth.py) but with a half-width of 50 %:
same LHS construction, n=2000, seed 20260816, same dynamic outputs (peak Ih,
peak day, 3-year cumulative infection), same PRCC estimator.  The uniform
variates u remain paired with tasks A/D/E.

One documented deviation from the literal +/-50 % definition: beta_vh is a
probability, and 0.70 * 1.5 = 1.05 is not valid, so its upper bound is capped
at 1.0 (flagged as beta_vh_capped_at_1).  All other bounds follow +/-50 %.
"""

import json
import os
import sys
import time
import multiprocessing as mp

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

import dv_common as C
import model as M
import sensitivity as S

SEED = 20260816
N_SAMPLES = 2000
HALF_WIDTH = 0.50

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


def equal_width_ranges():
    rng = {k: (v * (1.0 - HALF_WIDTH), v * (1.0 + HALF_WIDTH))
           for k, v in BASE_VALUES.items()}
    capped = False
    if rng['beta_vh'][1] > 1.0:
        rng['beta_vh'] = (rng['beta_vh'][0], 1.0)
        capped = True
    return rng, capped


def build_jobs(samples):
    jobs = []
    for row in samples:
        overrides = {
            'nv_ratio': float(row[0]),
            'beta_hv': float(row[1]),
            'beta_vh': float(row[2]),
            'mu_v': float(row[3]),
            'gamma_h': float(row[4]),
            'sigma_h': float(row[7]),
        }
        b_scale = float(row[5]) / C.B_REF
        sv_scale = (1.0 / float(row[6])) / C.SV_REF
        jobs.append((b_scale, sv_scale, overrides))
    return jobs


def five(x):
    x = x[np.isfinite(x)]
    return {'min': float(np.min(x)), 'p05': float(np.percentile(x, 5)),
            'median': float(np.median(x)), 'p95': float(np.percentile(x, 95)),
            'max': float(np.max(x)), 'mean': float(np.mean(x))}


def main():
    t_start = time.perf_counter()
    C.reset_rates()

    t0 = time.perf_counter()
    C.run_case(M.baseline_params())
    per_run = time.perf_counter() - t0

    ranges, capped = equal_width_ranges()
    samples = C.sample_lhs_ranges(S.PARAM_NAMES, ranges,
                                  n_samples=N_SAMPLES, seed=SEED)

    # paired-design check against task A's variates
    a_samples = S.sample_lhs(n_samples=N_SAMPLES, seed=SEED)
    lo_a = np.array([S.PARAM_RANGES[n][0] for n in S.PARAM_NAMES])
    wid_a = np.array([S.PARAM_RANGES[n][1] - S.PARAM_RANGES[n][0] for n in S.PARAM_NAMES])
    u_a = (a_samples - lo_a) / wid_a
    lo_g = np.array([ranges[n][0] for n in S.PARAM_NAMES])
    wid_g = np.array([ranges[n][1] - ranges[n][0] for n in S.PARAM_NAMES])
    u_g = (samples - lo_g) / wid_g
    paired = bool(np.allclose(u_a, u_g, atol=1e-12))
    print(f"[check] LHS variates paired with task A: {paired}")
    print(f"[check] beta_vh capped at 1.0: {capped}")

    jobs = build_jobs(samples)
    n_proc = max(1, min(7, mp.cpu_count() - 1))
    print(f"[run] {len(jobs)} solves on {n_proc} processes ...")
    if n_proc > 1:
        with mp.Pool(n_proc) as pool:
            results = pool.map(C.worker, jobs, chunksize=16)
    else:
        results = [C.worker(j) for j in jobs]
    wall = time.perf_counter() - t_start
    n_fail = sum(0 if r['ok'] else 1 for r in results)
    print(f"[run] done in {wall:.1f} s ({n_fail} failed)")

    peak_ih = np.array([r['peak_ih'] for r in results])
    peak_day = np.array([r['peak_day'] for r in results])
    cum = np.array([r['cum_infected'] for r in results])

    blocks = {}
    for out_name, arr in (('peak_ih', peak_ih), ('peak_day', peak_day),
                          ('cum_infected', cum)):
        rows, dropped = C.prcc_table(samples, arr, S.PARAM_NAMES)
        blocks[out_name] = {'prcc_table': rows, 'n_dropped_nonfinite': dropped}
        print(f"\n--- +/-50 % dynamic PRCC vs {out_name} (n={rows[0]['n_eff']}) ---")
        for r in rows:
            sig = '***' if r['p_value'] < 0.001 else (
                '**' if r['p_value'] < 0.01 else ('*' if r['p_value'] < 0.05 else 'n.s.'))
            print(f"{r['rank']:<3}{r['param']:<12}{r['prcc']:>9.3f}"
                  f"{C.fmt_p(r['p_value']):>12}  {sig}")

    # comparison with task A (published widths), same format as task D
    comparison = {}
    a_path = os.path.join(HERE, 'dynamic_prcc.json')
    if os.path.exists(a_path):
        with open(a_path, encoding='utf-8') as fh:
            a_json = json.load(fh)
        for out_name in blocks:
            a_map = {r['param']: r for r in a_json['outputs'][out_name]['prcc_table']}
            g_map = {r['param']: r for r in blocks[out_name]['prcc_table']}
            names = [r['param'] for r in blocks[out_name]['prcc_table']]
            a_rank = [a_map[n]['rank'] for n in names]
            g_rank = [g_map[n]['rank'] for n in names]
            comparison[out_name] = {
                'per_param': [
                    {'param': n,
                     'prcc_published_widths': a_map[n]['prcc'],
                     'prcc_half50': g_map[n]['prcc'],
                     'rank_published_widths': a_map[n]['rank'],
                     'rank_half50': g_map[n]['rank'],
                     'delta_prcc': g_map[n]['prcc'] - a_map[n]['prcc']}
                    for n in names
                ],
                'spearman_rank_correlation_of_prcc': float(
                    stats.spearmanr([a_map[n]['prcc'] for n in names],
                                    [g_map[n]['prcc'] for n in names]).statistic),
                'identical_ranking': bool(a_rank == g_rank),
                'max_abs_rank_change': int(max(abs(x - y) for x, y in zip(a_rank, g_rank))),
            }

    payload = {
        'task': 'G - dynamic PRCC with +/-50 % equal-width ranges',
        'method': (
            'Same pipeline as task A (LHS n=2000, seed 20260816, baseline '
            'scenario, t_span=(0,1095), daily output, Marino PRCC with '
            'df = n - 9) but every parameter is sampled on '
            '[0.5*theta_base, 1.5*theta_base].'
        ),
        'beta_vh_capped_at_1': capped,
        'beta_vh_note': (
            '0.70 * 1.5 = 1.05 is not a valid probability; the upper bound is '
            'capped at 1.0. The other seven ranges follow the literal +/-50 %.'
        ),
        'ranges': {k: {'lo': v[0], 'hi': v[1], 'baseline': BASE_VALUES[k]}
                   for k, v in ranges.items()},
        'lhs_variates_paired_with_task_A': paired,
        'outputs': blocks,
        'output_distributions': {'peak_ih': five(peak_ih),
                                 'peak_day': five(peak_day),
                                 'cum_infected': five(cum)},
        'comparison_vs_task_A': comparison,
        'n_failed_solves': n_fail,
        'meta': C.make_meta(
            script='run_equalwidth50.py', seed=SEED, n_samples=N_SAMPLES,
            wall_seconds=wall, per_run_seconds=per_run,
            extra={'parallel_processes': n_proc, 'half_width': HALF_WIDTH}),
    }
    C.write_json(os.path.join(HERE, 'equalwidth50_prcc.json'), payload)


if __name__ == '__main__':
    main()
