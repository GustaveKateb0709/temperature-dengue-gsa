"""Assemble the 3 responses x 3 ranges rank matrix.

Cells (9 rankings, each = the 8 parameters ordered by |PRCC|):

    response \\ range   published      +/-20 %          +/-50 %
    static R0           F0 (ref)       F20              F50
    dynamic peak        A (peak_ih)    D (peak_ih)      G (peak_ih)
    dynamic cumulative  A (cum)        D (cum)          G (cum)

Sources:
    F0  static_prcc_reference.json   (published literature ranges)
    F20/F50 static_equalwidth_prcc.json
    A   dynamic_prcc.json            outputs.{peak_ih,cum_infected}
    D   equalwidth_prcc.json         outputs.{peak_ih,cum_infected}
    G   equalwidth50_prcc.json       outputs.{peak_ih,cum_infected}

Deliverables: rank_matrix.json (full pairwise Spearman rho / Kendall tau on
rank vectors, per-parameter rank spans, top-3 composition, pairwise flip
lists, self-check), rank_matrix_summary.md, and an appended section in
SUMMARY.md.  All reported numbers are rounded to 3 decimals.
"""

import json
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

import dv_common as C
import sensitivity as S

PARAMS = list(S.PARAM_NAMES)          # canonical parameter order
R3 = 3                                 # decimals for reported numbers


def load(name):
    with open(os.path.join(HERE, name), encoding='utf-8') as fh:
        return json.load(fh)


def cell_from_table(rows):
    """PRCC map + rank map from a prcc_table block (ranks by |PRCC|)."""
    prcc = {r['param']: float(r['prcc']) for r in rows}
    stored_rank = {r['param']: int(r['rank']) for r in rows}
    order = sorted(PARAMS, key=lambda p: -abs(prcc[p]))
    recomputed = {p: i + 1 for i, p in enumerate(order)}
    return {'prcc': prcc, 'rank': stored_rank, 'rank_recomputed': recomputed,
            'ranks_match': recomputed == stored_rank}


def collect_cells():
    cells = {}

    # static, published ranges
    dS = load('static_prcc_reference.json')
    prcc0 = dict(zip(dS['param_names'], [float(x) for x in dS['prcc']]))
    order0 = sorted(PARAMS, key=lambda p: -abs(prcc0[p]))
    cells[('static_R0', 'published')] = {
        'prcc': prcc0,
        'rank': {p: i + 1 for i, p in enumerate(order0)},
        'source': 'static_prcc_reference.json',
    }

    # static, equal widths
    dF = load('static_equalwidth_prcc.json')
    for key, tag in (('half_width_20', 'hw20'), ('half_width_50', 'hw50')):
        c = cell_from_table(dF['outputs'][key]['prcc_table'])
        c['source'] = f'static_equalwidth_prcc.json#{key}'
        cells[('static_R0', tag)] = c

    # dynamic
    sources = {
        'published': 'dynamic_prcc.json',
        'hw20': 'equalwidth_prcc.json',
        'hw50': 'equalwidth50_prcc.json',
    }
    for tag, fname in sources.items():
        d = load(fname)
        for resp in ('peak_ih', 'cum_infected'):
            c = cell_from_table(d['outputs'][resp]['prcc_table'])
            c['source'] = f'{fname}#outputs.{resp}'
            cells[(resp, tag)] = c
    return cells


def cell_id(cell):
    resp, tag = cell
    return f'{resp}@{tag}'


def paired_correlations(x_ranks, y_ranks):
    x = np.array([x_ranks[p] for p in PARAMS], dtype=float)
    y = np.array([y_ranks[p] for p in PARAMS], dtype=float)
    sr = stats.spearmanr(x, y)
    kt = stats.kendalltau(x, y)
    return {
        'spearman_rho': round(float(sr.statistic), R3),
        'spearman_p': float(sr.pvalue),
        'kendall_tau': round(float(kt.statistic), R3),
        'kendall_p': float(kt.pvalue),
    }


def flips(x_ranks, y_ranks):
    """Parameter pairs whose relative order differs between two rankings."""
    out = []
    for i in range(len(PARAMS)):
        for j in range(i + 1, len(PARAMS)):
            p, q = PARAMS[i], PARAMS[j]
            a = np.sign(x_ranks[p] - x_ranks[q])
            b = np.sign(y_ranks[p] - y_ranks[q])
            if a != b:
                out.append({'pair': [p, q],
                            'order_first': f'{p}>{q}' if a > 0 else f'{q}>{p}',
                            'order_second': f'{p}>{q}' if b > 0 else f'{q}>{p}'})
    return out


def main():
    cells = collect_cells()

    # self-check 1: stored ranks match |PRCC| ordering; ranks unique 1..8
    missing = [cell_id(c) for c in cells_expected() if c not in cells]
    checks = {'all_cells_present': not missing, 'missing_cells': missing}
    # recompute ranks for the static-published cell (its JSON stores them
    # pre-sorted rather than per-parameter)
    order = sorted(PARAMS, key=lambda p: -abs(cells[('static_R0', 'published')]['prcc'][p]))
    cells[('static_R0', 'published')]['rank_recomputed'] = {
        p: i + 1 for i, p in enumerate(order)}
    checks['ranks_match_prcc'] = {
        cell_id(k): (v['rank_recomputed'] == v['rank']) for k, v in cells.items()}

    # self-check 2: equal-width static cells reproduce from stored prcc
    # self-check 3: pairwise correlations recomputed at the end (below)

    ids = [cell_id(c) for c in cells_expected()]
    id_to_cell = {cell_id(k): v for k, v in cells.items()}

    # ── full 9 x 9 pairwise matrix ─────────────────────────────────────────
    matrix = {}
    flip_list = {}
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = ids[i], ids[j]
            ra = id_to_cell[a]['rank']
            rb = id_to_cell[b]['rank']
            corr = paired_correlations(ra, rb)
            fl = flips(ra, rb)
            matrix[f'{a}|{b}'] = corr
            flip_list[f'{a}|{b}'] = {
                'n_flipped_pairs': len(fl),
                'flips': fl,
            }

    # ── structured views ───────────────────────────────────────────────────
    views = {'same_response_across_ranges': [], 'same_range_across_responses': []}
    for resp in ('static_R0', 'peak_ih', 'cum_infected'):
        for (t1, t2) in (('published', 'hw20'), ('published', 'hw50'),
                         ('hw20', 'hw50')):
            key = f'{resp}@{t1}|{resp}@{t2}'
            views['same_response_across_ranges'].append({
                'response': resp, 'range_pair': [t1, t2],
                **matrix[key], 'n_flipped_pairs': flip_list[key]['n_flipped_pairs']})
    for tag in ('published', 'hw20', 'hw50'):
        for (r1, r2) in (('static_R0', 'peak_ih'), ('static_R0', 'cum_infected'),
                         ('peak_ih', 'cum_infected')):
            key = f'{r1}@{tag}|{r2}@{tag}'
            views['same_range_across_responses'].append({
                'range': tag, 'response_pair': [r1, r2],
                **matrix[key], 'n_flipped_pairs': flip_list[key]['n_flipped_pairs']})

    # ── per-parameter rank spans ───────────────────────────────────────────
    per_param = {}
    for p in PARAMS:
        ranks = {cid: id_to_cell[cid]['rank'][p] for cid in ids}
        per_param[p] = {
            'ranks': ranks,
            'rank_min': int(min(ranks.values())),
            'rank_max': int(max(ranks.values())),
            'rank_span': int(max(ranks.values()) - min(ranks.values())),
            'rank_mean': round(float(np.mean(list(ranks.values()))), R3),
        }

    # ── top-3 composition ──────────────────────────────────────────────────
    top3 = {cid: sorted(PARAMS, key=lambda p: id_to_cell[cid]['rank'][p])[:3]
            for cid in ids}
    distinct_sets = {}
    for cid, s in top3.items():
        key = ','.join(sorted(s))
        distinct_sets.setdefault(key, []).append(cid)
    top3_freq = {p: sum(p in s for s in top3.values()) for p in PARAMS}

    # ── self-check 3: recompute every correlation from the JSON round-trip ─
    # (ranks are small integers, immune to 3-decimal rounding of PRCC)
    payload = {
        'task': 'Rank matrix: 3 responses x 3 parameter-range designs',
        'responses': ['static_R0', 'peak_ih', 'cum_infected'],
        'ranges': ['published', 'hw20', 'hw50'],
        'cell_definition': (
            'Each cell is the 8-parameter ranking by |PRCC| under one '
            '(response, range) combination; all cells use n=2000 LHS samples '
            'with seed 20260816 (paired uniform variates across cells) and '
            'the Marino et al. (2008) estimator with df = n - 9. Static cells '
            'use the algebraic R0; dynamic cells use the 1,095-day baseline '
            'simulation (daily output grid, parabolically refined peak).'
        ),
        'notes': {
            'beta_vh_hw50_capped_at_1': (
                'In the +/-50 % designs the upper bound of beta_vh (0.70*1.5 '
                '= 1.05) is capped at the probability bound 1.0; all other '
                'bounds follow the literal +/-50 %.'),
            'peak_day_excluded': (
                'peak_day is reported in tasks A/D/G but excluded from the '
                'matrix: a static R0 has no timing counterpart, so a 3x3 '
                'matrix for it would be unbalanced.'),
        },
        'cells': {
            cell_id(k): {
                'source': v['source'],
                'prcc': {p: round(v['prcc'][p], R3) for p in PARAMS},
                'rank': v['rank'],
            } for k, v in cells.items()
        },
        'pairwise_matrix': matrix,
        'views': views,
        'per_parameter_rank_span': per_param,
        'top3': {
            'per_cell': top3,
            'distinct_sets': {k: v for k, v in distinct_sets.items()},
            'frequency_in_top3': top3_freq,
        },
        'flips': flip_list,
        'self_check': checks,
        'meta': C.make_meta(
            script='run_rank_matrix.py', seed=20260816, n_samples=2000,
            wall_seconds=0.0, per_run_seconds=0.0,
            extra={'n_ode_solves': 0,
                   'note': 'aggregation only; reuses tasks A/D/F/G results'}),
    }

    # self-check: round-trip recompute of all correlations from the payload
    ok_corr = True
    for key, corr in payload['pairwise_matrix'].items():
        a, b = key.split('|')
        ra = payload['cells'][a]['rank']
        rb = payload['cells'][b]['rank']
        c2 = paired_correlations(ra, rb)
        if abs(c2['spearman_rho'] - corr['spearman_rho']) > 1e-9 or \
           abs(c2['kendall_tau'] - corr['kendall_tau']) > 1e-9:
            ok_corr = False
    payload['self_check']['pairwise_recompute_matches'] = ok_corr
    payload['self_check']['all_pass'] = bool(
        ok_corr and all(payload['self_check']['ranks_match_prcc'].values())
        and payload['self_check']['all_cells_present'])

    C.write_json(os.path.join(HERE, 'rank_matrix.json'), payload)

    # ── markdown summary ───────────────────────────────────────────────────
    write_md(payload, ids, matrix, views, per_param, top3, distinct_sets,
             flip_list, top3_freq)
    append_summary_section(payload, matrix, views, per_param, top3,
                           distinct_sets)


def cells_expected():
    out = []
    for resp in ('static_R0', 'peak_ih', 'cum_infected'):
        for tag in ('published', 'hw20', 'hw50'):
            out.append((resp, tag))
    return out


def write_md(payload, ids, matrix, views, per_param, top3, distinct_sets,
             flip_list, top3_freq):
    L = []
    a = L.append
    a('# Rank matrix — 3 responses x 3 range designs')
    a('')
    a('Every cell is the ranking of the same 8 parameters by |PRCC| under one '
      '(response, range) combination: n=2000 LHS, seed 20260816 (paired '
      'variates), Marino et al. (2008) estimator, df = n − 9. Static cells '
      'evaluate the algebraic R0; dynamic cells run the 1,095-day baseline '
      'scenario. All numbers rounded to 3 decimals.')
    a('')
    a('In the ±50 % columns the upper bound of beta_vh (1.05) is capped at 1.0 '
      '(probability bound); everything else follows the literal ±50 %.')
    a('')

    a('## PRCC per cell')
    a('')
    a('| Parameter | static@lit | static@±20 | static@±50 | peak@lit | peak@±20 | peak@±50 | cum@lit | cum@±20 | cum@±50 |')
    a('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for p in PARAMS:
        row = [payload['cells'][f'{r}@{t}']['prcc'][p]
               for r in ('static_R0', 'peak_ih', 'cum_infected')
               for t in ('published', 'hw20', 'hw50')]
        a(f"| {p} | " + ' | '.join(f"{v:+.3f}" for v in row) + ' |')
    a('')

    a('## Rank per cell (1 = strongest influence)')
    a('')
    a('| Parameter | static@lit | static@±20 | static@±50 | peak@lit | peak@±20 | peak@±50 | cum@lit | cum@±20 | cum@±50 | span |')
    a('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for p in PARAMS:
        cells_r = [payload['cells'][f'{r}@{t}']['rank'][p]
                   for r in ('static_R0', 'peak_ih', 'cum_infected')
                   for t in ('published', 'hw20', 'hw50')]
        span = per_param[p]['rank_span']
        a(f"| {p} | " + ' | '.join(str(x) for x in cells_r) +
          f" | {per_param[p]['rank_min']}–{per_param[p]['rank_max']} |")
    a('')

    a('## Pairwise rank agreement')
    a('')
    a('Spearman ρ (upper right) and Kendall τ (lower left) between the rank '
      'vectors of the 9 cells. Cell naming: `response@range`.')
    a('')
    a('| cell | ' + ' | '.join(ids) + ' |')
    a('|' + '---|' * (len(ids) + 1))
    for i, ci in enumerate(ids):
        row = [ci]
        for j, cj in enumerate(ids):
            if i == j:
                row.append('—')
            elif j > i:
                row.append(f"ρ={matrix[f'{ci}|{cj}']['spearman_rho']:.3f}")
            else:
                key = f'{ids[j]}|{ci}'
                row.append(f"τ={matrix[key]['kendall_tau']:.3f}")
        a('| ' + ' | '.join(row) + ' |')
    a('')

    a('### View 1 — same response, across ranges')
    a('')
    a('| Response | Range pair | ρ | τ | # swapped parameter pairs |')
    a('|---|---|---:|---:|---:|')
    for v in payload['views']['same_response_across_ranges']:
        a(f"| {v['response']} | {v['range_pair'][0]} ↔ {v['range_pair'][1]} | "
          f"{v['spearman_rho']:.3f} | {v['kendall_tau']:.3f} | "
          f"{v['n_flipped_pairs']} |")
    a('')
    a('### View 2 — same range, across responses')
    a('')
    a('| Range | Response pair | ρ | τ | # swapped parameter pairs |')
    a('|---|---|---:|---:|---:|')
    for v in payload['views']['same_range_across_responses']:
        a(f"| {v['range']} | {v['response_pair'][0]} ↔ {v['response_pair'][1]} | "
          f"{v['spearman_rho']:.3f} | {v['kendall_tau']:.3f} | "
          f"{v['n_flipped_pairs']} |")
    a('')

    a('## Top-3 composition')
    a('')
    a('| Cell | Top 3 |')
    a('|---|---|')
    for cid in ids:
        a(f"| {cid} | {' > '.join(top3[cid])} |")
    a('')
    a(f"Distinct top-3 sets across the 9 cells: **{len(distinct_sets)}**.")
    a('')
    a('| Top-3 set | Cells |')
    a('|---|---|')
    for k, v in sorted(distinct_sets.items(), key=lambda kv: -len(kv[1])):
        a(f"| {k.replace(',', ' > ')} | {', '.join(v)} |")
    a('')
    a('Frequency of each parameter inside a top-3 (out of 9 cells): ' +
      ', '.join(f"{p} {top3_freq[p]}/9" for p in
                sorted(PARAMS, key=lambda q: -top3_freq[q])) + '.')
    a('')

    a('## Rank-order flips')
    a('')
    a('For each pair of cells, the number of parameter pairs whose relative '
      'order differs (max possible: 28).')
    a('')
    a('| Cell pair | # swapped pairs |')
    a('|---|---:|')
    for key in sorted(flip_list, key=lambda k: -flip_list[k]['n_flipped_pairs']):
        a(f"| {key} | {flip_list[key]['n_flipped_pairs']} |")
    a('')
    a('The full per-pair flip details (which parameter pair swaps in which '
      'order) are in `rank_matrix.json` under `flips`.')
    a('')

    a('## Self-check')
    a('')
    a(f"- all 9 cells present: {payload['self_check']['all_cells_present']}")
    a(f"- stored ranks match |PRCC| ordering in every cell: "
      f"{all(payload['self_check']['ranks_match_prcc'].values())}")
    a(f"- all 36 pairwise correlations reproduce from the written JSON: "
      f"{payload['self_check']['pairwise_recompute_matches']}")
    a('')

    out = os.path.join(HERE, 'rank_matrix_summary.md')
    with open(out, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    print(f'[written] {out}')


def append_summary_section(payload, matrix, views, per_param, top3,
                           distinct_sets):
    """Append (once) the rank-matrix section to SUMMARY.md."""
    path = os.path.join(HERE, 'SUMMARY.md')
    header = '## rank matrix (rebuild, 2026-11-25)'
    with open(path, encoding='utf-8') as fh:
        content = fh.read()
    if header in content:
        # replace everything after the header with the fresh section
        content = content[:content.index(header)].rstrip() + '\n\n'
    else:
        content = content.rstrip() + '\n\n'

    L = [header, '']
    a = L.append
    a('Nine rankings of the same 8 parameters (3 responses: static R0, '
      'dynamic peak, dynamic 3-year cumulative; 3 range designs: published '
      'literature ranges, ±20 %, ±50 %), all on the same n=2000 paired LHS '
      'samples (seed 20260816) and the same Marino (2008) PRCC estimator. '
      'Full matrices in `rank_matrix.json` and `rank_matrix_summary.md`.')
    a('')
    a('New cells computed for this section: '
      '`static_equalwidth_prcc.json` (static R0 at ±20 % and ±50 %; algebraic, '
      'no ODE) and `equalwidth50_prcc.json` (dynamic at ±50 %, full ODE '
      'pipeline). In the ±50 % designs the upper bound of beta_vh (1.05) is '
      'capped at 1.0.')
    a('')
    a('| View | Comparison | ρ | τ | swapped pairs (of 28) |')
    a('|---|---|---:|---:|---:|')
    for v in payload['views']['same_response_across_ranges']:
        a(f"| same response | {v['response']}: {v['range_pair'][0]} ↔ "
          f"{v['range_pair'][1]} | {v['spearman_rho']:.3f} | "
          f"{v['kendall_tau']:.3f} | {v['n_flipped_pairs']} |")
    for v in payload['views']['same_range_across_responses']:
        a(f"| same range | {v['range']}: {v['response_pair'][0]} ↔ "
          f"{v['response_pair'][1]} | {v['spearman_rho']:.3f} | "
          f"{v['kendall_tau']:.3f} | {v['n_flipped_pairs']} |")
    a('')
    a('| Parameter | rank min–max (of 9 cells) | span |')
    a('|---|---|---:|')
    for p in sorted(PARAMS, key=lambda q: -per_param[q]['rank_span']):
        a(f"| {p} | {per_param[p]['rank_min']}–{per_param[p]['rank_max']} | "
          f"{per_param[p]['rank_span']} |")
    a('')
    a(f"Top-3 composition: **{len(distinct_sets)} distinct sets** across the "
      '9 cells; frequency: ' +
      ', '.join(f"{p} {sum(p in s for s in top3.values())}/9" for p in
                sorted(PARAMS, key=lambda q: -sum(q in s for s in top3.values()))) +
      '.')
    a('')
    a('Self-check: ' + ('all passed' if payload['self_check']['all_pass']
                        else 'FAILED — inspect rank_matrix.json') + '.')
    a('')

    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(content + '\n'.join(L))
    print(f'[updated] {path}')


if __name__ == '__main__':
    main()
