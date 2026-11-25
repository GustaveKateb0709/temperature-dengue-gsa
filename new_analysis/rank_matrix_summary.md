# Rank matrix — 3 responses x 3 range designs

Every cell is the ranking of the same 8 parameters by |PRCC| under one (response, range) combination: n=2000 LHS, seed 20260816 (paired variates), Marino et al. (2008) estimator, df = n − 9. Static cells evaluate the algebraic R0; dynamic cells run the 1,095-day baseline scenario. All numbers rounded to 3 decimals.

In the ±50 % columns the upper bound of beta_vh (1.05) is capped at 1.0 (probability bound); everything else follows the literal ±50 %.

## PRCC per cell

| Parameter | static@lit | static@±20 | static@±50 | peak@lit | peak@±20 | peak@±50 | cum@lit | cum@±20 | cum@±50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nv_ratio | +0.954 | +0.858 | +0.816 | +0.950 | +0.827 | +0.716 | +0.955 | +0.786 | +0.761 |
| beta_hv | +0.762 | +0.861 | +0.821 | +0.737 | +0.811 | +0.699 | +0.758 | +0.771 | +0.757 |
| beta_vh | +0.775 | +0.866 | +0.815 | +0.757 | +0.829 | +0.705 | +0.777 | +0.792 | +0.757 |
| mu_v | -0.722 | -0.904 | -0.869 | -0.476 | -0.623 | -0.504 | -0.645 | -0.866 | -0.821 |
| gamma_h | -0.639 | -0.864 | -0.825 | -0.729 | -0.943 | -0.870 | -0.559 | -0.773 | -0.735 |
| b | +0.573 | +0.962 | +0.947 | +0.553 | +0.946 | +0.901 | +0.571 | +0.934 | +0.926 |
| eip_days | -0.354 | -0.428 | -0.346 | -0.394 | -0.454 | -0.314 | -0.410 | -0.322 | -0.323 |
| sigma_h | +0.000 | +0.021 | +0.012 | +0.263 | +0.637 | +0.461 | +0.117 | +0.010 | +0.090 |

## Rank per cell (1 = strongest influence)

| Parameter | static@lit | static@±20 | static@±50 | peak@lit | peak@±20 | peak@±50 | cum@lit | cum@±20 | cum@±50 | span |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nv_ratio | 1 | 6 | 5 | 1 | 4 | 3 | 1 | 4 | 3 | 1–6 |
| beta_hv | 3 | 5 | 4 | 3 | 5 | 5 | 3 | 6 | 4 | 3–6 |
| beta_vh | 2 | 3 | 6 | 2 | 3 | 4 | 2 | 3 | 5 | 2–6 |
| mu_v | 4 | 2 | 2 | 6 | 7 | 6 | 4 | 2 | 2 | 2–7 |
| gamma_h | 5 | 4 | 3 | 4 | 2 | 2 | 6 | 5 | 6 | 2–6 |
| b | 6 | 1 | 1 | 5 | 1 | 1 | 5 | 1 | 1 | 1–6 |
| eip_days | 7 | 7 | 7 | 7 | 8 | 8 | 7 | 7 | 7 | 7–8 |
| sigma_h | 8 | 8 | 8 | 8 | 6 | 7 | 8 | 8 | 8 | 6–8 |

## Pairwise rank agreement

Spearman ρ (upper right) and Kendall τ (lower left) between the rank vectors of the 9 cells. Cell naming: `response@range`.

| cell | static_R0@published | static_R0@hw20 | static_R0@hw50 | peak_ih@published | peak_ih@hw20 | peak_ih@hw50 | cum_infected@published | cum_infected@hw20 | cum_infected@hw50 |
|---|---|---|---|---|---|---|---|---|---|
| static_R0@published | — | ρ=0.286 | ρ=0.214 | ρ=0.929 | ρ=0.262 | ρ=0.381 | ρ=0.976 | ρ=0.429 | ρ=0.476 |
| static_R0@hw20 | τ=0.143 | — | ρ=0.857 | ρ=0.262 | ρ=0.548 | ρ=0.619 | ρ=0.357 | ρ=0.929 | ρ=0.786 |
| static_R0@hw50 | τ=0.071 | τ=0.786 | — | ρ=0.214 | ρ=0.500 | ρ=0.667 | ρ=0.262 | ρ=0.786 | ρ=0.833 |
| peak_ih@published | τ=0.857 | τ=0.143 | τ=0.071 | — | ρ=0.524 | ρ=0.595 | ρ=0.905 | ρ=0.381 | ρ=0.405 |
| peak_ih@hw20 | τ=0.143 | τ=0.429 | τ=0.357 | τ=0.286 | — | ρ=0.952 | ρ=0.286 | ρ=0.524 | ρ=0.381 |
| peak_ih@hw50 | τ=0.286 | τ=0.429 | τ=0.500 | τ=0.429 | τ=0.857 | — | ρ=0.405 | ρ=0.643 | ρ=0.571 |
| cum_infected@published | τ=0.929 | τ=0.214 | τ=0.143 | τ=0.786 | τ=0.214 | τ=0.357 | — | ρ=0.524 | ρ=0.595 |
| cum_infected@hw20 | τ=0.286 | τ=0.857 | τ=0.643 | τ=0.286 | τ=0.429 | τ=0.429 | τ=0.357 | — | ρ=0.881 |
| cum_infected@hw50 | τ=0.357 | τ=0.643 | τ=0.714 | τ=0.357 | τ=0.214 | τ=0.357 | τ=0.429 | τ=0.786 | — |

### View 1 — same response, across ranges

| Response | Range pair | ρ | τ | # swapped parameter pairs |
|---|---|---:|---:|---:|
| static_R0 | published ↔ hw20 | 0.286 | 0.143 | 12 |
| static_R0 | published ↔ hw50 | 0.214 | 0.071 | 13 |
| static_R0 | hw20 ↔ hw50 | 0.857 | 0.786 | 3 |
| peak_ih | published ↔ hw20 | 0.524 | 0.286 | 10 |
| peak_ih | published ↔ hw50 | 0.595 | 0.429 | 8 |
| peak_ih | hw20 ↔ hw50 | 0.952 | 0.857 | 2 |
| cum_infected | published ↔ hw20 | 0.524 | 0.357 | 9 |
| cum_infected | published ↔ hw50 | 0.595 | 0.429 | 8 |
| cum_infected | hw20 ↔ hw50 | 0.881 | 0.786 | 3 |

### View 2 — same range, across responses

| Range | Response pair | ρ | τ | # swapped parameter pairs |
|---|---|---:|---:|---:|
| published | static_R0 ↔ peak_ih | 0.929 | 0.857 | 2 |
| published | static_R0 ↔ cum_infected | 0.976 | 0.929 | 1 |
| published | peak_ih ↔ cum_infected | 0.905 | 0.786 | 3 |
| hw20 | static_R0 ↔ peak_ih | 0.548 | 0.429 | 8 |
| hw20 | static_R0 ↔ cum_infected | 0.929 | 0.857 | 2 |
| hw20 | peak_ih ↔ cum_infected | 0.524 | 0.429 | 8 |
| hw50 | static_R0 ↔ peak_ih | 0.667 | 0.500 | 7 |
| hw50 | static_R0 ↔ cum_infected | 0.833 | 0.714 | 4 |
| hw50 | peak_ih ↔ cum_infected | 0.571 | 0.357 | 9 |

## Top-3 composition

| Cell | Top 3 |
|---|---|
| static_R0@published | nv_ratio > beta_vh > beta_hv |
| static_R0@hw20 | b > mu_v > beta_vh |
| static_R0@hw50 | b > mu_v > gamma_h |
| peak_ih@published | nv_ratio > beta_vh > beta_hv |
| peak_ih@hw20 | b > gamma_h > beta_vh |
| peak_ih@hw50 | b > gamma_h > nv_ratio |
| cum_infected@published | nv_ratio > beta_vh > beta_hv |
| cum_infected@hw20 | b > mu_v > beta_vh |
| cum_infected@hw50 | b > mu_v > nv_ratio |

Distinct top-3 sets across the 9 cells: **6**.

| Top-3 set | Cells |
|---|---|
| beta_hv > beta_vh > nv_ratio | static_R0@published, peak_ih@published, cum_infected@published |
| b > beta_vh > mu_v | static_R0@hw20, cum_infected@hw20 |
| b > gamma_h > mu_v | static_R0@hw50 |
| b > beta_vh > gamma_h | peak_ih@hw20 |
| b > gamma_h > nv_ratio | peak_ih@hw50 |
| b > mu_v > nv_ratio | cum_infected@hw50 |

Frequency of each parameter inside a top-3 (out of 9 cells): beta_vh 6/9, b 6/9, nv_ratio 5/9, mu_v 4/9, beta_hv 3/9, gamma_h 3/9, eip_days 0/9, sigma_h 0/9.

## Rank-order flips

For each pair of cells, the number of parameter pairs whose relative order differs (max possible: 28).

| Cell pair | # swapped pairs |
|---|---:|
| static_R0@published|static_R0@hw50 | 13 |
| static_R0@hw50|peak_ih@published | 13 |
| static_R0@published|static_R0@hw20 | 12 |
| static_R0@published|peak_ih@hw20 | 12 |
| static_R0@hw20|peak_ih@published | 12 |
| static_R0@hw50|cum_infected@published | 12 |
| static_R0@hw20|cum_infected@published | 11 |
| peak_ih@hw20|cum_infected@published | 11 |
| peak_ih@hw20|cum_infected@hw50 | 11 |
| static_R0@published|peak_ih@hw50 | 10 |
| static_R0@published|cum_infected@hw20 | 10 |
| peak_ih@published|peak_ih@hw20 | 10 |
| peak_ih@published|cum_infected@hw20 | 10 |
| static_R0@published|cum_infected@hw50 | 9 |
| static_R0@hw50|peak_ih@hw20 | 9 |
| peak_ih@published|cum_infected@hw50 | 9 |
| peak_ih@hw50|cum_infected@published | 9 |
| peak_ih@hw50|cum_infected@hw50 | 9 |
| cum_infected@published|cum_infected@hw20 | 9 |
| static_R0@hw20|peak_ih@hw20 | 8 |
| static_R0@hw20|peak_ih@hw50 | 8 |
| peak_ih@published|peak_ih@hw50 | 8 |
| peak_ih@hw20|cum_infected@hw20 | 8 |
| peak_ih@hw50|cum_infected@hw20 | 8 |
| cum_infected@published|cum_infected@hw50 | 8 |
| static_R0@hw50|peak_ih@hw50 | 7 |
| static_R0@hw20|cum_infected@hw50 | 5 |
| static_R0@hw50|cum_infected@hw20 | 5 |
| static_R0@hw50|cum_infected@hw50 | 4 |
| static_R0@hw20|static_R0@hw50 | 3 |
| peak_ih@published|cum_infected@published | 3 |
| cum_infected@hw20|cum_infected@hw50 | 3 |
| static_R0@published|peak_ih@published | 2 |
| static_R0@hw20|cum_infected@hw20 | 2 |
| peak_ih@hw20|peak_ih@hw50 | 2 |
| static_R0@published|cum_infected@published | 1 |

The full per-pair flip details (which parameter pair swaps in which order) are in `rank_matrix.json` under `flips`.

## Self-check

- all 9 cells present: True
- stored ranks match |PRCC| ordering in every cell: True
- all 36 pairwise correlations reproduce from the written JSON: True

