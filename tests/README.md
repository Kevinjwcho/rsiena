# tests/ — Three-way SAOM recovery tests

Organized to match the four statistic groups in
[`notes/three_way_statistics_catalog.tex`](../notes/three_way_statistics_catalog.tex).

## Layout

```
tests/
├── group1_basic/                     # Within-network + attribute effects
│   ├── recovery_within_Yi.R          # density/recip/triadic on Y[shared]
│   ├── recovery_within_Yself.R       # density/recip on Y[self]  (dead path)
│   ├── recovery_attributes_standard.R   # egoX/altX/simX  (STUB)
│   └── recovery_attributes_perceiver.R  # percX/percAltSame/...  (STUB)
│
├── group2_perception_to_self/        # Y[i] → Y[self]  (currently dead path)
│   ├── recovery_crprod_dir2.R
│   └── test_specification_crprod_dir2.R
│
├── group3_self_to_perception/        # Y[self] → Y[i]  (regression — broken)
│   ├── recovery_crprod_dir1.R
│   ├── identifiability_crprod_dir1.R
│   ├── test_specification_crprod_dir1.R
│   └── diagnostics/                  # Isolation tests from May 19 session
│       ├── baseline_crprodRecip.R
│       ├── no_share_crprodRecip.R
│       ├── estimation_mode.R
│       └── *.log                     # diagnostic run outputs
│
├── group4_perceiver_restricted/      # percRecip family  (blocked on G3)
│   ├── recovery_percRecip.R
│   ├── diagnose_percRecip.R
│   ├── diagnose_percRecip_simulator.R
│   └── artifacts/
│       ├── recovery_percRecip_summary.csv
│       ├── recovery_percRecip_raw.rds
│       ├── *_n10.{rds,csv}           # earlier runs with n=10
│       └── *.log                     # full-run debug logs
│
├── exploratory/                      # Pre-SAOM raw R simulations
│   ├── threeway-simulation.R
│   ├── threeway-simulation_prob.R
│   └── threeway-simulation_prob_separate.R
│
└── (root — RSiena standard test fixtures, leave alone)
    ├── parallel.R, parallel.Rout.save
    ├── s50.csv, s50e.csv, s50paj.csv
    ├── s50-network[1-3].dat, s50_d[1-3].net, s50e.dat
    ├── sienaDependent_test.r, sienaDependent_test_symmetric.r
    └── test_chains.r
```

## Status at a glance

See `notes/three_way_statistics_catalog.tex` Summary table for full
state. Quick version:

| Group | Status | Notes |
|---|---|---|
| 1a Within Y[i] | not re-verified | recommended first diagnostic |
| 1b Within Y[self] | dead | Y[self] non-evolution |
| 1c std attrs | works | standard SAOM |
| 1c perc attrs | not re-verified | recommended second diagnostic |
| 2 Y[i] → Y[self] | dead | Y[self] non-evolution |
| 3 Y[self] → Y[i] | broken (regressed) | next-session top priority |
| 4 Perceiver-restricted | blocked on G3 fix | wiring identical to G3 |

## Next-session priorities

Run in order to localize the broken layer:

1. `group1_basic/recovery_within_Yi.R` — does the basic Y[i] objective
   machinery still work? Recovers density and recip on Y[shared].
2. `group1_basic/recovery_attributes_perceiver.R` (after filling in the
   stub) — does the perceiver-role attribute path work?
3. `group3_self_to_perception/diagnostics/baseline_crprodRecip.R` —
   does standard `crprodRecip` recover under any configuration?
4. If 1 and 2 pass and 3 fails: bisect commits between `41c4d03` (last
   commit before this session) and current HEAD to find the regression.
   See `notes/three_way_statistics_catalog.tex` § Regression note.

## Conventions

* All recovery scripts write outputs to their own
  `artifacts/` subfolder via absolute paths so they work from any cwd.
* All scripts use `load_all("/Users/jinwoocho/Desktop/rsiena")` —
  update if you clone the repo to a different path.
* `*.log` files inside `diagnostics/` and `artifacts/` are debug run
  captures from this session; not regenerated automatically.

## Provenance

* Files like `simulation_crprod_recovery.R` from `tests/` root were
  renamed during this reorganization. Use `git log --follow` to track
  history across moves:
  ```sh
  git log --follow tests/group3_self_to_perception/recovery_crprod_dir1.R
  ```
