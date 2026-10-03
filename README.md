# TreeAge models in THALASSA

A THALASSA model whose logic lives in a TreeAge Pro file. Instead of writing the
model in R, as [cea-model-test](https://github.com/david-gomez-guillen/cea-model-test)
does, `thalassa_interface.R` reads a `.trex` file and runs it with a small TreeAge
engine written in R (`treeage.R`). The strategies, parameters, strata, states and
rewards all come from the file, so the same interface works for any decision tree
or Markov cohort model saved by TreeAge.

| File | What it is |
|---|---|
| `screening.trex` | The TreeAge model loaded by default. |
| `treeage.R` | Reads a `.trex` file and evaluates it as TreeAge's expected value (cohort) analysis does. |
| `thalassa_interface.R` | The THALASSA interface on top of it, plus an example calibration scheme for `screening.trex`. |
| `overview.md` | Hand-written description of `screening.trex`, shown before the description generated from the file. |
| `tests/run_tests.R` | Tests: expressions, a hand-written version of `screening.trex`, and TreeAge's own results for a reference model. |

## Using another TreeAge model

1. Save the model from TreeAge Pro as a `.trex` file and put it in the repository.
2. Set `TREX.FILE` at the top of `thalassa_interface.R` to its name, and
   `CYCLES.PER.STRATUM` to the size of the strata you want.
3. Replace `overview.md` with a description of the model, or delete it: the
   description generated from the file is shown either way.
4. Write a calibration scheme for it if you need one. The one in
   `thalassa_interface.R` is specific to `screening.trex` and is not offered for
   another model.

How the interface maps the TreeAge model to THALASSA:

- **Strategies** are the branches of the root decision node (or the root itself
  when it is not one), named after their labels.
- **Parameters** are the variables defined at the root as constants, such as `0.05`,
  `6358/2` or `DistSamp(1)` (a distribution is used at its mean). Formulas, and
  variables defined at other nodes, stay in the model and are recomputed from them.
  They are shown by the label of their TreeAge variable and grouped by its category,
  or by name prefix (`p`, `c`, `u`) without one.
- **Strata** are groups of `CYCLES.PER.STRATUM` Markov cycles. A parameter split
  into strata takes, at each cycle, the value of its stratum.
- **C** and **E** are the cost and effectiveness reward sets of the model's
  cost-effectiveness settings (the first two sets otherwise). `run.simulation()`
  also returns every reward set, the cohort trace and the flows into each state,
  per cycle and per stratum.
- **State diagrams** are drawn for every distinct Markov process, and the Code
  panel shows the tree as read from the file.

## What the engine supports

Decision, chance, label, terminal and Markov nodes; clones; variables defined at
any node (the nearest definition towards the root, evaluated where the variable is
used, as in TreeAge); traditional Markov rewards (initial, incremental, final) and
within-cycle correction (startup, cycle, event, trapezoidal rule); transition
rewards; global discounting; distributions at their mean; `#` probabilities; the
keywords `_stage` and `_stage_mid`; and the functions `if`, `Modulo`, `Min`, `Max`,
`Exp`, `Ln`, `Log`, `Round`, `Discount`, `RateToProb`, `ProbToRate` and a few more
(see `TREX.FUNCTIONS` in `treeage.R`).

Not supported, with an error saying so: tables, tunnel states (`_tunnel`), tracker
variables in expressions (tracker modifications are ignored, as they only matter
to microsimulation) and Markov nodes inside a Markov process. Within-cycle
correction always uses the trapezoidal rule, even in a model set up for one of
Simpson's rules.

Where TreeAge's documentation leaves the semantics open, the engine assumes that,
with within-cycle correction, the end-of-cycle half of a cycle reward uses the
expression evaluated at that cycle and is discounted at the next one. This only
matters for cycle rewards that depend on `_stage` or are discounted.

## Tests

```sh
Rscript tests/run_tests.R
```

The check against TreeAge's results uses a model that is not in this repository.
Point `TREEAGE_REFERENCE_DIR` at the folder holding
`MammaPrint-ICO-REAL_cyclelengthsOK-incr6months.trex` to run it; the engine
reproduces the costs and QALYs TreeAge reported for it, over three horizons, to
machine precision.

The model needs the `xml2` package, which THALASSA installs from conda-forge on its
own when it builds the model's environment.
