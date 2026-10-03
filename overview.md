# Screening for a chronic disease (TreeAge model)

This model is not written in R: it is the TreeAge Pro model `screening.trex`, read
and run as it is by `treeage.R`. The structure, the formulas and the base values
below all come from that file, which does not document where its values come from,
so treat them as illustrative unless their source is known.

## Model structure

A decision between two strategies, each a Markov cohort model of 50 cycles
(`_stage = totalCycles`):

- **No screening**: the disease is only found once it gives symptoms.
- **Screening Imperfect**: a test, with a sensitivity and a specificity below
  one, is offered every `screen_int` cycles, starting at the first one (cycles
  0, 4, 8, ... at the base value).

The states are:

- **Healthy**: free of the disease. 90% of the cohort starts here.
- **Stage 1 undiagnosed**: early disease nobody knows about. The other 10% of the
  cohort (`prev`) starts here. It costs nothing, since it is not treated.
- **Stage 1 diagnosed - Tx** (screening only): early disease found by a screen,
  and treated at `cStage1` per cycle.
- **Stage 2 symptomatic**: advanced disease, treated at `cStage2` per cycle.
- **Dead**: absorbing.

Each cycle the healthy develop the disease with probability `pInc_ann`,
undiagnosed stage 1 progresses to stage 2 with probability `pProg2_UnDx` and
treated stage 1 with `pProg2_Dx`, and stage 2 is fatal with probability `pDie`.
There is no recovery and no death from other causes.

In a screening cycle, everyone healthy or in undiagnosed stage 1 is screened at
`cScreen`. A healthy person tests positive with probability `1 - test_spec`, and
the work-up of that false positive costs `cStage1 / 4`. Someone in undiagnosed
stage 1 tests positive with probability `test_sens`, and moves to treated stage 1
unless the disease progresses to stage 2 that same cycle.

A year in good health is worth 1 QALY, a year in stage 1 (diagnosed or not)
`uStage1` and a year in stage 2 `uStage2`. Nothing is discounted.

## Rewards and outcomes

The model tracks eight reward sets: cost, QALYs, life years, time healthy, time in
stage 1, time in stage 2, deaths and number of screens. The cost-effectiveness
analysis compares the expected cost (`C`) and QALYs (`E`) of each strategy; the
others are reported in `rewards`.

The model uses TreeAge's within-cycle correction: the cost and utility of a
state in a cycle are those of the average of the cohort in it at the start and at
the end of the cycle. Screening, false positives and deaths are transition
rewards, accrued by the cohort taking that path during the cycle. The file also
holds traditional (initial, incremental, final) rewards for every state, among
them a `cStage1` cost for undiagnosed stage 1, which are not used when
within-cycle correction is on.

## Calibration

The `stage1_incidence` scheme is an example of how a scheme is set up for a
TreeAge model, with invented targets. It calibrates `pInc_ann`, one value per
stratum of 5 cycles, so that the no-screening strategy reproduces a target
incidence of stage 1: the mean, over the cycles of the stratum, of the share of
the cohort alive at the start of a cycle that enters undiagnosed stage 1 during
it.

The targets fall from 0.035 in the first stratum to 0.016 in the last, and are
reached exactly by values rising from about 0.042 to 0.133. Since nobody returns
to *Healthy*, the healthy pool shrinks as the disease builds up, so each stratum
needs a higher probability than the last to keep up even a falling incidence, and
it depends on every earlier one through the size of that pool: the strata cannot
be fitted one at a time. Incidence targets that kept rising with time can be out
of reach in the later strata, whatever the probability.

## Checks

`tests/run_tests.R` compares every reward set of both strategies with a
hand-written cohort model of the same structure, at the base values and at
others. Given the reference file, it also checks the engine against results
TreeAge Pro itself reported for a larger Markov model with traditional rewards,
discounting and clones.
