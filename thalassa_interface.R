# The THALASSA interface for screening.trex.
#
# This file holds everything that belongs to one TreeAge model; the interface
# itself, which reads the strategies, parameters, states and rewards from the
# .trex file, is in treeage_interface.R and is the same for every model. To load
# another model, change this file: point TREX.FILE at it, choose the strata,
# give its parameters display names and write calibration schemes for it (or
# delete the ones below).

# ==== Model =====================================================================

TREX.FILE <- 'screening.trex'

# The strata are groups of this many Markov cycles: 0-4, 5-9, ... A parameter
# split into strata takes, at each cycle, the value of the stratum of that cycle,
# and the cohort outputs are also reported by stratum.
CYCLES.PER.STRATUM <- 5

# The names the parameters are shown by, by the name of their TreeAge variable.
# They take the place of the labels of the variables in TreeAge, which
# screening.trex leaves empty. A parameter left out is shown by its TreeAge
# label, if it has one. The cycles of this model are years.
PARAMETER.DISPLAY.NAMES <- list(
  prev = 'Prevalence of stage 1 disease at the start',
  pInc_ann = 'Annual probability of developing the disease',
  pProg2_UnDx = 'Annual probability of progression to stage 2, undiagnosed',
  pProg2_Dx = 'Annual probability of progression to stage 2, treated',
  pDie = 'Annual probability of death in stage 2',
  test_sens = 'Sensitivity of the screening test',
  test_spec = 'Specificity of the screening test',
  screen_int = 'Screening interval (years)',
  cScreen = 'Cost of a screening test',
  cStage1 = 'Annual cost of treating stage 1',
  cStage2 = 'Annual cost of treating stage 2',
  uStage1 = 'Utility of stage 1',
  uStage2 = 'Utility of stage 2',
  totalCycles = 'Time horizon (years)'
)

# The parameters not shown in THALASSA, which always run at their base value.
HIDDEN.PARAMETERS <- c('totalCycles')

source('treeage_interface.R')


# ==== Calibration ===============================================================
#
# An example scheme: it calibrates the annual probability of developing the
# disease, one value per stratum, so that the no-screening strategy reproduces a
# target incidence of stage 1 disease. It is only offered when the model has
# what it needs, so pointing TREX.FILE at another model loads it without a
# Calibration tab until a scheme is written for it.

CALIBRATION <- list(
  parameter = 'pInc_ann',
  strategy = 'no_screening',
  state = 'Stage 1 undiagnosed - no Tx',
  dead = 'Dead',
  # Invented targets, one per stratum of 5 cycles, to show how a scheme is set
  # up. They are reached exactly by values of pInc_ann rising from about 0.042
  # in the first stratum to 0.133 in the last.
  target = c(`0-4` = 0.035, `5-9` = 0.031, `10-14` = 0.028, `15-19` = 0.025, `20-24` = 0.023,
             `25-29` = 0.021, `30-34` = 0.019, `35-39` = 0.018, `40-44` = 0.017, `45-49` = 0.016)
)

# Incidence of the state over a stratum: the mean, over the cycles of the
# stratum, of the share of the cohort alive at the start of the cycle that enters
# the state during it.
calibration.incidence <- function(results) {
  trace <- results$trace[[CALIBRATION$strategy]]
  inflow <- results$inflow[[CALIBRATION$strategy]]
  alive <- 1 - trace[[CALIBRATION$dead]][seq_len(nrow(inflow))]
  yearly <- inflow[[CALIBRATION$state]] / alive
  stratum <- STRATA$stratum.of(inflow$stage)
  setNames(vapply(seq_along(STRATA$names), function(k) mean(yearly[stratum == k]), numeric(1)), STRATA$names)
}

calibration.error <- function(pars, target) {
  target.incidence <- unlist(target$`Stage 1 incidence`)
  tryCatch({
    incidence <- calibration.incidence(run.simulation(CALIBRATION$strategy, pars))
    # Only the strata present in the target contribute to the error.
    list(error = sum((incidence[names(target.incidence)] - target.incidence)^2),
         output = list(`Stage 1 incidence` = incidence))
  }, error = function(e) {
    list(error = Inf, output = list(`Stage 1 incidence` = setNames(rep(NA_real_, length(STRATA$names)), STRATA$names)))
  })
}

generate.training.dataset <- function(initial_guess, n, ...) {
  variation <- list(...)$variation
  dataset <- t(replicate(n, pmin(1, initial_guess * runif(length(initial_guess), 1 - variation, 1 + variation))))
  dataset[sample(nrow(dataset)), , drop = FALSE]
}

if (CALIBRATION$parameter %in% names(PARAMETERS) &&
    CALIBRATION$strategy %in% STRATEGIES$name &&
    any(vapply(MODEL$states, function(s) all(c(CALIBRATION$state, CALIBRATION$dead) %in% s$labels), logical(1))) &&
    identical(STRATA$names, names(CALIBRATION$target))) {
  get.calibration.schemes <- function() {
    list(
      stage1_incidence = list(
        description = 'Incidence of stage 1 disease without screening',
        parameters = CALIBRATION$parameter,
        strata = get.strata(),
        target = list(`Stage 1 incidence` = as.list(CALIBRATION$target)),
        initial_guess = rep(PARAMETERS[[CALIBRATION$parameter]]$base.value, length(STRATA$names)),
        error_function = calibration.error,
        latent_space_training_set = generate.training.dataset,
        other.plots = NULL
      )
    )
  }
}
