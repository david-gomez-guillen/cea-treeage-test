# Tests of the TreeAge engine. Run from the root of the repository:
#
#   Rscript tests/run_tests.R
#
# The check against TreeAge's own results needs a model that is not in the
# repository; set TREEAGE_REFERENCE_DIR to the folder holding
# MammaPrint-ICO-REAL_cyclelengthsOK-incr6months.trex to run it.

source('thalassa_interface.R')

failures <- 0
check <- function(what, ok) {
  cat(if (isTRUE(ok)) 'ok  ' else 'FAIL', what, '\n')
  if (!isTRUE(ok)) failures <<- failures + 1
}
near <- function(a, b, tolerance = 1e-9) isTRUE(all(abs(a - b) <= tolerance * pmax(1, abs(b))))


# ==== Expressions ===============================================================

value.of <- function(text, variables = list(), stage = 0, decimal.comma = FALSE) {
  env <- new.env(parent = TREX.BASE.ENV)
  env$.stage <- stage
  env$.v <- function(name) variables[[name]]
  eval(trex.compile(text, decimal.comma), env)
}
check('arithmetic and precedence', value.of('1 + 2 * 3 ^ 2 / 6') == 4)
check('unary minus binds looser than ^', value.of('-2 ^ 2') == -4)
check('comparison with = and <>', value.of('(3 = 3) + (3 <> 3)') == 1)
check('if() with ; separators', value.of('if(_stage < 10; 1; 2)', stage = 12) == 2)
check('nested if()', value.of('if(_stage=1; 645; if(_stage=2; 597; 258))', stage = 2) == 597)
check('Modulo() of a variable', value.of('Modulo(_stage; screen_int) = 0', list(screen_int = 4), stage = 8) == 1)
check('logical operators', value.of('(1 & 0) | !0') == 1)
check('case-insensitive functions', value.of('MAX(1; 3) + min(2; 5)') == 5)
check('decimal comma', value.of('0,3 + 0,2', decimal.comma = TRUE) == 0.5)
check('constants are recognised', trex.is.constant(trex.compile('6358/2')) && !trex.is.constant(trex.compile('pRec * 2')))
check('unsupported functions are refused', inherits(try(trex.compile('Foo(1)'), silent = TRUE), 'try-error'))


# ==== screening.trex against a hand-written cohort model =========================
#
# The same model written directly as transition matrices, from its description
# in overview.md, with within-cycle correction by the trapezoidal rule.

hand.model <- function(screening, p) {
  states <- if (screening) c('H', 'S1u', 'S1d', 'S2', 'D') else c('H', 'S1u', 'S2', 'D')
  n <- length(states)
  cohort <- setNames(numeric(n), states)
  cohort[c('H', 'S1u')] <- c(1 - p$prev, p$prev)
  # Rewards accrued per cycle by each state, one column per reward set.
  per.cycle <- matrix(0, n, 8, dimnames = list(states, c('Cost', 'QALY', 'LYs', 'Time Healthy', 'Time Stage 1', 'Time Stage 2', 'Death', 'Screens')))
  per.cycle['H', c('QALY', 'LYs', 'Time Healthy')] <- 1
  per.cycle['S1u', c('QALY', 'LYs', 'Time Stage 1')] <- c(p$uStage1, 1, 1)
  if (screening) per.cycle['S1d', c('Cost', 'QALY', 'LYs', 'Time Stage 1')] <- c(p$cStage1, p$uStage1, 1, 1)
  per.cycle['S2', c('Cost', 'QALY', 'LYs', 'Time Stage 2')] <- c(p$cStage2, p$uStage2, 1, 1)

  total <- setNames(numeric(8), colnames(per.cycle))
  for (stage in 0:(p$totalCycles - 1)) {
    screen <- screening && stage %% p$screen_int == 0
    P <- matrix(0, n, n, dimnames = list(states, states))
    P['H', 'S1u'] <- p$pInc_ann
    P['S1u', 'S2'] <- p$pProg2_UnDx
    if (screen) P['S1u', 'S1d'] <- p$test_sens * (1 - p$pProg2_UnDx)
    if (screening) P['S1d', 'S2'] <- p$pProg2_Dx
    P['S2', 'D'] <- p$pDie
    diag(P) <- 1 - rowSums(P)

    # Transition rewards: the screen and its false positives, and the deaths.
    if (screen) {
      total['Cost'] <- total['Cost'] + (cohort['H'] + cohort['S1u']) * p$cScreen +
        cohort['H'] * (1 - p$test_spec) * p$cStage1 / 4
      total['Screens'] <- total['Screens'] + cohort['H'] + cohort['S1u']
    }
    total['Death'] <- total['Death'] + cohort['S2'] * p$pDie

    following <- drop(cohort %*% P)
    total <- total + (colSums(cohort * per.cycle) + colSums(following * per.cycle)) / 2
    cohort <- following
  }
  total
}

base <- setNames(lapply(PARAMETERS, function(p) p$base.value), names(PARAMETERS))
results <- run.simulation(c('no_screening', 'screening_imperfect'), base)
for (strategy in c('no_screening', 'screening_imperfect')) {
  expected <- hand.model(strategy == 'screening_imperfect', base)
  got <- unlist(results$rewards[results$rewards$strategy == strategy, names(expected)])
  check(paste('screening.trex,', strategy, 'matches the hand-written model in every reward set'), near(got, expected))
}

# The same with other parameter values, screening every 3 cycles.
changed <- modifyList(base, list(pInc_ann = 0.08, screen_int = 3, test_sens = 0.7, cScreen = 900, totalCycles = 30))
results <- run.simulation('screening_imperfect', changed)
check('screening.trex matches the hand-written model with changed parameters',
      near(unlist(results$rewards[1, -1]), hand.model(TRUE, changed)[names(results$rewards)[-1]]))

# A parameter split into strata takes, at each cycle, the value of its stratum.
stratified <- base
stratified$pInc_ann <- setNames(as.list(rep(0.05, length(STRATA$names))), STRATA$names)
check('a stratified parameter with one value everywhere changes nothing',
      near(run.simulation('no_screening', stratified)$summary$C, run.simulation('no_screening', base)$summary$C))
stratified$pInc_ann[['0-4']] <- 0
trace <- run.simulation('no_screening', stratified)$trace$no_screening
check('a stratified parameter applies to the cycles of its stratum only',
      near(trace$Healthy[1:6], rep(0.9, 6)) && near(trace$Healthy[7], 0.9 * 0.95))

check('invalid probabilities make the run fail', inherits(try(run.simulation('no_screening', modifyList(base, list(pInc_ann = 1.5))), silent = TRUE), 'try-error'))


# ==== Interface =================================================================

check('every strategy has a name and a display name', all(vapply(get.strategies(), function(s) !is.null(s$name) && !is.null(s$display.name), logical(1))))
check('hidden parameters are not shown', !any(HIDDEN.PARAMETERS %in% vapply(get.parameters(), function(p) p$name, character(1))))
check('every parameter has a display name', all(vapply(get.parameters(), function(p) nzchar(p$display.name %or% ''), logical(1))))
check('parameter display names are unique', !any(duplicated(Filter(nzchar, vapply(get.parameters(), function(p) p$display.name %or% '', character(1))))))
check('the summary has strategy, C and E', all(c('strategy', 'C', 'E') %in% names(results$summary)))
check('the overview renders', is.character(get.overview()) && nchar(get.overview()) > 0)
check('the model states have nodes and edges', all(vapply(get.model.states(), function(d) nrow(d$nodes) > 0 && nrow(d$edges) > 0, logical(1))))
scheme <- get.calibration.schemes()[[1]]
initial <- base
initial$pInc_ann <- setNames(as.list(scheme$initial_guess), STRATA$names)
evaluation <- scheme$error_function(initial, scheme$target)
check('the calibration error is a finite number with one output per stratum',
      is.finite(evaluation$error) && length(evaluation$output$`Stage 1 incidence`) == length(STRATA$names))


# ==== Against TreeAge's own results ==============================================
#
# Results TreeAge Pro reported for this model (sensitivity analysis on costChem,
# at its value of 5751.4667), for three horizons of 6-month cycles. They check the
# traditional rewards, global discounting, clones, label nodes and node-level
# definitions.

reference.dir <- Sys.getenv('TREEAGE_REFERENCE_DIR')
reference.file <- file.path(reference.dir, 'MammaPrint-ICO-REAL_cyclelengthsOK-incr6months.trex')
if (reference.dir != '' && file.exists(reference.file)) {
  reference <- read.trex(reference.file)
  expected <- list(
    `10` = c(7544.3686498167735, 3.5231134116055034, 8742.144929185446, 3.2910250724589614),
    `20` = c(11484.746519883789, 6.2817908530330415, 13627.904135678999, 5.739819191247468),
    `80` = c(16226.977385548078, 14.143041022822741, 19683.54652340052, 12.175718246312783)
  )
  for (cycles in names(expected)) {
    model <- reference
    for (id in model$markov.nodes) model$nodes[[id]]$termination.code <- call('==', as.name('.stage'), as.numeric(cycles))
    summary <- trex.simulate(model, parameters = list(costChem = 5751.466666666667))$summary
    check(paste('MammaPrint ICO-REAL over', cycles, 'cycles matches TreeAge'),
          near(c(summary$C[1], summary$E[1], summary$C[2], summary$E[2]), expected[[cycles]]))
  }
} else {
  cat('skip TreeAge reference results (TREEAGE_REFERENCE_DIR not set)\n')
}

cat('\n', if (failures == 0) 'All tests passed.' else paste(failures, 'test(s) failed.'), '\n')
quit(status = if (failures == 0) 0 else 1)
