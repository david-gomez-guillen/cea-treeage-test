source('treeage.R')

# ==== Configuration =============================================================
#
# The model is whatever TreeAge model TREX.FILE names: its strategies, parameters,
# states and rewards are all read from it when this file is sourced. Pointing it
# at another .trex file is all it takes to load another decision tree or Markov
# cohort model (the calibration scheme at the end of this file is the one part
# written for screening.trex, and is left out for any other model).

TREX.FILE <- 'screening.trex'

# The strata are groups of this many Markov cycles: 0-4, 5-9, ... A parameter
# split into strata takes, at each cycle, the value of the stratum of that cycle,
# and the cohort outputs are also reported by stratum.
CYCLES.PER.STRATUM <- 5

MODEL <- read.trex(TREX.FILE)
PARAMETERS <- trex.parameters(MODEL)
STRATEGIES <- trex.strategies(MODEL)

# The strata cover the cycles of the longest Markov process of the model with
# every parameter at its base value. A run made longer by its parameters counts
# the cycles past the last stratum as part of it. A model with no Markov process
# has a single stratum.
make.strata <- function(horizon, size) {
  if (horizon == 0) return(list(names = 'all', stratum.of = function(stage) 1L))
  starts <- seq(0, horizon - 1, by = size)
  ends <- pmin(starts + size - 1, horizon - 1)
  names <- ifelse(starts == ends, as.character(starts), paste0(starts, '-', ends))
  list(names = names, stratum.of = function(stage) pmin(stage %/% size + 1L, length(names)))
}
STRATA <- make.strata(trex.horizon(MODEL), CYCLES.PER.STRATUM)


# ==== Interface =================================================================

get.overview <- function() {
  # Markdown shown in the Overview tab: the text of overview.md, if there is one,
  # followed by a description of the model generated from the .trex file, which
  # is always in sync with what the model runs.
  text <- if (file.exists('overview.md')) paste(readLines('overview.md'), collapse = '\n') else ''
  paste(c(text, describe.model()), collapse = '\n\n')
}

get.model.settings <- function() {
  ce <- trex.ce.sets(MODEL)
  list(
    cost.unit = MODEL$preferences$currencySymbol %or% '€',
    effectiveness.unit = trex.reward.names(MODEL)[match(ce[['effect']], MODEL$reward.sets)],
    run.initial.guess = TRUE
  )
}

get.strategies <- function() {
  lapply(seq_len(nrow(STRATEGIES)), function(i) list(name = STRATEGIES$name[i], display.name = STRATEGIES$display.name[i]))
}

get.parameters <- function() {
  # The parameters are the constants defined at the root of the tree. They are
  # shown by the label of their TreeAge variable, and grouped by its category or,
  # without one, by what the usual prefix of their name says they are.
  labels <- vapply(PARAMETERS, function(p) if (p$label != '') p$label else p$comment, character(1))
  shared <- labels != '' & (duplicated(labels) | duplicated(labels, fromLast = TRUE))
  labels[shared] <- paste0(labels[shared], ' (', names(PARAMETERS)[shared], ')')

  lapply(seq_along(PARAMETERS), function(i) {
    p <- PARAMETERS[[i]]
    parameter <- list(name = p$name, base.value = p$base.value, class = parameter.class(p))
    if (labels[i] != '') parameter$display.name <- labels[i]
    parameter
  })
}

parameter.class <- function(p) {
  if (p$category != '') return(p$category)
  if (grepl('^(p|prob)([A-Z_]|$)', p$name)) return('Probabilities')
  if (grepl('^(c|cost)([A-Z_]|$)', p$name)) return('Costs')
  if (grepl('^(u|util)([A-Z_]|$)', p$name)) return('Utilities')
  'Other'
}

get.strata <- function() {
  STRATA$names
}

run.simulation <- function(strategies, pars) {
  # The parameters arrive as a list with one value each, or a list of values by
  # stratum for a parameter split into strata, which is what trex.simulate()
  # takes. The cohort outputs are added up by stratum on top of what it returns.
  pars <- pars[intersect(names(pars), names(PARAMETERS))]
  results <- trex.simulate(MODEL, strategies, pars, STRATA)

  by.stratum <- function(table, summarise) {
    if (is.null(table)) return(NULL)
    stratum <- STRATA$stratum.of(table$stage)
    out <- t(vapply(seq_along(STRATA$names), function(k) {
      rows <- table[stratum == k, -1, drop = FALSE]
      if (nrow(rows) == 0) rep(NA_real_, ncol(rows)) else summarise(rows)
    }, numeric(ncol(table) - 1)))
    if (ncol(table) == 2) out <- t(out)
    dimnames(out) <- list(STRATA$names, colnames(table)[-1])
    out
  }
  # The cohort in each state over a stratum, as its mean over the cycles of the
  # stratum (at their start), and the share of it entering each state.
  results$occupancy <- lapply(results$trace, function(trace) by.stratum(trace[trace$stage < max(trace$stage), ], colMeans))
  results$entries <- lapply(results$inflow, function(inflow) by.stratum(inflow, colSums))
  # The details of every Markov process are left out: they are large and the
  # app copies the results of every run between processes.
  results$markov <- NULL
  results
}

get.model.states <- function() {
  # One diagram per distinct Markov process: processes with the same states and
  # transitions, such as the copies of a clone, share one.
  diagrams <- list()
  for (id in MODEL$markov.nodes) {
    structure <- MODEL$states[[as.character(id)]]
    transitions <- trex.markov.transitions(MODEL, id)
    transitions <- transitions[transitions$from != transitions$to, ]
    signature <- paste(c(structure$labels, transitions$path), collapse = '|')
    path <- markov.path(id)
    if (!is.null(diagrams[[signature]])) {
      diagrams[[signature]]$also <- c(diagrams[[signature]]$also, path)
      next
    }

    nodes <- MODEL$nodes
    describe.state <- function(i) {
      state <- nodes[[structure$ids[i]]]
      lines <- paste0('Initial probability: ', state$source$prob %or% '#')
      fields <- if (trex.uses.wcc(MODEL)) c(Startup = 'startup', Rate = 'per cycle', Fixed = 'event') else c(Init = 'initial', Incr = 'per cycle', Final = 'final')
      reward.names <- trex.reward.names(MODEL)
      for (set in names(state$markov$state.rewards)) {
        values <- unlist(state$markov$state.rewards[[set]][names(fields)])
        used <- !is.na(values) & trimws(values) != '0'
        if (!any(used)) next
        name <- reward.names[match(as.integer(set), MODEL$reward.sets)]
        lines <- c(lines, paste0(name, ': ', paste0(values[used], ' (', fields[used], ')', collapse = ', ')))
      }
      paste(lines, collapse = '\n')
    }
    pairs <- unique(transitions[, c('from', 'to')])
    edges <- data.frame(
      from = pairs$from, to = pairs$to,
      label = vapply(seq_len(nrow(pairs)), function(k) {
        rows <- transitions[transitions$from == pairs$from[k] & transitions$to == pairs$to[k], ]
        paste(unique(sub('^.* > ', '', rows$path)), collapse = ' / ')
      }, character(1)),
      description = vapply(seq_len(nrow(pairs)), function(k) {
        rows <- transitions[transitions$from == pairs$from[k] & transitions$to == pairs$to[k], ]
        paste0(rows$path, ': ', rows$prob, collapse = '\n')
      }, character(1)),
      stringsAsFactors = FALSE
    )
    diagrams[[signature]] <- list(
      title = path,
      also = character(0),
      nodes = data.frame(id = structure$labels, label = structure$labels,
                         description = vapply(seq_along(structure$ids), describe.state, character(1)),
                         stringsAsFactors = FALSE),
      edges = edges
    )
  }
  diagrams <- lapply(diagrams, function(d) {
    d$description <- paste0('The states of the Markov process **', d$title, '** and the transitions between them, ',
                            'each labelled by the last branch of its path. Hover a state for its rewards, ',
                            'or a transition for every path it stands for and its probability.')
    if (length(d$also) > 0) d$description <- paste0(d$description, '\n\nAlso used by: ', paste(d$also, collapse = '; '), '.')
    d$also <- NULL
    d
  })
  names(diagrams) <- vapply(diagrams, function(d) d$title, character(1))
  diagrams
}

get.code.sample <- function() {
  # The tree as read from the .trex file, and the code that runs it.
  list(
    `Model tree` = list(code = trex.tree.text(MODEL), language = 'text'),
    `treeage.R` = readLines('treeage.R')
  )
}


# ==== Description of the model ==================================================

# The labels from the root to a Markov node, the root left out.
markov.path <- function(id) {
  labels <- character(0)
  while (!is.na(id) && id != 1L) {
    labels <- c(MODEL$nodes[[id]]$label, labels)
    id <- MODEL$nodes[[id]]$parent
  }
  paste(labels, collapse = ' > ')
}

describe.model <- function() {
  prefs <- MODEL$preferences
  ce <- trex.ce.sets(MODEL)
  reward.names <- trex.reward.names(MODEL)
  enabled <- MODEL$reward.sets <= as.integer(prefs$numberOfEnabledPayoffs %or% length(MODEL$reward.sets))
  md.escape <- function(x) gsub('([|*_])', '\\\\\\1', x)

  lines <- c(
    paste0('## Model read from `', basename(TREX.FILE), '`'),
    '',
    paste0('Everything below is read from the TreeAge file when the model is loaded, and is what `run.simulation()` runs. ',
           'The root of the tree is **', md.escape(MODEL$nodes[[1]]$label), '**.'),
    '',
    '### Strategies',
    '',
    paste0('- **', md.escape(STRATEGIES$display.name), '** (`', STRATEGIES$name, '`)'),
    '',
    '### Outcomes',
    '',
    paste0('The cost of a strategy (`C`) is its expected **', reward.names[match(ce[['cost']], MODEL$reward.sets)],
           '** and its effect (`E`) its expected **', reward.names[match(ce[['effect']], MODEL$reward.sets)],
           '**. Every reward set is reported in `rewards`: ', paste0(reward.names[enabled], collapse = ', '), '.')
  )

  if (length(MODEL$markov.nodes) > 0) {
    discounting <- if (identical(prefs$useGlobalDiscounting, 'true')) {
      rates <- vapply(MODEL$reward.sets[enabled], function(s) prefs[[paste0('globalDiscountingDiscountRate', s)]] %or% '0', character(1))
      paste0('Rewards are discounted (global discounting) at a yearly rate of ', paste0(rates, ' for ', reward.names[enabled], collapse = ', '),
             ', with cycles of ', prefs$globalDiscountingMarkovCycleLength %or% '1', ' years.')
    } else 'Rewards are not discounted by the model (global discounting is off), unless its own expressions do it.'
    correction <- if (trex.uses.wcc(MODEL)) {
      'The model uses **within-cycle correction**: startup rewards are accrued in the first cycle, event rewards by the cohort at the start of each cycle, and cycle rewards by the average of the cohort at the start and at the end of each cycle (trapezoidal rule).'
    } else {
      'The model uses the **traditional** rewards: initial rewards in the first cycle (`_stage` 0), incremental rewards in every other cycle and final rewards on the cohort once the process ends. Any half-cycle correction is the one the modeller wrote into them.'
    }
    # Each process with its states, which are only listed the first time a set of
    # them comes up (copies of a clone share them).
    seen <- character(0)
    processes <- vapply(MODEL$markov.nodes, function(id) {
      structure <- MODEL$states[[as.character(id)]]
      states <- paste(md.escape(structure$labels), collapse = ', ')
      listed <- if (states %in% seen) 'the same states as above' else states
      seen <<- c(seen, states)
      paste0('- **', md.escape(markov.path(id)), '**: ', length(structure$ids), ' states (', listed,
             '), until `', MODEL$nodes[[id]]$termination, '`.')
    }, character(1))
    lines <- c(lines, '', '### Markov processes', '',
               paste0('The model goes through ', length(MODEL$markov.nodes), ' Markov cohort process',
                      if (length(MODEL$markov.nodes) > 1) 'es' else '', ', of ', trex.horizon(MODEL),
                      ' cycles at most with the base values. ', correction, ' ', discounting),
               '', processes)
  }

  lines <- c(lines, '', '### Parameters', '',
             'The constants defined at the root of the tree. A distribution is used at its mean, as in a TreeAge cohort analysis.',
             '', '| Parameter | Base value | Definition in TreeAge | Description |', '|---|---|---|---|',
             vapply(PARAMETERS, function(p) paste0('| `', p$name, '` | ', signif(p$base.value, 6), ' | `', p$expression, '` | ',
                                                    md.escape(if (p$label != '') p$label else p$comment), ' |'), character(1)))

  formulas <- trex.formulas(MODEL)
  if (length(formulas) > 0)
    lines <- c(lines, '', '### Formulas', '',
               'Variables defined at the root as formulas, recomputed from the parameters at every node and cycle they are used in.',
               '', '| Variable | Definition |', '|---|---|',
               paste0('| `', names(formulas), '` | `', gsub('|', '\\|', formulas, fixed = TRUE), '` |'))

  lines <- c(lines, '', '### Strata', '',
             if (length(STRATA$names) == 1 && STRATA$names == 'all') 'The model has no Markov process, so it has a single stratum, `all`.'
             else paste0('The strata are groups of ', CYCLES.PER.STRATUM, ' Markov cycles, counted by `_stage` from 0: ',
                         paste0('`', STRATA$names, '`', collapse = ', '),
                         '. A parameter split into strata takes, at each cycle, the value of the stratum of that cycle; ',
                         'cycles past the last stratum count as part of it. `occupancy` gives the mean share of the cohort in each state over the cycles of a stratum, ',
                         'and `entries` the share entering each state from another one during them.'))
  paste(lines, collapse = '\n')
}


# ==== Calibration ===============================================================
#
# An example scheme, written for screening.trex: it calibrates the annual
# probability of developing the disease, one value per stratum, so that the
# no-screening strategy reproduces a target incidence of stage 1 disease. It is
# only offered when the model has what it needs, so another model loads without
# a Calibration tab until a scheme is written for it.

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
