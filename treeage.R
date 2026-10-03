# Reading and running TreeAge Pro models (.trex files) in R.
#
# A .trex file is the XML document TreeAge Pro saves a tree in. This file reads
# one into an R structure, and evaluates it the way TreeAge's expected value
# (cohort) analysis does, so that a decision tree or a Markov cohort model built
# in TreeAge can be run without TreeAge.
#
#   read.trex(path)                      the model: nodes, variables, preferences
#   trex.parameters(model)               the constants a user can change
#   trex.strategies(model)               the branches of the root decision node
#   trex.horizon(model)                  the number of Markov cycles at base values
#   trex.simulate(model, strategies, parameters, strata)
#                                        expected costs, effects and cohort traces
#
# What is supported:
#   - Decision, chance, label, terminal and Markov nodes, and clones.
#   - Variables defined at any node. As in TreeAge, a variable used at a node
#     takes the definition nearest to it on the way to the root, and that
#     definition is evaluated at the node where the variable is used.
#   - Markov cohort models with the traditional rewards (initial, incremental and
#     final, the half-cycle correction being whatever the modeller put in them)
#     or with within-cycle correction (startup, cycle and event rewards, the cycle
#     rewards integrated with the trapezoidal rule).
#   - Transition rewards, at any node of a transition subtree.
#   - Global discounting, with its cycle length, one rate per reward set.
#   - Distributions, through their mean (DistSamp(n) in a cohort analysis).
#   - The TreeAge functions listed in TREX.FUNCTIONS below.
#
# What is not: tables, tunnel states (_tunnel), trackers (they only matter to
# microsimulation, so their modifications are ignored and using one in an
# expression is an error) and Markov nodes inside a Markov process. A model that
# uses one of them fails to load or to run with a message saying which.
# Within-cycle correction always uses the trapezoidal rule, even in a model set
# up for one of Simpson's rules.

library(xml2)

# `a`, or `b` when `a` is missing, empty or blank.
`%or%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, '')) b else a


# ==== Expressions ===============================================================
#
# TreeAge expressions, such as `if(_stage < 10; pRec * 2; 0)`, are translated once
# into R calls, which are then evaluated as many times as the model needs them:
#
#   - a variable `x` becomes `.v("x")`, which looks the variable up in the tree,
#   - `_stage` becomes `.stage`, the current Markov cycle,
#   - a function `Name(a; b)` becomes `.fn_name(a, b)`, `if()` an R `if`.
#
# A probability given as `#`, the complement of its siblings, is not an
# expression and is kept as the string "#".

# Functions available to expressions, by lower-case TreeAge name.
TREX.FUNCTIONS <- list(
  abs = function(x) abs(x),
  average = function(...) mean(c(...)),
  ceil = function(x) ceiling(x),
  ceiling = function(x) ceiling(x),
  choose = function(i, ...) list(...)[[i]],
  discount = function(value, rate, time) value / (1 + rate)^time,
  exp = function(x) exp(x),
  floor = function(x) floor(x),
  int = function(x) trunc(x),
  ln = function(x) log(x),
  log = function(x, base = 10) log(x, base),
  log10 = function(x) log10(x),
  max = function(...) max(...),
  mean = function(...) mean(c(...)),
  min = function(...) min(...),
  mod = function(x, y) x %% y,
  modulo = function(x, y) x %% y,
  probtorate = function(p, time = 1) -log(1 - p) / time,
  ratetoprob = function(rate, time = 1) 1 - exp(-rate * time),
  round = function(x, digits = 0) round(x, digits),
  sqrt = function(x) sqrt(x),
  sum = function(...) sum(...),
  trunc = function(x) trunc(x)
)

# Functions resolved against the model being run rather than fixed.
TREX.MODEL.FUNCTIONS <- c('distsamp', 'distmean')

# Environment every compiled expression is evaluated in (through a child that
# holds the state of the run).
TREX.BASE.ENV <- local({
  env <- new.env(parent = baseenv())
  for (name in names(TREX.FUNCTIONS)) assign(paste0('.fn_', name), TREX.FUNCTIONS[[name]], envir = env)
  # Truth of a value: TreeAge has no booleans, any number other than 0 is true.
  env$.lgl <- function(x) !is.na(x) && x != 0
  env
})

# Splits an expression into tokens. `decimal.comma` reads `0,5` as a number, as
# TreeAge does in locales whose decimal separator is the comma; arguments are
# then only separated by `;`.
trex.tokenize <- function(text, decimal.comma = FALSE) {
  number <- if (decimal.comma) '^([0-9]+(,[0-9]*)?|,[0-9]+)([eE][-+]?[0-9]+)?' else '^([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][-+]?[0-9]+)?'
  patterns <- c(
    num = number,
    id = '^[A-Za-z_][A-Za-z0-9_.]*',
    op = '^(<=|>=|<>|!=|==|&&|\\|\\||[-+*/^=<>&|!])',
    lparen = '^\\(', rparen = '^\\)', sep = if (decimal.comma) '^;' else '^[;,]',
    hash = '^#', ws = '^\\s+'
  )
  tokens <- list()
  rest <- text
  while (nchar(rest) > 0) {
    matched <- FALSE
    for (type in names(patterns)) {
      m <- regmatches(rest, regexpr(patterns[[type]], rest, perl = TRUE))
      if (length(m) == 1 && nchar(m) > 0) {
        if (type != 'ws') tokens[[length(tokens) + 1]] <- list(type = type, value = m)
        rest <- substring(rest, nchar(m) + 1)
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      if (substr(rest, 1, 1) == '[') stop('TreeAge tables are not supported, in the expression: ', text, call. = FALSE)
      stop('Cannot read "', substr(rest, 1, 1), '" in the expression: ', text, call. = FALSE)
    }
  }
  tokens
}

# Translates a TreeAge expression into an R call (or a number, or "#").
trex.compile <- function(text, decimal.comma = FALSE) {
  if (is.null(text) || is.na(text)) return(NULL)
  text <- trimws(text)
  if (text == '') return(NULL)
  if (text == '#') return('#')

  tokens <- trex.tokenize(text, decimal.comma)
  pos <- 1
  peek <- function() if (pos <= length(tokens)) tokens[[pos]] else list(type = 'end', value = '')
  advance <- function() { token <- peek(); pos <<- pos + 1; token }
  fail <- function(what) stop(what, ' in the expression: ', text, call. = FALSE)
  expect <- function(type) {
    token <- advance()
    if (token$type != type) fail(paste0('Unexpected "', token$value, '"'))
    token
  }
  is.op <- function(token, ops) token$type == 'op' && token$value %in% ops
  is.word <- function(token, words) token$type == 'id' && toupper(token$value) %in% words

  parse.or <- function() {
    left <- parse.and()
    while (is.op(peek(), c('|', '||')) || is.word(peek(), 'OR')) {
      advance()
      left <- call('||', call('.lgl', left), call('.lgl', parse.and()))
    }
    left
  }
  parse.and <- function() {
    left <- parse.not()
    while (is.op(peek(), c('&', '&&')) || is.word(peek(), 'AND')) {
      advance()
      left <- call('&&', call('.lgl', left), call('.lgl', parse.not()))
    }
    left
  }
  parse.not <- function() {
    if (is.op(peek(), '!') || is.word(peek(), 'NOT')) {
      advance()
      return(call('!', call('.lgl', parse.not())))
    }
    parse.comparison()
  }
  parse.comparison <- function() {
    left <- parse.sum()
    ops <- c('=' = '==', '==' = '==', '<>' = '!=', '!=' = '!=', '<' = '<', '>' = '>', '<=' = '<=', '>=' = '>=')
    if (is.op(peek(), names(ops))) {
      op <- ops[[advance()$value]]
      left <- call(op, left, parse.sum())
    }
    left
  }
  parse.sum <- function() {
    left <- parse.product()
    while (is.op(peek(), c('+', '-'))) {
      op <- advance()$value
      left <- call(op, left, parse.product())
    }
    left
  }
  parse.product <- function() {
    left <- parse.unary()
    while (is.op(peek(), c('*', '/'))) {
      op <- advance()$value
      left <- call(op, left, parse.unary())
    }
    left
  }
  parse.unary <- function() {
    if (is.op(peek(), c('-', '+'))) {
      op <- advance()$value
      return(call(op, parse.unary()))
    }
    parse.power()
  }
  parse.power <- function() {
    base <- parse.primary()
    if (is.op(peek(), '^')) {
      advance()
      return(call('^', base, parse.unary()))
    }
    base
  }
  parse.primary <- function() {
    token <- advance()
    if (token$type == 'num') return(as.numeric(sub(',', '.', token$value, fixed = TRUE)))
    if (token$type == 'lparen') {
      inner <- parse.or()
      expect('rparen')
      return(inner)
    }
    if (token$type == 'hash') fail('"#" can only stand alone as a probability')
    if (token$type != 'id') fail(paste0('Unexpected "', token$value, '"'))

    name <- token$value
    if (peek()$type == 'lparen') {
      advance()
      args <- list()
      if (peek()$type != 'rparen') {
        repeat {
          args[[length(args) + 1]] <- parse.or()
          if (peek()$type != 'sep') break
          advance()
        }
      }
      expect('rparen')
      return(compile.function(name, args))
    }
    compile.identifier(name)
  }
  compile.function <- function(name, args) {
    key <- tolower(name)
    if (key == 'if') {
      if (length(args) != 3) fail('if() needs three arguments')
      return(call('if', call('.lgl', args[[1]]), args[[2]], args[[3]]))
    }
    if (!key %in% c(names(TREX.FUNCTIONS), TREX.MODEL.FUNCTIONS))
      fail(paste0('Unsupported TreeAge function ', name, '()'))
    as.call(c(as.name(paste0('.fn_', key)), args))
  }
  compile.identifier <- function(name) {
    key <- tolower(name)
    if (key == '_stage') return(as.name('.stage'))
    if (key == '_stage_mid') return(call('+', as.name('.stage'), 0.5))
    if (key %in% c('true', 'false')) return(as.numeric(key == 'true'))
    if (startsWith(name, '_')) fail(paste0('Unsupported TreeAge keyword ', name))
    call('.v', name)
  }

  result <- parse.or()
  if (peek()$type != 'end') fail(paste0('Unexpected "', peek()$value, '"'))
  result
}

# Whether a compiled expression is a constant: it uses no variable and no keyword.
# Distributions count as constants, since a cohort analysis uses their mean.
trex.is.constant <- function(code) {
  if (is.null(code) || identical(code, '#')) return(FALSE)
  !any(c('.v', '.stage') %in% all.names(code))
}

# Whether a compiled expression is the number 0, which lets a reward that is
# never anything else be skipped.
trex.is.zero <- function(code) is.null(code) || (is.numeric(code) && code == 0)


# ==== Reading the file ==========================================================

# The text of a TreeAge label, with its line breaks folded into spaces, which is
# also how a Markov transition names the state it jumps to.
trex.clean.label <- function(label) {
  if (is.na(label)) return('')
  gsub('\\s+', ' ', trimws(label))
}

# Reads a .trex file into a model.
read.trex <- function(path) {
  doc <- read_xml(path)
  ns <- xml_ns(doc)
  tree <- xml_find_first(doc, '//tree:Tree', ns)
  if (inherits(tree, 'xml_missing')) stop('No TreeAge tree in ', path)

  attr.of <- function(x, name) {
    value <- xml_attr(x, name)
    if (is.na(value)) NULL else value
  }

  # ---- Preferences --------------------------------------------------------------
  preferences <- list()
  for (p in xml_find_all(tree, './PreferenceSet/Preference')) preferences[[xml_attr(p, 'Name')]] <- xml_attr(p, 'Value')

  # A model saved in a locale with a decimal comma writes `0,5`. Without the
  # preference, it is recognised by a number written that way.
  decimal.comma <- if (!is.null(preferences$decimalSeparator)) {
    preferences$decimalSeparator == ','
  } else {
    values <- xml_attr(xml_find_all(tree, './/*[@Value]'), 'Value')
    any(grepl('^\\s*-?[0-9]+,[0-9]+\\s*$', values))
  }
  compile <- function(text) trex.compile(text, decimal.comma)

  # ---- Nodes --------------------------------------------------------------------
  # Every node of the tree as read from the file, clone copies still empty.
  read.node <- function(x) {
    node <- list(
      name.id = attr.of(x, 'NameID'),
      label = trex.clean.label(xml_attr(x, 'Label')),
      type = attr.of(x, 'NodeType'),
      clone.master = attr.of(x, 'CloneMasterIndex'),
      clone.master.name = attr.of(x, 'CloneMasterName'),
      clone.of = attr.of(x, 'AttachToCloneMaster'),
      definitions = list(), prob = NULL, payoffs = list(), termination = NULL,
      markov = NULL, children = list()
    )
    for (d in xml_find_all(x, './Definition')) node$definitions[[xml_attr(d, 'Variable')]] <- xml_attr(d, 'Value')
    prob <- xml_find_first(x, './Prob')
    if (!inherits(prob, 'xml_missing')) node$prob <- xml_attr(prob, 'Value')
    for (p in xml_find_all(x, './Payoff')) node$payoffs[[xml_attr(p, 'Set')]] <- xml_attr(p, 'Value')

    # A Markov node may hold one termination condition per reward set, which are
    # in practice the same; the one without a set, or the first, is used.
    terminations <- xml_find_all(x, './Termination')
    if (length(terminations) > 0) {
      sets <- xml_attr(terminations, 'Set')
      node$termination <- xml_attr(terminations[[if (any(is.na(sets))) which(is.na(sets))[1] else 1]], 'Value')
    }

    markov <- xml_find_first(x, './MarkovData')
    if (!inherits(markov, 'xml_missing')) {
      kind <- xml_attr(markov, 'xsi:type', ns)
      node$markov <- list(
        kind = if (is.na(kind)) NA else sub('^.*:', '', kind),
        jump = trex.clean.label(xml_attr(markov, 'MarkovJumpState')),
        state.rewards = list(),
        transition.rewards = list()
      )
      for (r in xml_find_all(markov, './StateReward')) {
        fields <- list()
        for (f in c('Init', 'Incr', 'Final', 'Startup', 'Fixed', 'Rate')) {
          value <- xml_find_first(r, paste0('./', f))
          fields[[f]] <- if (inherits(value, 'xml_missing')) NA else xml_attr(value, 'Value')
        }
        node$markov$state.rewards[[xml_attr(r, 'Set')]] <- fields
      }
      for (r in xml_find_all(markov, './TransitionReward')) {
        value <- xml_attr(r, 'Value')
        if (!is.na(value)) node$markov$transition.rewards[[xml_attr(r, 'Set')]] <- value
      }
    }

    node$children <- lapply(xml_find_all(x, './Node'), read.node)
    node
  }
  roots <- lapply(xml_find_all(tree, './Node'), read.node)
  if (length(roots) != 1) stop('Expected one root node in ', path, ', found ', length(roots))

  # ---- Clones -------------------------------------------------------------------
  # A clone copy is saved without children: they are those of its master, which
  # it shares while keeping its own label, probability and definitions.
  masters <- list()
  collect.masters <- function(node) {
    if (!is.null(node$clone.master)) masters[[node$clone.master]] <<- node
    for (child in node$children) collect.masters(child)
  }
  collect.masters(roots[[1]])
  expand.clones <- function(node, depth = 0) {
    if (depth > 100) stop('Clones nested too deep (a clone inside its own master?)')
    if (!is.null(node$clone.of) && length(node$children) == 0) {
      master <- masters[[node$clone.of]]
      if (is.null(master)) stop('Node "', node$label, '" is a copy of clone master ', node$clone.of, ', which is not in the file')
      node$children <- master$children
    }
    node$children <- lapply(node$children, expand.clones, depth = depth + 1)
    node
  }
  root <- expand.clones(roots[[1]])

  # ---- Flattening ---------------------------------------------------------------
  # The nodes are numbered depth first, and each one knows its parent and
  # children by number, which is what the evaluation walks.
  nodes <- list()
  flatten <- function(node, parent) {
    id <- length(nodes) + 1
    nodes[[id]] <<- NULL
    entry <- node[c('name.id', 'label', 'type', 'clone.master.name', 'clone.of', 'termination', 'markov')]
    entry$id <- id
    entry$parent <- parent
    entry$source <- list(definitions = node$definitions, prob = node$prob, payoffs = node$payoffs)
    entry$definitions <- lapply(node$definitions, compile)
    entry$prob <- compile(node$prob)
    entry$payoffs <- lapply(node$payoffs, compile)
    entry$termination.code <- compile(node$termination)
    if (!is.null(node$markov)) {
      entry$markov$state.code <- lapply(node$markov$state.rewards, function(fields) lapply(fields, compile))
      entry$markov$transition.code <- lapply(node$markov$transition.rewards, compile)
    }
    nodes[[id]] <<- entry
    children <- integer(0)
    for (child in node$children) children <- c(children, flatten(child, id))
    nodes[[id]]$children <<- children
    id
  }
  flatten(root, NA_integer_)

  # ---- Variables and distributions ----------------------------------------------
  # A variable names its category by the XMI id of the category, as "#<id>".
  categories <- list()
  for (c in xml_find_all(tree, './CategoriesRoot//Category')) categories[[xml_attr(c, 'xmi:id', ns)]] <- xml_attr(c, 'Name')
  variables <- list()
  for (v in xml_find_all(tree, './Variable')) {
    category <- sub('^#', '', xml_attr(v, 'Category'))
    variables[[xml_attr(v, 'NameID')]] <- list(
      label = trex.clean.label(xml_attr(v, 'Label')),
      comment = trex.clean.label(xml_attr(v, 'Comment')),
      category = if (is.na(category)) NULL else categories[[category]]
    )
  }
  trackers <- xml_attr(xml_find_all(tree, './Tracker'), 'NameID')
  distributions <- list()
  for (d in xml_find_all(tree, './Distribution')) {
    parameters <- list()
    for (p in xml_find_all(d, './Parameter')) parameters[[xml_attr(p, 'Name')]] <- xml_attr(p, 'Value')
    distributions[[length(distributions) + 1]] <- list(
      name = xml_attr(d, 'NameID'),
      label = trex.clean.label(xml_attr(d, 'Label')),
      index = xml_attr(d, 'Index'),
      type = xml_attr(d, 'Type'),
      parameters = parameters
    )
  }

  model <- list(
    path = path,
    nodes = nodes,
    variables = variables,
    trackers = trackers,
    distributions = distributions,
    preferences = preferences,
    decimal.comma = decimal.comma
  )
  model$markov.nodes <- which(vapply(nodes, function(n) identical(n$type, 'MarkovNode'), logical(1)))
  model$states <- lapply(setNames(model$markov.nodes, model$markov.nodes), function(id) trex.markov.structure(model, id))
  model$reward.sets <- trex.reward.sets(model)
  model$distribution.means <- vapply(distributions, trex.distribution.mean, numeric(1), decimal.comma = decimal.comma)
  model
}

# The states of a Markov node and, for every terminal node of their transition
# subtrees, the state it jumps to.
trex.markov.structure <- function(model, markov.id) {
  nodes <- model$nodes
  markov <- nodes[[markov.id]]
  state.ids <- markov$children
  labels <- vapply(state.ids, function(id) nodes[[id]]$label, character(1))
  if (any(duplicated(labels)))
    stop('Markov node "', markov$label, '" has two states named "', labels[duplicated(labels)][1], '"')

  find.state <- function(name, from) {
    index <- match(name, labels)
    if (is.na(index)) index <- match(tolower(name), tolower(labels))
    if (is.na(index)) stop('In Markov node "', markov$label, '", "', nodes[[from]]$label,
                           '" jumps to "', name, '", which is not one of its states')
    index
  }

  jumps <- integer(0)
  check.subtree <- function(id) {
    node <- nodes[[id]]
    if (identical(node$type, 'MarkovNode')) stop('Markov node "', node$label, '" inside Markov node "', markov$label, '" is not supported')
    if (length(node$children) == 0) {
      if (is.null(node$markov) || node$markov$jump == '')
        stop('In Markov node "', markov$label, '", the transition "', node$label, '" does not jump to any state')
      jumps[as.character(id)] <<- find.state(node$markov$jump, id)
    }
    for (child in node$children) check.subtree(child)
  }
  for (id in state.ids) for (child in nodes[[id]]$children) check.subtree(child)

  list(ids = state.ids, labels = labels, jumps = jumps,
       absorbing = vapply(state.ids, function(id) length(nodes[[id]]$children) == 0, logical(1)))
}

# The reward sets the model uses, which are at least those TreeAge has enabled.
trex.reward.sets <- function(model) {
  sets <- seq_len(max(1, as.integer(model$preferences$numberOfEnabledPayoffs %or% 1)))
  for (node in model$nodes) {
    sets <- c(sets, as.integer(names(node$payoffs)))
    if (!is.null(node$markov)) {
      used <- function(code) !all(vapply(code, trex.is.zero, logical(1)))
      sets <- c(sets, as.integer(names(Filter(used, node$markov$state.code))),
                as.integer(names(Filter(function(code) !trex.is.zero(code), node$markov$transition.code))))
    }
  }
  sort(unique(sets))
}

# The mean of a distribution, which is what a cohort analysis samples it at.
trex.distribution.mean <- function(distribution, decimal.comma = FALSE) {
  p <- lapply(distribution$parameters, function(value) {
    number <- suppressWarnings(as.numeric(if (decimal.comma) sub(',', '.', value, fixed = TRUE) else value))
    if (is.na(number)) value else number
  })
  names(p) <- tolower(names(p))
  get <- function(...) {
    for (name in c(...)) if (!is.null(p[[name]])) return(p[[name]])
    NULL
  }
  type <- tolower(distribution$type)
  mean <- switch(type,
    uniform = (get('low', 'min') + get('high', 'max')) / 2,
    normal = get('mean', 'mu'),
    triangular = (get('min', 'low') + get('likeliest', 'mode') + get('max', 'high')) / 3,
    beta = {
      alpha <- get('alpha', 'a'); beta <- get('beta', 'b')
      if (!is.null(alpha) && !is.null(beta)) alpha / (alpha + beta) else get('mean')
    },
    gamma = {
      alpha <- get('alpha'); lambda <- get('lambda', 'rate')
      if (!is.null(alpha) && !is.null(lambda)) alpha / lambda else get('mean')
    },
    lognormal = {
      mu <- get('mu'); sigma <- get('sigma')
      if (!is.null(mu) && !is.null(sigma)) exp(mu + sigma^2 / 2) else get('mean')
    },
    constant = get('value', 'mean'),
    get('mean')
  )
  if (is.null(mean) || !is.numeric(mean))
    stop('Cannot work out the mean of distribution "', distribution$name, '" (', distribution$type, ')')
  mean
}


# ==== What the model exposes ====================================================

# The parameters of the model: every variable defined at the root node as a
# constant, such as `0.05`, `6358/2` or `DistSamp(1)`. Definitions that are
# formulas of other variables or of the cycle stay in the model and are
# recomputed from these.
trex.parameters <- function(model) {
  root <- model$nodes[[1]]
  run <- trex.evaluator(model, list())
  parameters <- list()
  for (name in names(root$definitions)) {
    code <- root$definitions[[name]]
    if (!trex.is.constant(code)) next
    variable <- model$variables[[name]]
    parameters[[name]] <- list(
      name = name,
      base.value = run$eval(code, 1, 0),
      label = variable$label %or% '',
      comment = variable$comment %or% '',
      category = variable$category %or% '',
      expression = root$source$definitions[[name]]
    )
  }
  parameters
}

# The root-level formulas, the definitions that are not parameters.
trex.formulas <- function(model) {
  root <- model$nodes[[1]]
  keep <- !vapply(root$definitions, trex.is.constant, logical(1))
  unlist(root$source$definitions[keep])
}

# The strategies of the model: the branches of the root decision node, or the
# root itself when it is not a decision node.
trex.strategies <- function(model) {
  root <- model$nodes[[1]]
  ids <- if (identical(root$type, 'DecisionNode')) root$children else 1L
  labels <- vapply(ids, function(id) model$nodes[[id]]$label, character(1))
  names <- tolower(gsub('[^A-Za-z0-9]+', '_', iconv(labels, to = 'ASCII//TRANSLIT', sub = '')))
  names <- gsub('^_|_$', '', names)
  names[names == ''] <- paste0('strategy_', which(names == ''))
  names <- make.unique(names, sep = '_')
  data.frame(name = names, display.name = labels, node = ids, stringsAsFactors = FALSE)
}

# The name of each reward set, as the model's preferences call it.
trex.reward.names <- function(model) {
  prefs <- model$preferences
  vapply(model$reward.sets, function(set) {
    custom <- prefs[[paste0('customPayoffName', set)]]
    if (!identical(prefs$useCustomPayoffNames, 'false') && !is.null(custom) && custom != '') custom else paste('Payoff', set)
  }, character(1))
}

# The reward sets that are the cost and the effect of a cost-effectiveness
# analysis. A model not set up as one has its first set as the cost and its
# second, if there is one, as the effect.
trex.ce.sets <- function(model) {
  prefs <- model$preferences
  sets <- model$reward.sets
  if (identical(prefs$calcType, 'ct_costEff') && !is.null(prefs$ceCostPayoff) && !is.null(prefs$ceEffPayoff))
    return(c(cost = as.integer(prefs$ceCostPayoff), effect = as.integer(prefs$ceEffPayoff)))
  c(cost = sets[1], effect = if (length(sets) > 1) sets[2] else sets[1])
}

# Whether the Markov rewards follow within-cycle correction (startup, cycle and
# event rewards) rather than the traditional initial, incremental and final ones.
trex.uses.wcc <- function(model) identical(model$preferences$wccEnabled, 'true')

# The number of cycles of the longest Markov process of the model, with every
# parameter at its base value. 0 when the model has no Markov node.
trex.horizon <- function(model) {
  if (length(model$markov.nodes) == 0) return(0)
  run <- trex.evaluator(model, list())
  max.stages <- as.numeric(model$preferences$maxMarkovStages %or% 10000)
  max(vapply(model$markov.nodes, function(id) {
    code <- model$nodes[[id]]$termination.code
    if (is.null(code)) stop('Markov node "', model$nodes[[id]]$label, '" has no termination condition')
    stage <- 0
    while (!run$truth(code, id, stage)) {
      stage <- stage + 1
      if (stage > max.stages) stop('The termination condition of Markov node "', model$nodes[[id]]$label, '" is never met')
    }
    stage
  }, numeric(1)))
}


# ==== Evaluation ================================================================

# What evaluates expressions in one run of the model. `parameters` overrides the
# root definitions of the same name: a number, or a list with one value per
# stratum, the stratum of a cycle being given by `stratum.of(stage)`.
trex.evaluator <- function(model, parameters, stratum.of = function(stage) NA_character_) {
  nodes <- model$nodes
  env <- new.env(parent = TREX.BASE.ENV)
  env$.stage <- 0
  node <- 1L                  # the node the expression being evaluated belongs to
  shadowed <- list()          # variables being defined: name -> defining node
  cache <- new.env(hash = TRUE)

  parameter.value <- function(name, stage) {
    value <- parameters[[name]]
    if (is.list(value) || length(value) > 1) {
      stratum <- stratum.of(stage)
      value <- if (!is.na(stratum) && !is.null(value[[stratum]])) value[[stratum]] else NULL
      # A stratum the parameter has no value for keeps the one of the model.
      if (is.null(value)) return(NULL)
    }
    as.numeric(value)
  }

  # The value of a variable at the current node: the definition nearest to it on
  # the way to the root, evaluated at the current node. Inside a definition of
  # `x`, `x` refers to the definition above it.
  env$.v <- function(name) {
    start <- if (is.null(shadowed[[name]])) node else nodes[[shadowed[[name]]]]$parent
    key <- paste(start, node, name, env$.stage)
    if (!is.null(cache[[key]])) return(cache[[key]])

    at <- start
    while (!is.na(at) && is.null(nodes[[at]]$definitions[[name]])) at <- nodes[[at]]$parent
    value <- NULL
    if (!is.na(at) && at == 1L && name %in% names(parameters)) value <- parameter.value(name, env$.stage)
    if (is.null(value)) {
      if (is.na(at)) {
        index <- match(name, vapply(model$distributions, function(d) d$name, character(1)))
        if (is.na(index) && name %in% model$trackers)
          stop('Tracker "', name, '" is used at node "', nodes[[node]]$label, '", but trackers only exist in microsimulation', call. = FALSE)
        if (is.na(index)) stop('Variable "', name, '" is not defined at node "', nodes[[node]]$label, '"', call. = FALSE)
        value <- model$distribution.means[[index]]
      } else {
        previous <- shadowed[[name]]
        shadowed[[name]] <<- at
        value <- tryCatch(eval(nodes[[at]]$definitions[[name]], env), finally = { shadowed[[name]] <<- previous })
      }
    }
    cache[[key]] <- value
    value
  }
  dist.mean <- function(which) {
    index <- match(as.character(which), vapply(model$distributions, function(d) d$index %or% '', character(1)))
    if (is.na(index)) index <- match(as.character(which), vapply(model$distributions, function(d) d$name, character(1)))
    if (is.na(index)) stop('Distribution ', which, ' is not in the model', call. = FALSE)
    model$distribution.means[[index]]
  }
  env$.fn_distsamp <- function(which, ...) dist.mean(which)
  env$.fn_distmean <- function(which, ...) dist.mean(which)

  evaluate <- function(code, at, stage) {
    if (is.numeric(code)) return(code)
    node <<- at
    env$.stage <- stage
    value <- eval(code, env)
    if (is.logical(value)) value <- as.numeric(value)
    if (length(value) != 1 || is.na(value)) stop('The expression ', deparse1(code), ' at node "', nodes[[at]]$label, '" is not a number', call. = FALSE)
    value
  }

  list(
    eval = evaluate,
    truth = function(code, at, stage) evaluate(code, at, stage) != 0,
    # The probabilities of the children of a node, `#` (or no probability at
    # all) being what the others leave.
    child.probs = function(id, stage) {
      children <- nodes[[id]]$children
      probs <- numeric(length(children))
      complement <- integer(0)
      for (i in seq_along(children)) {
        code <- nodes[[children[i]]]$prob
        if (is.null(code) || identical(code, '#')) complement <- c(complement, i)
        else probs[i] <- evaluate(code, children[i], stage)
      }
      if (length(complement) > 1) stop('Node "', nodes[[id]]$label, '" has more than one branch with probability #')
      if (length(complement) == 1) probs[complement] <- 1 - sum(probs)
      if (any(probs < -1e-9 | probs > 1 + 1e-9) ||
          (!identical(model$preferences$allowProbabilitiesNotSumTo1, 'true') && abs(sum(probs) - 1) > 1e-6))
        stop('The probabilities of the branches of "', nodes[[id]]$label, '" at cycle ', stage, ' are not valid: ',
             paste(signif(probs, 4), collapse = ', '), call. = FALSE)
      probs
    }
  )
}

# The discount factor of each reward set at a cycle, under the model's global
# discounting. 1 for every set when it is off.
trex.discounter <- function(model, run) {
  prefs <- model$preferences
  sets <- model$reward.sets
  if (!identical(prefs$useGlobalDiscounting, 'true')) return(function(stage) rep(1, length(sets)))
  value <- function(text) {
    code <- trex.compile(text %or% '', model$decimal.comma)
    if (is.null(code)) NA else run$eval(code, 1L, 0)
  }
  cycle.length <- value(prefs$globalDiscountingMarkovCycleLength)
  if (is.na(cycle.length)) cycle.length <- 1
  rates <- vapply(sets, function(set) value(prefs[[paste0('globalDiscountingDiscountRate', set)]]), numeric(1))
  rates[is.na(rates)] <- 0
  function(stage) 1 / (1 + rates)^(stage * cycle.length)
}

# Runs one Markov node: the cohort enters its states with their initial
# probabilities and moves through their transition subtrees until the
# termination condition is met.
#
# Each cycle `_stage`, TreeAge first checks the termination condition; when it is
# met the process ends, after the final rewards of the traditional method. Then
# the cohort accrues the rewards of the states it is in, and goes down the
# transition subtrees, accruing the transition rewards on the way, into the
# states of the next cycle.
#
#   Traditional:  initial rewards at cycle 0, incremental ones at every other
#                 cycle, final rewards on the cohort at the end.
#   WCC:          startup rewards at cycle 0; event rewards on the cohort at the
#                 start of the cycle; cycle rewards on the average of the cohort
#                 at the start and at the end of the cycle (trapezoidal rule).
#
# Returns the rewards accrued (discounted), the trace of the cohort through the
# states, one row per cycle start, and the flows between states of each cycle.
trex.run.markov <- function(model, run, markov.id, discount) {
  nodes <- model$nodes
  structure <- model$states[[as.character(markov.id)]]
  states <- structure$ids
  n <- length(states)
  sets <- as.character(model$reward.sets)
  wcc <- trex.uses.wcc(model)
  termination <- nodes[[markov.id]]$termination.code
  if (is.null(termination)) stop('Markov node "', nodes[[markov.id]]$label, '" has no termination condition')
  max.stages <- as.numeric(model$preferences$maxMarkovStages %or% 10000)

  # The rewards a field (e.g. 'Incr') of the states gives at a cycle, one row per
  # state, one column per reward set.
  state.rewards <- function(field, stage) {
    out <- matrix(0, n, length(sets))
    for (i in seq_len(n)) {
      code <- nodes[[states[i]]]$markov$state.code
      for (j in seq_along(sets)) {
        expr <- code[[sets[j]]][[field]]
        if (!trex.is.zero(expr)) out[i, j] <- run$eval(expr, states[i], stage)
      }
    }
    out
  }
  transition.reward <- function(id, stage) {
    code <- nodes[[id]]$markov$transition.code
    out <- numeric(length(sets))
    if (is.null(code)) return(out)
    for (j in seq_along(sets)) if (!trex.is.zero(code[[sets[j]]])) out[j] <- run$eval(code[[sets[j]]], id, stage)
    out
  }

  cohort <- run$child.probs(markov.id, 0)
  rewards <- numeric(length(sets))
  trace <- list(cohort)
  flows <- list()
  stage <- 0
  repeat {
    if (run$truth(termination, markov.id, stage)) {
      if (!wcc) rewards <- rewards + colSums(cohort * state.rewards('Final', stage)) * discount(stage)
      break
    }
    if (stage >= max.stages) stop('Markov node "', nodes[[markov.id]]$label, '" did not terminate within ', max.stages, ' cycles')

    # ---- State rewards at the start of the cycle ----
    d <- discount(stage)
    if (wcc) {
      if (stage == 0) rewards <- rewards + colSums(cohort * state.rewards('Startup', stage)) * d
      rewards <- rewards + colSums(cohort * state.rewards('Fixed', stage)) * d
      cycle.rewards <- state.rewards('Rate', stage)
      rewards <- rewards + colSums(cohort * cycle.rewards) * d / 2
    } else {
      rewards <- rewards + colSums(cohort * state.rewards(if (stage == 0) 'Init' else 'Incr', stage)) * d
    }

    # ---- Transitions ----
    flow <- matrix(0, n, n, dimnames = list(structure$labels, structure$labels))
    for (i in seq_len(n)) {
      if (structure$absorbing[i]) {
        flow[i, i] <- cohort[i]
        next
      }
      walk <- function(id, p) {
        probs <- run$child.probs(id, stage)
        for (k in seq_along(probs)) {
          child <- nodes[[id]]$children[k]
          q <- p * probs[k]
          rewards <<- rewards + cohort[i] * q * transition.reward(child, stage) * d
          if (length(nodes[[child]]$children) == 0) {
            to <- structure$jumps[[as.character(child)]]
            flow[i, to] <<- flow[i, to] + cohort[i] * q
          } else {
            walk(child, q)
          }
        }
      }
      walk(states[i], 1)
    }
    cohort <- colSums(flow)

    # The end of the cycle's integral of the cycle rewards.
    if (wcc) rewards <- rewards + colSums(cohort * cycle.rewards) * discount(stage + 1) / 2

    trace[[length(trace) + 1]] <- cohort
    flows[[length(flows) + 1]] <- flow
    stage <- stage + 1
  }

  names(rewards) <- sets
  list(
    rewards = rewards,
    cycles = stage,
    trace = do.call(rbind, lapply(trace, function(x) setNames(x, structure$labels))),
    flows = flows
  )
}

# The expected rewards of a subtree, and the Markov processes it goes through
# with the probability of reaching each.
trex.rollback <- function(model, run, id, discount, path = character(0)) {
  nodes <- model$nodes
  node <- nodes[[id]]
  sets <- as.character(model$reward.sets)
  path <- c(path, node$label)

  if (identical(node$type, 'MarkovNode')) {
    result <- trex.run.markov(model, run, id, discount)
    return(list(rewards = result$rewards,
                markov = list(c(list(node = id, path = path, probability = 1), result))))
  }

  if (length(node$children) == 0) {
    rewards <- vapply(sets, function(set) {
      code <- node$payoffs[[set]]
      if (trex.is.zero(code)) 0 else run$eval(code, id, 0)
    }, numeric(1))
    return(list(rewards = rewards, markov = list()))
  }

  branches <- lapply(node$children, function(child) trex.rollback(model, run, child, discount, path))

  if (identical(node$type, 'DecisionNode')) {
    # A decision inside a strategy takes its best branch: in a cost-effectiveness
    # model the one of highest net monetary benefit, otherwise the best value of
    # the main reward set.
    prefs <- model$preferences
    ce <- trex.ce.sets(model)
    score <- vapply(branches, function(b) {
      if (identical(prefs$calcType, 'ct_costEff'))
        b$rewards[[as.character(ce[['effect']])]] * as.numeric(prefs$willingnessToPay %or% 0) - b$rewards[[as.character(ce[['cost']])]]
      else {
        main <- prefs$mainPayoff %or% '1'
        sign <- if (identical(prefs[[paste0('optType', main)]], 'opt_low')) -1 else 1
        sign * b$rewards[[main]]
      }
    }, numeric(1))
    return(branches[[which.max(score)]])
  }

  probs <- run$child.probs(id, 0)
  rewards <- numeric(length(sets))
  markov <- list()
  for (k in seq_along(branches)) {
    rewards <- rewards + probs[k] * branches[[k]]$rewards
    markov <- c(markov, lapply(branches[[k]]$markov, function(m) { m$probability <- m$probability * probs[k]; m }))
  }
  names(rewards) <- sets
  list(rewards = rewards, markov = markov)
}


# ==== Simulation ================================================================

# Runs the strategies of a model with a set of parameters.
#
# `parameters` is a named list of values for parameters of trex.parameters();
# those left out keep the value of the model. A parameter given as a list keyed
# by stratum takes, at each cycle, the value of the stratum of that cycle, which
# `strata` says: a list with the names of the strata and `stratum.of(stage)`,
# the index of the stratum a cycle belongs to.
#
# Returns, for every strategy:
#   summary     data frame of strategy, C (cost) and E (effect)
#   rewards     data frame of strategy and the expected value of every enabled
#               reward set, by the name TreeAge gives it
#   trace       the cohort in each state at the start of every cycle, summed over
#               the Markov processes of the strategy weighted by the probability
#               of reaching them (states with the same name are added together)
#   inflow      per cycle, the share of the cohort entering each state from
#               another one
#   markov      every Markov process the strategy reaches, on its own
trex.simulate <- function(model, strategies = NULL, parameters = list(), strata = NULL) {
  all.strategies <- trex.strategies(model)
  if (is.null(strategies)) strategies <- all.strategies$name
  unknown <- setdiff(strategies, all.strategies$name)
  if (length(unknown) > 0) stop('Unknown strategy: ', paste(unknown, collapse = ', '))
  unknown <- setdiff(names(parameters), names(trex.parameters.cached(model)))
  if (length(unknown) > 0) stop('Unknown parameter: ', paste(unknown, collapse = ', '))

  stratum.of <- if (is.null(strata)) function(stage) NA_character_ else function(stage) strata$names[strata$stratum.of(stage)]
  ce <- as.character(trex.ce.sets(model))
  # Reported are the reward sets TreeAge has enabled, and the cost and effect.
  enabled <- as.integer(model$preferences$numberOfEnabledPayoffs %or% length(model$reward.sets))
  reported <- model$reward.sets <= enabled | as.character(model$reward.sets) %in% ce
  reward.names <- trex.reward.names(model)[reported]

  summary <- data.frame()
  rewards <- data.frame()
  trace <- list()
  inflow <- list()
  markov <- list()
  for (strategy in strategies) {
    node <- all.strategies$node[all.strategies$name == strategy]
    run <- trex.evaluator(model, parameters, stratum.of)
    result <- trex.rollback(model, run, node, trex.discounter(model, run))

    summary <- rbind(summary, data.frame(strategy = strategy, C = result$rewards[[ce[1]]], E = result$rewards[[ce[2]]]))
    rewards <- rbind(rewards, cbind(data.frame(strategy = strategy), as.data.frame(as.list(setNames(result$rewards[reported], reward.names)), check.names = FALSE)))

    combined <- trex.combine.markov(result$markov)
    trace[[strategy]] <- combined$trace
    inflow[[strategy]] <- combined$inflow
    markov[[strategy]] <- result$markov
  }
  list(summary = summary, rewards = rewards, trace = trace, inflow = inflow, markov = markov)
}

# trex.parameters() evaluated once per model, since trex.simulate() checks the
# names it is given against it on every run.
trex.parameters.cached <- local({
  cache <- list()
  function(model) {
    key <- model$path
    if (is.null(cache[[key]])) cache[[key]] <<- trex.parameters(model)
    cache[[key]]
  }
})

# The Markov processes of a strategy added together, each weighted by the
# probability of reaching it. A process that ends earlier than another keeps its
# cohort where it ended.
trex.combine.markov <- function(processes) {
  if (length(processes) == 0) return(list(trace = NULL, inflow = NULL))
  states <- unique(unlist(lapply(processes, function(p) colnames(p$trace))))
  cycles <- max(vapply(processes, function(p) p$cycles, numeric(1)))
  trace <- matrix(0, cycles + 1, length(states), dimnames = list(NULL, states))
  inflow <- matrix(0, cycles, length(states), dimnames = list(NULL, states))
  for (p in processes) {
    rows <- c(seq_len(nrow(p$trace)), rep(nrow(p$trace), cycles + 1 - nrow(p$trace)))
    trace[, colnames(p$trace)] <- trace[, colnames(p$trace)] + p$probability * p$trace[rows, , drop = FALSE]
    for (s in seq_along(p$flows)) {
      f <- p$flows[[s]]
      entering <- colSums(f) - diag(f)
      inflow[s, colnames(f)] <- inflow[s, colnames(f)] + p$probability * entering
    }
  }
  list(trace = cbind(data.frame(stage = 0:cycles), as.data.frame(trace, check.names = FALSE)),
       inflow = cbind(data.frame(stage = seq_len(cycles) - 1), as.data.frame(inflow, check.names = FALSE)))
}


# ==== Description ===============================================================

# The tree as indented text, for whoever wants to read the model as TreeAge
# shows it.
trex.tree.text <- function(model) {
  nodes <- model$nodes
  reward.names <- setNames(trex.reward.names(model), model$reward.sets)
  # Only the state rewards the model's method uses are shown: the file keeps the
  # others, which TreeAge ignores.
  fields <- if (trex.uses.wcc(model)) c(Startup = 'startup', Rate = 'cycle', Fixed = 'event') else c(Init = 'initial', Incr = 'incremental', Final = 'final')
  lines <- character(0)
  add <- function(depth, text) lines <<- c(lines, paste0(strrep('  ', depth), text))
  describe <- function(id, depth) {
    node <- nodes[[id]]
    kind <- sub('Node$', '', node$type %or% 'Node')
    prob <- if (!is.null(node$source$prob)) paste0('  [p = ', node$source$prob, ']') else ''
    clone <- if (!is.null(node$clone.of)) '  (clone)' else ''
    add(depth, paste0(kind, ': ', node$label, prob, clone))
    for (name in names(node$source$definitions)) add(depth + 2, paste0(name, ' = ', node$source$definitions[[name]]))
    if (!is.null(node$termination)) add(depth + 2, paste0('terminates when ', node$termination))
    if (!is.null(node$markov)) {
      if (node$markov$jump != '') add(depth + 2, paste0('-> ', node$markov$jump))
      for (set in sort(names(node$markov$state.rewards))) {
        values <- unlist(node$markov$state.rewards[[set]][names(fields)])
        used <- !is.na(values) & trimws(values) != '0'
        if (any(used)) add(depth + 2, paste0(reward.names[[set]], ': ', paste(fields[used], values[used], sep = ' = ', collapse = ', ')))
      }
      for (set in sort(names(node$markov$transition.rewards)))
        if (trimws(node$markov$transition.rewards[[set]]) != '0')
          add(depth + 2, paste0(reward.names[[set]], ' on transition: ', node$markov$transition.rewards[[set]]))
    }
    for (set in names(node$source$payoffs)) add(depth + 2, paste0(reward.names[[set]], ' payoff: ', node$source$payoffs[[set]]))
    for (child in node$children) describe(child, depth + 1)
  }
  describe(1L, 0)
  paste(lines, collapse = '\n')
}

# The transitions of a Markov node between its states, each with the branches it
# goes through and the product of their probabilities.
trex.markov.transitions <- function(model, markov.id) {
  nodes <- model$nodes
  structure <- model$states[[as.character(markov.id)]]
  out <- data.frame(from = character(0), to = character(0), path = character(0), prob = character(0), stringsAsFactors = FALSE)
  walk <- function(id, from, labels, probs) {
    for (child in nodes[[id]]$children) {
      node <- nodes[[child]]
      l <- c(labels, node$label)
      p <- c(probs, node$source$prob %or% '#')
      if (length(node$children) == 0) {
        out[nrow(out) + 1, ] <<- list(from, structure$labels[structure$jumps[[as.character(child)]]],
                                      paste(l, collapse = ' > '), paste(p, collapse = ' * '))
      } else walk(child, from, l, p)
    }
  }
  for (i in seq_along(structure$ids)) if (!structure$absorbing[i]) walk(structure$ids[i], structure$labels[i], character(0), character(0))
  out
}
