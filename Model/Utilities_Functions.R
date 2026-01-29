##*****************************************************************************************************
## UTILITIES FUNCTIONS
##***************************************************************************************************** 
  
##*********************************************************************
## Function to automatically install (if missing) and load libraries
##*********************************************************************
load_libraries <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
    } else {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
    }
  }
}

##*********************************************************************
## safe rbind
##*********************************************************************
safe_rbind <- function(dfs) {
  out <- do.call(rbind, dfs)
  out
}

##*********************************************************************
## journal fill color
##*********************************************************************
get_journal_fill_scale <- function(journal) {
  journal <- tolower(journal)

  if (journal == "lancet") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for Lancet palette.")
    return(ggsci::scale_fill_lancet())
  }

  if (journal == "nejm") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for NEJM palette.")
    return(ggsci::scale_fill_nejm())
  }

  if (journal == "jama") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for JAMA palette.")
    return(ggsci::scale_fill_jama())
  }

  if (journal == "bmj") {
    # BMJ doesn’t have a dedicated palette; use a clean blue palette
    return(scale_fill_brewer(palette = "Blues"))
  }

  if (journal == "nature") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for Nature palette.")
    return(ggsci::scale_fill_nature())
  }

  # default safe palette
  return(scale_fill_brewer(palette = "Set2"))
}

##*********************************************************************
## journal colors
##*********************************************************************
get_journal_color_scale <- function(journal) {
  journal <- tolower(journal)

  if (journal == "lancet") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for Lancet palette.")
    return(ggsci::scale_color_lancet())
  }

  if (journal == "nejm") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for NEJM palette.")
    return(ggsci::scale_color_nejm())
  }

  if (journal == "jama") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for JAMA palette.")
    return(ggsci::scale_color_jama())
  }

  if (journal == "bmj") {
    # BMJ doesn’t have a dedicated palette; use a clean blue palette
    return(scale_color_brewer(palette = "Blues"))
  }

  if (journal == "nature") {
    if (!requireNamespace("ggsci", quietly = TRUE)) stop("Install 'ggsci' for Nature palette.")
    return(ggsci::scale_color_nature())
  }

  # default safe palette
  return(scale_color_brewer(palette = "Set2"))
}

##*********************************************************************
## number formatter
##*********************************************************************
fmt_ci_journal <- function(
  mean, 
  lower = NULL, 
  upper = NULL,
  journal = c(
    "value_in_health",
    "health_economics",
    "health_affairs",
    "jama",
    "nejm",
    "bmj",
    "lancet",
    "science",
    "pnas",
    "plos",
    "ejhe",
    "custom"
  ),
  digits = 0,
  scale = 1,
  linebreak = "\n",
  big.mark = NULL,
  decimal.mark = NULL,
  ci_sep = NULL
) {

  journal <- match.arg(journal)

  presets <- list(
    value_in_health = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    health_economics = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    health_affairs = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    jama = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    nejm = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    bmj = list(big.mark = ",", decimal.mark = ".", ci_sep = " - "),
    lancet = list(big.mark = ",", decimal.mark = ".", ci_sep = " - "),
    science = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    pnas = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    plos = list(big.mark = ",", decimal.mark = ".", ci_sep = ", "),
    ejhe = list(big.mark = ".", decimal.mark = ",", ci_sep = "–")
  )

  if (journal != "custom") {
    big.mark <- presets[[journal]]$big.mark
    decimal.mark <- presets[[journal]]$decimal.mark
    if (is.null(ci_sep)) ci_sep <- presets[[journal]]$ci_sep
  } else {
    if (is.null(big.mark) || is.null(decimal.mark)) {
      stop("For journal = 'custom', big.mark and decimal.mark must be specified.")
    }
    if (is.null(ci_sep)) ci_sep <- ", "
  }

  #if(is.null(lower) | is.null(upper)){ci_sep <- ""}

  f <- function(x) {
    formatC(
      x / scale,
      format = "f",
      digits = digits,
      big.mark = big.mark,
      decimal.mark = decimal.mark
    )
  }

  estimate <- f(mean)

  if (!is.null(lower) && !is.null(upper)) {
    paste0(
      estimate, linebreak,
      "(", f(lower), ci_sep, f(upper), ")"
    )
  } else {
    estimate
  }
}

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 