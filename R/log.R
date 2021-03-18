#!/usr/bin/env Rscript

#
#  This file is part of the `intercell` R package
#
#  Copyright
#  2021
#  Saez Lab, Uniklinik RWTH Aachen, Heidelberg University
#
#  Author(s): Daniel Dimitrov
#             Charlotte Boys
#             Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the MIT (Expat) License.
#  See accompanying file `LICENSE` or find a copy at
#      https://directory.fsf.org/wiki/License:Expat
#
#  Git repo: https://github.com/saezlab/Cell_Cell_Investigation
#


#' Dispatch log message to the OmnipathR logger
#'
#' @importFrom logger log_level
#' @importFrom magrittr %<>%
logg <- function(level, ...){

    level %<>% (OmnipathR:::ensure_loglevel)
    log_level(level, ..., namespace = 'OmnipathR')

}


log_trace <- function(...){logg(level = 'trace', ...)}
log_info <- function(...){logg(level = 'info', ...)}
log_success <- function(...){logg(level = 'success', ...)}
log_warn <- function(...){logg(level = 'warn', ...)}
log_error <- function(...){logg(level = 'error', ...)}
log_fatal <- function(...){logg(level = 'fatal', ...)}