#!/usr/bin/env Rscript


CLEAN_SAMPLES_DECONTAM_conc <- function(samples,
                                        concentration,
                                        batch,
                                        method='frequency'
){
  ps <- samples
  prac_conc = concentration
  out <- isContaminant(ps,
                       conc=prac_conc,
                       batch=batch,
                       detailed = T,
                       method=method
  )
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == T] <- 0
      return(x)
    }  ) %>% t()
  }
  return(decontaminated_samples = remove_out_samps(samples, out$contaminant) )
}

## LOAD LIBRARIES ###
suppressMessages({library(tidyverse, quietly=TRUE)})
suppressMessages({library(stringr, quietly=TRUE)})
suppressMessages({library(dplyr, quietly=TRUE)})
suppressMessages({library(decontam, quietly=TRUE)})
suppressMessages({library(data.table, quietly=TRUE)})


set.seed(42)


data <- read.csv('tmp-dataset/data.csv', row.names=1)
md <- read.csv('tmp-dataset/md.csv', row.names=1)


CLEAN_SAMPLES_DECONTAM_conc(samples=as.matrix(data),
                            concentration = md$mu,
                            batch = md$bt,
                            method='frequency'
) %>% write.csv('tmp-dataset/decontam_raw_out.csv')

### QUIT ###
q(status=0)








