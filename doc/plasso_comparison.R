rm(list=ls())
library(plasso)
library(glmnet)
library(parallel)
library(doParallel)
library(microbenchmark)
library(tidyverse)

simulate = function(n=100,p=100,parallel=FALSE,kf=10){
  
  ##### preparation #####
  
  # set seed
  set.seed(42)
  
  x = matrix(runif(n*p,-pi,pi),ncol=p)
  y = runif(n)
  
  
  
  ##### benchmarking #####
  
  m = microbenchmark(
    "plasso" = plasso::plasso(x,y,kf=kf,parallel=parallel),
    times = 3,
    unit = "s"
  )
  m = summary(m)
  
  ##### output #####
  
  output <- data.frame(
    "N" = n,
    "P" = p,
    "Parallel" = parallel,
    "N_Folds" = kf,
    "Median_Time" = m[, "median"],
    "Mean_Time" = m[, "mean"],
    "Iterations" = m[, "neval"]
  )
  
  return(output)
  
}

n_obs = c(1000,10000)
n_variables = c(5,20)
para = c(FALSE,TRUE)
kfolds = c(2,10)

comparison = NULL

for (i in 1:length(n_obs)) {
  for (j in 1:length(n_variables)){
    for (p in 1:length(para)){
      for (k in 1:length(kfolds)){
        print(paste0("--- Simulation for n = ", n_obs[i], ", p = ", n_variables[j], ", parallel = ", para[p], " & kf = " ,kfolds[k]), " ---")
        comparison = rbind(comparison, simulate(n=n_obs[i],p=n_variables[j],parallel=para[p],kf=kfolds[k]))
      }
    }
  }
}

comparison$N_P = paste0(as.character(comparison$N), "_", as.character(comparison$P))

ggplot(comparison, aes(x = K,
                       y = Median_Time,
                       color = Parallel)) +
  geom_line(size = 0.8, alpha = 0.8) +
  facet_wrap(. ~ N_P) +
  theme_gray() +
  scale_colour_brewer(palette = "Paired") +
  labs(x = "Number of folds", y = "Time (in seconds)", color = "Parallelization") +
  scale_x_continuous(breaks = seq(2, max(kfolds), 2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, NA))
