#======================================#
#==========Hypothesis Testing==========#
#======================================#

#====================================#
#==========Draw Illustrator==========#
#====================================#
draw.outcome.path <- function(realized.assignment.path, realized.outcome.path)
{
  all.ones.outcome.path = generate_outcomes(assignment.path_ = all.ones)
  all.zeros.outcome.path = generate_outcomes(assignment.path_ = all.zeros)
  
  min.ylab = min(c(all.ones.outcome.path, all.zeros.outcome.path, realized.outcome.path))
  max.ylab = max(c(all.ones.outcome.path, all.zeros.outcome.path, realized.outcome.path))
  
  assignment.ones = which(realized.assignment.path == 1)
  assignment.zeros = which(realized.assignment.path == 0)
  realization.at.ones = realized.outcome.path[assignment.ones]
  realization.at.zeros = realized.outcome.path[assignment.zeros]
  
  my.png.name = as.character(paste(gsub(":", "", as.character(Sys.time())), "Observed_Outcomes.png"))
  png(filename = my.png.name, width = 1000, height = 600)

  plot(1:time.horizon, realized.outcome.path, ylim = c(min.ylab, max.ylab), type = 'l', lwd = 2, lty = 1, xlab = "Time horizon", ylab = "Observed Outcomes", main = "Observed Outcomes Following One Realization of Assignments")
  lines(1:time.horizon, all.ones.outcome.path, col = rgb(0.4, 0.7, 1), lwd = 2, lty = 2)
  lines(1:time.horizon, all.zeros.outcome.path, col = rgb(1, 0.4, 0.7), lwd = 2, lty = 2)
  
  lines(assignment.ones, realization.at.ones, col = rgb(0, 0.5, 1), lwd = 2, type = 'p')
  lines(assignment.zeros, realization.at.zeros, col = rgb(1, 0, 0.5), lwd = 2, type = 'p', pch = 2)
  
  legend("bottomright", legend = c("treatment", "control", "potential outcomes under consecutive treatments", "potential outcomes under consecutive controls"),
         col = c(rgb(0, 0.5, 1), rgb(1, 0, 0.5), rgb(0.4, 0.7, 1), rgb(1, 0.4, 0.7)), pch = c(1, 2, NA, NA), lty = c(NA, NA, 2, 2), cex = 1)
  dev.off()
}

#================================================#
#==========m misspecified causal effect==========#
#================================================#
generate_misspecified_estimand <- function(assignment.path_ = assignment.path)
{
  causal.effect.path = rep(0, p.lag.length)
  
  for (digit.temp in (p.lag.length+1):time.horizon)
  {
    counterfactual.ones.assignment.path = assignment.path_
    counterfactual.zeros.assignment.path = assignment.path_
    
    counterfactual.ones.assignment.path[(digit.temp-p.lag.length) : digit.temp] = rep(1, (p.lag.length+1))
    counterfactual.zeros.assignment.path[(digit.temp-p.lag.length) : digit.temp] = rep(0, (p.lag.length+1))
    
    counterfactual.ones.outcome.path = generate_outcomes(counterfactual.ones.assignment.path)
    counterfactual.zeros.outcome.path = generate_outcomes(counterfactual.zeros.assignment.path)
    
    this.digit.causal.effect = counterfactual.ones.outcome.path[digit.temp] - counterfactual.zeros.outcome.path[digit.temp]
    
    causal.effect.path = c(causal.effect.path, this.digit.causal.effect)
  }
  return(sum(causal.effect.path)/(time.horizon - p.lag.length))
}

#======================================#
#==========Variance Estimator==========#
#======================================#
generate_variance_UB_OPTDesign <- function(randomization.points_ = randomization.points,
                                           assignment.path_,
                                           potential.outcome.path_,
                                           p.lag.length_ = p.lag.length,
                                           which.upper.bound = 2)
{
  K.many.randomizations = length(randomization.points_) 
  #==Note: 1. n in the paper corresponds to (K.many.randomizations+2)
  #        2. the paper starts with k=0; the codes starts with K.many.randomizations=1
  observed.chunks = c()
  for(this.randomization in 1:(K.many.randomizations+1))
  {
    temp.sum = sum(potential.outcome.path_[(this.randomization*p.lag.length_ + 1):((this.randomization+1)*p.lag.length_)])
    observed.chunks = c(observed.chunks, temp.sum)
  }
  if(which.upper.bound == 1)
  {
    variance.estimator = 6 * (observed.chunks[1])^2 + 6 * (observed.chunks[K.many.randomizations+1])^2
    for(this.randomization in 2:(K.many.randomizations))
    {
      if(assignment.path_[(this.randomization-1)*p.lag.length_+1] == assignment.path_[this.randomization*p.lag.length_+1])
      {
        variance.estimator = variance.estimator + 24 * (observed.chunks[this.randomization])^2
      }
    }
    for(this.randomization in 1:(K.many.randomizations))
    {
      if(assignment.path_[(this.randomization-1)*p.lag.length_+1] == assignment.path_[this.randomization*p.lag.length_+1] &&
         assignment.path_[(this.randomization+1)*p.lag.length_+1] == assignment.path_[this.randomization*p.lag.length_+1])
      {
        variance.estimator = variance.estimator + 16 * observed.chunks[this.randomization] * observed.chunks[this.randomization+1]
      }
    }
  }
  if(which.upper.bound == 2)
  {
    variance.estimator = 8 * (observed.chunks[1])^2 + 8 * (observed.chunks[K.many.randomizations+1])^2
    for(this.randomization in 2:(K.many.randomizations))
    {
      if(assignment.path_[(this.randomization-1)*p.lag.length_+1] == assignment.path_[this.randomization*p.lag.length_+1])
      {
        variance.estimator = variance.estimator + 32 * (observed.chunks[this.randomization])^2
      }
    }
  }
  estimator.return = variance.estimator / ((time.horizon - p.lag.length_)^2)
  return(estimator.return)
}

#===============================================#
#==========Fisher's exact p-value test==========#
#===============================================#
my_fisher_test <- function(randomization.points_ = randomization.points,
                           realized.outcome.path_ = realized.outcome.path,
                           realized.HT.estimator_ = realized.HT.estimator,
                           SAMPLE.TIMES.fisher.test_ = SAMPLE.TIMES.fisher.test)
{
  indicator.more.extreme = c()
  for (TIMES.temp in 1:SAMPLE.TIMES.fisher.test_)
  {
    simulated.assignment.path = generate_assignments(time.horizon, randomization.points_)
    simulated.HT.estimator = generate_estimator(randomization.points_, simulated.assignment.path, realized.outcome.path_)
    
    if(abs(simulated.HT.estimator) > abs(realized.HT.estimator_))
    {
      indicator.more.extreme = c(indicator.more.extreme, 1)
    }
    else
    {
      indicator.more.extreme = c(indicator.more.extreme, 0)
    }
  }
  return(sum(indicator.more.extreme)/SAMPLE.TIMES.fisher.test_)
}

#============================================#
#==========Neyman's asymptotic test==========#
#============================================#
my_neyman_test <- function(randomization.points_ = randomization.points,
                           assignment.path_,
                           potential.outcome.path_,
                           realized.HT.estimator_ = realized.HT.estimator,
                           which.upper.bound = 2)
{
  estimated.variance = generate_variance_UB_OPTDesign(randomization.points_, assignment.path_, potential.outcome.path_, p.lag.length, which.upper.bound)
  place.of.zero = abs(realized.HT.estimator_) / sqrt(estimated.variance)
  
  return(2 - 2 * pnorm(place.of.zero))
}


#==========INPUTs==========#
set.seed(111111)

time.horizon = 120
SAMPLE.TIMES.fisher.test = 100000 #100000

#=====Other fixed inputs=====#
all.ones = rep(1, time.horizon)
all.zeros = rep(0, time.horizon)

experiment6.output.matrix = c()

#===Correct m===#
no.m.carryover = 2
p.lag.length = 2

two.tests.output.matrix = c()
for (SPEC.temp in 1:3)
{
  set.seed(111111)
  delta.coef.only = SPEC.temp
  
  #=====Outcome Model=====#
  mu.fixed.effect = 0
  alpha.fixed.effects = log(1:time.horizon) #runif(time.horizon, min = 0, max = 1) #rnorm(time.horizon, mean = 0, sd = 1)
  epsilon.noises = rnorm(time.horizon, mean = 0, sd = 1)
  
  delta.coef.1 = delta.coef.only
  delta.coef.2 = delta.coef.only
  delta.coef.3 = delta.coef.only
  
  estimand = generate_estimand()
  
  #===Optimal design===#
  randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-2*p.lag.length+1, by = p.lag.length))
  #---One Realization---#
  realized.assignment.path = generate_assignments(time.horizon, randomization.points)
  realized.outcome.path = generate_outcomes(assignment.path_ = realized.assignment.path)
  realized.HT.estimator = generate_estimator(randomization.points, realized.assignment.path, realized.outcome.path)
  draw.outcome.path(realized.assignment.path, realized.outcome.path)
  
  #=====Exact Inference=====#
  exact.inference.PValue = my_fisher_test(randomization.points, realized.outcome.path, realized.HT.estimator, SAMPLE.TIMES.fisher.test)

  #=====Asymptotic Inference=====#
  realized.estimated.variance = generate_variance_UB_OPTDesign(randomization.points, realized.assignment.path, realized.outcome.path, p.lag.length, 2)
  asymptotic.inference.PValue = my_neyman_test(randomization.points, realized.assignment.path, realized.outcome.path, realized.HT.estimator, 2)

  hypothesis.testing.output = rbind(realized.assignment.path, realized.outcome.path)
  my.file.name = as.character(paste(gsub(":", "", as.character(Sys.time())), "Hypothesis_Testing_correct.txt"))
  write.table(exact.inference.output, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  two.tests.output.matrix = rbind(two.tests.output.matrix, c(estimand, NA, realized.HT.estimator, realized.estimated.variance, exact.inference.PValue, asymptotic.inference.PValue))
}
experiment6.output.matrix = rbind(experiment6.output.matrix, two.tests.output.matrix)


#===Misspecified m: Case 1: when we overestimate m===#
no.m.carryover = 2
p.lag.length = 3

two.tests.output.matrix = c()
for (SPEC.temp in 1:3)
{
  set.seed(111111)
  delta.coef.only = SPEC.temp
  
  #=====Outcome Model=====#
  mu.fixed.effect = 0
  alpha.fixed.effects = log(1:time.horizon) #runif(time.horizon, min = 0, max = 1) #rnorm(time.horizon, mean = 0, sd = 1)
  epsilon.noises = rnorm(time.horizon, mean = 0, sd = 1)
  
  delta.coef.1 = delta.coef.only
  delta.coef.2 = delta.coef.only
  delta.coef.3 = delta.coef.only
  
  estimand = generate_estimand()
  
  #===Optimal design===#
  randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-2*p.lag.length+1, by = p.lag.length))
  #---One Realization---#
  realized.assignment.path = generate_assignments(time.horizon, randomization.points)
  realized.outcome.path = generate_outcomes(assignment.path_ = realized.assignment.path)
  realized.HT.estimator = generate_estimator(randomization.points, realized.assignment.path, realized.outcome.path)
  draw.outcome.path(realized.assignment.path, realized.outcome.path)
  
  #=====Exact Inference=====#
  exact.inference.PValue = my_fisher_test(randomization.points, realized.outcome.path, realized.HT.estimator, SAMPLE.TIMES.fisher.test)

  #=====Asymptotic Inference=====#
  realized.estimated.variance = generate_variance_UB_OPTDesign(randomization.points, realized.assignment.path, realized.outcome.path, p.lag.length, 2)
  asymptotic.inference.PValue = my_neyman_test(randomization.points, realized.assignment.path, realized.outcome.path, realized.HT.estimator, 2)
  
  hypothesis.testing.output = rbind(realized.assignment.path, realized.outcome.path)
  my.file.name = as.character(paste(gsub(":", "", as.character(Sys.time())), "Hypothesis_Testing_over.txt"))
  write.table(exact.inference.output, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  two.tests.output.matrix = rbind(two.tests.output.matrix, c(estimand, NA, realized.HT.estimator, realized.estimated.variance, exact.inference.PValue, asymptotic.inference.PValue))
}
experiment6.output.matrix = rbind(experiment6.output.matrix, two.tests.output.matrix)


#===Misspecified m: Case 2: when we underestimate m===#
no.m.carryover = 2
p.lag.length = 1

two.tests.output.matrix = c()
for (SPEC.temp in 1:3)
{
  set.seed(111111)
  delta.coef.only = SPEC.temp
  
  #=====Outcome Model=====#
  mu.fixed.effect = 0
  alpha.fixed.effects = log(1:time.horizon) #runif(time.horizon, min = 0, max = 1) #rnorm(time.horizon, mean = 0, sd = 1)
  epsilon.noises = rnorm(time.horizon, mean = 0, sd = 1)
  
  delta.coef.1 = delta.coef.only
  delta.coef.2 = delta.coef.only
  delta.coef.3 = delta.coef.only
  
  #===Optimal design===#
  randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-2*p.lag.length+1, by = p.lag.length))
  #---One Realization---#
  realized.assignment.path = generate_assignments(time.horizon, randomization.points)
  misspecified.estimand = generate_misspecified_estimand(realized.assignment.path)
  realized.outcome.path = generate_outcomes(assignment.path_ = realized.assignment.path)
  realized.HT.estimator = generate_estimator(randomization.points, realized.assignment.path, realized.outcome.path)
  draw.outcome.path(realized.assignment.path, realized.outcome.path)
  
  #=====Exact Inference=====#
  exact.inference.PValue = my_fisher_test(randomization.points, realized.outcome.path, realized.HT.estimator, SAMPLE.TIMES.fisher.test)

  #=====Asymptotic Inference=====#
  realized.estimated.variance = generate_variance_UB_OPTDesign(randomization.points, realized.assignment.path, realized.outcome.path, p.lag.length, 2)
  asymptotic.inference.PValue = my_neyman_test(randomization.points, realized.assignment.path, realized.outcome.path, realized.HT.estimator, 2)
  
  hypothesis.testing.output = rbind(realized.assignment.path, realized.outcome.path)
  my.file.name = as.character(paste(gsub(":", "", as.character(Sys.time())), "Hypothesis_Testing_under.txt"))
  write.table(exact.inference.output, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  two.tests.output.matrix = rbind(two.tests.output.matrix, c(NA, misspecified.estimand, realized.HT.estimator, realized.estimated.variance, exact.inference.PValue, asymptotic.inference.PValue))
}
experiment6.output.matrix = rbind(experiment6.output.matrix, two.tests.output.matrix)

write.table(experiment6.output.matrix, file = paste(Sys.Date(), "experiment6.txt"), append = FALSE, sep = "\t", col.names=FALSE)

experiment6.output.matrix.copy = experiment6.output.matrix
experiment6.output.matrix.copy[,c(3,4)] = round(experiment6.output.matrix.copy[,c(3,4)], digits = 2)
experiment6.output.matrix.copy[,c(5,6)] = round(experiment6.output.matrix.copy[,c(5,6)], digits = 3)
write.table(experiment6.output.matrix.copy, file = paste(Sys.Date(), "experiment6_rounded.txt"), append = FALSE, sep = "\t", col.names=FALSE)
