#===================================================#
#==========Identify Carryover Effect Order==========#
#===================================================#


#==========INPUTs==========#
set.seed(111111)

time.horizon = 1020

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
  realized.assignment.path = generate_assignments(time.horizon, randomization.points, re.randomization = FALSE)
  realized.outcome.path = generate_outcomes(assignment.path_ = realized.assignment.path)
  realized.HT.estimator = generate_estimator(randomization.points, realized.assignment.path, realized.outcome.path)
  
  #=====Exact Inference=====#
  exact.inference.PValue = -1
  
  #=====Asymptotic Inference=====#
  realized.estimated.variance = generate_variance_UB_OPTDesign(randomization.points, realized.assignment.path, realized.outcome.path, p.lag.length, 2)
  asymptotic.inference.PValue = my_neyman_test(randomization.points, realized.assignment.path, realized.outcome.path, realized.HT.estimator, 2)
  
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
  realized.assignment.path = generate_assignments(time.horizon, randomization.points, re.randomization = FALSE)
  realized.outcome.path = generate_outcomes(assignment.path_ = realized.assignment.path)
  realized.HT.estimator = generate_estimator(randomization.points, realized.assignment.path, realized.outcome.path)
  
  #=====Exact Inference=====#
  exact.inference.PValue = -1
  
  #=====Asymptotic Inference=====#
  realized.estimated.variance = generate_variance_UB_OPTDesign(randomization.points, realized.assignment.path, realized.outcome.path, p.lag.length, 2)
  asymptotic.inference.PValue = my_neyman_test(randomization.points, realized.assignment.path, realized.outcome.path, realized.HT.estimator, 2)
  
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
  realized.assignment.path = generate_assignments(time.horizon, randomization.points, re.randomization = FALSE)
  misspecified.estimand = generate_misspecified_estimand(realized.assignment.path)
  realized.outcome.path = generate_outcomes(assignment.path_ = realized.assignment.path)
  realized.HT.estimator = generate_estimator(randomization.points, realized.assignment.path, realized.outcome.path)
  
  #=====Exact Inference=====#
  exact.inference.PValue = -1
  
  #=====Asymptotic Inference=====#
  realized.estimated.variance = generate_variance_UB_OPTDesign(randomization.points, realized.assignment.path, realized.outcome.path, p.lag.length, 2)
  asymptotic.inference.PValue = my_neyman_test(randomization.points, realized.assignment.path, realized.outcome.path, realized.HT.estimator, 2)
  
  two.tests.output.matrix = rbind(two.tests.output.matrix, c(NA, misspecified.estimand, realized.HT.estimator, realized.estimated.variance, exact.inference.PValue, asymptotic.inference.PValue))
}
experiment6.output.matrix = rbind(experiment6.output.matrix, two.tests.output.matrix)


#=====INPUTs=====#
three.experiments = experiment6.output.matrix

#=====Test 1: p_1=2 vs. p_2=3 =====#
difference.in.means = experiment6.output.matrix[1,3]-experiment6.output.matrix[4,3]
sum.variances = experiment6.output.matrix[1,4]+experiment6.output.matrix[4,4]
2*(1-pnorm(abs(difference.in.means / sqrt(sum.variances))))

#=====Test 2: p_1=1 vs. p_2=2 =====#
difference.in.means = experiment6.output.matrix[1,3]-experiment6.output.matrix[7,3]
sum.variances = experiment6.output.matrix[1,4]+experiment6.output.matrix[7,4]
2*(1-pnorm(abs(difference.in.means / sqrt(sum.variances))))



#=====Test 1: p_1=2 vs. p_2=3 =====#
difference.in.means = experiment6.output.matrix[2,3]-experiment6.output.matrix[5,3]
sum.variances = experiment6.output.matrix[2,4]+experiment6.output.matrix[5,4]
2*(1-pnorm(abs(difference.in.means / sqrt(sum.variances))))

#=====Test 2: p_1=1 vs. p_2=2 =====#
difference.in.means = experiment6.output.matrix[2,3]-experiment6.output.matrix[8,3]
sum.variances = experiment6.output.matrix[2,4]+experiment6.output.matrix[8,4]
2*(1-pnorm(abs(difference.in.means / sqrt(sum.variances))))



#=====Test 1: p_1=2 vs. p_2=3 =====#
difference.in.means = experiment6.output.matrix[3,3]-experiment6.output.matrix[6,3]
sum.variances = experiment6.output.matrix[3,4]+experiment6.output.matrix[6,4]
2*(1-pnorm(abs(difference.in.means / sqrt(sum.variances))))

#=====Test 2: p_1=1 vs. p_2=2 =====#
difference.in.means = experiment6.output.matrix[3,3]-experiment6.output.matrix[9,3]
sum.variances = experiment6.output.matrix[3,4]+experiment6.output.matrix[9,4]
2*(1-pnorm(abs(difference.in.means / sqrt(sum.variances))))
