#================================================#
#==========Testing Different Mechanisms==========#
#==This experiment evaluates the Risk Functions==#
#==for different randomization mechanisms, all===#
#==under the second linear model in the paper====#
#================================================#



#==================================#
#==========Run experiment==========#
#==================================#
run_experiment1_different_mechanisms <- function(randomization.points_ = randomization.points, SAMPLE.TIMES.expected.risk_ = SAMPLE.TIMES.expected.risk)
{
  HT.estimator.vec = c()
  for (times.temp in 1:SAMPLE.TIMES.expected.risk_)
  {
    assignment.path = generate_assignments(time.horizon, randomization.points_)
    
    potential.outcome.path = generate_outcomes(assignment.path_ = assignment.path)
    
    HT.estimator = generate_estimator(randomization.points_, assignment.path, potential.outcome.path, p.lag.length)
    
    HT.estimator.vec = c(HT.estimator.vec, HT.estimator)
  }
  return(HT.estimator.vec)
}

run_experiment1_deterministic_assignment <- function(assignment.path_ = assignment.path) #Using the Difference-in-means estimator
{
  potential.outcome.path = generate_outcomes(assignment.path_)
  
  indicator.ones = c()
  indicator.zeros = c()
  for(t in (p.lag.length+1):time.horizon)
  {
    if(sum(assignment.path_[(t-p.lag.length) : t] == 1) == (p.lag.length+1))
    {
      indicator.ones = c(indicator.ones, t)
    }
    if(sum(assignment.path_[(t-p.lag.length) : t] == 1) == 0)
    {
      indicator.zeros = c(indicator.zeros, t)
    }
  }
  
  HT.estimator = mean(potential.outcome.path[indicator.ones]) - mean(potential.outcome.path[indicator.zeros])
  
  return(HT.estimator)
}

#=====INPUTs=====#
set.seed(111111)
time.horizon = 120
no.m.carryover = 2
p.lag.length = 2
SAMPLE.TIMES.expected.risk = 100000 #100000

#=====Other fixed inputs=====#
all.ones = rep(1, time.horizon)
all.zeros = rep(0, time.horizon)

tttt = Sys.time()

experiment1.output.matrix = c()
for (SPEC.temp in 1:8)
{
  #=====Outcome Model=====#
  set.seed(111111)
  mu.fixed.effect = 0
  alpha.fixed.effects = log(1:time.horizon) #runif(time.horizon, min = 0, max = 1) #rnorm(time.horizon, mean = 0, sd = 1)
  epsilon.noises = rnorm(time.horizon, mean = 0, sd = 1)
  
  SPEC.temp = SPEC.temp - 1
  delta.coef.1 = (SPEC.temp >= 4) * 2 + (SPEC.temp < 4) * 1
  delta.coef.2 = ((SPEC.temp %% 4) >=2) * 2 + ((SPEC.temp %% 4) < 2) * 1
  delta.coef.3 = ((SPEC.temp %% 2) >=1) * 2 + ((SPEC.temp %% 2) < 1) * 1
  
  estimand = generate_estimand(p.lag.length)
  
  #===Optimal design===#
  randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-2*p.lag.length+1, by = p.lag.length))
  HT.estimator.vec = run_experiment1_different_mechanisms(randomization.points_ = randomization.points, SAMPLE.TIMES.expected.risk)
  expected.risk = var((HT.estimator.vec - estimand)) * (SAMPLE.TIMES.expected.risk - 1) / SAMPLE.TIMES.expected.risk #var() is the sample variance! Not the population variance...
  
  #===Naive Design 1 ===#
  randomization.points.H1 = c(1:time.horizon)
  HT.estimator.vec.H1 = run_experiment1_different_mechanisms(randomization.points_ = randomization.points.H1, SAMPLE.TIMES.expected.risk)
  expected.risk.H1 = var((HT.estimator.vec.H1 - estimand)) * (SAMPLE.TIMES.expected.risk - 1) / SAMPLE.TIMES.expected.risk
  
  #===Naive Design 2 ===#
  randomization.points.H2 = seq(1, time.horizon, by = (p.lag.length+1))
  HT.estimator.vec.H2 = run_experiment1_different_mechanisms(randomization.points_ = randomization.points.H2, SAMPLE.TIMES.expected.risk)
  expected.risk.H2 = var((HT.estimator.vec.H2 - estimand)) * (SAMPLE.TIMES.expected.risk - 1) / SAMPLE.TIMES.expected.risk
  
  # #===Deterministic Design 1 ===#
  # deterministic.assignment.path = c(rep(1,30), rep(0,30))
  # HT.estimator.D1 = run_experiment1_deterministic_assignment(assignment.path_ = deterministic.assignment.path)
  # risk.D1 = (HT.estimator.D1 - estimand)^2
  # 
  # #===Deterministic Design 2 ===#
  # deterministic.assignment.path = rep(c(1,1,1,0,0,0),10)
  # HT.estimator.D2 = run_experiment1_deterministic_assignment(assignment.path_ = deterministic.assignment.path)
  # risk.D2 = (HT.estimator.D1 - estimand)^2
  
  
  output.line = c(deltacoef1 = delta.coef.1,
                  deltacoef2 = delta.coef.2,
                  deltacoef3 = delta.coef.3,
                  causal.estimand = estimand, 
                  avg.estimator = summary(HT.estimator.vec)[4],
                  avg.estimator.H1 = summary(HT.estimator.vec.H1)[4],
                  avg.estimator.H2 = summary(HT.estimator.vec.H2)[4],
                  avg.risk = expected.risk,
                  avg.risk.H1 = expected.risk.H1,
                  avg.risk.H2 = expected.risk.H2)
  experiment1.output.matrix = rbind(experiment1.output.matrix, t(output.line))
}

Sys.time() - tttt

my.file.name = as.character(paste(Sys.Date(), "Experiment1.txt"))
write.table(experiment1.output.matrix, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

experiment1.output.matrix.copy = experiment1.output.matrix
experiment1.output.matrix.copy[,c(8,9,10)] = round(experiment1.output.matrix.copy[,c(8,9,10)], digits = 2)
experiment1.output.matrix.copy[,c(5,6,7)] = round(experiment1.output.matrix.copy[,c(5,6,7)], digits = 3)
my.file.name.copy = as.character(paste(Sys.Date(), "Experiment1_rounded.txt"))
write.table(experiment1.output.matrix.copy, file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)
