#================================================#
#==========Testing Different Mechanisms==========#
#================================================#
#==This experiment evaluates the Risk Functions==#
#==for different randomization mechanisms, all===#
#==under the worst case potential outcomes=======#
#================================================#


#==================================#
#==========Run experiment==========#
#==================================#
run_experiment2_worst_case <- function(randomization.points_ = randomization.points, SAMPLE.TIMES.expected.risk_ = SAMPLE.TIMES.expected.risk)
{
  HT.estimator.vec = c()
  for (times.temp in 1:SAMPLE.TIMES.expected.risk_)
  {
    assignment.path = generate_assignments(time.horizon, randomization.points_)
    
    potential.outcome.path = generate_outcomes_worst_case(assignment.path_ = assignment.path)
    
    HT.estimator = generate_estimator(randomization.points_, assignment.path, potential.outcome.path, p.lag.length)
    
    HT.estimator.vec = c(HT.estimator.vec, HT.estimator)
  }
  return(HT.estimator.vec)
}


#=====INPUTs=====#
set.seed(111111)
time.horizon = 120
no.m.carryover = 2
p.lag.length = 2
SAMPLE.TIMES.expected.risk = 100000

#=====Other fixed inputs=====#
all.ones = rep(1, time.horizon)
all.zeros = rep(0, time.horizon)

#=====Outcome Model=====#
delta.coef.1 = 0
delta.coef.2 = 0
delta.coef.3 = 0

estimand = generate_estimand_worst_case(p.lag.length)

#===Optimal design===#
randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-2*p.lag.length+1, by = p.lag.length))

HT.estimator.vec = run_experiment2_worst_case(randomization.points_ = randomization.points)
expected.risk = var((HT.estimator.vec - estimand)) * (SAMPLE.TIMES.expected.risk - 1) / SAMPLE.TIMES.expected.risk #var() is the sample variance! Not the population variance...

#===Naive Design 1 ===#
randomization.points.H1 = c(1:time.horizon)

HT.estimator.vec.H1 = run_experiment2_worst_case(randomization.points_ = randomization.points.H1)
expected.risk.H1 = var((HT.estimator.vec.H1 - estimand)) * (SAMPLE.TIMES.expected.risk - 1) / SAMPLE.TIMES.expected.risk

#===Naive Design 2 ===#
randomization.points.H2 = seq(1, time.horizon, by = (p.lag.length+1))

HT.estimator.vec.H2 = run_experiment2_worst_case(randomization.points_ = randomization.points.H2)
expected.risk.H2 = var((HT.estimator.vec.H2 - estimand)) * (SAMPLE.TIMES.expected.risk - 1) / SAMPLE.TIMES.expected.risk

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

my.file.name = as.character(paste(Sys.Date(), "Experiment2.txt"))
write.table(t(output.line), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)





# hist(HT.estimator.vec,
#      main = "Distribution of the estimator",
#      xlab = "Value of Estimator",
#      xlim = c(min(HT.estimator.vec), max(HT.estimator.vec)),
#      freq = FALSE)
# abline(v = mean(HT.estimator.vec), col="red", lwd=2, lty=2)
