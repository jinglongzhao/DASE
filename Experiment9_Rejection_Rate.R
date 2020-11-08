#==============================================#
#==========Randomization Distribution==========#
#==============================================#


#=====INPUTs=====#
set.seed(111111)

SAMPLE.TIMES.fisher.test = 1000
SAMPLE.TIMES.rejection.rate = 1000 #100000

#===Correct m===#
no.m.carryover = 2
p.lag.length = 2

rejected.exact.matrix = c()
rejected.asymptotic.matrix = c()

for(SPEC.temp in 1:3)
{
  set.seed(111111)
  delta.coef.only = SPEC.temp
  
  rejected.exact.vec  = c()
  rejected.asymptotic.vec = c()
  for(temp.horizon in seq(60,600,60))
  {
    time.horizon = temp.horizon
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
    
    #=====Randomization Distribution=====#
    rej.rate.rejected.exact = 0
    rej.rate.rejected.asymptotic = 0
    for (TIMES.temp in 1:SAMPLE.TIMES.rejection.rate)
    {
      rej.rate.assignment.path = generate_assignments(time.horizon, randomization.points)
      rej.rate.outcome.path = generate_outcomes(assignment.path_ = rej.rate.assignment.path)
      rej.rate.HT.estimator = generate_estimator(randomization.points, rej.rate.assignment.path, rej.rate.outcome.path)
      
      #=====Exact Inference=====#
      exact.inference.PValue = my_fisher_test(randomization.points, rej.rate.outcome.path, rej.rate.HT.estimator, SAMPLE.TIMES.fisher.test)
      
      #=====Asymptotic Inference=====#
      asymptotic.inference.PValue = my_neyman_test(randomization.points, rej.rate.assignment.path, rej.rate.outcome.path, rej.rate.HT.estimator, 2)
      
      if(exact.inference.PValue < 0.1)
      {
        rej.rate.rejected.exact = rej.rate.rejected.exact + 1
      }
      if(asymptotic.inference.PValue < 0.1)
      {
        rej.rate.rejected.asymptotic = rej.rate.rejected.asymptotic + 1
      }
    }
    
    rejected.exact.vec = c(rejected.exact.vec, rej.rate.rejected.exact / SAMPLE.TIMES.rejection.rate)
    rejected.asymptotic.vec = c(rejected.asymptotic.vec, rej.rate.rejected.asymptotic / SAMPLE.TIMES.rejection.rate)
  }
  
  rejected.exact.matrix = rbind(rejected.exact.matrix, rejected.exact.vec)
  rejected.asymptotic.matrix = rbind(rejected.asymptotic.matrix, rejected.asymptotic.vec)
}

for(delta.temp in 1:3)
{
  my.png.name = as.character(paste("Time_dependence,delta=",delta.temp,".png"))
  png(filename = my.png.name, width = 600, height = 600)
  par(mar=c(4.5,5,0,0.5))
  plot(seq(30,300,30), rejected.exact.matrix[delta.temp,], type = 'p', pch = 16, col = "blue", ylim = c(0,1), xlab = "T/m", ylab = "Rejection Rate", cex = 3, cex.axis = 2.5, cex.lab = 2.5)
  lines(seq(30,300,30), rejected.asymptotic.matrix[delta.temp,], type = 'p', pch = 15, col = "red", cex = 2.5)
  legend("bottomright", legend = c("exact inference", "asymptotic inference"), col = c("blue", "red"), pch = c(16, 15), cex = 2.5)
  dev.off()
}

write.table(rejected.exact.matrix, file = paste(Sys.Date(), "experiment9_exact.txt"), append = FALSE, sep = "\t", col.names=FALSE)
write.table(rejected.asymptotic.matrix, file = paste(Sys.Date(), "experiment9_asymptotic.txt"), append = FALSE, sep = "\t", col.names=FALSE)
