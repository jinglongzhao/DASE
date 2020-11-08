#==============================================#
#==========Randomization Distribution==========#
#==============================================#


#=====INPUTs=====#
set.seed(111111)

time.horizon = 120
SAMPLE.TIMES.find.distribution = 100000 #100000

#===Correct m===#
no.m.carryover = 2
p.lag.length = 2

dist.HT.estimator.matrix.mCorrect = c()
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
  
  #=====Randomization Distribution=====#
  dist.HT.estimator.vec = c()
  dist.variance.UB1.vec = c()
  dist.variance.UB2.vec = c()
  for (TIMES.temp in 1:SAMPLE.TIMES.find.distribution)
  {
    dist.assignment.path = generate_assignments(time.horizon, randomization.points)
    dist.outcome.path = generate_outcomes(assignment.path_ = dist.assignment.path)
    dist.HT.estimator = generate_estimator(randomization.points, dist.assignment.path, dist.outcome.path)
    dist.variance.UB1 = generate_variance_UB_OPTDesign(randomization.points, dist.assignment.path, dist.outcome.path, p.lag.length, 1)
    dist.variance.UB2 = generate_variance_UB_OPTDesign(randomization.points, dist.assignment.path, dist.outcome.path, p.lag.length, 2)
    
    dist.HT.estimator.vec = c(dist.HT.estimator.vec, dist.HT.estimator)
    dist.variance.UB1.vec = c(dist.variance.UB1.vec, dist.variance.UB1)
    dist.variance.UB2.vec = c(dist.variance.UB2.vec, dist.variance.UB2)
  }
  dist.HT.estimator.matrix.mCorrect = rbind(dist.HT.estimator.matrix.mCorrect, 
                                            c(estimand, NA,
                                              mean(dist.HT.estimator.vec), var(dist.HT.estimator.vec)*(SAMPLE.TIMES.find.distribution-1)/SAMPLE.TIMES.find.distribution,
                                              mean(dist.variance.UB1.vec), mean(dist.variance.UB2.vec)))
  
  if(delta.coef.only == 1)
  {
    my.png.name = as.character(paste("Correct_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
  if(delta.coef.only == 2)
  {
    my.png.name = as.character(paste("Correct_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
  if(delta.coef.only == 3)
  {
    my.png.name = as.character(paste("Correct_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
}



#===Misspecified m: Case 1: when we overestimate m===#
no.m.carryover = 2
p.lag.length = 3

dist.HT.estimator.matrix.mOver = c()
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
  randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-p.lag.length, by = p.lag.length))
  
  #=====Randomization Distribution=====#
  dist.HT.estimator.vec = c()
  dist.variance.UB1.vec = c()
  dist.variance.UB2.vec = c()
  for (TIMES.temp in 1:SAMPLE.TIMES.find.distribution)
  {
    dist.assignment.path = generate_assignments(time.horizon, randomization.points)
    dist.outcome.path = generate_outcomes(assignment.path_ = dist.assignment.path)
    dist.HT.estimator = generate_estimator(randomization.points, dist.assignment.path, dist.outcome.path)
    dist.variance.UB1 = generate_variance_UB_OPTDesign(randomization.points, dist.assignment.path, dist.outcome.path, p.lag.length, 1)
    dist.variance.UB2 = generate_variance_UB_OPTDesign(randomization.points, dist.assignment.path, dist.outcome.path, p.lag.length, 2)
    
    dist.HT.estimator.vec = c(dist.HT.estimator.vec, dist.HT.estimator)
    dist.variance.UB1.vec = c(dist.variance.UB1.vec, dist.variance.UB1)
    dist.variance.UB2.vec = c(dist.variance.UB2.vec, dist.variance.UB2)
  }
  dist.HT.estimator.matrix.mOver = rbind(dist.HT.estimator.matrix.mOver, 
                                         c(estimand, NA,
                                           mean(dist.HT.estimator.vec), var(dist.HT.estimator.vec)*(SAMPLE.TIMES.find.distribution-1)/SAMPLE.TIMES.find.distribution,
                                           mean(dist.variance.UB1.vec), mean(dist.variance.UB2.vec)))
  
  if(delta.coef.only == 1)
  {
    my.png.name = as.character(paste("Overestimated_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
  if(delta.coef.only == 2)
  {
    my.png.name = as.character(paste("Overestimated_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
  if(delta.coef.only == 3)
  {
    my.png.name = as.character(paste("Overestimated_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
}



#===Misspecified m: Case 2: when we overestimate m===#
no.m.carryover = 2
p.lag.length = 1

dist.HT.estimator.matrix.mUnder = c()
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
  randomization.points = c(1,seq(2*p.lag.length+1, time.horizon-p.lag.length, by = p.lag.length))
  
  #=====Randomization Distribution=====#
  dist.HT.estimator.vec = c()
  dist.variance.UB1.vec = c()
  dist.variance.UB2.vec = c()
  misspecified.estimand.vec = c()
  for (TIMES.temp in 1:SAMPLE.TIMES.find.distribution)
  {
    dist.assignment.path = generate_assignments(time.horizon, randomization.points)
    dist.outcome.path = generate_outcomes(assignment.path_ = dist.assignment.path)
    dist.HT.estimator = generate_estimator(randomization.points, dist.assignment.path, dist.outcome.path)
    dist.variance.UB1 = generate_variance_UB_OPTDesign(randomization.points, dist.assignment.path, dist.outcome.path, p.lag.length, 1)
    dist.variance.UB2 = generate_variance_UB_OPTDesign(randomization.points, dist.assignment.path, dist.outcome.path, p.lag.length, 2)
    misspecified.estimand = generate_misspecified_estimand(dist.assignment.path)
    
    dist.HT.estimator.vec = c(dist.HT.estimator.vec, dist.HT.estimator)
    dist.variance.UB1.vec = c(dist.variance.UB1.vec, dist.variance.UB1)
    dist.variance.UB2.vec = c(dist.variance.UB2.vec, dist.variance.UB2)
    misspecified.estimand.vec = c(misspecified.estimand.vec, misspecified.estimand)
  }
  dist.HT.estimator.matrix.mUnder = rbind(dist.HT.estimator.matrix.mUnder, 
                                          c(NA, mean(misspecified.estimand.vec),
                                            mean(dist.HT.estimator.vec), var(dist.HT.estimator.vec)*(SAMPLE.TIMES.find.distribution-1)/SAMPLE.TIMES.find.distribution,
                                            mean(dist.variance.UB1.vec), mean(dist.variance.UB2.vec)))
  
  if(delta.coef.only == 1)
  {
    my.png.name = as.character(paste("Underestimated_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
  if(delta.coef.only == 2)
  {
    my.png.name = as.character(paste("Underestimated_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col="red", lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
  if(delta.coef.only == 3)
  {
    my.png.name = as.character(paste("Underestimated_m,delta=",delta.coef.only,".png"))
    png(filename = my.png.name, width = 600, height = 600)
    par(mar=c(4.5,5,0,0.5))
    hist(dist.HT.estimator.vec,
         main = "",
         xlab = "Value of Estimator",
         xlim = c(min(dist.HT.estimator.vec), max(dist.HT.estimator.vec)),
         freq = FALSE, breaks = 21, cex = 3, cex.axis = 2.5, cex.lab = 2.5)
    abline(v = mean(dist.HT.estimator.vec), col=2, lwd=2, lty=2)
    text(mean(dist.HT.estimator.vec),0, paste(round(mean(dist.HT.estimator.vec), digits = 3)), col = 2, adj = c(-0.2, 0), cex = 3)
    dev.off()
  }
}


experiment3.output.matrix = rbind(dist.HT.estimator.matrix.mCorrect, dist.HT.estimator.matrix.mOver, dist.HT.estimator.matrix.mUnder)

write.table(experiment3.output.matrix, file = paste(Sys.Date(), "experiment3.txt"), append = FALSE, sep = "\t", col.names=FALSE)

experiment3.output.matrix.copy = experiment3.output.matrix
experiment3.output.matrix.copy[,c(4,5,6)] = round(experiment3.output.matrix.copy[,c(4,5,6)], digits = 2)
experiment3.output.matrix.copy[,c(3)] = round(experiment3.output.matrix.copy[,c(3)], digits = 3)
write.table(experiment3.output.matrix.copy, file = paste(Sys.Date(), "experiment3_rounded.txt"), append = FALSE, sep = "\t", col.names=FALSE)
