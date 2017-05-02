shotData<- c(1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0,
             1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1,
             0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0)

#figure 1 from blog, likelihood curve for 58/100 shots

x = seq(.001, .999, .001) ##Set up for creating the distributions
y2 = dbeta(x, 1 + 58, 1 + 42) # data for likelihood curve, plotted as the posterior from a beta(1,1)

plot(x, y2, xlim=c(0,1), ylim=c(0, 1.25 * max(y2,1.6)), type = "l", ylab= "Density", lty = 3,
     xlab= "Probability of success", las=1, main="Likelihood Curve for 3-pt Shots", sub= "(Binomial Data, 58/100)",lwd=2,
     cex.lab=1.5, cex.main=1.5, col = "darkorange", axes=FALSE)
axis(1, at = seq(0,1,.2)) #adds custom x axis
axis(2, las=1) # custom y axis

## Function for plotting priors, likelihoods, and posteriors for binomial data
## Output consists of a plot and various statistics
## PS and PF determine the shape of the prior distribution.
## PS = prior success, PF = prior failure for beta dist. 
## PS = 1, PF = 1 corresponds to uniform(0,1) and is default. If left at default, posterior will be equivalent to likelihood
## k = number of observed successes in the data, n = total trials. If left at 0 only plots the prior dist.
## null = Is there a point-null hypothesis? null = NULL leaves it out of plots and calcs
## CI = Is there a relevant X% credibility interval? .95 is recommended and standard

plot.beta <- function(PS = 1, PF = 1, k = 0, n = 0, null = NULL, CI = NULL, ymax = "auto", main = NULL) {
  
  x = seq(.001, .999, .001) ##Set up for creating the distributions
  y1 = dbeta(x, PS, PF) # data for prior curve
  y3 = dbeta(x, PS + k, PF + n - k) # data for posterior curve
  y2 = dbeta(x, 1 + k, 1 + n - k) # data for likelihood curve, plotted as the posterior from a beta(1,1)
  
  if(is.numeric(ymax) == T){ ##you can specify the y-axis maximum
    y.max = ymax
  }        
  else(
    y.max = 1.25 * max(y1,y2,y3,1.6) ##or you can let it auto-select
  )
  
  if(is.character(main) == T){
    Title = main
  }
  else(
    Title = "Prior-to-Posterior Transformation with Binomial Data"
  )
  
  
  plot(x, y1, xlim=c(0,1), ylim=c(0, y.max), type = "l", ylab= "Density", lty = 2,
       xlab= "Probability of success", las=1, main= Title,lwd=3,
       cex.lab=1.5, cex.main=1.5, col = "skyblue", axes=FALSE)
  
  axis(1, at = seq(0,1,.2)) #adds custom x axis
  axis(2, las=1) # custom y axis
  
  
  
  if(n != 0){
    #if there is new data, plot likelihood and posterior
    lines(x, y2, type = "l", col = "darkorange", lwd = 2, lty = 3)
    lines(x, y3, type = "l", col = "darkorchid1", lwd = 5)
    legend("topleft", c("Prior", "Posterior", "Likelihood"), col = c("skyblue", "darkorchid1", "darkorange"), 
           lty = c(2,1,3), lwd = c(3,5,2), bty = "n", y.intersp = .55, x.intersp = .1, seg.len=.7)
    
    ## adds null points on prior and posterior curve if null is specified and there is new data
    if(is.numeric(null) == T){
      ## Adds points on the distributions at the null value if there is one and if there is new data
      points(null, dbeta(null, PS, PF), pch = 21, bg = "blue", cex = 1.5)
      points(null, dbeta(null, PS + k, PF + n - k), pch = 21, bg = "darkorchid", cex = 1.5)
      abline(v=null, lty = 5, lwd = 1, col = "grey73")
      ##lines(c(null,null),c(0,1.11*max(y1,y3,1.6))) other option for null line
    }
  }
  
  ##Specified CI% but no null? Calc and report only CI
  if(is.numeric(CI) == T && is.numeric(null) == F){
    CI.low <- qbeta((1-CI)/2, PS + k, PF + n - k)
    CI.high <- qbeta(1-(1-CI)/2, PS + k, PF + n - k)
    
    SEQlow<-seq(0, CI.low, .001)
    SEQhigh <- seq(CI.high, 1, .001)
    ##Adds shaded area for x% Posterior CIs
    cord.x <- c(0, SEQlow, CI.low) ##set up for shading
    cord.y <- c(0,dbeta(SEQlow,PS + k, PF + n - k),0) ##set up for shading
    polygon(cord.x,cord.y,col='orchid', lty= 3) ##shade left tail
    cord.xx <- c(CI.high, SEQhigh,1) 
    cord.yy <- c(0,dbeta(SEQhigh,PS + k, PF + n - k), 0)
    polygon(cord.xx,cord.yy,col='orchid', lty=3) ##shade right tail
    
    return( list( "Posterior CI lower" = round(CI.low,3), "Posterior CI upper" = round(CI.high,3)))
  }
  
  ##Specified null but not CI%? Calculate and report BF only 
  if(is.numeric(null) == T && is.numeric(CI) == F){
    null.H0 <- dbeta(null, PS, PF)
    null.H1 <- dbeta(null, PS + k, PF + n - k)
    CI.low <- qbeta((1-CI)/2, PS + k, PF + n - k)
    CI.high <- qbeta(1-(1-CI)/2, PS + k, PF + n - k)
    return( list("BF01 (in favor of H0)" = round(null.H1/null.H0,3), "BF10 (in favor of H1)" = round(null.H0/null.H1,3)
    ))
  }
  
  ##Specified both null and CI%? Calculate and report both
  if(is.numeric(null) == T && is.numeric(CI) == T){
    null.H0 <- dbeta(null, PS, PF)
    null.H1 <- dbeta(null, PS + k, PF + n - k)
    CI.low <- qbeta((1-CI)/2, PS + k, PF + n - k)
    CI.high <- qbeta(1-(1-CI)/2, PS + k, PF + n - k)
    
    SEQlow<-seq(0, CI.low, .001)
    SEQhigh <- seq(CI.high, 1, .001)
    ##Adds shaded area for x% Posterior CIs
    cord.x <- c(0, SEQlow, CI.low) ##set up for shading
    cord.y <- c(0,dbeta(SEQlow,PS + k, PF + n - k),0) ##set up for shading
    polygon(cord.x,cord.y,col='orchid', lty= 3) ##shade left tail
    cord.xx <- c(CI.high, SEQhigh,1) 
    cord.yy <- c(0,dbeta(SEQhigh,PS + k, PF + n - k), 0)
    polygon(cord.xx,cord.yy,col='orchid', lty=3) ##shade right tail
    
    return( list("BF01 (in favor of H0)" = round(null.H1/null.H0,3), "BF10 (in favor of H1)" = round(null.H0/null.H1,3),
                 "Posterior CI lower" = round(CI.low,3), "Posterior CI upper" = round(CI.high,3)))
  }
  
}

#plot dimensions (415,550) for the blog figures

#Initial Priors
plot.beta(1,1,ymax=3.2,main="Uniform Prior, Beta(1,1)")
plot.beta(.5,.5,ymax=3.2,main="Jeffreys's Prior, Beta(1/2,1/2)")
plot.beta(4,9,ymax=3.2,main="Informed Prior, Beta(4,9)")

#Posteriors after Round 1
plot.beta(1,1,13,25,main="Beta(1,1) to Beta(14,13)",ymax=10)
plot.beta(.5,.5,13,25,main="Beta(1/2,1/2) to Beta(13.5,12.5)",ymax=10)
plot.beta(4,9,13,25,main="Beta(4,9) to Beta(17,21)",ymax=10)

#Posteriors after Round 2
plot.beta(14,13,12,25,ymax=10,main="Beta(14,13) to Beta(26,26)")
plot.beta(13.5,12.5,12,25,ymax=10,main="Beta(13.5,12.5) to Beta(25.5,25.5)")
plot.beta(17,21,12,25,ymax=10,main="Beta(17,21) to Beta(29,34)")

#Posteriors after Round 3
plot.beta(26,26,14,25,ymax=10,main="Beta(26,26) to Beta(40,37)")
plot.beta(25.5,25.5,14,25,ymax=10,main="Beta(25.5,25.5) to Beta(39.5,36.5)")
plot.beta(29,34,14,25,ymax=10,main="Beta(29,34) to Beta(43,45)")

#Posteriors after Round 4
plot.beta(40,37,19,25,ymax=10,main="Beta(40,37) to Beta(59,43)")
plot.beta(39.5,36.5,19,25,ymax=10,main="Beta(39.5,36.5) to Beta(58.5,42.5)")
plot.beta(43,45,19,25,ymax=10,main="Beta(43,45) to Beta(62,51)")

#Initial Priors and final Posteriors after all rounds at once
plot.beta(1,1,58,100,ymax=10,main="Beta(1,1) to Beta(59,43)")
plot.beta(.5,.5,58,100,ymax=10,main="Beta(1/2,1/2) to Beta(58.5,42.5)")
plot.beta(4,9,58,100,ymax=10,main="Beta(4,9) to Beta(62,51)")

