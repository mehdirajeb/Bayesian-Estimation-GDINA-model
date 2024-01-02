library(GDINA)
library(rjags)

####################################################################
#Data form GDINA package
###################################################################

Y <- Rep3_score[,-1]
Q <- Q
all.patterns <- GDINA::attributepattern(ncol(Q))


###################################################################################
#GDINA model information
# For 3 attribute
# n= individual
# i = item
# j = attribute
# alpha(n,k) = individual's attribute profile
###################################################################################

##############################################
#JAGS model
##############################################

jags.gdina <- function() {
  for (n in 1:N) {
    for (i in 1:I) { # for k= 3, we may estimate maximum 8 parameters for GDINA model. 
      
      eta1[n, i] <- lamda1[i]* alpha[n,1] + lamda2[i]* alpha[n,2]+ lamda3[i]* alpha[n,3]
      
      eta2[n,i] <- lamda12[i]* alpha[n,1]* alpha[n,2] + lamda13[i]* alpha[n,1]* alpha[n,3] + lamda23[i]* alpha[n,2]* alpha[n,3]
      
      eta3[n,i] <- lamda123[i]* alpha[n,1]* alpha[n,2]* alpha[n,3]
      
      logit(p[n, i]) <- lamda0[i]+ eta1[n,i] + eta2[n,i]+ eta3[n,i] ## used logit link
      
      Y[n, i] ~ dbern(p[n, i])
    }
    
    for (k in 1:K) {
      alpha[n, k] <- all.patterns[c[n], k]
    }
    c[n] ~ dcat(pai[1:C])
  }
  pai[1:C] ~ ddirch(delta[1:C])
  
  for (i in 1:I) {
    lamda0[i] ~ dnorm(-1.096, 0.25)
    xlamda1[i] ~ dnorm(0, 0.25)%_% T(0,)
    xlamda2[i] ~ dnorm(0, 0.25)%_% T(0,) 
    xlamda3[i] ~ dnorm(0, 0.25)%_% T(0,) 
    xlamda12[i] ~ dnorm(0, 0.25)
    xlamda13[i] ~ dnorm(0, 0.25)
    xlamda23[i] ~ dnorm(0, 0.25)
    xlamda123[i] ~ dnorm(0, 0.25)
    lamda1[i] <- xlamda1[i] * Q[i, 1]
    lamda2[i] <- xlamda2[i] * Q[i, 2]
    lamda3[i] <- xlamda3[i] * Q[i, 3]
    lamda12[i] <- xlamda12[i] * Q[i, 1] * Q[i, 2]
    lamda13[i] <- xlamda13[i] * Q[i, 1] * Q[i, 3]
    lamda23[i] <- xlamda23[i] * Q[i, 2] * Q[i, 3]
    lamda123[i] <- xlamda123[i] * Q[i, 1] * Q[i, 2] * Q[i, 3]}
  
  
  ## the posterior predictive model checking## I need to understand this part more thoroughly
  for (n in 1:N) {
    for (i in 1:I) {
      teststat[n, i] <- pow(Y[n, i] - p[n, i], 2)/(p[n, i] * (1 -
                                                                p[n, i]))
      Y_rep[n, i] ~ dbern(p[n, i])
      teststat_rep[n, i] <- pow(Y_rep[n, i] - p[n, i], 2)/(p[n, i] *
                                                             (1 - p[n, i]))
    }
  }
  teststatsum <- sum(teststat[1:N, 1:I])
  teststatsum_rep <- sum(teststat_rep[1:N, 1:I])
  ppp <- step(teststatsum_rep - teststatsum)
}
######################################################################################

## USEd 200 items for practice purpose
N <- 935 
I <- nrow(Q)
K <- ncol(Q)
C <- nrow(all.patterns)
delta <- rep(1, C)
jags.data <- list("N", "I", "K", "Y", "Q", "C", "all.patterns", "delta")
jags.parameters <- c("lamda0","lamda1","lamda2", "lamda3", "lamda12","lamda13","lamda23","lamda123", "c", "pai", "ppp")  
#c: estimated attribute pattern of each respondent 
#ppp: posterior predictive probability 
jags.inits <- NULL #Initial values are not specified.
time1 <- as.POSIXlt(Sys.time())

#####################################################################################
#RJAGs code
####################################################################################

sim <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.parameters,
            model.file = jags.gdina, n.chains = 2, n.iter = 20000, n.burnin = 10000,
            n.thin = 1, DIC = TRUE)
time2 <- as.POSIXlt(Sys.time())
use.time <- difftime(time2, time1, units = "secs")

######################################################################################
#Outputs
######################################################################################
sim1 <- sim$BUGSoutput
E.pattern <- cbind(sim1$median$c)
sim1$mean$ppp

write.csv(sim1$summary, "Rep3_summary_gdina_20k.csv")

E.itempar <- cbind(sim1$mean$lamda0, sim1$sd$lamda0, sim1$mean$lamda1,sim1$sd$lamda1,
                   sim1$mean$lamda2, sim1$sd$lamda2, sim1$mean$lamda3,sim1$sd$lamda3,
                   sim1$mean$lamda12, sim1$sd$lamda12, sim1$mean$lamda13,sim1$sd$lamda13,
                   sim1$mean$lamda23, sim1$sd$lamda23, sim1$mean$lamda123,sim1$sd$lamda123)
colnames(E.itempar) <- c("lamda0", "SE(lamda0)", "lamda1", "SE(lamda1)", "lamda2", "SE(lamda2)", 
                         "lamda3", "SE(lamda3)","lamda12", "SE(lamda12)", "lamda13", "SE(lamda13)",
                         "lamda23", "SE(lamda23)","lamda123", "SE(lamda123)")
E.itempar
