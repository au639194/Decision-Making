seed_id = 1982
set.seed(seed_id)

#install.packages("pacman")
#pacman::p_load(tidyverse, R2jags, parallel, polspline, ggplot2, glue, dplyr, coda, rjags, bayesplot)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

### looking at describing outcome in public goods from donation-coefficient 

########################
# tidy pre-processing
########################


##### info about game #####
groupSize <- 4
ntrials <- 10
pi <- 1.6 # multiplication factor in game from Josh's slides 
ntokens <- 20
vals <- seq(0,ntokens,1) #possible values to contribute - from 0 to 20 tokens

##### load data #####
raw <- read.csv('Conditional Cooperation/Module5/data/HerrmannThoeniGaechterDATA.csv', skip = 3, header = T)

### setting up cities with their corresponding nation id + getting their data 
city = c("Melbourne", "Copenhagen", "Bonn", "Athens", 
         "Seoul", "Zurich", "St. Gallen", "Istanbul", "Nottingham", "Boston")
nation = c(1,2,3,4,5,6,6,7,8,9) # 2, 3, 8, 12 removed - lack of data

id = seq(1,9,1)
donation = c(.22, .71, .76, .16, .16, .50, .96, .50, .20)
ind = c(90, 74, 67, 35, 18, 68, 37, 89, 91)

### first round of merge
data <- data.frame(id, ind, donation) %>% # gdp, trust, donation removed
  rename(nation = id)
cities = data.frame(city, nation)
data = merge(data, cities)

### second round
redDat <- left_join(raw, data , on = "city") %>% 
  filter(row_number() %% 3 == 1)


##############################################
# full_PGG_analysis - clean and arrange data
##############################################

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# remove data for where there are no subjects
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)
missing <- array(0,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}

# data for punishment condition #
singlechild <- array(0,c(groupSize,ngroups)) # ADDED THIS CUZ WAS MISSING
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    singlechild[s,g] <- redDat$singlechild[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][s+(s-1)*10]
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
}


# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish

Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish


Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

c_choice_index <- c

donation <- array(0, ngroups)
Nation <- array(0, ngroups)
for (g in 1:ngroups) {
  donation[g] <- mean(redDat$donation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
  Nation[g] <- mean(redDat$nation[redDat$groupid==group_names[g]&redDat$p=="N-experiment"])
}

######################
# removing NA's
######################

c_win <- c_no_punish[,,!is.na(donation)]
c_keep <- rep()-c_win

Gga_punish <- Gga_punish[,!is.na(donation)]
Gga_no_punish <- Gga_no_punish[,!is.na(donation)]

c <- c[,,!is.na(donation),]
Gga <- Gga[,!is.na(donation),]
Ggas <- Ggas[,,!is.na(donation),]
Gc <- Gc[,,!is.na(donation),]
Ga <- Ga[,,!is.na(donation),]
donation <- donation[!is.na(donation)]
Nation <- Nation[!is.na(Nation)]

# redefine number of groups after removing those without civic scores
ngroups <- length(donation)

# aggregate donation to just 1 number per Nation-index (using the mean here should be unproblematic since all groups within a given nation should have been given the same donation-coefficient)
donation <- aggregate(donation~Nation, FUN=mean)[,2]

nnations <- length(donation)

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

################################################################################
########################### Conditional cooperation model ######################
################################################################################

# JZS priors for partial correlation. Method described here
# https://link.springer.com/article/10.3758/s13423-012-0295-x
# Code available here
# https://github.com/MicheleNuijten/BayesMed/blob/master/R/jzs_corSD.R
# Paper where code is used here (mediation paper)
# https://link.springer.com/article/10.3758/s13428-014-0470-2

##############################################################
# Regresss donation on winnings - same as preprocessing script
##############################################################

# standardize and make arrays of X and Y
Y = array((winnings - mean(winnings)) / sd(winnings))
X = array((donation - mean(donation)) / sd(donation))
#X = array((donation-mean(donation))/sd(donation))

invSigma <- solve(t(X)%*%X) # required for JZS priors

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX", "nat_mu") 

# run jags
win.samples <- jags.parallel(data, inits=NULL, params,
                             model.file ="Conditional Cooperation/Module5/win_corr.txt",
                             n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)


#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  Regress donation on belief weights and slope of prefs in CC model ---------------

# standardise covariate
X <- donation
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

Ga_old <- Ga
Ga <- Ggas

data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("betaX_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                   model.file ="Conditional Cooperation/Module5/CC_corr.txt",
                   n.chains=3, n.iter=30000, n.burnin=8000, n.thin=20, n.cluster=4)
end_time = Sys.time()
end_time - start_time

# saving the model
#save(CC.samples, file= "Conditional Cooperation/saved_JAGS/CC.samples.Rdata")

# reading in the model
#load("Conditional Cooperation/saved_JAGS/CC.samples_burn2.Rdata")

#################################################################
#------------------ Plotting ---------------------------
#################################################################
# looking at convergence
#traceplot(CC.samples)

# diagnostics for convergence
#gelman.diag(as.mcmc(CC.samples))
#coda::effectiveSize(as.mcmc(CC.samples))

# ------ Create empirical and parameter arrays for plots -------
# empirical group winnings data - means and standard deviations
empirical.win <- array(0,c(3,length(donation)))
for (i in 1:length(donation)) {
  empirical.win[1,i] <- mean(winnings[Nation==i]) - sd(winnings[Nation==i]) 
  empirical.win[2,i] <- mean(winnings[Nation==i]) 
  empirical.win[3,i] <- mean(winnings[Nation==i]) + sd(winnings[Nation==i])
}

empirical.initial <- array(0,c(3,length(donation)))
for (i in 1:length(donation)) {
  empirical.initial[1,i] <- mean(colMeans(c[,1,,1])[Nation==i]) - sd(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[2,i] <- mean(colMeans(c[,1,,1])[Nation==i]) 
  empirical.initial[3,i] <- mean(colMeans(c[,1,,1])[Nation==i]) + sd(colMeans(c[,1,,1])[Nation==i])
}

#----- empirical correlation - donation and initial contribution --------------------
plot(c(.1,1), c(0,20), type = "n", main = "Data - Initial Contrib.", 
     xlab = "National donation Coefficient", ylab = "Initial Contribution",axes=FALSE)
for (i in 1:length(donation)) {
  lines(c(donation[i],donation[i]),c(empirical.initial[1,i],empirical.initial[3,i]))
  points(donation[i],empirical.initial[2,i])
}
axis(1)
axis(2)

#----- posterior of effect - donation and group winnings --------------------
plot(density(win.samples$BUGSoutput$sims.list$betaX),frame=FALSE,lwd=2,ylim=c(0,7),
     cex=2,xlab = "Standardised donation Estimate",ylab = " ",main="B: Contribution")
COL <- adjustcolor(c("red"))
lines(c(win.samples$BUGSoutput$summary[2,3],win.samples$BUGSoutput$summary[2,7]),c(-.1,-.1),col=COL,lwd=2)
points(win.samples$BUGSoutput$summary[2,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title


#---------- Initial Belief -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_alpha),frame=FALSE,lwd=2,ylim=c(0,8),
     cex=2,xlab = "Standardised donation Estimate",ylab = " ",main="D: Initial Belief")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[3,3],CC.samples$BUGSoutput$summary[3,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[3,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2)) # reset margin because no title


#---------- Belief Learning -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_omega),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardised donation Estimate",ylab = " ",main="E: Belief Learning Weight")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[4,3],CC.samples$BUGSoutput$summary[4,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[4,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))


#---------- Preference Slope -------------------------------------------------------

plot(density(CC.samples$BUGSoutput$sims.list$betaX_rho),frame=FALSE,lwd=2,ylim=c(0,10),
     cex=2,xlab = "Standardised donation Estimate",ylab = " ",main="F: Conditional Preferences")
COL <- adjustcolor(c("red"))
lines(c(CC.samples$BUGSoutput$summary[5,3],CC.samples$BUGSoutput$summary[5,7]),c(-.1,-.1),col=COL,lwd=2)
points(CC.samples$BUGSoutput$summary[5,1],c(-.1),pch=19,col=COL)
abline(h=0,v=0,col="gray60")

par(mar=c(5,3,1,2))
