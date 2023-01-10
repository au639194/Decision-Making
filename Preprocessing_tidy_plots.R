#install.packages("pacman")
#pacman::p_load(tidyverse, R2jags, parallel, polspline, dplyr)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

##################################
## pre-processing
##################################

##### info about game #####
groupSize <- 4
ntrials <- 10
pi <- 1.4 # multiplication factor in game
ntokens <- 20
vals <- seq(0,ntokens,1) #possible values to contribute - from 0 to 20 tokens

##### load data #####
raw <- read.csv('Conditional Cooperation/data/HerrmannThoeniGaechterDATA.csv', skip = 3, header = T)

### removing punishment condition - remove line if you want to retain 
raw <- filter(raw, p == "N-experiment")

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
df <- merge(raw, data , on = "city") %>% 
  filter(row_number() %% 3 == 1)

#### setting up other things ####

## group average contribution on each trial 
ngroups = length(unique(df$groupid))

Ga_np <- df %>% 
  group_by(groupid, period) %>% # add p if you are interested in punishment conditions
  summarise(mean_send = mean(senderscontribution))
# make it an array 
Ga <- array(data = Ga_np$mean_send, dim = c(ntrials, ngroups))

## winnings
winnings <- df %>% 
  group_by(groupid) %>% 
  summarise(win = sum(senderscontribution)*pi)
# make it an array
winnings <- array(data = unlist(winnings[2]))

## getting nation for every group 
Nation <- df %>% 
  group_by(groupid) %>% 
  summarise(nation = mean(nation)) %>% 
  select(nation)
# make it an array
Nation <- array(data = unlist(Nation))

################# JAGS - CORRELATION: WINNING, donation ###############

ngroups = length(unique(df$groupid))
nnations = length(unique(df$nation))

# standardize and make arrays of X and Y
Y = array((winnings - mean(winnings)) / sd(winnings))
X = array((donation - mean(donation)) / sd(donation))

invSigma <- solve(t(X)%*%X) # required for JZS priors

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX", "nat_mu") 


#######################################
# run jags code
#######################################
win.samples <- jags.parallel(data, inits=NULL, params,
                             model.file ="Conditional Cooperation/win_corr.txt",
                             n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

# With a precision prior on the slope of just 1 (rather than a JZS prior)
win_priors.samples <- jags.parallel(data, inits=NULL, params,
                                    model.file ="Conditional Cooperation/win_corr_priors.txt",
                                    n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)


###############################
### Plotting 
###############################

# Picking the maximum a posteori probability + 95% credible interval of the nat_mu
MAP_nat_mu <- array(NA,nnations)
nat_mu_cred <- array(NA,c(nnations,2))
for (n in 1:nnations) {
  MAP_nat_mu[n] <- MPD(win_priors.samples$BUGSoutput$sims.list$nat_mu[,n])
  nat_mu_cred[n,] <- quantile(win_priors.samples$BUGSoutput$sims.list$nat_mu[,n],c(0.025,0.975))
}

# getting an index of the sorted the donation-values
idx <- order(donation)

# plotting the estimated national means (nat_mu) of the z-scored winnings against the "raw" donation-coefficients
plot(donation[idx], MAP_nat_mu[idx], ylim=c(-0.6,0.8))
for (n in 1:nnations) {
  lines(rep(donation[idx][n],2), nat_mu_cred[idx,][n,])
}

# plotting the "raw" national means of the winnings against the "raw" donation-coefficients
mean_win <- array(NA, length(donation))
par(mfrow = c(1,1))
for (n in idx) {
  if (n==idx[1]){
    plot(donation[n], mean(winnings[Nation==n]), xlim=c(.15, 1), ylim=c(-3, 800))
  } else {
    points(donation[n], mean(winnings[Nation==n]))
  }
  mean_win[n] <- mean(winnings[Nation==n])
  lines(c(donation[n],donation[n]),
        c(mean(winnings[Nation==n])-sd(winnings[Nation==n]), 
          mean(winnings[Nation==n])+sd(winnings[Nation==n]))
  )
}

# plotting with ggplot and a smoothed linear fit
df <- data.frame(donation[idx], mean_win[idx])
pl <- ggplot(df, aes(x = donation.idx.,
                     y = mean_win.idx.)) + 
  geom_point() + 
  labs(title = "Correlation of donation and mean winnings",
       x = "Donation amount as % of gross national income",
       y = "Mean winnings") +
  geom_smooth(method = "lm", se = T, formula = "y ~ x")
pl