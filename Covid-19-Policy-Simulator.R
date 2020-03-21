## Simulation of management of Covid-19 outbreak in the UK.
library(nleqslv)

# Parameters (unit: days)
incubation_time <- 5.1
start_trans_symp <- incubation_time - 0.5
start_trans_asymp <- 4.6
prop_symptoms <- 0.5 
mean_generation_time <- 6.5
R0 <- 2.4
inf_symp_asymp_ratio <- 2
UK_population_size <- 66.44e6
pop_size <- UK_population_size

# stages:
stages <- c(
  "S", # susceptible
  "Is1", # symptomatic infected no trans
  "Is2", # symptomatic infected trans
  "Ia1", # symptomatic asymptomatic no trans
  "Ia2", # symptomatic asymptomatic trans
  "R", # recovered or deceased
  c() )
n_stages <- length(stages)

# Linear dynamics:
SS_calc <- function(Is2_inf,rec_rate){
  SS <- matrix(0,nrow = n_stages,ncol= n_stages,
               dimnames = list(Next=stages,Current=stages) )
  SS["S","S"] <- 1
  SS["Is2","Is1"] <- 1/start_trans_symp
  SS["Is1","Is1"] <- 1-1/start_trans_symp
  SS["Ia2","Ia1"] <- 1/start_trans_asymp
  SS["Ia1","Ia1"] <- 1-1/start_trans_asymp
  SS["R","Ia2"] <- rec_rate
  SS["Ia2","Ia2"] <- 1-rec_rate
  SS["R","Is2"] <- rec_rate
  SS["Is2","Is2"] <- 1-rec_rate
  SS["R","R"] <- 1
  return(SS)
}

FF0_calc <- function(Is2_inf,rec_rate){
  FF0 <- matrix(0,nrow = n_stages,ncol= n_stages,
                dimnames = list(Next=stages,Current=stages) )
  FF0["Is1","Is2"] <- Is2_inf * prop_symptoms 
  FF0["Ia1","Is2"] <- Is2_inf * (1-prop_symptoms)
  FF0["Is1","Ia2"] <- Is2_inf * prop_symptoms  / inf_symp_asymp_ratio
  FF0["Ia1","Ia2"] <- Is2_inf * (1-prop_symptoms)  / inf_symp_asymp_ratio
  return(FF0)
}

FF_calc <- function(Is2_inf,rec_rate){
  FF0 <- FF0_calc(Is2_inf,rec_rate)
  return(FF0 * pop_size)
}

get_mean_generation_time <- function(SS,FF){
  L <- SS+FF
  w <- eigen(L)$vectors[,1]
  v <- eigen(t(L))$vectors[,1]
  T_gen <- 1 + (t(v) %*% SS %*% w) / (t(v) %*% FF %*% w) 
  return(T_gen)
}

R_and_Tgen <- function(Is2_inf,rec_rate){
  SS <- SS_calc(Is2_inf,rec_rate)
  FF <- FF_calc(Is2_inf,rec_rate)
  L <- SS + FF
  Tgen <- get_mean_generation_time(SS,FF)
  r <- eigen(L)$values[1]
  return(c(r*Tgen,Tgen))
}

R_and_Tgen(1/pop_size,1/6)
R_and_Tgen_diff <- function(x){
  R_and_Tgen(x[1],x[2]) - c(2.4,6.5)
}

fit <- nleqslv(c(Is2_inf=1/pop_size,rec_rate=1/6.5),R_and_Tgen_diff);
fit$x

SS <- with(as.list(fit),SS_calc(Is2_inf,rec_rate))
FF <- with(as.list(fit),FF_calc(Is2_inf,rec_rate))
FF0 <- with(as.list(fit),FF0_calc(Is2_inf,rec_rate))
L <- SS + FF
  
IStages <- c("Is1", "Is2", "Ia1", "Ia2")
#get_mean_generation_time(SS[IStages,IStages],FF[IStages,IStages])
log(eigen(L)$values[1]) # linear growth rate r should be around 0.2, OK

# Non-linear dynamics
## Simulate policy interventions with uncertainty
n_steps <- 365
state <- matrix(0, nrow=n_steps+1, ncol=n_stages)
colnames(state) <- stages
S_index <- which(stages == "S")
state[1,"Ia1"] <- 100
state[1,"S"] <- pop_size - state[1,"Ia1"]
dfac <- rep(NA,n_steps+1) # social distancing factor, the controll parameter
dfac[1] <- 1 # social distancing factor, the controll parameter
for(i in 1:n_steps){
  dfac[i+1] <- pnorm(qnorm(dfac[i])+rnorm(1,1/30,1/30)) # decline of compliance
  if(i %% 7 == 0){# Every week, adjust policy:
    if(state[i,"Is2"]/pop_size > 0.005) # more measures
      dfac[i+1] <- dfac[i+1] * runif(1,0.2,1) # outcome is uncertain
    else if(state[i,"Is2"]/pop_size < 0.005) # less measures
      dfac[i+1] <- dfac[i+1]*min(1,runif(1,1,1.5)) # outcome uncertain
  }
  Lx <- SS + FF0 * state[i,"S"] * dfac[i+1]
  state[i+1,] <- Lx %*% state[i,]
  state[i+1,"S"] <- pop_size - sum(state[i+1,-S_index])
}
dev.off()
par(mfrow=c(3,1))
plot(state[,"R"]/pop_size,type="l",ylim=c(0,1),ylab="Proportion immunized",xlab="Days")
plot(dfac,col="blue",type="l",ylim=c(0,1),ylab="Social Distancing Factor",xlab="Days")
plot(rowSums(state[,IStages])/pop_size,col="red",type="l",ylim=c(0,1),ylab="Proportion infected",xlab="Days")

