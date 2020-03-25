## Simulation of management of Covid-19 outbreak in the UK.
library(nleqslv) # for nonlinear solver 'nleqslv'

#### Parameters (time unit: days)
## Epidemiological parameters:
incubation_time <- 5.1 # Ferguson et al. 
incubation_time <- incubation_time +1 # fudge to reproduce R0 and generation time 
start_trans_symp <- incubation_time - 0.5 # Ferguson et al. 
start_trans_asymp <- start_trans_symp # Ferguson et al. 
#prop_symptoms <- 1/3 # Ferguson et al. 
prop_symptoms <- 0.2 # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30567-5/fulltext
mean_generation_time <- 6.5 # Ferguson et al. 
#generation_time_alpha <- 3 # shape parameter of an assumed gamma distribution, https://academic.oup.com/biostatistics/article/12/2/303/280518
generation_time_alpha <- 2 # shape parameter of an assumed gamma distribution, fudge to reproduce r 
R0 <- 2.4 # Ferguson et al. 
inf_symp_asymp_ratio <- 2  # ratio of symptomatic to asymptomatic, Ferguson et al.  
UK_population_size <- 66.44e6
pop_size <- UK_population_size

## Parameters controlling mortality:
case_fatality_rate_1 <- 0.009 # hospitals below capacity
case_fatality_rate_2 <- 0.09 # hospitals over capacity
hospitalization_rate <- 0.044 # proportion of infected hospitlized
critical_care_rate <- 0.044 * 0.3 # proportion of infected needing critical care
critical_care_hospital_time <- 10
# number of infected over number of critical care cases at any moment:
hospital_load_factor <- 
  critical_care_rate * critical_care_hospital_time / 
  mean_generation_time
critical_care_beds <- 8175 # uk number of respirators 
## Threshold for proportion infected:
critical_care_beds / (hospital_load_factor * pop_size) 

## Initial linear growth rate of epidemy:
# Formula assuming gamma distribution for generation time
r0 <- (R0^(1/generation_time_alpha)-1) * generation_time_alpha/mean_generation_time;
r0
# doubling time
#exp(r*T_double)==2 -> T_double == log(2)/r
T_double <- log(2)/r0; T_double # consistency check


## Policy and behaviour parameters:
policy_intervention_rate <- 1/6 # https://lab.gedidigital.it/gedi-visual/2020/coronavirus-i-contagi-in-italia/


## The stage-structured model:
# Stages:
stages <- c(
  "S", # susceptible
  "Is1", # symptomatic infected no trans
  "Is2", # symptomatic infected trans
  "Ia1", # symptomatic asymptomatic no trans
  "Ia2", # symptomatic asymptomatic trans
  "R", # recovered or deceased
  c() )
n_stages <- length(stages)

# Infected stages
IStages <- c("Is1", "Is2", "Ia1", "Ia2")
# Infectious (test postive) stages
PStages <- c("Is2", "Ia2")


## Linear dynamics:
# Later we still need to adjust two parameters, which will be left open for now:
# Is2_inf: infectiousness of Is2 stage (which then scales other stages)
# rec_rate: recovery rate of transmitting individuals
# This adjustment is done using the full linear model such as to fit given values for R0 and mean_generation_time.

## The "survivial" part of the population matrix
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

## Density indepent factor of the "fertility" (= transmission) matrix
FF0_calc <- function(Is2_inf,rec_rate){
  FF0 <- matrix(0,nrow = n_stages,ncol= n_stages,
                dimnames = list(Next=stages,Current=stages) )
  FF0["Is1","Is2"] <- Is2_inf * prop_symptoms 
  FF0["Ia1","Is2"] <- Is2_inf * (1-prop_symptoms)
  FF0["Is1","Ia2"] <- Is2_inf * prop_symptoms  / inf_symp_asymp_ratio
  FF0["Ia1","Ia2"] <- Is2_inf * (1-prop_symptoms)  / inf_symp_asymp_ratio
  return(FF0)
}

## Calculate the "fertility" (= transmission) part of the population matrix
FF_calc <- function(Is2_inf,rec_rate){
  FF0 <- FF0_calc(Is2_inf,rec_rate)
  return(FF0 * pop_size)
}

## Decide which is the index of the eigenvalue determining linear growth rate:
eigenvalue_index <- function(eig_values){ # input is a list of eigenvalues
  w <- which(abs(eig_values)==eig_values & eig_values != 1)
  return(w[which.max(eig_values[w])])
}

get_mean_generation_time <- function(SS,FF){
  # See page 842 of Bienvenu, F.; Legendre, S. (2015). "A New Approach to the Generation Time in Matrix Population Models" (PDF). The American Naturalist. 185 (6): 834–843. doi:10.1086/681104. PMID 25996867.
  L <- SS+FF
  eig <- eigen(L)
  i <- eigenvalue_index(eig$values)
  w <- eig$vectors[,i]
  t_eig <- eigen(t(L))
  j <- min(which(abs(t_eig$values - eig$values[i])<1e-10))
  v <- t_eig$vectors[,j]
  if(abs(t(v) %*% w)<1e-10) {
    print(c(i,j))
    print(eig$values)
    print(t_eig$values)
    print(w)
    print(v)
    stop("Potential eigenvector index mixup")
  }
  T_gen <- 1 + (t(v) %*% SS %*% w) / (t(v) %*% FF %*% w) 
  return(T_gen)
}

## Computes both R0 and mean generation time for given SS and FF:
R0_and_Tgen <- function(SS,FF){
  L <- SS + FF
  Tgen <- get_mean_generation_time(SS,FF)
  eig <- eigen(L[IStages,IStages])
  i <- eigenvalue_index(eig$values)
  r <- log(eig$values[i])
  # https://beast.community/estimating_R0
  # https://www.sciencedirect.com/science/article/pii/S1755436518300847
  generation_time_beta <- generation_time_alpha / Tgen
  R0 <- (1+r/generation_time_beta)^generation_time_alpha
  return(c(R0,Tgen))
}

## Computes deviation of R0 and mean generation time from target values, for give for given Is2_inf, rec_rate, passed as vector in argument x:
R0_and_Tgen_diff <- function(x){
  as.numeric(
    with(as.list(x),
         R0_and_Tgen(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)) - 
           c(R0,mean_generation_time) ))
}
initialGuess <- c(Is2_inf=0.1/pop_size,rec_rate=0.01 )
# #Tests
# with(as.list(initialGuess),
#      SS_calc(Is2_inf,rec_rate)+FF_calc(Is2_inf,rec_rate))
# with(as.list(initialGuess),
#      get_mean_generation_time(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)))
# with(as.list(initialGuess),
#      R0_and_Tgen(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)))
R0_and_Tgen_diff(initialGuess)

## Find parameters to fit R0 and generation time
fit <- nleqslv(initialGuess,R0_and_Tgen_diff);
fit$message  # should be "Function criterion near zero" or "x-values within tolerance 'xtol'"
fit$x
## double check
with(as.list(fit$x),
     R0_and_Tgen(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)))

## Compute population matrices:
SS <- with(as.list(fit$x),SS_calc(Is2_inf,rec_rate))
FF <- with(as.list(fit$x),FF_calc(Is2_inf,rec_rate))
FF0 <- with(as.list(fit$x),FF0_calc(Is2_inf,rec_rate))
L <- SS + FF

  
linear_growth_rate <- # linear growth rate r should be around 0.2
  log(eigen(L[IStages,IStages])$values[1]) 
c(r0,linear_growth_rate) ## should be the same

## Now, simulate non-linear dynamics
n_steps <- 365 * 1.5
state <- matrix(0, nrow=n_steps+1, ncol=n_stages)
colnames(state) <- stages
S_index <- which(stages == "S")
state[1,"Ia1"] <- 100
state[1,"S"] <- pop_size - state[1,"Ia1"]
deaths <- rep(NA,n_steps+1) # deaths per step
deaths[1] <- 0

compute_deaths <- function(state_snippet){
  total_infected <- sum(state_snippet[1,IStages])
  respirator_demand <- hospital_load_factor * total_infected
  with_respirator <- min(critical_care_beds, respirator_demand)
  without_respirator <- respirator_demand - with_respirator
  fatality_rate <-
    (with_respirator * case_fatality_rate_1 + without_respirator * case_fatality_rate_2)/
    respirator_demand
  rec <- state_snippet[2,"R"] - state_snippet[1,"R"]
  deaths <- rec * fatality_rate
  return(deaths)
}

if(FALSE){ ## this passage is deactivated, used only for testing and model comparisons:
  ## Simulate without policy interventions
  plot(c(0),type="n",ylim=c(0,1),xlim=c(0,n_steps),xlab="Days",ylab="Proportion")
  for(i in 1:n_steps){
    Lx <- SS + FF0 * state[i,"S"]
    state[i+1,] <- Lx %*% state[i,]
    state[i+1,"S"] <- pop_size - sum(state[i+1,-S_index])
    deaths[i+1] <- compute_deaths(state[c(0,1)+i,])
  }
  lines((state[,"R"]-cumsum(deaths))/pop_size,col="black",type="l")
  lines(rowSums(state[,IStages])/pop_size,col="red",type="l")
  lines(cumsum(deaths)/pop_size*1e1,col="brown",type="l")
  legend(x="bottomright",
         legend = c("Recovered","Infected","cum. mortality x 10"),
         col=c("black","red","brown"),lty=1,bg=rgb(1,1,1,alpha = 0.4))
}

## Simulate policy interventions with uncertainty
dfac <- rep(NA,n_steps+1) # "social distancing factor", scales infection rate
dfac[1] <- 1 # value 1 corresponds to no distancing, the initial assumption
growth_rate <- linear_growth_rate # growth_rate is adjusted in simulations to inform policy
n_runs <- 200 # how often to repeat simulation
if(!exists("plotting")) plotting <- TRUE # set to FALSE if plotting takes too much time
other_runs_alpha <- 10/n_runs # intensity of plotting all but last simulation
final_mortality <- rep(NA,n_runs) # record total mortality at end of each simulation
final_immunization <- rep(NA,n_runs) # record total immunication rate at end of each simulation

# Prepare plot of results:
plot(c(0),type="n",ylim=c(0,1),xlim=c(0,n_steps),xlab="Days",ylab="Proportion")
# run n_runs simulations
for(run in 1:n_runs){
  policy_intervention_times <- c()
  for(i in 1:n_steps){
    dfac[i+1] <- 1-exp(min(0,(log(1-dfac[i])+rnorm(1,-1/(30*6),1/30)))) # decline of compliance
    if(runif(1)< policy_intervention_rate){# adjust policy:
      if(state[i,"Is2"]/pop_size > 0.0005 && growth_rate > 0){ # more measures
        dfac[i+1] <- dfac[i+1] * runif(1,0.2,1) # outcome is uncertain
        policy_intervention_times <-
          append(policy_intervention_times,i)
      }else if(dfac[i+1] < 1 && state[i,"Is2"]/pop_size < 0.0001 && growth_rate <  0){ # less measures
        dfac[i+1] <- min(1,dfac[i+1]*runif(1,1,2)) # outcome uncertain
        policy_intervention_times <-
          append(policy_intervention_times,i)
      }
    }
    Lx <- SS + FF0 * state[i,"S"] * dfac[i+1]
    state[i+1,] <- Lx %*% state[i,]
    state[i+1,"S"] <- pop_size - sum(state[i+1,-S_index])
    deaths[i+1] <- compute_deaths(state[c(0,1)+i,])
    growth_rate <- log(sum(state[i+1,PStages])/sum(state[i,PStages]))
  }
  if(plotting || run == n_runs){
    line_alpha <- ifelse(run > n_runs-1,1,other_runs_alpha)*255
    plotcol <- do.call(rgb,as.list(c(col2rgb("green"),alpha=line_alpha,max = 255)))
    lines((state[,"R"]-cumsum(deaths))/pop_size,col=plotcol,type="l")
    plotcol <- do.call(rgb,as.list(c(col2rgb("blue"),alpha=line_alpha,max = 255)))
    lines(1-dfac,col=plotcol,type="l")
    plotcol <- do.call(rgb,as.list(c(col2rgb("red"),alpha=line_alpha,max = 255)))
    lines(100*rowSums(state[,IStages])/pop_size,col=plotcol,type="l")
    plotcol <- do.call(rgb,as.list(c(col2rgb("black"),alpha=line_alpha,max = 255)))
    lines(cumsum(deaths)/pop_size * 1e1,col=plotcol,type="l")
  }
  final_mortality[run] <- sum(deaths)/pop_size
  final_immunization[run] <- 
    (state[n_steps,"R"]-sum(deaths))/(pop_size-sum(deaths))
}
for(t in policy_intervention_times){
  lines(c(t,t),c(0,1), col="black", lty=3)
}
legend(x="topright",
       legend = c("Recovered","Distancing","Infected x 100","cum. mort. x 10"),
       col=c("green","blue","red","black"),lty=1,bg=rgb(1,1,1,alpha = 0.8))
#plot(final_immunization,final_mortality,pch="+",cex=0.3)
hist(final_mortality)
