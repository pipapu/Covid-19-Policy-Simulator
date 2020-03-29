#### Simulation of management of Covid-19 outbreak in the UK. ####
library(nleqslv) # for nonlinear solver 'nleqslv'
#set.seed(9) #  3
# Parameters (unit: days) ------------------------------------------------------

## Parameters controlling disease dynamics:
incubation_time <- 5.1 # Period from infection to symptoms from Ferguson et al. 
start_trans_symp <- incubation_time - 0.5 # Period from infection to the start of infectiousness for symptomatic cases Ferguson et al. 
start_trans_asymp <- start_trans_symp # Period from infection to the start of infectiousness for asymptomatic casesFerguson et al. 
prop_symptoms <- 0.2 # Proportion of infected people developing symptoms https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30567-5/fulltext
#prop_symptoms <- 1/3 # Ferguson et al. 
prop_critical <- 0.06 # Proportion of cases becoming critical # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30567-5/fulltext
mean_generation_time <- 6.5 # Ferguson et al. 
#generation_time_alpha <- 3 # shape parameter of an assumed gamma distribution, https://academic.oup.com/biostatistics/article/12/2/303/280518
generation_time_alpha <- 2 # shape parameter of an assumed gamma distribution, fudge to reproduce r
R0 <- 2.4 # Basic reproductive number for the infection, Ferguson et al. 
inf_symp_asymp_ratio <- 2  # ratio infection rates by symptomatic to asymptomatic, Ferguson et al.  
UK_population_size <- 66e6
pop_size <- UK_population_size

## Parameters controlling mortality:
# We assume that case fataliy is much higher for cases beyond hospital capacity:
case_fatality_rate_1 <- 0.009 # case fatality rate when hospitals below capacity https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/ (as of 27 Mar 2020)
case_fatality_rate_2 <- 0.09 # case fatality rate when hospitals over capacity  https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/ (as of 27 Mar 2020)
hospitalization_rate <- 0.044 # proportion of cases needing hospitalisation Ferguson et al. 
critical_care_rate <- hospitalization_rate * 0.3 # proprtion of cases needing critical care Ferguson et al. 
critical_care_hospital_time <- 10 # Ferguson et al. 
# number of critical care cases over number of cases at any moment:
critical_care_beds <- 8175 # uk number of respirators available
## Threshold for proportion of cases p.p.:
critical_care_beds / (hospital_load_factor * pop_size) 

## Initial linear growth rate of epidemy:
# Formula assuming gamma distribution for generation time
r0 <- (R0^(1/generation_time_alpha)-1) * generation_time_alpha/mean_generation_time;
r0
# doubling time
#exp(r*T_double)==2 -> T_double == log(2)/r
T_double <- log(2)/r0; T_double # consistency check

## Policy and behaviour parameters:
policy_intervention_rate <- 1/3# 
policy_intervention_threshold <- 20*10/1e6 # case p.p. above/below which policy acts
policy_relaxation_growth_rate_threshold <- -0.05 # don't relax if growth is larger

compliance_sd <- 1/30 # Day to day standard deviation of compliance.
compliance_decline <- 1/(6*30) # Day to day decline of compliance.

## For plotting results:
date_matching_cumulative_fatality <- 55
date_matching_cumulative_date <- as.Date("2020-03-16")

# Stages: set up population categories for the different classes of infected and uninfected people 
stages <- c(
  "S", # susceptible
  "Is1", # symptomatic infected no transmission
  "Is2", # symptomatic infected transmission
  "Ia1", # symptomatic asymptomatic no transmission
  "Ia2", # symptomatic asymptomatic transmission
  "R", # recovered or deceased
  c() )
n_stages <- length(stages)

# Infected stages
IStages <- c("Is1", "Is2", "Ia1", "Ia2")

# Infectious (test postive) stages
PStages <- c("Is2", "Ia2")


# Linear dynamics -------------------------------------------------------------

# Later we still need to adjust two parameters, which will be left open for now: 
# Is2_inf: infectiousness of the Is2 stage (which then scales other stages)
# rec_rate: recovery rate of transmitting individuals
# This adjustment is done using the full linear model such as to fit given values for R0 and mean_generation_time.

# Generate the "survival" part of the population projection matrix

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

# Density independent factor of the "fertility" (= transmission) matrix
FF0_calc <- function(Is2_inf,rec_rate){
  FF0 <- matrix(0,nrow = n_stages,ncol= n_stages,
                dimnames = list(Next=stages,Current=stages) )
  FF0["Is1","Is2"] <- Is2_inf * prop_symptoms 
  FF0["Ia1","Is2"] <- Is2_inf * (1-prop_symptoms)
  FF0["S"  ,"Is2"] <- -Is2_inf
  FF0["Is1","Ia2"] <- Is2_inf * prop_symptoms  / inf_symp_asymp_ratio
  FF0["Ia1","Ia2"] <- Is2_inf * (1-prop_symptoms)  / inf_symp_asymp_ratio
  FF0["S"  ,"Ia2"] <- -Is2_inf / inf_symp_asymp_ratio
  return(FF0)
}

# The "fertility" (= transmission) part of the population projection matrix
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

## Calculates both R0 and mean generation time for given SS and FF:
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

# Calculate the difference between observed values for R0 and values calculated from the survival and transmission matrices
R0_and_Tgen_diff <- function(x){
  as.numeric(
    with(as.list(x),
         R0_and_Tgen(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)) - 
           c(R0,mean_generation_time) ))
}

# Initial values for Is2_inf and rec_rate
initialGuess <- c(Is2_inf=0.1/pop_size,rec_rate=0.01 )
# #Tests
# with(as.list(initialGuess),
#      SS_calc(Is2_inf,rec_rate)+FF_calc(Is2_inf,rec_rate))
# with(as.list(initialGuess),
#      get_mean_generation_time(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)))
# with(as.list(initialGuess),
#      R0_and_Tgen(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)))

# Differences between observed R0 and generation time calculated using the "initial guess" values and the observed values
R0_and_Tgen_diff(initialGuess)

## Generate values of Is2_inf and rec_rate which fit R0 and generation time
fit <- nleqslv(initialGuess,R0_and_Tgen_diff);
fit$message  # should be "Function criterion near zero" or "x-values within tolerance 'xtol'"
fit$x
## double check
with(as.list(fit$x),
     R0_and_Tgen(SS_calc(Is2_inf,rec_rate),FF_calc(Is2_inf,rec_rate)))
with(as.list(fit$x),
     R0_and_Tgen(SS_calc(Is2_inf,rec_rate)[IStages,IStages],FF_calc(Is2_inf,rec_rate)[IStages,IStages]))

## infection rate:
with(as.list(fit$x), pop_size * Is2_inf )

hospital_load_factor <- critical_care_rate * critical_care_hospital_time / 
  with(as.list(fit$x),incubation_time + 1/rec_rate)



## Compute population matrices:
SS <- with(as.list(fit$x),SS_calc(Is2_inf,rec_rate))
FF <- with(as.list(fit$x),FF_calc(Is2_inf,rec_rate))
FF0 <- with(as.list(fit$x),FF0_calc(Is2_inf,rec_rate))
L <- SS + FF # the complete population projection matrix


linear_growth_rate <- # linear growth rate r should be around 0.2
  log(eigen(L[IStages,IStages])$values[1]) 
c(r0,linear_growth_rate) ## should be the same

# Non-linear dynamics ----------------------------------------------------------

# Set up matrix for numbers in each state at each timepoint
n_steps <- 365 * 1.5 #How long to simulate for
state <- matrix(0, nrow=n_steps+1, ncol=n_stages) 
colnames(state) <- stages
S_index <- which(stages == "S") # index of S stage
state[1,"Ia1"] <- 10 # Seed with some asymptotic cases
state[1,"S"] <- pop_size - sum(state[1,-S_index]) # Initial number of uninfected individuals
deaths <- rep(NA,n_steps+1) # deaths per step
deaths[1] <- 0

# Function to compute the number of deaths
compute_deaths <- function(state_snippet){
  # We assume that symptomatic transmitting patients are "case"
  total_cases <- sum(state_snippet[1,IStages])
  respirator_demand <- hospital_load_factor * total_cases
  with_respirator <- min(critical_care_beds, respirator_demand)
  without_respirator <- respirator_demand - with_respirator
  if(total_cases > 0 ){
    fatality_rate <-
      (with_respirator * case_fatality_rate_1 + without_respirator * case_fatality_rate_2)/
      respirator_demand
  }else{
    fatality_rate <- 0
  }
  rec <- state_snippet[2,"R"]-state_snippet[1,"R"] # "recoveries" including deaths
  deaths <- rec * fatality_rate
  return(deaths)
}

# Simulate without policy interventions ----------------------------------------

for(i in 1:n_steps){
  Lx <- SS + FF0 * state[i,"S"] # nonlinear population projection matrix
  state[i+1,] <- Lx %*% state[i,]
  deaths[i+1] <- compute_deaths(state[c(0,1)+i,])
}
date_matching_i <- sum(cumsum(deaths) <= date_matching_cumulative_fatality)

if(FALSE){ ## this passage is deactivated, used only for testing and model comparisons:
  ## Simulate without policy interventions
  dates <- date_matching_cumulative_date-date_matching_i+seq(0,n_steps)
  plot(dates,rep(0,n_steps+1),type="n",ylim=c(0,1),xlab="Time",ylab="Proportion",xaxt="n")
  selectedDates <- dates[format(dates, "%d")=="01"]
  axis(1, selectedDates, labels = F , cex.axis = .7)
  text(selectedDates, par("usr")[3]-0.05, 
       srt = 40, adj= 1, xpd = TRUE,
       labels = format(selectedDates, "%d %b '%y"), cex=0.65)
  lines(dates,(state[,"R"]-cumsum(deaths))/pop_size,col="green",type="l")
  lines(dates,rowSums(state[,IStages])/pop_size,col="red",type="l")
  lines(dates,cumsum(deaths)/pop_size*1e1,col="black",type="l")
  legend(x="topright",
         legend = c("Recovered","Infected","cum. mort. x 10"),
         col=c("green","red","black"),lty=c(1,1,1),bg=rgb(1,1,1,alpha = 0.7))
}


# Simulate policy interventions with uncertain outcomes ---------------

dfac <- rep(NA,n_steps+1) # social distancing factor, scales infection rate
dfac[1] <- 1 # Value 1 corresponds to no distancing, the initial assumption
growth_rate <- linear_growth_rate # current growth_rate is used in simulations to inform policy
n_runs <- 200 # How often to repeat simulation
if(!exists("plotting")) plotting <- TRUE # set to FALSE if plotting takes too much time
other_runs_alpha <- 10/n_runs # Intensity of plotting all but last simulation
final_mortality <- rep(NA,n_runs) # Record total mortality at end of each simulation
final_immunization <- rep(NA,n_runs) # Record total immunization rate at end of each simulation
# Prepare plot of results:
date_shift <- date_matching_cumulative_date-date_matching_i
dates <- seq(0,n_steps) + date_shift
plot(dates,rep(0,n_steps+1),type="n",ylim=c(0,1),xlab="Time",ylab="Proportion",xaxt="n")
# Draw x axis:
selectedDates <- dates[format(dates, "%d")=="01"]  # where to draw ticks
axis(1, selectedDates, labels = F , cex.axis = .7)
text(selectedDates, par("usr")[3]-0.05, 
     srt = 40, adj= 1, xpd = TRUE,
     labels = format(selectedDates, "%d %b '%y"), cex=0.65)


# run n_runs simulations
for(run in 1:n_runs){
  policy_intervention_times <- c() # record dates when inverventions occured
  for(i in 1:n_steps){
    # Simulate variation of compliance through time
    dfac[i+1] <- 1-exp(min(0,(log(1-dfac[i])-(1-dfac[i])*rnorm(1,compliance_decline,compliance_sd)))) 
    # Simulate policy response and interventions affecting social distancing
    if(runif(1)< policy_intervention_rate){# adjust policy:
      if(state[i,"Is2"]/pop_size > policy_intervention_threshold && 
         growth_rate > 0){
        # more measures
        dfac[i+1] <- dfac[i+1] * runif(1,0.5,1) # outcome is uncertain
        policy_intervention_times <-
          append(policy_intervention_times,i)
      }else if(dfac[i+1] < 1 && 
               state[i,"Is2"]/pop_size < policy_intervention_threshold / 10 && 
               growth_rate < policy_relaxation_growth_rate_threshold){ 
        # less measures
        dfac[i+1] <- min(1,dfac[i+1]*runif(1,1,1.5)) # outcome uncertain
        policy_intervention_times <-
          append(policy_intervention_times,i)
      }
    }
    Lx <- SS + FF0 * state[i,"S"] * dfac[i+1]
    state[i+1,] <- Lx %*% state[i,]
    deaths[i+1] <- compute_deaths(state[c(0,1)+i,])
    growth_rate <- log(sum(state[i+1,"Is2"])/sum(state[i,"Is2"]))
  }
  if(plotting || run == n_runs){
    # Plot changes through time
    line_alpha <- ifelse(run > n_runs-1,1,other_runs_alpha)*255
    plotcol <- do.call(rgb,as.list(c(col2rgb("green"),alpha=line_alpha,max = 255)))
    lines(dates,(state[,"R"]-cumsum(deaths))/pop_size,col=plotcol,type="l")
    plotcol <- do.call(rgb,as.list(c(col2rgb("blue"),alpha=line_alpha,max = 255)))
    lines(dates,1-dfac,col=plotcol,type="l")
    plotcol <- do.call(rgb,as.list(c(col2rgb("red"),alpha=line_alpha,max = 255)))
    lines(dates,10*rowSums(state[,IStages])/pop_size,col=plotcol,type="l")
    plotcol <- do.call(rgb,as.list(c(col2rgb("black"),alpha=line_alpha,max = 255)))
    lines(dates,10*cumsum(deaths)/pop_size,col=plotcol,type="l")
  }
  # Record main outcomes at end of simulation:
  final_mortality[run] <- sum(deaths)/pop_size
  final_immunization[run] <- 
    state[n_steps,"R"]/pop_size
}
for(t in policy_intervention_times + date_shift){
  lines(c(t,t),c(0,1), col="black", lty=3)
}
legend(x="topright",
       legend = c("Recovered","Distancing","Infected x 10","cum. mort. x 10","policy change"),
       col=c("green","blue","red","black","black"),lty=c(1,1,1,1,3),bg=rgb(1,1,1,alpha = 0.7))
# plot(final_immunization,final_mortality,cex=min(1,200/n_runs))
# hist(final_mortality*pop_size)
median(final_mortality)

