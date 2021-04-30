leadership_simul <- function(
# PARAMETERS 
N = 100,									# Size of population
GPSZ = 5,									# Group size
RNDS = 10,									# Number of rounds played within a generation
GENS = 10000, 								# Number of generations
V0 = 1	,									# Baseline fitness
dataheavy = F	,							# Set this to TRUE if you want to save all generations' data, or to FALSE if you just want the end result
B = 1.3,										# Collective production multiplier
candidate_cost = 0.1,						# Cost of volunteering to lead
noncandidate_punishment = 0.1,				# Punishment incurred for not volunteering, when punished
noncandidate_punishment_cost = 0.1,			# Cost of punishing non-volunteers, when punishing

hiL_fraction = 0.5,				        	# The fraction of the population with high leadership ability
hiL = 2,						        	# The leadership ability of high-ability individuals
loL = 1,							        # The leadership ability of low-ability individuals
hiP_fraction = 0.5,					        # The fraction of the population with high leadership ability
hiP = 10,							        # The productivity of high-productivity individuals
loP = 3,							        # The productivity of low-productivity individuals
PLcor = 0,                                  # The correlation between leadership and productivity

leaderlessL = 0.8,							# Effective leadership ability multiplier in groups without leaders

InvF = function (x) exp(x)/4,					# Returns to investment function										
prop_invested = 0.8,     # the percent of individual production that group members give to the group
extraction_coefficient = 0.9,					# If less than 1, resources are devalued when they're extracted
overturn_cost = 0.5,						# The cost of overturning a leader 
overturn_punishment = 2,                  # The extra punishment leaders get when they are overturned
overturn_crit_mass = 0.25,					# The critical mass of the group needed to overturn a leader

mutant_sigma = .05,				    	    # The standard deviation of the differences between mutants and parents
mutation = 0.01,                    # The probability of mutation

# Here create list of values that each trait can take
# For continuous traits this is a range, for categorical a vector of possible options that mutants will pick from randomly
trait_options = list (
  V = c(0,1),      # The probability of volunteering
  I = c(0.4,0.7,1),      # The amount invested by leaders
  E = c(0,20),		 # The amount extracted by leaders
  O = c(0.5,1.5),  # The ratio of leader to groupmember returns at which you vote to overturn leader
  A = c(0.5,1.5),  # The ratio of leader to groupmember returns at which you choose to abdicate
  C = c("L","P")
),

# Set the starting values
StartingVals = c(1,  # everyone volunteers
                 mean(trait_options$I), #everyone invests the mean possible value
                 mean(trait_options$E), #everyone extracts the mean possible value
                 mean(trait_options$O),
                 mean(trait_options$A),
                 "L"),
continuous_traits = c("V","I","E","O","A") # List of traits that vary continuously rather than categorically
) {

## ERROR MESSAGES ###
if("C" %in% continuous_traits) stop ("ERROR: VOTE CHOICE TRAIT CANNOT BE CONTINUOUS")

##################################################
# DEFINE TRAITS AND ENDOWMENTS

# List the evolving decisions variables possessed by individuals here
# V=volunteer I=invest E=extract O=overturn A=abdicate C=choice (high leadership or high prod) P= Punish? R=Reward?
Traits <- c("V","I","E","O","A","C")
# List the endowments (conditions) that are randomly assigned within each generation, 
Endows <- c("L", "P")
EndowProb <- c(hiL_fraction,hiP_fraction) # Specify the probability of having the endowment

##################################################
##################################################
# SIMULATION

##### Prep variables and other data objects

## Prep the trait matrix ts
# The ts matrix gives the trait values for every individual in the population
# What traits does the population start with?
ts <- matrix(StartingVals,ncol=length(Traits),nrow=N,byrow=T)
colnames(ts) <- Traits

ts <- rep(list(as.data.frame(ts)),4)
names(ts) <- c("HiL_HiP","HiL_LoP","LoL_HiP","LoL_LoP")

## Prep data objects for saving all generations' data
report <- c("DATA NOT SAVED FOR ALL GENERATIONS")
if (dataheavy==T) {
  report <- list(
     endows = list(),
     mutants = list(),
     mutant_value = list()
                 )
	# ...
}

## Tidy things up
# Modify N if it not divisible by the groupsize GPSZ
N <- GPSZ * floor(N / GPSZ)
Ngps <- floor(N/GPSZ)

##################################################

### Loop GENS generations
gen <- 1
for (gen in 1:GENS) {

    ## Prep the endowment matrix es
    # The es matrix gives the endowment values for every individual in the population
    # First generate random normal variables with preset correlation
    es <- as.data.frame(MASS::mvrnorm(n=N, mu=c(0, 0), Sigma=matrix(c(1, PLcor, PLcor, 1), nrow=2), empirical=TRUE))
    colnames(es) <- Endows
    # Now set the endowment according to preset probabilities
    for (i in 1:length(Endows)) {es[,i] <- as.numeric(es[,i]> qnorm(EndowProb[i]))} 

    # Create vector of which of the 4 endowment quadrants ppl fall into
    endow_quads <- c()
    endow_quads[es$L==1 & es$P==1] <- "HiL_HiP"
    endow_quads[es$L==1 & es$P==0] <- "HiL_LoP"
    endow_quads[es$L==0 & es$P==1] <- "LoL_HiP"
    endow_quads[es$L==0 & es$P==0] <- "LoL_LoP"


    ## Restrict expressed traits according to endowment
    exp_ts <- as.data.frame(ts[[1]])
    for(i in 1:N){
    exp_ts[i,] <- ts[[endow_quads[i]]][i,]
    }
    
	# gen represent the current generation, from 1 to GENS
	
	# Create the fitness or payoff vector Vs and fill it with zeros
	Vs <- rep(0,N)
	
	#### Create groups
	
	# The gps vector lists each individual's group.
	# Randomly assign each individual to a group.
	gps <- sample(rep(1:Ngps,GPSZ),N, replace=F) # (we calculated Ngps when prepping variables above)
	
	# The gpLdr vector list each group's leader.
	# Each group's leader is listed in order, or is NA if there is no leader. 
	# Start with a null vector of leaders for these groups
	gpLdr <- rep(NA,Ngps)	
	
	#### Repeat RNDS rounds
	for (rnd in 1:RNDS) {
		# rnd represent the current generation, from 1 to RNDS

		########## Election phase ##########
		## Is there a leader already in place from a past round? If not, then do election phase
		
		# Create a list of leaderless groups
		leaderlessgps <- which(is.na(gpLdr)) 
		gp <- 1
		
		# Establish who is volunteering 
		volunteers <- rbinom(N,1,as.numeric(exp_ts$V))
		
		## Repeat for each leaderless groups
		for (gp in leaderlessgps) { # gp represents the current group having an election
			
			# For now, 2 voting strategies: choose randomly between high P or high L folks
			voting_strategy <- sample(names(which.max(table(exp_ts$C[which(gps==gp)]))),1)
			gp_volunteers <- which(gps==gp & volunteers==1)
            suitable_volunteers <- which(gps==gp & volunteers==1 & es[,voting_strategy]==1)
			
			if (length(suitable_volunteers)>0) {
				gpLdr[gp] <- sample(suitable_volunteers,1)		
			}else if(length(gp_volunteers)>0){
			gpLdr[gp] <- sample(gp_volunteers,1)}	
		}
		
		# Now, each group either no leader, with gpLdr=NA, or a leader whose index in the population, somewhere in the range 1:N
		leaderedgps <- which(!is.na(gpLdr))
		leads <- gpLdr[leaderedgps]
		
		# Create the vectors of Ls (leadership endowments) for each group
		k <- which(Endows=="L")
		leadLs <- es[leads,k]
		leadLs[leadLs==1] <- hiL
		leadLs[leadLs==0] <- loL
		grpL <- rep(leaderlessL,Ngps)
		grpL[leaderedgps] <- leadLs
		
		# Create the vectors of Is (leaders' investment) for each group
		k <- which(Traits=="I")
		leadIs <- as.numeric(exp_ts[leads,k])
		grpI <- rep(0,Ngps)
		grpI[leaderedgps] <- leadIs
	
		# Create the vectors of Es (leaders' extraction) for each group
		k <- which(Traits=="E")
		leadEs <- as.numeric(exp_ts[leads,k])
		grpE <- rep(0,Ngps)
		grpE[leaderedgps] <- leadEs

		# Create the vectors of Ps (leaders' productivity) for each group
		k <- which(Endows=="P")
		leadPs <- as.numeric(es[leads,k])
		leadPs[leadPs==1] <- hiP
		leadPs[leadPs==0] <- loP
		grpP <- rep(0,Ngps)
		grpE[leaderedgps] <- leadPs

		########## Production phase ##########
		#grpProd <- GPSZ * B * grpL * grpI 
		#grpMemberV <- grpProd / GPSZ - grpE / (GPSZ - 1)
		#grpLeaderbonus <- grpE + grpE / (GPSZ - 1) # this can't be right?
		
		grpProd <- GPSZ * B * grpL * InvF(grpI*grpP) - grpE
		grpMemberV <- grpProd / GPSZ
		grpLeaderbonus <- grpE*extraction_coefficient - grpI*grpP
		
		# Subtract cost of volunteering
		Vs <- Vs - volunteers * candidate_cost
		
		# Add individual production
		Vs <- Vs + loP + es$P*(hiP-loP)
		
		#Update Vs
		Vs <- Vs + grpMemberV[gps]
		if(length(leads)>0){
		Vs[leads] <- Vs[leads] + grpLeaderbonus[leaderedgps]
		}

		
		########## End-of-round phase ##########

		##### Abdicate?
		
		# Calculate leader to follower return ratio
		groupVs <- aggregate(Vs, list(gps),mean)
		l2f <- Vs[leads]/groupVs$x[leaderedgps]
		
		# For groups with leaders whose return ratio is too small, set gpLdr to NA
		abdicatingdgps <- which(exp_ts$A[leads]>=l2f)
		gpLdr[ abdicatingdgps ] <- NA
		##### Overturn?
		leaderedgps <- which(!is.na(gpLdr))
		gp <- 1
		## Repeat for each leaderless groups
		for (gp in leaderedgps) { # gp represents the current group with a leader
			# right now, we see if there is critical mass of group members to overturn
			overturners <- which(gps==gp & exp_ts$O<=l2f[gp])
			if (length(overturners) >= overturn_crit_mass * GPSZ) {
				gpLdr[gp] <- NA
				Vs[gpLdr[gp]] <- Vs[gpLdr[gp]] - overturn_punishment
				Vs[overturners] <- Vs[overturners] - overturn_cost
			}
		}
	} # End of round loop
	
	#### Reproduce in proportion to payoffs to create the next generation
	# exp_ts has the vectors of traits that are getting passed on, and Vs has the fitness
	Vtemp <- V0 + Vs
	Vfinal <- Vtemp / mean(Vtemp)
	
	parents <- sample(1:N, N, replace=T, prob=Vfinal)
	exp_ts_old <- exp_ts
	exp_ts <- exp_ts_old[parents,]
	colnames(exp_ts) <- Traits
	
	# Mutate trait values
	mutant1s <- matrix(rbinom(n=N*ncol(exp_ts),size=1,prob=mutation),nrow=nrow(exp_ts),ncol=ncol(exp_ts))
	colnames(mutant1s) <- Traits
	
	for(trait in Traits){
	  if (trait %in% continuous_traits){
	    exp_ts[,trait] <- as.numeric(exp_ts[,trait]) + mutant1s[,trait]*rnorm(N,0,mutant_sigma)*(trait_options[[trait]][2]-trait_options[[trait]][1])
	  
	 # Fix traits that go beyond the range  
	    exp_ts[exp_ts[,trait]>trait_options[[trait]][2],trait] <- trait_options[[trait]][2]
	    exp_ts[exp_ts[,trait]<trait_options[[trait]][1],trait] <- trait_options[[trait]][1]
	    
	  }else {
	    exp_ts[mutant1s[,trait]==1,trait] <- sapply(exp_ts[mutant1s[,trait]==1,trait], function(x)sample(setdiff(trait_options[[trait]], x)))
	  } 
	}

	# Redistribute expressed genes to their respective locations
	for(i in 1:N){
	  ts[[endow_quads[i]]][i,] <- exp_ts[i,] 
	}

	########## Record keeping
	if (dataheavy==T) {
	report$endows[[gen]] <- endow_quads
	report$mutants[[gen]] <- which(mutant1s==1)
	report$mutant_value[[gen]] <- as.matrix(exp_ts)[which(mutant1s==1)]
	report$final_rd_l2f[[gen]] <- l2f
	}
	print(gen)
} # End of generation loop
return(list(
  final_traits = ts,
  records = report,
  params = c(N = N, GPSZ = GPSZ, RNDS = RNDS,GENS = GENS, V0=V0,
mutation=mutation,
  B = B,candidate_cost = candidate_cost,noncandidate_punishment = noncandidate_punishment, noncandidate_punishment_cost = noncandidate_punishment_cost,hiL_fraction = hiL_fraction,
hiL = hiL,	loL =	loL, hiP_fraction=hiP_fraction,hiP=hiP,loP=loP,PLcor=PLcor,leaderlessL = leaderlessL,	leaderlessI=leaderlessI,	InvF=InvF,leader_reward=leader_reward,extraction_coefficient = extraction_coefficient,overturn_cost=overturn_cost,
overturn_punishment=overturn_punishment,overturn_crit_mass=overturn_crit_mass,mutant_sigma=mutant_sigma,mutation=mutation),

continuous_traits = continuous_traits,
trait_options = trait_options,
StartingVals=StartingVals
))
}




