########################################################################################################
########################################################################################################
#
# An illustration of the stacked survival models using the COVID of south Korea dataset. 
########################################################################################################
########################################################################################################


# 'root' is the directory where the files are stored.  This is computer specific
root <- "/Users/Lenovo/Downloads/covid/"



library(MASS)
library(survival)
library(randomForestSRC)


source(paste(root, "coxSurvEst_function.R", sep = ""))
source(paste(root, "ParametricEst_function.R", sep = ""))
source(paste(root, "stacking_code.R", sep = ""))


#############################################################################################
#############################################################################################
#############################################################################################

# this section imports the data and makes appropriate changes

covid_data <- read.csv(paste(root, "patient.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
dat1 <- covid_data



dat1$time  <- dat1$time
dat1$Event <- as.numeric(dat1$status)

dat1 <- dat1[order(dat1$time), ]	# ordering of the data set by observed time


potVars <- c("age", "sex")


dat2 <- dat1[complete.cases(dat1[, c("time", "Event", potVars)]), ]

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# This section estimates the conditional survival functions for each model.


# Setting the seed so the out-of-bag estimates are reproducible
set.seed(1234)


# covariate matrix
potX <- as.matrix(dat2[, potVars])



# fitting parametric and semi-parametric models
mod.ln <- survreg(Surv(time, Event) ~ age +sex, 
				data = dat2, dist = "lognormal")

mod.exp <- survreg(Surv(time, Event) ~ age +sex, 
				data = dat2, dist = "exponential")

mod.weib <- survreg(Surv(time, Event) ~ age +sex, 
				data = dat2, dist = "weib")

mod.cox <- coxph(Surv(time, Event) ~ age + sex, 
				data = dat2)



# Estimating the survival matricies for each patient for the parametric models: 
logNormSurv <- sapply(as.vector(t(mod.ln$coefficients) %*% 
		t(cbind(1,dat2[,dimnames(attr(mod.ln$terms, "factors"))[[2]]]))),
		function(x){ plnorm(dat2$time[dat2$Event == 1], meanlog = x, sdlog = mod.ln$scale, lower.tail = FALSE) } )
	



# Estimating the survival matricies for each patient for the semi-parametric Cox PH model requires a baseline
#	hazard estimate first.  'survfit' will do this automatically with the 'newdata' option.  Otherwise, care
#	should be taken as the 'survival' package centers covariates
survAll.X<- survfit(mod.cox, newdata = dat2[, dimnames(attr(mod.cox$terms, "factors"))[[2]]])
coxSurvAll <- matrix(unlist(sapply({1:length(survAll.X$time)}[survAll.X$n.event > 0], function(x){
				if(survAll.X$n.event[x] == 1){ tmp1 <- survAll.X$surv[x, ] }
				if(survAll.X$n.event[x] >  1){ tmp1 <- matrix(rep(survAll.X$surv[x, ], 2), ncol = 2) }
				return(tmp1) })), ncol = sum(dat2$Event))


# previous code required the construction of a list in this format; should only need the 'Surv_predict' part now
logNormAll  <- list(Surv_predict = t(logNormSurv), timeInterest = dat2$time[dat2$Event == 1], 
		cens = dat2$Event, time = dat2$time, predictors = dat2[, dimnames(attr(mod.ln$terms, "factors"))[[2]]])	

coxModAll <- list(Surv_predict = coxSurvAll, cens = dat2$Event, time = dat2$time, 
				predictors = dat2[, dimnames(attr(mod.cox$terms, "factors"))[[2]]])

	
	

# Now I'm going to fit the out-of-bag parametric and semi-parametric models.  The functions are in the sourced
#	files.
coxMod  <- cox.oob(dat2$time, dat2[, dimnames(attr(mod.cox$terms, "factors"))[[2]]], dat2$Event)
logNMod <- lognormal.oob(dat2$time, dat2[, dimnames(attr(mod.ln$terms, "factors"))[[2]]], dat2$Event)


#######################################################################################################################
#######################################################################################################################

# Time to fit the non-parametric random survival forests.  This requires somewhat more care than the 
#	parmametric and semi-parametric models.  Also note that 'rsf' estimates the out-of-bag estimate
#	at the same time as the estimate with all of the observations.


mod.rsf <- rfsrc(Surv(time, Event) ~ age + sex, 
	data = dat2, nodesize = 64)



# getting the appropriate setup (i.e., dimensions) for the 'rsf' matrix with all the data
rsfSurvAll <- sapply(dat2$time[dat2$Event == 1], function(x){ exp(-mod.rsf$chf[, mod.rsf$time.interest == x]) })



rsfAll <- list(Surv_predict = rsfSurvAll, timeInterest = dat2$time[dat2$Event == 1], 
		cens = dat2$Event, time = dat2$time, predictors = potX)
		


# getting the appropriate setup (i.e., dimensions) for the 'rsf' matrix with out-of-bag data
rsfSurv <- sapply(dat2$time[dat2$Event == 1], function(x){ exp(-mod.rsf$chf.oob[, mod.rsf$time.interest == x]) })
	
rsfPred.oob <- list(Surv_predict = rsfSurv, timeInterest = dat2$time[dat2$Event == 1], 
		cens = dat2$Event, time = dat2$time, predictors = potX)




# We have to make a list for the survival function matricies: 'surv.lst1' is the list for the out-of-bag
#	estimates and 'surv.lst2' is the list for the estimate with all of the observations.  The two lists need
#	to have the same order of survival models.
surv.lst1 <- list(logNMod$Surv_predict, coxMod$Surv_predict, rsfPred.oob$Surv_predict)
surv.lst2 <- list(logNormAll$Surv_predict, coxMod$Surv_predict, rsfAll$Surv_predict)



# 'stacking_ensemble' minimizes the Brier Score; the last argument determines over what time points
stacked.est  <- stacking_ensemble(dat2$time, dat2$Event, potX, surv.lst1, surv.lst2, 
			sapply(1:9/10, function(x){ quantile(dat2$time[dat2$Event == 1], prob = x) }))

# Note that the weight estimates are slightly different from the paper due to differences between the
#	'randomSurvivalForest' package and the 'randomForestSRC' package







