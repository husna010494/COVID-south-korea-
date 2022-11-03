########################################################################################################
########################################################################################################
#
# An illustration of the stacked survival models using the German Breast Cancer study.  See the main
#	paper for a more detailed description of the analysis.
#
########################################################################################################
########################################################################################################


# 'root' is the directory where the files are stored.  This is computor specific
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
coba=sapply({1:length(survAll.X$time)}[survAll.X$n.event > 0], function(x){
				if(survAll.X$n.event[x] == 1){ tmp1 <- survAll.X$surv[x, ] }
				if(survAll.X$n.event[x] >  1){ tmp1 <- matrix(rep(survAll.X$surv[x, ], 2), ncol = 2) }
				return(tmp1) })
coba2=matrix(unlist(coba),ncol = sum(dat2$Event))
#coba2=ldply(coba, fill = TRUE)


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





# the GBCS analysis investigates the effect of covariates on five-year survival.  This selects the column that
#	corresponds to five-year survival

pred5.TI  <- stacked.est$Surv_predict[, sum(stacked.est$timeInterest <= 1826.225)]
pred5.Cox <- coxModAll$Surv_predict[, sum(stacked.est$timeInterest <= 1826.225)]
pred5.lgn <- logNormAll$Surv_predict[, sum(stacked.est$timeInterest <= 1826.225)]
pred5.wbl <- weibullAll$Surv_predict[, sum(stacked.est$timeInterest <= 1826.225)]
pred5.rsf <- rsfAll$Surv_predict[, sum(stacked.est$timeInterest <= 1826.225)]







library(mgcv)	# for fitting b-splines


pred.dat <- data.frame(prob5.TI = pred5.TI, prob5.Cox = pred5.Cox, prob5.lgn = pred5.lgn, 
	prob5.wbl = pred5.wbl, prob5.rsf = pred5.rsf, x1 = stacked.est$predictors[,1], x2 = stacked.est$predictors[,2], 
	x3 = stacked.est$predictors[,3], x4 = stacked.est$predictors[,4], x5 = stacked.est$predictors[,5], 
	x6 = stacked.est$predictors[,6], x7 = stacked.est$predictors[,7], x8 = stacked.est$predictors[,8])



mod.TI  <- gam(prob5.TI ~ s(x1, bs = "cr") + s(x2, bs = "cr") + x3 + s(x4, bs = "cr") + x5 + 
			s(x6, bs = "cr") + s(x7, bs = "cr") + x8, data = pred.dat)
mod.Cox <- gam(prob5.Cox ~ s(x1, bs = "cr") + s(x2, bs = "cr") + x3 + s(x4, bs = "cr") + x5 + 
			s(x6, bs = "cr") + s(x7, bs = "cr") + x8, data = pred.dat)
mod.lgn <- gam(prob5.lgn ~ s(x1, bs = "cr") + s(x2, bs = "cr") + x3 + s(x4, bs = "cr") + x5 + 
			s(x6, bs = "cr") + s(x7, bs = "cr") + x8, data = pred.dat)
mod.wbl <- gam(prob5.wbl ~ s(x1, bs = "cr") + s(x2, bs = "cr") + x3 + s(x4, bs = "cr") + x5 + 
			s(x6, bs = "cr") + s(x7, bs = "cr") + x8, data = pred.dat)
mod.rsf <- gam(prob5.rsf ~ s(x1, bs = "cr") + s(x2, bs = "cr") + x3 + s(x4, bs = "cr") + x5 + 
			s(x6, bs = "cr") + s(x7, bs = "cr") + x8, data = pred.dat)



pred.DATA <- data.frame(x1 = 0.8 * 1:100, x2 = 1.2 * 1:100, x3 = rep(mean(stacked.est$predictors[,3]), 100), 
	x4 = 1:100 / 2, x5 = rep(mean(stacked.est$predictors[,5]), 100), x6 = 24 * 1:100, x7 = 11 * 1:100,
	x8 = rep(mean(stacked.est$predictors[,8]), 100))


pred.DATA.t1 <- data.frame(x1 = rep(median(stacked.est$predictors[,1]), 100), x2 = 1.2 * 1:100, 
	x3 = rep(median(stacked.est$predictors[,3]), 100), x4 = rep(median(stacked.est$predictors[,4]), 100),
	x5 = rep(median(stacked.est$predictors[,5]), 100), x6 = rep(median(stacked.est$predictors[,6]), 100), 
	x7 = rep(median(stacked.est$predictors[,7]), 100), x8 = rep(median(stacked.est$predictors[,8]), 100))


pred.DATA.t2 <- data.frame(x1 = rep(median(stacked.est$predictors[,1]), 100), 
	x2 = rep(median(stacked.est$predictors[,2]), 100), x3 = rep(median(stacked.est$predictors[,3]), 100), 
	x4 = 1:100 / 2, x5 = rep(median(stacked.est$predictors[,5]), 100), 
	x6 = rep(median(stacked.est$predictors[,6]), 100), x7 = rep(median(stacked.est$predictors[,7]), 100),
	x8 = rep(median(stacked.est$predictors[,8]), 100))







TI.pred.5ys1  <- predict(mod.TI, newdata = pred.DATA.t1, type = "link")
Cox.pred.5ys1 <- predict(mod.Cox, newdata = pred.DATA.t1, type = "link")
lgn.pred.5ys1 <- predict(mod.lgn, newdata = pred.DATA.t1, type = "link")
wbl.pred.5ys1 <- predict(mod.wbl, newdata = pred.DATA.t1, type = "link")
rsf.pred.5ys1 <- predict(mod.rsf, newdata = pred.DATA.t1, type = "link")

TI.pred.5ys2  <- predict(mod.TI, newdata = pred.DATA.t2, type = "link")
Cox.pred.5ys2 <- predict(mod.Cox, newdata = pred.DATA.t2, type = "link")
lgn.pred.5ys2 <- predict(mod.lgn, newdata = pred.DATA.t2, type = "link")
wbl.pred.5ys2 <- predict(mod.wbl, newdata = pred.DATA.t2, type = "link")
rsf.pred.5ys2 <- predict(mod.rsf, newdata = pred.DATA.t2, type = "link")






par(mfrow = c(1, 2))


plot(TI.pred.5ys1 ~ I(1.2 * 1:100), type = "l", ylim = c(0, 1), xlab = "Tumor Size (mm)", 
	ylab = "Five-Year Survival")
lines(Cox.pred.5ys1 ~ I(1.2 * 1:100), col = 1, lty = 3)
lines(wbl.pred.5ys1 ~ I(1.2 * 1:100), col = 1, lty = 2)
lines(rsf.pred.5ys1 ~ I(1.2 * 1:100), col = 1, lty = 4)
lines(lgn.pred.5ys1 ~ I(1.2 * 1:100), col = 1, lty = 5)
rug(jitter(stacked.est$predictors[,2]), tick = 0.03)
	
legend(x = 2, y = 0.475,c("RSF", "Stacked", "log-Normal", "Weibull", "Cox"), 
	lty = c(4, 1, 5, 2, 3), col = rep(1,5), cex = 0.65)
	
	
plot(TI.pred.5ys2 ~ I(0.5 * 1:100), type = "l", ylim = c(0, 1), xlab = "Number of Nodes", 
	ylab = "Five-Year Survival")
lines(Cox.pred.5ys2 ~ I(0.5 * 1:100), col = 1, lty = 3)
lines(wbl.pred.5ys2 ~ I(0.5 * 1:100), col = 1, lty = 2)
lines(rsf.pred.5ys2 ~ I(0.5 * 1:100), col = 1, lty = 4)
lines(lgn.pred.5ys2 ~ I(0.5 * 1:100), col = 1, lty = 5)
rug(stacked.est $predictors[,4], tick = 0.03)






