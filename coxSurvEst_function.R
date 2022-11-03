### This will fit the Cox model using five fold cross-validation.  That is the data will be split up into five
### almost equal parts, then five Cox models will be fit to estimate the survival curves by not using the data 
### the in the group.



library(survival)

cox.oob<- function(y, x, delta){


	
	first  <- order(y)

	X <- as.matrix(x)[first, , drop = FALSE]
	y <- y[first]
	delta <- delta[first]

	zdat   <- 1:length(y)


#########################################################################################################
#########################################################################################################
	

# The following randomly splits the data into five parts.

# This creates a variable for the remainder (need to know whether the length is divisible by 5)
	aL <- format(round(length(y) / 5, 1), nsmall = 1)
	aR <- as.numeric(substr(aL, nchar(aL), nchar(aL)))


# Note: The 'ind' vector is going to define the sub-sample of an observation
# If the length is not divisible by 5, then we need to add on the length to the end...
	if(aR != 0){
		ind <- c(rep(1:5, floor(length(y) / 5)), 1:{aR / 2})
	}
# If the length is divisible by 5, then everything is easy
	if(aR == 0){
		ind <- rep(1:5, floor(length(y) / 5))
	}

	grps <- sample(ind, length(ind), replace = FALSE)


# Finally, we can create five vectors indicating which observations are considering oob for that
# particular iteration

	d1 <- zdat[grps == 1]
	d2 <- zdat[grps == 2]
	d3 <- zdat[grps == 3]
	d4 <- zdat[grps == 4]
	d5 <- zdat[grps == 5]
	
	
	subgrp1 <- ! zdat %in% d1
	subgrp2 <- ! zdat %in% d2
	subgrp3 <- ! zdat %in% d3
	subgrp4 <- ! zdat %in% d4
	subgrp5 <- ! zdat %in% d5

#########################################################################################################
#########################################################################################################

# Now we can fit the cox models

	CoxMod1 <- coxph(Surv(y, delta) ~ X, subset = subgrp1)
	CoxMod2 <- coxph(Surv(y, delta) ~ X, subset = subgrp2)
	CoxMod3 <- coxph(Surv(y, delta) ~ X, subset = subgrp3)
	CoxMod4 <- coxph(Surv(y, delta) ~ X, subset = subgrp4)
	CoxMod5 <- coxph(Surv(y, delta) ~ X, subset = subgrp5)


# This will retrive the baseline survival curves
	bl.surv1 <- survfit(CoxMod1)
	bl.surv2 <- survfit(CoxMod2)
	bl.surv3 <- survfit(CoxMod3)
	bl.surv4 <- survfit(CoxMod4)
	bl.surv5 <- survfit(CoxMod5)


# These baseline survival curves are not on the original dimensions, unfortunately.  So that needs fixing...
	bl.survf1 <- sapply(y[delta == 1], function(x){ 
						if(mean(bl.surv1$time < x) == 1){ tmp1 <- min(bl.surv1$surv) }
						if(mean(bl.surv1$time < x) != 1){ tmp1 <- bl.surv1$surv[bl.surv1$time >= x][1] }
						return(tmp1) })
	bl.survf2 <- sapply(y[delta == 1], function(x){ 
						if(mean(bl.surv2$time < x) == 1){ tmp1 <- min(bl.surv2$surv) }
						if(mean(bl.surv2$time < x) != 1){ tmp1 <- bl.surv2$surv[bl.surv2$time >= x][1] }
						return(tmp1) })
	bl.survf3 <- sapply(y[delta == 1], function(x){ 
						if(mean(bl.surv3$time < x) == 1){ tmp1 <- min(bl.surv3$surv) }
						if(mean(bl.surv3$time < x) != 1){ tmp1 <- bl.surv3$surv[bl.surv3$time >= x][1] }
						return(tmp1) })
	bl.survf4 <- sapply(y[delta == 1], function(x){ 
						if(mean(bl.surv4$time < x) == 1){ tmp1 <- min(bl.surv4$surv) }
						if(mean(bl.surv4$time < x) != 1){ tmp1 <- bl.surv4$surv[bl.surv4$time >= x][1] }
						return(tmp1) })
	bl.survf5 <- sapply(y[delta == 1], function(x){ 
						if(mean(bl.surv5$time < x) == 1){ tmp1 <- min(bl.surv5$surv) }
						if(mean(bl.surv5$time < x) != 1){ tmp1 <- bl.surv5$surv[bl.surv5$time >= x][1] }
						return(tmp1) })



	coxSurv1 <- t(sapply(1:length(d1), function(x){ bl.survf1 ^ exp(sum((X[d1[x], ] - CoxMod1$means) * 
					CoxMod1$coefficients)) }))
	coxSurv2 <- t(sapply(1:length(d2), function(x){ bl.survf2 ^ exp(sum((X[d2[x], ] - CoxMod2$means) * 
					CoxMod2$coefficients)) }))
	coxSurv3 <- t(sapply(1:length(d3), function(x){ bl.survf3 ^ exp(sum((X[d3[x], ] - CoxMod3$means) * 
					CoxMod3$coefficients)) }))
	coxSurv4 <- t(sapply(1:length(d4), function(x){ bl.survf4 ^ exp(sum((X[d4[x], ] - CoxMod4$means) * 
					CoxMod4$coefficients)) }))
	coxSurv5 <- t(sapply(1:length(d5), function(x){ bl.survf5 ^ exp(sum((X[d5[x], ] - CoxMod5$means) * 
					CoxMod5$coefficients)) }))
	
	
# Now I'm going to combine all of the estimates into one matrix,  woohoo!
	coxSurv <- matrix(0, nrow = length(y), ncol = sum(delta))
	
	coxSurv[d1, ] <- coxSurv1
	coxSurv[d2, ] <- coxSurv2
	coxSurv[d3, ] <- coxSurv3
	coxSurv[d4, ] <- coxSurv4
	coxSurv[d5, ] <- coxSurv5
		

	coxMod <- list(Surv_predict = coxSurv, timeInterest = y[delta == 1], cens = delta, time = y, predictors = X)
}


