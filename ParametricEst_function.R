### This will fit parametric models with five-cross validation to obtain out of bag estimates.


library(survival)

lognormal.oob <- function(y, x, delta, Inits = NULL){
	
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

	mod1 <- survreg(Surv(y, delta) ~ X, subset = subgrp1, dist = "lognormal", init = Inits)
	mod2 <- survreg(Surv(y, delta) ~ X, subset = subgrp2, dist = "lognormal", init = Inits)
	mod3 <- survreg(Surv(y, delta) ~ X, subset = subgrp3, dist = "lognormal", init = Inits)
	mod4 <- survreg(Surv(y, delta) ~ X, subset = subgrp4, dist = "lognormal", init = Inits)
	mod5 <- survreg(Surv(y, delta) ~ X, subset = subgrp5, dist = "lognormal", init = Inits)


# This will retrive the baseline survival curves
	surv1 <- sapply(as.vector(t(mod1$coefficients) %*% t(cbind(1, X))), function(x){
				plnorm(y[delta == 1], meanlog = x, sdlog = mod1$scale, lower.tail = FALSE) } )[, d1]
	surv2 <- sapply(as.vector(t(mod2$coefficients) %*% t(cbind(1, X))), function(x){
				plnorm(y[delta == 1], meanlog = x, sdlog = mod2$scale, lower.tail = FALSE) } )[, d2]
	surv3 <- sapply(as.vector(t(mod3$coefficients) %*% t(cbind(1, X))), function(x){
				plnorm(y[delta == 1], meanlog = x, sdlog = mod3$scale, lower.tail = FALSE) } )[, d3]
	surv4 <- sapply(as.vector(t(mod4$coefficients) %*% t(cbind(1, X))), function(x){
				plnorm(y[delta == 1], meanlog = x, sdlog = mod4$scale, lower.tail = FALSE) } )[, d4]
	surv5 <- sapply(as.vector(t(mod5$coefficients) %*% t(cbind(1, X))), function(x){
				plnorm(y[delta == 1], meanlog = x, sdlog = mod5$scale, lower.tail = FALSE) } )[, d5]



# Now I'm going to combine all of the estimates into one matrix,  woohoo!
	fnlSurv <- matrix(0, nrow = length(y), ncol = sum(delta))
	
	fnlSurv[d1, ] <- t(surv1)
	fnlSurv[d2, ] <- t(surv2)
	fnlSurv[d3, ] <- t(surv3)
	fnlSurv[d4, ] <- t(surv4)
	fnlSurv[d5, ] <- t(surv5)
		

	fnlMod <- list(Surv_predict = fnlSurv, timeInterest = y[delta == 1], cens = delta, time = y, predictors = X)

	return(fnlMod)
}







weibull.oob <- function(y, x, delta, Inits = NULL){
	
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

	mod1 <- survreg(Surv(y, delta) ~ X, subset = subgrp1, dist = "weibull", init = Inits)
	mod2 <- survreg(Surv(y, delta) ~ X, subset = subgrp2, dist = "weibull", init = Inits)
	mod3 <- survreg(Surv(y, delta) ~ X, subset = subgrp3, dist = "weibull", init = Inits)
	mod4 <- survreg(Surv(y, delta) ~ X, subset = subgrp4, dist = "weibull", init = Inits)
	mod5 <- survreg(Surv(y, delta) ~ X, subset = subgrp5, dist = "weibull", init = Inits)

	logNormSurv <- 
# This will retrive the baseline survival curves
	surv1 <- sapply(as.vector(t(mod1$coefficients) %*% t(cbind(1, X))), function(x){
				pweibull(y[delta == 1], scale = exp(x), shape = I(1 / mod1$scale), lower.tail = FALSE) } )[, d1]
	surv2 <- sapply(as.vector(t(mod2$coefficients) %*% t(cbind(1, X))), function(x){
				pweibull(y[delta == 1], scale = exp(x), shape = I(1 / mod2$scale), lower.tail = FALSE) } )[, d2]
	surv3 <- sapply(as.vector(t(mod3$coefficients) %*% t(cbind(1, X))), function(x){
				pweibull(y[delta == 1], scale = exp(x), shape = I(1 / mod3$scale), lower.tail = FALSE) } )[, d3]
	surv4 <- sapply(as.vector(t(mod4$coefficients) %*% t(cbind(1, X))), function(x){
				pweibull(y[delta == 1], scale = exp(x), shape = I(1 / mod4$scale), lower.tail = FALSE) } )[, d4]
	surv5 <- sapply(as.vector(t(mod5$coefficients) %*% t(cbind(1, X))), function(x){
				pweibull(y[delta == 1], scale = exp(x), shape = I(1 / mod5$scale), lower.tail = FALSE) } )[, d5]


# Now I'm going to combine all of the estimates into one matrix,  woohoo!
	fnlSurv <- matrix(0, nrow = length(y), ncol = sum(delta))
	
	fnlSurv[d1, ] <- t(surv1)
	fnlSurv[d2, ] <- t(surv2)
	fnlSurv[d3, ] <- t(surv3)
	fnlSurv[d4, ] <- t(surv4)
	fnlSurv[d5, ] <- t(surv5)
		

	fnlMod <- list(Surv_predict = fnlSurv, timeInterest = y[delta == 1], cens = delta, time = y, predictors = X)

	return(fnlMod)
}









