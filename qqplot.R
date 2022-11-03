library(car)
par(mfrow = c(1, 3), cex = 0.6)

lin.pred <- predict(mod.ln, type = "lp")[dat2$status == 1]
lin.resid <- log(dat2$status[dat1$status == 1]) - lin.pred
weib.pred <- predict(mod.weib, type = "lp")[dat2$status == 1]
weib.resid <- log(dat2$status[dat1$status == 1]) - weib.pred
exp.pred <- predict(mod.exp, type = "lp")[dat2$status == 1]
exp.resid <- log(dat2$status[dat2$status == 1]) - exp.pred



qqPlot(exp(weib.resid), dist = "weibull",shape = 1/mod.ln$scale,
main = "Q-Q plot weibull", xlab = "Theor. quantiles", ylab = "Emp. quantiles")

qqPlot(lin.resid, dist = "norm",sd = mod.ln$scale,main = "Q-Q plot lognormal", 
xlab = "Theor. quantiles (normal)", ylab = "Emp. quantiles")

qqPlot(exp.resid, dist = "exp",rate = mean(abs(exp.resid)),main = "Q-Q plot exponential", 
xlab = "Theor. quantiles (normal)", ylab = "Emp. quantiles")