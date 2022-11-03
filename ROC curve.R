ln.predict=predict(mod.ln)
cox.predict=predict(mod.cox,type="lp")
rsf.predict=mod.rsf$predicted
stacked.predict=stacked.est$alphas[1]*ln.predict+
stacked.est$alphas[2]*cox.predict+
stacked.est$alphas[3]*rsf.predict

library(survivalROC)
library(tydiverse)
## Define a helper function to evaluate at various t
survivalROC_helper <- function(t) {
    survivalROC(Stime        = dat2$time,
                status       = dat2$status,
                marker       = ln.predict,
                predict.time = t,
                method       = "NNE",
                span = 0.25 * nrow(ovarian)^(-0.20))
}
## Evaluate every 14 days
survivalROC_data <- data_frame(t = 14 * c(1,2,3,4)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
               as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
## Plot
survivalROC_data %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_point() +
    geom_line() +
    geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
               mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
    facet_wrap( ~ t) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
