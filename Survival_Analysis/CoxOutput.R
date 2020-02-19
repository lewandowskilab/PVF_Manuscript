##### Multivariate analyses with Cutoff groups #####
# The following script is used for creating the final three coxmodels for the SPP1 and COL6a1 analyses. 
# Should be working if ran all at once, for any questions contact me through Sebastian.
# Last edit: 2020-01-27
# I removed other frailty testing and cohort specific models from this script.
# After every cox figure is finished Schoenfeld test results are printed in the console.
# F.S.Sanders

#### Making the dataframe with cutoffs. #####
setwd("~/Documents/ALS-research/")

allcoxdf <- combinedALS %>% 
  transmute(survobj,
            SPP1 = RawB.53,
            COL6A1 = RawB.123,
            NEFL = RawB.193,
            Gender = gender,
            Onset = factor(onset_bulbar, labels = c("Thorasic/Spinal", "Bulbar")),
            Sampling_Age = sampling_age,
            Sampling_delay = gap,
            cohort)

#Making cutoffs
single.cut <- surv_cutpoint(combinedALS, time = "survival_months", event = "died", variables = c("RawB.53"),minprop = 0.1)
max53 <- single.cut[["RawB.53"]]
cut53 <- max53$estimate
single.cut <- surv_cutpoint(combinedALS, time = "survival_months", event = "died", variables = c("RawB.123"),minprop = 0.1)
max123 <- single.cut[["RawB.123"]]
cut123 <- max123$estimate
single.cut <- surv_cutpoint(combinedALS, time = "survival_months", event = "died", variables = c("RawB.193"),minprop = 0.1)
max193 <- single.cut[["RawB.193"]]
cut193 <- max193$estimate

# Stratifying on earlier establish cutoffs. 
allcoxdf <- allcoxdf %>% mutate(SPP1_Level = case_when(
  allcoxdf$SPP1 <= cut53 ~ "Lower", 
  allcoxdf$SPP1 > cut53 ~ "Upper"))
allcoxdf <- allcoxdf %>% mutate(COL6A1_Level = case_when(
  allcoxdf$COL6A1 <= cut123 ~ "Lower", 
  allcoxdf$COL6A1 > cut123 ~ "Upper"))
allcoxdf <- allcoxdf %>% mutate(NEFL_Level = case_when(
  allcoxdf$NEFL <= cut193 ~ "Lower", 
  allcoxdf$NEFL > cut193 ~ "Upper"))

#### Helper function ####
library(broom)
.get_data <- function(fit, data = NULL, complain = TRUE) {
  if(is.null(data)){
    if (complain)
      warning ("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
    if (is.null(data))
      stop("The `data` argument should be provided either to ggsurvfit or survfit.")
  }
  data
}


##### Continuous Univariate Cox ####

main <-  "Univariate Cox Models: Hazard ratio"
cpositions <- c(0.02, 0.22, 0.4)
fontsize <-  0.7
noDigits <- 2

covariates <- c("SPP1", "COL6A1", "NEFL", "Gender", "Onset", "Sampling_Age")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('survobj ~', x,'+ strata(Sampling_delay) + cluster(cohort)')))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = allcoxdf)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(y){ 
                         coef <- as.data.frame(tidy(y))
                         conf.low <- signif(coef$conf.low, digits = 3)
                         conf.high <- signif(coef$conf.high, digits = 3)
                         x <- summary(y)
                         N <- x$n
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=2)
                         estimate <-signif(x$coef[1], digits=3);#coeficient beta
                         estimate.1 <-signif(x$coef[2], digits=3);#exp(beta)
                         conf.low.1 <- signif(x$conf.int[,"lower .95"], 3)
                         conf.high.1 <- signif(x$conf.int[,"upper .95"],3)
                         ci <- paste0("(", 
                                      conf.low.1, "-", conf.high.1, ")")
                         res<-c(N, p.value, estimate, estimate.1, conf.low,
                                conf.high, conf.low.1, conf.high.1, ci)
                         names(res)<-c( "N","p.value", "estimate", "estimate.1", "conf.low",
                                        "conf.high","conf.low.1", "conf.high.1", "ci")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- tibble::rownames_to_column(as.data.frame(res), "var")

#Correct N numbers
datalist <- list()
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset")
res$N <- as.numeric(as.character(res$N))
for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(allcoxdf[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[4:5] <- n[4:5,]
res$N <- paste0("(N=",res$N,")")
# Change datatypes
res$N <- as.character(res$N)
res$p.value <- as.numeric(as.character(res$p.value))
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(res$estimate.1)
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$conf.low.1 <- as.character(res$conf.low.1)
res$conf.high.1 <- as.character(res$conf.high.1)
res$ci <- as.character(res$ci)
# Add stars p.value to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"

#flip order
res <- res[nrow(res):1, ]

#Setting plot sizes and scales
rangeb <- range(res$conf.low, res$conf.high, na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# increase white space on right for p-vals:
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# y-coordinates for labels:
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

res$var[res$var == "Onset"] <- "Onset Bulbar"
res$var[res$var == "Gender"] <- "Gender Male"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

library(grid)
annot_size_mm <- fontsize *
  as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))

p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Univariate Cox Models. Continuous data."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# switch off clipping for p-vals, bottom annotation:
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
# invisible(p)
pdf(file = "2020_Univ_Cox_Continuous.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()
cat("Continuous Coxmodel finished, proportionality for all univariate models are tested seperately. See below.")
unitest <- coxph(survobj ~ SPP1 + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest)
unitest2 <- coxph(survobj ~ COL6A1+ strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest2)
unitest3 <- coxph(survobj ~ NEFL + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest3)

###### Continuous Multivariate Cox ####

# Final cox model
res.clu <- coxph(survobj ~ SPP1 + COL6A1 + NEFL +
                   Gender + Onset + Sampling_Age + strata(Sampling_delay) +
                   cluster(cohort), data = allcoxdf)

# Extracting all essential stats for forestplot
x <- summary(res.clu)
N <- x$n
coef <- as.data.frame(x$coefficients)
conf <- as.data.frame(x$conf.int)
p.value <- signif(coef["Pr(>|z|)"], digits = 3)
estimate <-signif(coef[1], digits=3);#coeficient beta
estimate.1 <-signif(coef[2], digits=3);#exp(beta)
conf.low <- signif(conf[3], 3)
conf.high <- signif(conf[4],3)
# Binding together and renaming
res<-cbind(N, p.value, estimate, estimate.1, conf.low, conf.high)
res <- tibble::rownames_to_column(as.data.frame(res), "var")
res <- res %>% transmute(var,
                         N,
                         p.value = res[,"Pr(>|z|)"],
                         estimate = coef,
                         estimate.1 = exp(coef),
                         conf.low = res[,"lower .95"],
                         conf.high = res[,"upper .95"],
                         ci = paste0("(", 
                                     conf.low, "-", conf.high, ")"))
names(res)<-c( "var","N","p.value", "estimate", "estimate.1", 
               "conf.low", "conf.high", "ci")
# Correct N for categorical variables and replacing the N column in res
datalist <- list()
covariates <- c("SPP1", "COL6A1", "NEFL", "Gender", "Onset")
data  <- .get_data(res.clu, data = allcoxdf)

for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(data[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[4:5] <- n[4:5,]
res$N <- paste0("(N=",res$N,")")
# Add stars p.value to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"
# Change datatypes
res$N <- as.character(res$N)
res$p.value <- as.character(res$p.value)
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(round(res$estimate.1,2))
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$ci <- as.character(res$ci)
# #flip order
res <- res[nrow(res):1, ]

# Labels and parameters for plotting
main <-  "Hazard Ratio, Continuous Protein Data"
cpositions <- c(0.02, 0.22, 0.4)
fontsize <-  0.7
noDigits <- 2

#Setting plot sizes and scales
rangeb <- range(log(res$conf.low), log(res$conf.high), na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# increase white space on right for p-vals:
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# y-coordinates for labels:
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

res$var[res$var == "OnsetBulbar"] <- "Onset Bulbar"
res$var[res$var == "GenderM"] <- "Gender Male"
res$var[res$var == "SPP1"] <- "SPP1"
res$var[res$var == "COL6A1"] <- "COL6A1"
res$var[res$var == "NEFL"] <- "NEFL"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

#Actual plotting
library(grid)
annot_size_mm <- fontsize *
  as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))

p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Multivariate coxmodels, includes corrections for sampling delay and cohort identity."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# switch off clipping for p-vals, bottom annotation:
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
# invisible(p)
pdf(file = "2020_Multiv_Cox_Cont.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()
cat("Multivariate Cox model with continuous protein data is finished, proportionality test:")
cox.zph(res.clu)
##### Univariate with cutoff #####

main <-  "Univariate Cox Models: Hazard ratio"
cpositions <- c(0.02, 0.22, 0.4)
fontsize <-  0.7
noDigits <- 2


covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset", "Sampling_Age")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('survobj ~', x,'+ strata(Sampling_delay) + cluster(cohort)')))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = allcoxdf)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(y){ 
                         coef <- as.data.frame(tidy(y))
                         conf.low <- signif(coef$conf.low, digits = 3)
                         conf.high <- signif(coef$conf.high, digits = 3)
                         x <- summary(y)
                         N <- x$n
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=2)
                         estimate <-signif(x$coef[1], digits=3);#coeficient beta
                         estimate.1 <-signif(x$coef[2], digits=3);#exp(beta)
                         conf.low.1 <- signif(x$conf.int[,"lower .95"], 3)
                         conf.high.1 <- signif(x$conf.int[,"upper .95"],3)
                         ci <- paste0("(", 
                                      conf.low.1, "-", conf.high.1, ")")
                         res<-c(N, p.value, estimate, estimate.1, conf.low,
                                conf.high, conf.low.1, conf.high.1, ci)
                         names(res)<-c( "N","p.value", "estimate", "estimate.1", "conf.low",
                                        "conf.high","conf.low.1", "conf.high.1", "ci")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- tibble::rownames_to_column(as.data.frame(res), "var")
#Correct N numbers
datalist <- list()
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset")
res$N <- as.numeric(as.character(res$N))
for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(allcoxdf[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[1:5] <- n[1:5,]
res$N <- paste0("(N=",res$N,")")
# Change datatypes
res$p.value <- as.numeric(as.character(res$p.value))
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(res$estimate.1)
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$conf.low.1 <- as.character(res$conf.low.1)
res$conf.high.1 <- as.character(res$conf.high.1)
res$ci <- as.character(res$ci)
# Add stars p.value to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"

#flip order
res <- res[nrow(res):1, ]

#Setting plot sizes and scales
rangeb <- range(res$conf.low, res$conf.high, na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# increase white space on right for p-vals:
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# y-coordinates for labels:
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

res$var[res$var == "Onset"] <- "Onset Bulbar"
res$var[res$var == "Gender"] <- "Gender Male"
res$var[res$var == "SPP1_Level"] <- "High level SPP1"
res$var[res$var == "COL6A1_Level"] <- "High level COL6A1"
res$var[res$var == "NEFL_Level"] <- "High level NEFL"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

library(grid)
annot_size_mm <- fontsize *
  as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))

p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Univariate Cox Models.\nHigh protein level categories based on cutpoints."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# switch off clipping for p-vals, bottom annotation:
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
# invisible(p)
pdf(file = "2020_Univ_Thresh_Cox.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()
cat("Univariate Thresholded Cox model finished. Proportionalities are tested seperately for all variables using cox.zph(). See below.")
unitest <- coxph(survobj ~ SPP1_Level + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest)
unitest2 <- coxph(survobj ~ COL6A1_Level + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest2)
unitest3 <- coxph(survobj ~ NEFL_Level + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest3)

#### Multivariate with cutoff Cox  #####
res.clu <- coxph(survobj ~ SPP1_Level + COL6A1_Level + NEFL_Level +
                   Gender + Onset + Sampling_Age + strata(Sampling_delay) +
                   cluster(cohort), data = allcoxdf)

# Extracting all essential stats for forestplot
x <- summary(res.clu)
N <- x$n
coef <- as.data.frame(x$coefficients)
conf <- as.data.frame(x$conf.int)
p.value <- signif(coef["Pr(>|z|)"], digits = 3)
estimate <-signif(coef[1], digits=3);#coeficient beta
estimate.1 <-signif(coef[2], digits=3);#exp(beta)
conf.low <- signif(conf[3], 3)
conf.high <- signif(conf[4],3)
# Binding together and renaming
res<-cbind(N, p.value, estimate, estimate.1, conf.low, conf.high)
res <- tibble::rownames_to_column(as.data.frame(res), "var")
res <- res %>% transmute(var,
                         N,
                         p.value = res[,"Pr(>|z|)"],
                         estimate = coef,
                         estimate.1 = exp(coef),
                         conf.low = res[,"lower .95"],
                         conf.high = res[,"upper .95"],
                         ci = paste0("(", 
                                     conf.low, "-", conf.high, ")"))
names(res)<-c( "var","N","p.value", "estimate", "estimate.1", 
               "conf.low", "conf.high", "ci")
# Correct N for categorical variables and replacing the N column in res
datalist <- list()
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset")
data  <- .get_data(res.clu, data = allcoxdf)

for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(data[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[1:5] <- n[1:5,]
res$N <- paste0("(N=",res$N,")")
# Add stars p.value to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"
# Change datatypes
res$N <- as.character(res$N)
res$p.value <- as.character(res$p.value)
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(round(res$estimate.1,2))
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$ci <- as.character(res$ci)
# #flip order
res <- res[nrow(res):1, ]

# Labels and parameters for plotting
main <-  "Hazard ratio"
cpositions <- c(0.02, 0.22, 0.4)
fontsize <-  0.7
noDigits <- 2

#Setting plot sizes and scales
rangeb <- range(log(res$conf.low), log(res$conf.high), na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# increase white space on right for p-vals:
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# y-coordinates for labels:
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

res$var[res$var == "OnsetBulbar"] <- "Onset Bulbar"
res$var[res$var == "GenderM"] <- "Gender Male"
res$var[res$var == "SPP1_LevelUpper"] <- "High level SPP1"
res$var[res$var == "COL6A1_LevelUpper"] <- "High level COL6A1"
res$var[res$var == "NEFL_LevelUpper"] <- "High level NEFL"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

#Actual plotting
library(grid)
annot_size_mm <- fontsize *
  as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))

p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Multivariate coxmodels, includes corrections for sampling delay and cohort identity."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# switch off clipping for p-vals, bottom annotation:
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
# invisible(p)
pdf(file = "2020_Multiv_Cox_Thresholds.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()
cat("Categorical Multivariate Coxmodel done finished, cox.zph results:")
cox.zph(res.clu)
