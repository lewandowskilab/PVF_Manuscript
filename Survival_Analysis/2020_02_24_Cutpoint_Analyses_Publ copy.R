# Cutpoint analyses script for combined cohorts
# For publication
# Last Edit: 2020-02-24, F.S.Sanders

library(ggplot2)
library(gridExtra)
library(survminer)
library(ggbeeswarm)
library(survival)

setwd("~/Documents/ALS-research/")
# Run helper function created within the package
.dichotomize <- function(x, cutpoint, labels = c("low", "high")){
  grps <- x
  grps[x <= cutpoint] = labels[1]
  grps[x > cutpoint] = labels[2]
  
  grps
}

Cutpoint_analyses <- function(df1,df2,protein){
  # Cutpoint and density
  print(protein)
  single.cut <- surv_cutpoint(df1, time = "survival_months", event = "died", variables = protein ,minprop = 0.1)
  max_stat <- single.cut[[protein]]
  cutpoint <- as.numeric(max_stat$estimate)
  print(cutpoint)
  # Create the data frame to plot
  p_data <- data.frame(
    stats = max_stat$stats,
    cuts = max_stat$cuts,
    grps = .dichotomize(max_stat$cuts, cutpoint)
  )
  # variables for plotting
  vline_df <- data.frame(x1 = cutpoint, x2 = cutpoint,
                         y1 = 0, y2 = max(p_data$stats))
  down <- floor(min(df2[[protein]]))
  top <- ceiling(max(df2[[protein]]))
  limits <- c(down,top)
  posx <- top-0.3
  cutpoint_label <- paste0("Cutpoint: ", round(cutpoint,2))
  x1 <- y1 <- x2 <- y2 <- NULL
  posmidy <- 1
  
  # Plotting
  m <- ggplot(data = p_data, mapping=aes_string("cuts", "stats"))+
    geom_point(aes_string(color = "grps"), shape = 19, size = 0.5)+
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                 data = vline_df, linetype = "dashed", size = 0.5)+
    ggplot2::annotate("text", x = posx, y=posmidy,
                      label = cutpoint_label, size = 3)+
    labs(y = "Standardized Log-Rank Statistic", x = paste0(" MFI"))+
    scale_x_continuous(limits = limits)+ coord_flip()
  ## Apply cutpoints
  controls <- subset(df2, df2$class == "Control" | df2$class == "Neuro control")
  controls$class <- "Controls"
  controls$level <- ifelse(controls[[protein]] <= cutpoint, "Low", "High")
  df1$level <- ifelse(df1[[protein]] <= cutpoint, "Low", "High")
  df1$class <- "ALS"
  
  #Plotting Patients
  a <- ggboxplot(df1, x= "class", y= protein, ylab = paste0(protein," MFI"), xlab = "ALS Patients", outlier.shape = " ") +
    geom_hline(yintercept = cutpoint, linetype = 2)+
    geom_beeswarm(aes(color = df1[[protein]] < cutpoint), show.legend = FALSE)+
    scale_y_continuous(limits = limits)
  # Plotting controls
  b <- ggboxplot(controls, x= "cohort", y= protein, ylab = paste0(protein," MFI"), xlab = "Controls", outlier.shape = " ") +
    geom_hline(yintercept = cutpoint, linetype = 2)+
    geom_beeswarm(aes(color = controls[[protein]] < cutpoint), show.legend = FALSE)+
    scale_y_continuous(limits = limits)
  # Density
  e <- ggplot(data = df1, aes(x=df1[[protein]], fill=level))+
    geom_density( alpha=0.3)+
    geom_vline(xintercept = cutpoint, linetype = 2)+
    scale_x_continuous(limits = limits)+
    labs(x= paste0(protein," MFI"), y="Density")+
    coord_flip()
  #Custom legend showing number of censored can be checked with the following:
  namedf <- "Combined"
  t <- table(df1$level)
  high <- paste0("High: ",t[1])
  low <- paste0("Low: ",t[2])
  highlow <- c(high,low)
  fit <- survfit(Surv(survival_months, died) ~ level, df1)
  km <- ggsurvplot(fit, data = df1, palette = c("#ff3333","#0099ff" ), conf.int = FALSE,legend.title=paste0(protein," KM ",namedf),
                   legend=c(0.8,0.8), pval = TRUE, legend.labs = highlow, pval.size = 4)
  kmp <- km$plot
  risk <- km$table
  lay <- rbind(c(1,2,3,4,5,5),
               c(1,2,3,4,5,5))
  # c(NA,NA,NA,6,6,6))
  plotlist <<- list(b,a,m,e,kmp)
  n <- table(df1$died)
  alive <- n[1]
  dead <- n[2]
  print(paste0("Alive: ",alive,",Dead: ",dead))
  pdf(file = paste0(protein,"_cutpoint__testanalyses_",namedf,".pdf"), width = 16,height = 6)
  gridExtra::grid.arrange(grobs = plotlist, layout_matrix = lay)
  dev.off()
}

Cutpoint_analyses(combinedALS,cohorts, "RawB.53")
Cutpoint_analyses(combinedALS,cohorts, "RawB.123")

