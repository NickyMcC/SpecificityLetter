pacman::p_load(tidyverse, cowplot, grid, gridExtra, ggpubr, 
               here, rio, janitor, boot)


setwd("C:/Users/Nicky/Dropbox/standalone bits of work/specificity letter")
data<-read.table("data_community.csv",header=TRUE,sep=",")
data_all<-read.table("data_presumptive_risk.csv",header=TRUE,sep=",")
legend_data<-rbind(data,c("","extra",rep(NA,5)))


cor(data[data$Test!="UltraWithTrace","Prevalence"], data[data$Test!="UltraWithTrace","Specificity"])
cor(data[data$Test=="UltraWithTrace","Prevalence"], data[data$Test=="UltraWithTrace","Specificity"])
cor(data[data$Test!="UltraWithTrace","CalcPrevInTest"], data[data$Test!="UltraWithTrace","Specificity"])
cor(data[data$Test=="UltraWithTrace","CalcPrevInTest"], data[data$Test=="UltraWithTrace","Specificity"])


data_risk<-data_all[data_all$Setting!="Community",]
data_presump<-data_risk[data_risk$Type=="Presumptive",]

cor(data_risk[data_risk$Type=="Presumptive","Prevalence"], data_risk[data_risk$Type=="Presumptive","Specificity"])
cor(data_risk[data_risk$Type=="RiskGroup","Prevalence"], data_risk[data_risk$Type=="RiskGroup","Specificity"])
cor(data_presump[data_presump$Test=="Mtb/Rif","Prevalence"], data_presump[data_presump$Test=="Mtb/Rif","Specificity"])
cor(data_presump[data_presump$Test=="UltraWithTrace","Prevalence"], data_presump[data_presump$Test=="UltraWithTrace","Specificity"])


data_false_pos_rif<-as.data.frame(matrix(NA,nrow=6991,ncol=5))
colnames(data_false_pos_rif)<-c("Prevalence","true_pos","specificity","false_pos","rif_ratio")
data_false_pos_rif[,1]<-seq(10,7000)
data_false_pos_rif$true_pos<-0.618*data_false_pos_rif$Prevalence
data_false_pos_rif$specificity<- (-0.00000387) * data_false_pos_rif$Prevalence + 0.99929255
data_false_pos_rif$false_pos<-(100000 - data_false_pos_rif$Prevalence) * (1 - data_false_pos_rif$specificity)
data_false_pos_rif$ratio<-data_false_pos_rif$false_pos / data_false_pos_rif$true_pos

data_false_pos_trace<-as.data.frame(matrix(NA,nrow=6991,ncol=5))
colnames(data_false_pos_trace)<-c("Prevalence","true_pos","specificity","false_pos","trace_ratio")
data_false_pos_trace[,1]<-seq(10,7000)
data_false_pos_trace$true_pos<-0.69*data_false_pos_trace$Prevalence
data_false_pos_trace$specificity<- (-0.00000664) * data_false_pos_trace$Prevalence + 0.99732916
data_false_pos_trace$false_pos<-(100000 - data_false_pos_trace$Prevalence) * (1 - data_false_pos_trace$specificity)
data_false_pos_trace$ratio<-data_false_pos_trace$false_pos / data_false_pos_trace$true_pos

data_false_pos<-as.data.frame(cbind(data_false_pos_rif$Prevalence,data_false_pos_rif$ratio,data_false_pos_trace$ratio))
colnames(data_false_pos)<-c("Prevalence","rif","trace")


#Bootstrapping the ratio for rif and trace
glimpse(data_false_pos_rif)
glimpse(data_false_pos_trace)

# Function to bootstrap the ratio for a single row
bootstrap_ratio <- function(fp, tp, R = 1000) {
  ratios <- numeric(R)
  for (i in 1:R) {
    fp_sample <- rpois(1, lambda = fp)
    tp_sample <- rpois(1, lambda = tp)
    # Avoid division by zero
    if (tp_sample == 0) {
      ratios[i] <- NA
    } else {
      ratios[i] <- fp_sample / tp_sample
    }
  }
  # Remove NAs and compute percentiles
  ratios <- na.omit(ratios)
  quantile(ratios, probs = c(0.025, 0.975))
}

# Apply the bootstrap function row by row
cis_rif <- t(mapply(bootstrap_ratio, data_false_pos_rif$false_pos, data_false_pos_rif$true_pos))

# Add CI columns to dataframe
data_false_pos_rif$rif_ci_lower <- cis_rif[, 1]
data_false_pos_rif$rif_ci_upper <- cis_rif[, 2]

# Apply the bootstrap function row by row
cis_trace <- t(mapply(bootstrap_ratio, data_false_pos_trace$false_pos, data_false_pos_trace$true_pos))

# Add CI columns to dataframe
data_false_pos_trace$trace_ci_lower <- cis_trace[, 1]
data_false_pos_trace$trace_ci_upper <- cis_trace[, 2]

glimpse(data_false_pos_trace)
glimpse(data_false_pos_rif)

trace <- data_false_pos_trace |> select(Prevalence, trace_ratio, trace_ci_lower, trace_ci_upper)
rif <- data_false_pos_rif |> select(Prevalence, rif_ratio, rif_ci_lower, rif_ci_upper)


df <- trace |> inner_join(rif) |>   
  pivot_longer(
    cols = c(rif_ratio, trace_ratio),
    names_to = "ratio_type",
    values_to = "ratio"
  ) |> 
  mutate(
    ci_lower = case_when(
      ratio_type == "rif_ratio" ~ rif_ci_lower,
      ratio_type == "trace_ratio" ~ trace_ci_lower
    ),
    ci_upper = case_when(
      ratio_type == "rif_ratio" ~ rif_ci_upper,
      ratio_type == "trace_ratio" ~ trace_ci_upper
    )
  )
####################



data$which_plot<-"first"
data$which_plot[data$Test=="UltraWithTrace"]<-"second"
plot_raw_all<-ggplot(data=data) +
  theme_bw() +
  geom_smooth(aes(x=Prevalence,y=Specificity),method = "lm",show.legend = TRUE,colour="grey60") +
  geom_point(aes(x=Prevalence,y=Specificity,colour=Test),size=2.5) +
  geom_errorbar(aes(x=Prevalence,ymin=Min, ymax=Max), width=.2) +
  theme(legend.position="none") +
  ggtitle("a) Estimated Xpert specificity by community TB prevalence, community-wide screening") +
  scale_x_continuous(limits=c(275, 1605),name="Estimated community TB prevalence per 100,000") +
  scale_y_continuous(name="Estimated specificity") +
  coord_cartesian(ylim = c(0.92,1)) +
  scale_color_manual(values=c("blue", "red","green4"))
plot_raw_all<-plot_raw_all + facet_wrap(vars(which_plot),ncol=2, labeller = labeller(which_plot = 
                                                                                       c("first" = "Algorithm 1: Xpert MTB/RIF and Xpert Ultra (trace as negative)",
                                                                                         "second" = "Algorithm 2: Xpert Ultra (trace as positive)")))

dat_text <- data.frame(
  label = c("r = -0.75","r = -0.91"),
  which_plot   = c("first","second"),
  x     = c(600,600),
  y     = c(0.945,0.945)
)
plot_raw_all<-plot_raw_all + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label),size=5
)

plot_adj_all<-ggplot(data=data) +
  theme_bw() +
  geom_smooth(aes(x=CalcPrevInTest,y=Specificity),method = "lm",show.legend = TRUE,colour="grey60") +
  geom_point(aes(x=CalcPrevInTest,y=Specificity,colour=Test),size=2.5) +
  geom_errorbar(aes(x=CalcPrevInTest,ymin=Min, ymax=Max), width=.2) +
  theme(legend.position="none") +
  ggtitle("b) Estimated Xpert specificity by TB prevalence in people screen positive, community-wide screening") +
  scale_x_continuous(limits=c(275, 6810),name="Estimated TB prevalence in people screen positive per 100,000") +
  scale_y_continuous(name="Estimated specificity") +
  coord_cartesian(ylim = c(0.92,1)) +
  scale_color_manual(values=c("blue", "red","green4"))
plot_adj_all<-plot_adj_all + facet_wrap(vars(which_plot),ncol=2, labeller = labeller(which_plot = 
                                                                                       c("first" = "Algorithm 1: Xpert MTB/RIF and Xpert Ultra (trace as negative)",
                                                                                         "second" = "Algorithm 2: Xpert Ultra (trace as positive)")))
dat_text2 <- data.frame(
  label = c("r = -0.96","r = -0.98"),
  which_plot   = c("first","second"),
  x     = c(1500,1500),
  y     = c(0.945,0.945)
)
plot_adj_all<-plot_adj_all + geom_text(
  data    = dat_text2,
  mapping = aes(x = x, y = y, label = label),size=5
)

data_risk$which_plot<-"first"
plot_risk<-ggplot(data=data_risk[data_risk$Type=="RiskGroup",]) +
  theme_bw() +
  geom_smooth(aes(x=Prevalence,y=Specificity),method = "lm",show.legend = TRUE,colour="grey60") +
  geom_point(aes(x=Prevalence,y=Specificity,colour=Test),size=2.5) +
  geom_errorbar(aes(x=Prevalence,ymin=Min, ymax=Max), width=.2) +
  theme(legend.position="none") +
  ggtitle("c) Estimated Xpert specificity by TB prevalence\nin people tested, people in high-risk groups") +
  scale_x_continuous(name="Estimated TB prevalence in people tested per 100,000") +
  scale_y_continuous(name="Estimated specificity") +
  coord_cartesian(ylim = c(0.82,1)) +
  annotate("text", x=9000, y=0.865, label= "r = -0.72", size=5) +
  scale_color_manual(values=c("blue","green4"))
plot_risk<-plot_risk + facet_wrap(vars(which_plot),ncol=1, labeller = labeller(which_plot = 
                                                                                       c("first" = "Algorithm 1: Xpert MTB/RIF")))

rif$which_plot<-"first"
plot_false_pos<-ggplot() +
  theme_bw() +
  geom_ribbon(data=rif,aes(x=Prevalence,ymin=rif_ci_lower,ymax=rif_ci_upper),fill="purple",alpha=0.2) +
  geom_ribbon(data=trace,aes(x=Prevalence,ymin=trace_ci_lower,ymax=trace_ci_upper),fill="green4",alpha=0.2) +
  geom_line(data=data_false_pos,aes(x=Prevalence,y=rif),colour="purple",size=0.7) +
  geom_line(data=data_false_pos,aes(x=Prevalence,y=trace),colour="green4",size=0.7) +
  theme(legend.position="none") +
  ggtitle("d) Predicted ratio of false positive to true positive\n  diagnoses by prevalence in community settings") +
  scale_x_continuous(name="Estimated TB prevalence in people tested per 100,000") +
  scale_y_continuous(name="Estimated ratio of false+ to true+ diagnoses") +
  coord_cartesian(ylim = c(0.5,2))

plot_false_pos<-plot_false_pos + facet_wrap(vars(which_plot),ncol=1, labeller = labeller(which_plot = 
                                                                                 c("first" = "Algorithm 1 and Algorithm 2")))


plot_legend<-ggplot(data=legend_data) +
  theme_bw() +
  geom_point(aes(x=CalcPrevInTest,y=Specificity,colour=Test),size=4) +
  scale_color_manual(values=c("blue", "red", "green4","purple"),
                     labels = c("Xpert MTB/RIF","Xpert Ultra\n(trace as negative)","Xpert Ultra\n(trace as positive)",
                     "Xpert MTB/RIF\nand Xpert Ultra\n(trace as negative)")) +
  theme(legend.title=element_blank()) +
  theme(legend.key.spacing.y = unit(0.6, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))
legend <- cowplot::get_legend(plot_legend)


pdf(paste0("Figure_final.pdf"), width=11,height=12,onefile=FALSE)
ggarrange(
  ggarrange(plot_raw_all,legend,ncol=2,widths=c(2.06,0.3)),
  ggarrange(NULL,ncol=1,widths=c(1)),
  ggarrange(plot_adj_all,NULL,ncol=2,widths=c(2.06,0.3)),
  ggarrange(NULL,ncol=1,widths=c(1)),
  ggarrange(plot_risk,NULL,plot_false_pos,ncol=4,widths=c(1,0.06,1,0.3)),
  nrow=5,heights=c(1,0.08,1,0.08,1))
dev.off()

 