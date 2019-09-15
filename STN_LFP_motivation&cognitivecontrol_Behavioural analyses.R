# The following code has been used to analyze the data presented in the paper:
# Subthalamic nucleus local field potentials recordings reveal subtle effects of promised reward during conflict resolution in Parkinson's disease' (Duprez et al., XXXX)
# This script preprocesses behavioral data and computes RT and accuracy averages, corresponding boxplots and an example of the statistical analyses applied to the data;

# This scripts needs the {ggplot2}, {nlme}, {lme4}, and {car} packages
require(ggplot2) # for the boxplots
require(nlme) # for linear mixed models
require(lme4) # for non-linear mixed models
require(car) # to get the anova table for non-linear mixed models


# Import data
data<-read.table("data2.txt", header=TRUE)
ndata = length(data$rt) # used to compute % of removed data

# Apply behavioral preprocessing
data<-droplevels(data[which(data$rt>200 & data$rt<1500),]) # removes trials with RT<100 & >1500
thresh<-aggregate(data$rt,list(data$n),function(x) mean(x)+3*sd(x)) # creates a RT threshold of mean + 3sd
data<-merge(data,thresh,by.x="n",by.y="Group.1", all=TRUE) # merges the thresh column with data
data<-droplevels(data[data$rt<data$x,]) # removes trials above threshold
pct_outliers<-(((ndata-length(data$rt))/(ndata)))*100 # Computes the % of removed trials

# Create dataframe with  correct data only
dataBR=droplevels(data[which(data$acc!="0"),])


# Average RT and Accuracy -------------------------------------------------

# Mean RT
mean_RT_ncm<-with(dataBR, aggregate(rt, list(dataBR$n, dataBR$cong, dataBR$gain),mean)) # Mean RT according to n and congruence and motivation cue
colnames(mean_RT_ncm)<-c("n", "Congruence","Mot","RT")

# Mean accuracy
mean_ACC_ncm<-with(data, aggregate(acc, list(data$n, data$cong, data$gain),mean)) # Mean ACC according to n and congruence and motivation cue
colnames(mean_ACC_ncm)<-c("n", "Congruence","Mot","ACC")

## Boxplots
# RT

# Rename conditions
mean_RT_ncm$Mot[which(mean_RT_ncm$Mot=="0")]="Fake"
mean_RT_ncm$Mot[which(mean_RT_ncm$Mot=="1")]="Cent"
mean_RT_ncm$Mot[which(mean_RT_ncm$Mot=="100")]="Euro"

mean_RT_ncm$Congruence[which(mean_RT_ncm$Congruence=="1")]="Congruent"
mean_RT_ncm$Congruence[which(mean_RT_ncm$Congruence=="0")]="Incongruent"

mean_RT_ncm$Mot = as.factor(mean_RT_ncm$Mot)
mean_RT_ncm$Congruence = as.factor(mean_RT_ncm$Congruence)

# Relevel so that "Fake" is the first condition
mean_RT_ncm$Mot <- with(mean_RT_ncm, relevel(mean_RT_ncm$Mot, "Fake"))


# 570*408 for size of the image

RTbox <- ggplot(mean_RT_ncm, aes(x=mean_RT_ncm$Congruence, group=mean_RT_ncm$Congruence, y=mean_RT_ncm$RT, color=mean_RT_ncm$Congruence)) + 
  geom_boxplot(position=position_dodge(width=0.9), lwd=2)+
  geom_line(aes(group=n), color = "darkgrey")+
  geom_point(shape = 21, color = "dimgrey", fill = "white", size = 1.5)+
  facet_grid(. ~ Mot, switch = "x")+
  scale_y_continuous(name = "RT (ms)",
                     limits=c(300, 800), breaks = c (300,500,700))
 RTbox + scale_color_manual(values=c("steelblue4", "steelblue2")) + theme_classic() + 
   theme(strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"))+
   theme(legend.position= c(0.85, 0.9), legend.title=element_blank()) +theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=14)) + theme(text = element_text(size=20))

# Accuracy

# Rename conditions so that ggplot understand they're not numbers

mean_ACC_ncm$Mot[which(mean_ACC_ncm$Mot=="0")]="Fake"
mean_ACC_ncm$Mot[which(mean_ACC_ncm$Mot=="1")]="Cent"
mean_ACC_ncm$Mot[which(mean_ACC_ncm$Mot=="100")]="Euro"

mean_ACC_ncm$Congruence[which(mean_ACC_ncm$Congruence=="1")]="Congruent"
mean_ACC_ncm$Congruence[which(mean_ACC_ncm$Congruence=="0")]="Incongruent"

mean_ACC_ncm$Mot = as.factor(mean_ACC_ncm$Mot)
mean_ACC_ncm$Congruence = as.factor(mean_ACC_ncm$Congruence)

# Relevel so that "Fake" is the first condition
mean_ACC_ncm$Mot <- with(mean_ACC_ncm, relevel(mean_ACC_ncm$Mot, "Fake"))


# 570*408 for size of the image

Accbox <- ggplot(mean_ACC_ncm, aes(x=mean_ACC_ncm$Congruence, group=mean_ACC_ncm$Congruence, y=mean_ACC_ncm$ACC, color=mean_ACC_ncm$Congruence)) + 
  geom_boxplot(position=position_dodge(width=0.9), lwd=2)+
  geom_line(aes(group=n), color = "darkgrey")+
  geom_point(shape = 21, color = "dimgrey", fill = "white", size = 1.5)+
  facet_grid(. ~ Mot, switch = "x")+
  scale_y_continuous(name = "Accuracy",
                     limits=c(0.75, 1.05), breaks = c (0.8, 0.9, 1))
Accbox + scale_color_manual(values=c("steelblue4", "steelblue2")) + theme_classic() + 
  theme(strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"))+
  theme(legend.position= c(0.2, 0.1), legend.title=element_blank()) +theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size=14)) + theme(text = element_text(size=20))


# Statistical analyses

# Example on RT and accuracy data 
# the same kind of mixed models were used to test condition differences in power/ITPCz extracted from the different time-frequency windows

data$congruence<-as.factor(data$congruence)
data$motivation<-as.factor(data$motivation)

# RT - linear mixed model

mod_rt<-lme(RT~Congruence*Mot, data=dataBR, random= ~ 1 |n)
anova(mod_rt)

# graphical inspection of the residuals
residus<-residuals(mod_rt)
qqnorm(mod_rt)
qqline(mod_rt)

plot(fitted(mod_rt), residuals(mod_rt),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(mod_rt), residuals(mod_rt)))

# Accuracy - non-linear mixed model

mod_acc<-glmer(ACC~Congruence*Mot+(1 | n), family=binomial(link=logit),data=data, glmerControl(optimizer="bobyqa"))
Anova(mod_acc)