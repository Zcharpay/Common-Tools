require(ggplot2)
require(dplyr);require(tidyr)

rats <- read.csv("ratsweights.csv",header = TRUE)
rats.summary <- data.frame(sd=with(rats,tapply(weight,diet,sd)),n=with(rats,tapply(weight,diet,length)),
                           mean=with(rats,tapply(weight,diet,mean)))
group.sizes <- setNames(rats.summary$n,rownames(rats.summary))
group.means <- setNames(rats.summary$mean,rownames(rats.summary))
group.labels <- as.character(unique(rats$diet))

# visualise data. Diet A looks different to the others
ggplot(rats,aes(diet,weight))+geom_boxplot()+ggtitle("Weight gain of rats, by diet")

# ANOVA test. Are any of the group means significantly different?
# REMEMBER: class variable must be a factor, or you will get incorrect output
rats.aov <- aov(weight ~ diet, data=rats)
# Check assumptions: independence. Can't verify this without knowing how the data was collected, go back and check
# Check assumptions: normally distributed residuals. Q-Q plot looks OK, they're approximately normally distributed
par(mfrow=c(1,2));plot(rats.aov,which = c(1:2))
# Check assumptions: equal variance within each treatment (group). Residuals vs fitted plot looks approximately even,
# also check ratio of variances, <2 so looks OK
max(rats.summary$sd)/min(rats.summary$sd)
# could also do Bartlett test to check whether the variances are close enough to equal
# p-value is very large, so do not reject H0 (all variances are equal)
bartlett.test(weight ~ diet, data=rats)
# Assumptions are satisfied! Our ANOVA is valid, so let's review the output
# p-value <0.05, so there is evidence at least one of the means is significantly different from the others.
print(summary(rats.aov))

# Which of the groups is different from the others?
# 1. Is this a balanced study? A: No, number of observations in each group is different
print(group.sizes)
# 2. If balanced study would use TukeyHSD as it is much easier than Bonferroni, use TukeyHSD()
#TukeyHSD(rats.aov)
# 3. Not a balanced study, use Bonferroni
g <- length(group.sizes) # the number of groups (treatments)
k <- choose(g,2)        # the number of pairwise comparisons between groups
a <- 0.05               # alpha, the type 1 error (significance level)
astar = a/k          # alpha star, the bonferroni adjusted significance level
N=length(rats$weight)   # the total number of observations in the dataset
dof=N-g                 # degrees of freedom
# calculate the t-statistic critical value for significance cutoff. Observed t-value (tobs) must be higher than
# this critical value to be considered significant, at significance level = alpha
tquant <- qt(1-(a/(2*k)),df=dof)
# Calculate Bonferroni results. Can see A-C difference in means > critical value, and A-C tobs > tquant (two different ways of reaching same conclusion)
# Group A is significantly different from group C. A-B and B-C aren't significantly different at alpha sig level
SS.resid <- 0
for (grp in group.labels){
        mean.g <- mean((rats %>% filter(diet==grp))$weight)
        SS.resid <- SS.resid + sum(((rats %>% filter(diet==grp))$weight-mean.g)^2)
}
MS.resid <- SS.resid/(N-g); Sp <- sqrt(MS.resid)
rats.bonferroni <- data.frame(t(combn(group.labels,2)),stringsAsFactors = FALSE) %>% rename(g1=X1,g2=X2) %>% mutate(n1=group.sizes[g1],n2=group.sizes[g2],
                                                                                          crit.val=tquant*Sp*sqrt((1/n1)+(1/n2)),
                                                                                          diff.mean1.mean2=abs(group.means[g1]-group.means[g2]),
                                                                                          is.above.crit.val=diff.mean1.mean2>crit.val,
                                                                                          tquant=tquant,tobs=diff.mean1.mean2/(Sp*sqrt((1/n1)+(1/n2))))