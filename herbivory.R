# ------------------------------------------
# Herbivory
# Ewout Knoester, Veerle Plug, Susan Okoth, Tinka Murk, Ronald Osinga
# Created 03-May-2021
# ------------------------------------------

# Set R and packages
rm(list=ls()) # Clear workspace
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set directory at current directory

library(lme4)
library(lattice) #QQ-plots
library(plyr)
library(tidyverse) # Data manipulation and plotting
library(jtools) #for transformaing model summaries
library(effects) #for plotting parameter effects
library(betareg) # Beta regression
library(lmtest) # likelihood testing
library(glmmTMB) # Beta regression
library(boot) # Inverse logit
library(emmeans) # Pairwise comparisons
library(mvtnorm) # Simulations
library(nlme)
library(rstan)

# Load & organize data
### data info:
### Data collected by Veerle plug
### Grazing (%) calculated as: Weight(out)/(Weight(in)*(1-Weight_control))
### Weight = wet weight of macrophyte (shaken 10 times)
### Weigth_control = grazing on control (caged) structures, calculated for each location x macrophyte combination. On average 19%
### Weight_control was 1 of the 11 structures, randomly chosen. Not included in this dataset
### Grazing > 100 was set to 100 and grazing < 0 was set to 0
### ~1 day interval between measurement Weight(out) and Weigh(in), standardized to 24hrs exactly

Assays <- read.csv("Assays.csv", check.names = FALSE, header = TRUE)
Assays$ID <- factor(Assays$ID)
Assays$Assay <- factor(Assays$Assay)
Assays$Structure <- factor(Assays$Structure)
Assays$Location <- factor(Assays$Location)
Assays$Date <- factor(Assays$Date)
Assays$Species <- factor(Assays$Species)

# Insert new columns for Protection and fill based on study site
Assays <- as.data.frame(append(Assays,list(Protection = ""),after = 4))

Assays$Protection <- ifelse(
  Assays$Location %in% c("Firefly", "Pilli Pipa"),"Fishing",
  ifelse(Assays$Location %in% c("Dolphin Point", "Lower Mpunguti"), "Reserve", "No-take"))

# Transform percentages to fractions
Assays <- as.data.frame(append(Assays,list(Grazing.fraction = ""),after = 9))
Assays$Grazing.fraction <- Assays$Grazing/100

str(Assays) # Check variables
summary(Assays) # Check data

Assays %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather() # Check missing data

# Plotting Grazing across Assay nested within Location and Protection
## Data bounded between 0 and 100: not suited for normal linear models
## Variance seems equal

ggplot() +
  geom_point(aes(x = Assay, y = Grazing), data = Assays) +
  facet_grid(. ~ Protection + Location, switch = "x") +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey50"),
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    panel.grid.minor = element_line(colour = "grey98", size = 0.5),
    strip.background = element_rect(fill = "grey80", colour = "grey50"),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 11)
  ) +
  ylab("Grazing (%)") +
  xlab("Survey")

# Null model
lme.null <- lme(Grazing ~ 1, random = ~ 1 | Location, data = Assays)

# ANOVA taking nesting of Location into account, three different methods:
aov.1 <- aov(Grazing ~ Protection + Error(Location), data = Assays)
summary(aov.1)

lme.1 <- lme(Grazing ~ Protection, random = ~ 1 | Location, data = Assays)
anova(lme.1)

lmm.1 <- lmer(Grazing ~ Protection + (1|Location), data = Assays)
anova(lmm.1)

# Comparing ANOVA versus Null model
## Despite Protection not being significant in ANOVA, it is a significant improvement versus null model
## But can't fulfill assumptions

lrtest(lme.1, lme.null)

# Simple analysis without nested structure, by averaging Structures and get single value per Assay
Assays.avg <- ddply(Assays, ~Assay+Species, summarise, Grazing.mean = mean(Grazing) / 100, Protection = Protection[1], Location = Location[1])

# Transform data, because Beta regression doesn't accept 0s and 1s
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}
Assays.avg$Grazing.scaled <- transform01(Assays.avg$Grazing.mean)

bm.null <- betareg(Grazing.scaled ~ 1, data = Assays.avg)

bm.1 <- betareg(Grazing.scaled ~ Protection, data = Assays.avg) 
summary(bm.1)

# Likelihood testing: bm.1 significantly improves over bm.null model, and both are much better than linear model
lrtest(bm.1, bm.null)
AIC(bm.1, bm.null, lme.1)

# Compare predictions versus observed data
dbeta2 <- function(X, mu, phi, ...) {
  dbeta(X, shape1 = mu * phi, shape2 = (1 - mu) * phi, ...)
}

rbeta2 <- function(N, mu, phi, ...) {
  rbeta(N, shape1 = mu * phi, shape2 = (1 - mu) * phi, ...)
}

# extract coefficients of beta regression model
coefs.bm.1 <- coef(bm.1)

# create vector spanning the transformed 0-1 interval
n.bm.1 <- length(fitted(bm.1))
x.range <- seq(0.5/n.bm.1 , 1-0.5/n.bm.1 , length.out = 200)
x.range.bt <- (x.range*n.bm.1 - 0.5)/(n.bm.1-1)

# control
plot(x.range.bt, dbeta2(x.range, inv.logit(coefs.bm.1["(Intercept)"]), coefs.bm.1["(phi)"]), col = "blue",
     type = "l",  lwd = 2,
     ylab = "Probability density", xlab = "Grazing",
     ylim=c(0,5)
)

# No-take
lines(x.range.bt, dbeta2(x.range, inv.logit(coefs.bm.1["(Intercept)"] + coefs.bm.1[2]), coefs.bm.1["(phi)"]),lwd = 2, col = "red")

# Reserve
lines(x.range.bt, dbeta2(x.range, inv.logit(coefs.bm.1["(Intercept)"] + coefs.bm.1[3]), coefs.bm.1["(phi)"]), col = "orange", lwd = 2)

rug(Assays.avg$Grazing.mean[Assays.avg$Protection=="Fishing"],col="blue", lwd=1.5, pos=5)
rug(Assays.avg$Grazing.mean[Assays.avg$Protection=="Reserve"], col="orange", pos = 4.75, side = 3, lwd=1.5)
rug(Assays.avg$Grazing.mean[Assays.avg$Protection=="No-take"],col="red", pos=4.5, side = 3,lwd=1.5)

legend("topright", lwd = 2, lty = c(1, 1, 1), col = c("blue", "orange", "red"), legend = c("Fishing", "Reserve", "No-take"), bty = "n")

# Check residuals
plot(resid(bm.1) ~ fitted(bm.1))

# Post hoc
test(pairs(emmeans(bm.1, ~ Protection, mode = "link")))

# Beta distribution and nested structure
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}
Assays$Grazing.scaled <- transform01(Assays$Grazing.fraction)

gmm.1 <- glmmTMB(Grazing.scaled ~ Protection + (1 | Location), data = Assays, family = list(family = "beta", link = "logit"))
summary(gmm.1)






# Linear
lmm1 <- lmer(Grazing ~ Protection*Species + (1|Location/Assay/Structure), Assays) 
summary(lmm1)
plot(lmm1)
plot(Assays$Grazing, resid(lmm1))
plot(Assays$Grazing, fitted(lmm1))
qqmath(lmm1)

# Generalized
glmm1 <- glmer(formula = Grazing.fraction ~ Protection * Species + (1|Location/Assay/Structure),
                    family = binomial(link = "logit"),
                    data = Assays)

glmm1 <- glmer(formula = Grazing.fraction ~ Protection * Species + (1 + Protection|ID),
               family = binomial(link = "logit"),
               data = Assays)

glmm1 <- glmer(formula = Grazing.fraction ~ 1 + (1|Location),
                               family = binomial(logit),
                               data = Assays)



summary(glmm1)
summ(glmm1)
plot(glmm1)
plot(allEffects(glmm1))

#TO DO: resolve biased pattern residuals: try fraction optiom GLLM

Assays %>%
  ggplot(aes(x = Location, y = Grazing.fraction, color = as.factor(Species))) +
  geom_point(alpha = .1, position = "jitter")+
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "binomial")) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(0, 1))

Assays %>%
  ggplot(aes(x = Species, y = Grazing.fraction, color = as.factor(Location))) +
  geom_point(alpha = .1, position = "jitter")+
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "quasibinomial")) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(0, 1))



