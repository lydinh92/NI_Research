#install.packages('ergm.ego')
#library(ergm.ego)
#install.packages('devtools',repos = "http://cran.us.r-project.org")
library(devtools)
#devtools::install_github("timonelmer/netglm")
library(netglm)
#install.packages('lme4',repos = "http://cran.us.r-project.org")
library(lme4)
#install.packages("fitdistrplus",repos = "http://cran.us.r-project.org") ## check distribution
library(fitdistrplus) 
#install.packages("goft",repos = "http://cran.us.r-project.org")
library(goft)
#install.packages("car",repos = "http://cran.us.r-project.org")
library(car)
#install.packages("actuar",repos = "http://cran.us.r-project.org")
library(actuar)
#install.packages("DHARMa")
library(DHARMa)
#install.packages("sjPlot")
library(sjPlot)
#install.packages("glmmTMB")
library(glmmTMB)
#install.packages("MuMIn")
require(MuMIn)
#install.packages("performance")
library(performance)
#install.packages("partR2")
library(partR2) # for part R2 values
#install.packages("ggpubr")
library(ggpubr)


## LINEAR MIXED EFFECT MODELS for NI data 
# Ly Dinh

set.seed(5)
setwd("/Users/lydinh/Documents/USF/NI")

NI = read.csv("NI_glmer_data.csv")
##rename columns
colnames(NI)[4] ="ni"
colnames(NI)[5] ="time"
#NI

##set variables
#y= c(NI$before_papers, NI$after_papers)
papers = NI$papers
coauthors = NI$coauthors
time = NI$time
target = NI$target
NI = NI$ni

################# 1. Check distribution of data ###########
## Shapiro-Wilks test 
mvshapiro_test(papers)  ## p-value<0.05 = NOT normal
mvshapiro_test(coauthors) ## p-value<0.05 = NOT normal
hist(papers, xlab="Richness of papers of a given ego", ylab="Frequency of egos with a given paper richness")
hist(coauthors,xlab="Richness of co-authors of a given ego", ylab="Frequency of egos with a given co-author richness")

##########GLMER ######### (for non-normal continuous DV)

#Papers model
## estimated using poisson distribution
glmer_papers_poisson = glmer(papers ~ NI * time + (1 | target), family="poisson"(link='log'))
summary(glmer_papers_poisson)

### Residual Diagnostics
# 1. Examine Residual vs Fitted values
plot(glmer_papers_poisson)

# 2. Examine Residual vs. predictor variables
plot(papers,resid(glmer_papers_poisson))

# 3. QQ Plot: quantile-quantile plot of the standardized residuals against the expected quantiles from the assumed distribution
qqnorm(residuals(glmer_papers_poisson)) ## no issues with normality of the weighted residuals

#5. Random effects plot
plot_model(glmer_papers_poisson,type = "re",facet.grid = FALSE) 
## some variability of papers here in the "1" of NI & time

###6. R2 of model (proportion of variation in DV explained by IVs [fixed effects])
r.squaredGLMM(glmer_papers_poisson) ## compare with other distributions
model_performance(glmer_papers_poisson)
#Marginal R2 is concerned with variance explained by fixed factors
#Conditional R2 is concerned with variance explained by both fixed and random factors

## odds ratio 
cc = confint(glmer_papers_poisson,parm="beta_")
ctab = cbind(est=fixef(glmer_papers_poisson),cc)
rtab = exp(ctab)
print(rtab,digits=3)
## get confidence interval
cc

### PLOT MODEL
set_theme(base = theme_classic(), #To remove the background color & grids
          #   theme.font = 'arial',   #To change the font type
          axis.title.size = 2.0,  #To change axis title size
          axis.textsize.x = 1.3,  #To change x axis text size
          axis.textsize.y = 1.3)  #To change y axis text size

model = plot_model(
  glmer_papers_poisson,
  colors = "Paired",  ## change colors here
  show.values = TRUE,
  value.offset = .45,
  value.size = 6.5,
  dot.size = 10,
  line.size = 5,
  vline.color = "red",
  width = 0.7,
  title="",
  # title ="Regression Estimates for GLMER_Papers",
  axis.title = c("Incidence Rate Ratio","hi", "hi")) + scale_x_discrete(labels=c("Membership x Time","Time","NI Membership")) + scale_y_continuous(limits = c(0.3, 5.5), breaks = c(0.5,1.0,1.5,2.0,2.5,3,3.5,4,4.5,5,5.5))
model

#Co-Authors model
## estimated using poisson distribution
glmer_coauthors_poisson = glmer(coauthors ~ NI * time + (1 | target), family="poisson"(link='log'))
summary(glmer_coauthors_poisson)

### Residual Diagnostics
# 1. Examine Residual vs Fitted values
plot(glmer_coauthors_poisson)

# 2. Examine Residual vs. predictor variables
plot(papers,resid(glmer_coauthors_poisson))

# 3. QQ Plot: quantile-quantile plot of the standardized residuals against the expected quantiles from the assumed distribution
qqnorm(residuals(glmer_papers_poisson)) ## no issues with normality of the weighted residuals

#5. Random effects plot
plot_model(glmer_coauthors_poisson,type = "re",facet.grid = FALSE) 
## some variability of papers here in the "1" of NI & time

###6. R2 of model (proportion of variation in DV explained by IVs [fixed effects])
r.squaredGLMM(glmer_coauthors_poisson) ## compare with other distributions
model_performance(glmer_coauthors_poisson)
#Marginal R2 is concerned with variance explained by fixed factors
#Conditional R2 is concerned with variance explained by both fixed and random factors

## odds ratio 
cc = confint(glmer_coauthors_poisson,parm="beta_")
ctab = cbind(est=fixef(glmer_coauthors_poisson),cc)
rtab = exp(ctab)
print(rtab,digits=3)
## get confidence interval
cc

### PLOT MODEL
set_theme(base = theme_classic(), #To remove the background color & grids
          #   theme.font = 'arial',   #To change the font type
          axis.title.size = 2.0,  #To change axis title size
          axis.textsize.x = 1.3,  #To change x axis text size
          axis.textsize.y = 1.3)  #To change y axis text size

model = plot_model(
  glmer_coauthors_poisson,
  colors = "Paired",  ## change colors here
  show.values = TRUE,
  value.offset = .45,
  value.size = 6.5,
  dot.size = 10,
  line.size = 5,
  vline.color = "red",
  width = 0.7,
  title="",
  # title ="Regression Estimates for GLMER_Papers",
  axis.title = c("Incidence Rate Ratio","hi", "hi")) + scale_x_discrete(labels=c("Membership x Time","Time","NI Membership")) + scale_y_continuous(limits = c(0.3, 5.5), breaks = c(0.5,1.0,1.5,2.0,2.5,3,3.5,4,4.5,5,5.5))
model