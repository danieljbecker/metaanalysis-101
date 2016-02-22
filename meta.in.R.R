## run a meta-analysis in R
## Daniel Becker
## dbecker@uga.edu
## last updated 12/8/2015

## see what I did up there?
## before we start, always annotate your code!

## clear workspace
rm(list=ls()) 
graphics.off()

## load packages
library(metafor)
library(multcomp)
library(plotrix)
library(MuMIn)

## load in data
## we'll use data within the metafor package
## in your real dataset, you can chose to convert test stats to effect size
## either in your excel sheet or within R

## for now, here's one way to convert stats to standardized effect size
## CO2 and tree data (mean differences)
## http://www.jstor.org/stable/4221855?seq=1#page_scan_tab_contents
data1=dat.curtis1998
## more free data is provided with the data() function

## plants grown in ambient or elevated CO2 conditions
## so the effect size will be a standardized difference between treatment & control
## calculation of this requires the mean, standard deviation, and sample size of each
## means are m1i and m2i
## sd are sd1i and sd2i
## sample size is n1i and n2i

## we use the escalc function here
## define the data frame (data=whatever)
## define the appropriate measures (here m1i=m1i, etc)
## tell R what measure to calculate
## here we use standardized mean difference (Hedges' g) using measure="SMD"
## see ?escalc for examples and help
data1=escalc(measure="SMD",data=data1,m1i=m1i,m2i=m2i,sd1i=sd1i,sd2i=sd2i,
       n1i=n1i,n2i=n2i,append=T)
## you're telling R which variables are the means for group 1, group 2, and resp sds

## append=T automatically attaches two variables to your dataset
## yi = effect size (SMD)
## vi = sampling variance 

## another common effect size measure for ecology/evolution is for 2 continuous vars
## i.e., the correlation between degree of urbanization and oh, say, melanin in birds
## studies usually do an OK job reporting correlation coefficients (r) or R2

## as an example, the relationship between conscientiousness and medication adherence
## do people who think more before acting take their medicine?
## correlations between 5 point conscientiousness scale and how often they take meds
## http://link.springer.com/article/10.1007/s12160-013-9524-4
data2=dat.molloy2014

## usually you wouldn't be this lucky and have all your rs reported for you
## have to manually convert different effect sizes into r first
## extract r statistics and sample size
## correlations are in ri, sample size
## we typically convert the correlation coefficient to Fisher's Z (ZCOR)
data2=escalc(measure="ZCOR",data=data2,ri=ri,ni=ni,append=T)

## again, effect size in yi and sampling variance in vi
## we use Fisher's Z because this transformation helps stabilize variance
## its also a good normalizing transformation for correlations

## a quick note on conversions
## Hedges g and Fisher's Z can be transformed from one to the other
## these are pretty good base effect sizes to use in a meta-analysis
## particularly for ecology/evolution/physiology, where study designs arent standardized
## mainly these are flexible for transforming different raw effect sizes
## for example, Chisq stats, F stats, and odds ratios can all be converted to r

## you can write out several functions that convert from X to r
## and then use an ifelse function to automate effect size conversion within r
## for example, let's say you've collected raw effect sizes from your studies
## they consist of correlation coefficients, odds ratios, and X2 statistics
## transforming everything to r is your best bet (and the most interpretable!)

## function for odds ratio conversion
## Pearson 1900 (though some debate on when to use this vs other conversions)
ORc=function(x){
  r=cos(pi/(1+(x^.5)))
  return(r)
}

## function for Chisquare stat conversion
## Rosenberg 2010 PLOS ONE
X2c=function(x,n){
  r=sqrt(x/n)
  return(r)
}

## ideally you will have a column of raw test stats and of their labels
## e.g., a column called stat with the name (OR,X2,R)
## and a column called effect with the values (actual stats) and n (sample size)
## ifelse function to convert everything at once
## odds ratios and risk ratios
## data$r=ifelse(data$stat=="OR",ORc(data$effect),
##              ifelse(data$stat=="X2",X2c(data$effect,data$n),data$effect))
## if you don't understand what's going on with ifelse, it's a really handy tool

## OK so now we have effect sizes...what can we do?
## i would say there are three main classes of models for a basic meta-analysis
## all of these include some form of weighting, where weights = 1/vi 
## i.e., small sampling variance (large sample size ) = large weight
## also assume the sampling variances are known in the form of sample size
## this contrasts with assumptions of standard stats (lms or glms)
## so there's some value to using meta-analysis specific models

## the first is to assess signal for an overall trend (mean effect)
## this is known as a fixed- or random-effects model
## fixed effects are basically = you have all the studies out there
## conditional inference = what is the average effect among my included studies
## random effects = your studies are but a sample of a larger "population"
## unconditional inference = true effect from larger set of studies
## people used fixed effects models a lot in biomedicine, but not as much in EEB

## 1) random effects models
## studies are going to inherently differ in methods, sample characeristics, etc
## this introduces heterogeneity among the effect sizes
## you can treat heterogeneity as just random, effects normally distributed
## then estimated the mean effect (mu) and the variance (tau^2)
## if tau^2 = 0, homogeneity among true effect sizes

## we use the rma() function for this for both our datasets
model1=rma(yi,vi,data=data1)
model2=rma(yi,vi,data=data2)
## models are fit a bunch of different ways
## rma automatically chooses restricted maximum likleihood (REML)
## this is good except for if you want to compare models
## then you'll need to use just maximum likelihood 
model1=rma(yi,vi,data=data1,method="ML")
model2=rma(yi,vi,data=data2,method="ML")
## see ?rma for some other options under method (i.e., DL, HE, etc)

## summary works as with lm or glm or lmer or whatever model fits you use
summary(model1)
summary(model2)

## OK so let's walk through this here with model2 since the dataset is smaller
summary(model2)
## we get first basic fit statistics (AIC, AICc, etc)
## estimate of tau2 (heterogeneity 0 or not)
## I2 and H2 = kinda wonky (ignore H2 for now)
## I2 is like the % variability due to underlying heterogeneity rather than error
## http://psych.colorado.edu/~willcutt/pdfs/Higgins_2002.pdf
## Q statistic = test if tau2 = 0 (more or less)
## estimate is mu, our mean true effect

## so there's lots of underlying heterogeneity (I2 = 58%)
## we reject tau2 = 0 (Q = 38, p < 0.001), so significant heterogeneity
## mean effect is a Z of 1.28 (z = 5, p < 0.0001)
## you could back transform this into r if you wanted for interpretation
## or use never convert r into Z in the first place

## is this different from just calculating the mean?
mean(data2$yi)
model2$b
## by about 0.01, so not greatly

## what about from running the same kind of linear model?
coef(lm(yi~1,data=data2))
## in this case, the mean is not really different at all
## you could argue that the mean and lm overestimate the true effect (mu)
## sometimes this is the case, other times not
## but the standard error and p values will be different
## recall that lm() assumes no weights or that variance is unknown
## http://www.metafor-project.org/doku.php/tips:rma_vs_lm_and_lme

## there are two other ways we can fit the random effects model
## the first is called the leave one out method
## this repeatedly fits the model, leaving out one study at a time
## sort of a sensitivity analysis of the estimate of mu and tau2
leave1out(model2)

## the other is called a cumulative meta-analysis
## here the estimate of mu is obtained sequentially
## studies are added to the model in chronological order
model3=cumul(model2,order=order(data2$year))
model3

## the most common way to showcase the random effect model is with a forest plot
## this will show the effect sizes and the estimated true effect, in order
## really easy in R, just use forest()
par(mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0))
forest(model2)

## this is ugly though: no study label, points are everywhere
## it's common to have the study name on the left
## you can combine the authors and year columns in r
data2$study=paste(data2$authors,data2$year)

## assign study label with slab and order observations by magintude
## notice that points are scaled by the weights (inverse of variance)
forest(model2,slab=data2$study,order="obs",annotate=T,pch=21,bg="black")

## we can also plot the cumulative meta-analysis to see how mu changes over time
forest(model3,pch=21,bg="black")
## the last data point is then the same as the REM estimate above

## so if we go back to the study (association between conscientious and adherence)
## it would seem that across all studies we find a positive mean correlation
## in that people that think more before they act take meds on schedule
## surprise!

## but we also observed significant heterogeneity among the true effects
## that brings us back to basic model 2 in meta-analysis

## 2) moderator analysis (or mixed effect models or other names)
## studying wildlife for example, we know that study conditions are insanely different
## study species, habitat conditions, study design, etc are all over the place
## so we might be less interested in estimating the mean
## and more interested in explaining variation in effect sizes

## lets look back at our forest plot
forest(model2,slab=data2$study,order="obs",annotate=T,pch=21,bg="black")
## most points are to the right (positive correlation)
## but one is negative, any many are on 0 or very small
## are there aspects of the studies that could influence that variation?

## the mixed effect models / moderator analysis is just a multiple regression
## like the random effect model, has sample size weights + known variance

## the model form is the same with rma()
## alter slightly with standard outcome~var+var+var formula in R
## let's look at covariates
names(data2)

## logical predictors might be study design and mean age of the cohort
## could studies with older cohorts have a higher chance of finding a positive r?
## what about a snapshot study versus one with comparison groups?

## how you handle these is probably based on your modeling philosophy
## for now, we have a small sample size, so lets independently look at each

## make sure early that factor variables are factors!
class(data2$design)
data2$design=factor(data2$design)

## mixed effects models
mem1=rma(yi,vi,mods=~design,method="ML",data=data2)
mem2=rma(yi,vi,mods=~meanage,method="ML",data=data2)

## lets look at mem1 with design
summary(mem1)
## similar to basic rma
## tau2, I2
## now have R2, which should be familiar
## this is a pseudoR2, essentially compares moderator model with random effect model
## i.e., including study design explains 10% more variation than intercept only
## can use the anova(full,reduced) formula to compare and get the same thing
anova(mem1,model2)

## we also get Q again, here is the residual heterogeneity (ignore)
## also QM, which is the test of significance for the moderator
## in the model results, we get the estimate, z val, and p val
## notice the p val here is the same as the moderator test
## which makes sense because we only have 1 moderator in the model

## so study design doesn't influence the effect size -- that's probably good
## what about the age of the cohort?
summary(mem2)
## QM = 1.9, p < 0.16, R2 = 35%
## so not a significant result, but tells us a bit more
## r2 of 35% isn't small, especially for social science results

## we can also look into the model for basic diagnostics
## diag=influence(model), then plot(diag)
diag=influence(mem2)
plot(diag)
## points in red illustrate potential outliers
## good to look at DFFITS and Cook's distance
out=diag$inf[which(diag$inf$cook.d>1),]; print(out)

## save the row number of the outlier
out=as.numeric(rownames(out))
## you could remove this and see how the results change 
data3=data2[-out,]

## redo the random effects model
remod=rma(yi,vi,data=data3,method="ML"); summary(remod)
## stronger true effect, smaller tau2 (between study heterogeneity)

## redo the age moderator analysis
remem2=rma(yi,vi,mods=~meanage,method="ML",data=data3); summary(remem2)
## looks like the outlier was heavily influential (p went from 0.15 to 0.8)
## removing it completely alters the marginal sigificant negative slope

## other model diagnostics include plotting the residuals agains the predict values
## or histograms of the residuals to check normality
plot(predict(mem2)$pred,resid(mem2),pch=21,bg="grey",
     xlab="predicted values (effect size ~ mean age)",ylab="model residuals")
abline(h=0,lty=2)
hist(resid(mem2))
shapiro.test(resid(mem2)) ## if p < alpha, reject normality

## can also look at Q-Q normal plots of the REM and MEM
## plot theoretical quantiles of normal dist against observed quantiles of residuals
qqnorm(model2)
qqnorm(mem2)

## ignoring the outlier, we can plot the results of the age model easily in standard R
plot(yi~meanage,data=data2,pch=21,bg="grey",xlab="mean age of cohort",
     ylab="Fisher's Z transformed correlation coefficient",las=1)

## add the model fit with abline(model), just like with univariate lm()
abline(mem2)

## for meta-analysis results, best to include sampling variance in the plot
## i.e., scale the circles by the sample size or inverse of variance (1/vi = weight)
plot(yi~meanage,data=data2,pch=21,bg="grey",xlab="mean age of cohort",
     ylab="Fisher's Z transformed correlation coefficient",las=1,
     cex=rescale(1/data2$vi,c(0.75,3)))
abline(mem2)
## you can see how the trendline most matches the points with largest n or smallest vi
## add fit for basic lm to compare
abline(lm(yi~meanage,data=data2),lty=2)
## so again, in this example the meta-analysis models don't differ that much

## we can already tell from these two models tha age is a better predictor than design
## you could compare with AICc if you wanted that for sure 
AICc(mem2,mem1)

## look at model weights
Weights(AICc(mem2,mem1))
## 30% more support for age vs design

## so again, what does this mean?
## well age of the population explains marginally significant heterogeneity in r
## marginal evidence that studies with an older sample were more likely to find
## positive correlations between conscientiousness and med adherence

## this isn't the most exciting dataset, but what have we done so far?
## quantified the average true effect among the pop of studies (random effects)
## examined drivers of heterogeneity among effects (mixed effects / moderator)
## the last basic technique is to test for publication bias

## 3) methods for testing evidence of publication bias
## publication bias is the publication of one type of result over another
## this usually means more publishing of significant findings in the expected direction
## can arise from intentional selective publishing
## but also can arise from the more subtle "file drawer" effect (Rosenthal 1979)
## non-significant or controversial findings get "stored" away

## one way to display publication bias is through the funnel plot
## the funnel plot shows the relationship between effect size and precision
## unbiased sets of studies should show a triangle formation around the true effect
## studies with high error should be on either side of the average
## with decreasing variation in effect size with more precision
## to show this, use funnel(random effects model name)
funnel(model2); box()

## in this case, effect sizes look distributed
## white bounds are CIs, +/- 1.96 SE
## use addtau2=T to make CI = +/- 1.96*sqrt(SE^2 + tau2)
funnel(model2,addtau2=T); box()

## we can test for bias by correlating precision measures with effect size
## regression test is regtest.rma(model), specifiy predictor (sei, vi, ni, etc)
regtest.rma(model2,model="rma",predictor="sei")

## same with the sample size or the sampling variance
regtest.rma(model2,model="rma",predictor="ni")
regtest.rma(model2,model="rma",predictor="vi")
## find no evidence for publication bias
## this makes sense, as points are spread evenly in our plot

## we can also test publication bias using the moderator analysis method
bias=rma(yi,vi,mods=~ni,method="ML",data=data2)
summary(bias)
## this is the same p value as in the "ni" example above

## the next step is to test if any studies are missing from the analysis
## if there is pub bias from "file drawer" effect, the funnel will be skewed
## this method figures where "missing" studies should go to correct assymetry
## known as the trim and fill analysis (Duval & Tweedie 2000)
## two estimator methods (R0 and LO)
## L0 estimates the # of studies missing
## R0 actually tests if that number is different from 0
## can observe how the random effects model estimate would differ with missing
trimL0=trimfill(model2,estimator="L0"); summary(trimL0)
trimR0=trimfill(model2,estimator="R0"); summary(trimR0)
## again, nothing missing, no evidence for pub bias

## then just use the funnel command to plot the new funnel with missing studies added
funnel(trimR0,addtau2=T)
funnel(trimL0,addtau2=T)
## can see two studies added to correct for assymetry with L0
## but the R0 estimator tells us 2 isnt different from 0
## in general it is good practice to check with R0 and L0

## so those are the three most basic models for meta-analysis
## random effects model, moderator analysis, publication bias tests
## AKA weighted average, weighted regression, impute missing studies
## there's some related ways to test for outliers and the like

## now what are some assumptions in these models?
## well, lots, but the main one is that effect sizes are independent
## much like in normal linear models, this isn't always the case
## especially in ecology/evolution, studies of wildlife, etc
## one way effect sizes can be non-independent arises from study clustering
## e.g., one study could report multiple effect sizes for your analysis
## e.g., an urbanization/melanin study could quantify the effect for several species
## so effect sizes could be more similar within studies
## likewise, effect sizes could be similar per species
## so you could consider random effects of study and species
## plenty of other options
## ignoring random effects could skew estimates of mu & tau2

## luckily the metafor package has a nice way of handling this
## the rma.mv() function, multilevel meta-analysis
## works a lot like rma
## except as rma.mv(yi=yi,V=vi,random=list(~1|thing,~1|thing2),mods=~1,data=data)

## in the data2, let's make a random effect factor (or two)
set.seed(5)
data2$re1=factor(round(runif(nrow(data2),1,4)))
data2$re2=factor(round(runif(nrow(data2),1,7)))

## now we can redo the random effect model and moderator analysis
## in multilevel meta-analysis form
base=rma.mv(yi=yi,V=vi,random=list(~1|re1,~1|re2),method="ML",mods=~1,data=data2)
summary(base)

## how different is the estimate of mu?
summary(model2)$b
summary(base)$b
## different by a little under 0.01
## you can do plot the forestplot again (forest())
## but you can't do the publication bias tests (oh well)
## so if you have data that lends itself to the multilevel approach
## do the random effect and moderator analyses with rma.mv but pub bias test with rma

## one different aspect of multilevel models is you can examine the variance components
summary(base)$sigma
## both are greater than 0, so it may be important to include them in your model
## though you can also argue to include them even if they're 0 on the basis of design
## since effect sizes from the same study aren't independent data

## you can also check the variance components to ensure the MLE is identifiable
par(mfrow=c(1,2),oma=c(0,0,1,0),mar=c(4,4,1,1))
profile(base,sigma2=1)
profile(base,sigma2=2)
## by identifiable, i mean that the loglike is maximized at the MLE estimates
## if the likelihood surface is flat, you have a problem
## i.e., your model is probably overparameterized
## and here sigma 2's estimate is actually a little lower than the max loglik

## next do the moderator analysis with multilevel approach
mlmem1=rma.mv(yi=yi,V=vi,random=list(~1|re1,~1|re2),method="ML",mods=~design,data=data2)
mlmem2=rma.mv(yi=yi,V=vi,random=list(~1|re1,~1|re2),method="ML",mods=~meanage,data=data2)
summary(mlmem1)
summary(mlmem2)

## the design effect is the same, but the age effect is now significant
## so after accounting for clustering of effect size
## cohort age explains substantial variation in effect size
## lets plot again and add model fit with rma and rma.mv
## have confidence intervals for rma.mv

## make new vector of age data based on the range
lim=range(data2$meanage)
lim[1]=lim[1]-5
lim[2]=lim[2]+5
nage=seq(lim[1],lim[2],length=50)

## predict function gives you CIs and prediction intervals
preds=predict(mlmem2,newmods=nage)

## now plot
par(mfrow=c(1,1),mar=c(4.5,4.5,0.5,0.5))
plot(yi~meanage,data=data2,pch=21,bg="grey",xlab="mean age of study cohort",
     ylab="Fisher's Z transformed correlation coefficient",las=1,
     cex=rescale(1/data2$vi,c(0.75,3)))

## add shading for CIs and trend
polygon(c(rev(nage),nage),c(rev(preds$ci.ub),preds$ci.lb),col='grey95',border=NA)
lines(nage,preds$ci.ub)
lines(nage,preds$ci.lb)
abline(mlmem2)

## add points back in again
points(yi~meanage,data=data2,pch=21,bg="grey",cex=rescale(1/data2$vi,c(0.75,3)))
box()

## add rma line
abline(mem2,lty=2)

## add legend
legend("topright",c("rma","rma.mv"),lty=c(2,1),bty="n")

## slope does indeed look stronger
## so after accounting for some clusting in the dataset
## studies with older populations are more likely to find a strong correlation
## between conscientiousness and medical adherence
## put another way, kids dont show much of a relation between adherence & thinking

## recall that with rma, we got a value for R2
summary(mem2)$R2

## but we dont have that built into rma.mv
summary(mlmem2)$R2

## instead we can calculate a pseudo R2
## proportional reduction in the total variance of the MEM compared to REM
## essentially comparing variance in the multilevel random effect to moderator
## variants of this given in Nakagawa & Schielzeth 2013 for models w/ random effects
pr2=function(full,reduced){
  r2=(sum(reduced$sigma2) - sum(full$sigma2)) / sum(reduced$sigma2)
  print(r2)}
pr2(mlmem2,base)
## so after accounting for two random effects in the data
## age explains up to 70% of r variation compared to the base random effects model
## seems a little high, could be due to small sample size

## we can again check and plot the variance components
## par(mfrow=c(1,2),oma=c(0,0,1,0),mar=c(4,4,1,1))
summary(mlmem2)$sigma
## profile(mlmem2,sigma2=1)
## profile(mlmem2,sigma2=2)
## both close to zero but still relevant
## you would probably never include these plots in a paper
## but you would want to state the values of sigma2
## also would want to note the MLEs were identifiable

## those are essentially the basics
## for more information, metafor has a great website
## http://www.metafor-project.org/
## code examples can be found in "Analysis Examples" and "Tips and Notes"
## more complex analyses can be done with rma.mv