########################################################################
########  	Population demographics of the Black River Turtle   ########  	
########  	(Rhinoclemmys funerea) at La Selva, Costa Rica      ########  	
########  	              B Folt, April 2020                    ########  	
########################################################################

# Clear the workspace
rm(list=ls())

# Set the working directory
setwd("/Users/Brian/Dropbox/La Selva Turtles/Rhinoclemmys funerea project/GitHub")
setwd("C:/Users/bpf0006/Desktop/river-turtle-demography")


# Load the table
datum = read.csv("captures.csv", header = TRUE)
head(datum)
colnames(datum)[1] = "Species"
str(datum$Species)

# Remove Chelydra acutirostris and juvenile R. funerea dataset
datum = droplevels(subset(datum, datum$Species == "Rhinoclemmys funerea"))#remove Chelydra
datum = droplevels(subset(datum, datum$Juvenile != 1))#remove J
str(datum)
head(datum)


###### Part I --
###### Summarize the capture-mark recapture data
###### and build mark-recapture histories for the study site
studyarea = droplevels(subset(datum, datum$Survey == "Yes"))

# Number of captures in study area
(caps = length(studyarea[,1]))

# Number of individuals captured in study area
(inds = length(summary(as.factor(studyarea$Number))))

# Number of individuals recaptured in study area
recaptured = summary(as.factor(studyarea$Number)) > 1
(recaptured = length(recaptured[recaptured==TRUE]))

# Percentage of animals recaptured (recapture rate)
(recapRate = recaptured/inds)

# Mean number of captures per individual
(mrecaps = mean(table(studyarea$Number)))

# Range of captures per individuals
range(summary(as.factor(studyarea$Number)))

### Look at behavior of turtles
summary(droplevels(subset(datum, datum$Behavior != ""))$Behavior)
table(droplevels(subset(datum, datum$Behavior != ""))$Behavior)


## Create mark-recapture histories for individuals in the study area

# Set up a matrix to save capture histories
inds = sort(as.numeric(as.character(levels(factor(studyarea$Number)))))
surveys = sort(as.numeric(as.character(levels(factor(studyarea$Julian)))))
  # This does miss out on one survey, the last one, 
  # because no individuals were capture
caphist = matrix(NA, length(inds), length(surveys))
rownames(caphist) = inds
colnames(caphist) = surveys

# Mark-recapture histories will be saved state-specific,
# where individuals observed each survey are specified as:
# juveniles = 1, females = 2, males = 3
for (j in 1:length(inds)){
  ind = subset(studyarea, studyarea$Number == inds[j])
  for (k in 1:length(surveys)){
    survey = subset(ind, ind$Julian == surveys[k])
    # Specify J = 1, F = 2, M = 3, or NA for unsampled years
    if(surveys[k] %in% levels(as.factor(studyarea$Julian)) == FALSE) {
      caphist[j,k] = NA} else {
        if(length(survey$Julian) > 0 && survey$AgeSex[1] == "J") {
          caphist[j,k] = 1} else {
            if(length(survey$Julian) > 0 && survey$AgeSex[1] == "F") {
              caphist[j,k] = 2} else {
                if(length(survey$Julian) > 0 && survey$AgeSex[1] == "M") {
                  caphist[j,k] = 3} else {caphist[j,k] = 0}}}}
  }
}

caphist # Capture histories by age-sex state!
caphist[caphist > 0] = 1  # Change all states to 1

# Create a simple dataframe for group categories
group = table(studyarea$Number, studyarea$AgeSex)
group[group > 0] = 1

sex = group[,1]
sex[sex == 1] = "F"
sex[sex == 0] = "M"


###### Part II -- 
###### Mark-recapture analysis
###### Build simple mark-recapture models to understand
###### demographic features of the population
library(assertr)

df = unname(col_concat(caphist, ""))
ch = data.frame(cbind(df,sex), stringsAsFactors=FALSE)
colnames(ch) = c("ch","sex")
ch[,2] = as.factor(sex)

# Use package RMark to build open population models
library(RMark)

# Specify the time interval between samples
t.int = c(2,1,1,2,3,1,1,2,2,3,1,5,2,10,1,2,1,4,84,3,93,9,3,1,1,19,1,1,43,210,1,1)
t.int = t.int/365

# Process the data
turtles.popan = process.data(data=ch, model="POPAN", groups=c("sex"), time.intervals=t.int)

# Run RELEASE goodness-of-fit tests
release.gof(turtles.popan)

# Fit a complete model set of all parameters varying (or not) by sex
mod1 = mark(turtles.popan) # null model
mod2 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex)))
mod3 = mark(turtles.popan, model.parameters=list(p=list(formula=~sex)))
mod4 = mark(turtles.popan, model.parameters=list(pent=list(formula=~sex)))
mod5 = mark(turtles.popan, model.parameters=list(N=list(formula=~sex)))
mod6 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),p=list(formula=~sex)))
mod7 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),pent=list(formula=~sex)))
mod8 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),N=list(formula=~sex)))
mod9 = mark(turtles.popan, model.parameters=list(p=list(formula=~sex),pent=list(formula=~sex)))
mod10 = mark(turtles.popan, model.parameters=list(p=list(formula=~sex),N=list(formula=~sex)))
mod11 = mark(turtles.popan, model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)))
mod12 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),p=list(formula=~sex),pent=list(formula=~sex)))
mod13 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),p=list(formula=~sex),N=list(formula=~sex)))
mod14 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),pent=list(formula=~sex),N=list(formula=~sex)))
mod15 = mark(turtles.popan, model.parameters=list(p=list(formula=~sex),pent=list(formula=~sex),N=list(formula=~sex)))
mod16 = mark(turtles.popan, model.parameters=list(Phi=list(formula=~sex),p=list(formula=~sex),pent=list(formula=~sex),N=list(formula=~sex)))

(models = collect.models(type="POPAN"))
  ## In looking at the data and the results, female survival is strikingly low. 
  ## I think this reflects a stochastic ~outlier result where when I sampled
  ## in June 2015, zero females were captured in the study area. 

summary(mod2) # Top model
summary(mod6) # Second-best model
summary(mod8) # Second-best model
summary(mod7) # Second-best model

# Top model results
summary(mod2)
mod2$results$real[1:4,]
mod2$results$derived$`N Population Size`
mod2$results$derived$`Gross N* Population Size`


popsize = popan.derived(turtles.popan, models)
str(popsize)
surveys2 = c(surveys,surveys)
abundance = cbind(popsize$N,surveys2)
write.csv(abundance, "DerivedAbundance.csv")

### Save top-models and model table to GitHub page
saveRDS(mod2, "POPAN-model-top.rds")
saveRDS(mod6, "POPAN-model-second.rds")
saveRDS(mod8, "POPAN-model-third.rds")
saveRDS(mod7, "POPAN-model-fourth.rds")

topmod = readRDS("POPAN-model-top.rds")
secondmod = readRDS("POPAN-model-second.rds")
thirdmod = readRDS("POPAN-model-third.rds")
fourthmod = readRDS("POPAN-model-fourth.rds")

# Look at secondary model parameter estimates
secondmod$results$real[1:4,]



###### Part III -- 
###### Simple home range size analysis
###### Compare observed home ranges along the stream 
###### between the sexes

# Load the home range data
hr = read.csv("space-use.csv", header=TRUE)
head(hr)

# Remove a juvenile
hr = droplevels(subset(hr, hr$ReproSex != "J"))

# Summarize maximum observed home range size by sex
model = glm(hr$Distance_m ~ hr$ReproSex, family="poisson")
summary(model)
confint(model)

# Translate estimates back to reality from log-link function
exp(summary(model)$coefficients[1,1]) # Female mean
exp(summary(model)$coefficients[1,2]) # Female SE
exp(confint(model)[1,1]) # Female lower CI
exp(confint(model)[1,2]) # Female lower CI
length(subset(hr, hr$ReproSex=="F")$Distance_m)

exp(summary(model)$coefficients[1,1] + summary(model)$coefficients[2,1]) # Male mean
exp(summary(model)$coefficients[1,2] + summary(model)$coefficients[2,2]) # Male SE
exp(confint(model)[1,1] + confint(model)[2,1]) # Male lower CI
exp(confint(model)[1,2] + confint(model)[2,2]) # Male upper CI
length(subset(hr, hr$ReproSex=="M")$Distance_m)


# Compare home-range sizes by sex using non-parametrix Wilcoxon tests
res = wilcox.test(subset(hr, hr$ReproSex=="F")$Distance_m, subset(hr, hr$ReproSex=="M")$Distance_m, alternative="less")
res
  # Males have larger homerange sizes than females, but this 
  # result was not statistically significant P = 0.156

# remove zeroes
hr0 = droplevels(subset(hr, hr$Distance_m != 0))

# compare home-range sizes by sex again, this time without zeroes
res = wilcox.test(subset(hr0, hr$ReproSex=="F")$Distance_m, subset(hr0, hr$ReproSex=="M")$Distance_m, alternative="less")
res
  # With zeroes removed, this analysis supports males having
  # significantly larger home range sizes than females (P = 0.002)


###### Part IV -- 
###### Morphometric analysis
###### Perform a series of morphometric analyses to explore 
###### morphological variation between the sexes

# Log transform length and mass measurements
logCL <- log(datum$CL)
logPL <- log(datum$PL)
logMass <- log(datum$Mass)
TailRatio <- datum$PV/datum$TL
datum <- cbind(datum,logCL,logPL,logMass,TailRatio)

library(lawstat)

morph <- subset(datum, datum$Status!="recap") #remove recaptures
head(morph)
morph

# Create objects for each sex
males <- subset(morph, morph$AgeSex=="M"); females <- subset(morph, morph$AgeSex=="F")

##### Do males and females differ in body size?

# Box-and-whiskers plots of adult males and females
plot(morph$AgeSex, morph$CL, ylab="Carapace length (cm)", xlab="Sex")

# Test for assumptions of parametric statistics
shapiro.test(males$CL)		#male carapace length is normal
shapiro.test(females$CL)	#female carapace length is normal
levene.test(morph$CL, morph$AgeSex)	#homoscedasticity is OK

# Compare CL of adult males & females with a t-test
adults = droplevels(subset(morph, morph$AgeSex != "J"))
adults = droplevels(subset(adults, adults$AgeSex != "JF"))
adults = droplevels(subset(adults, adults$AgeSex != "JM"))

(res <- lm(CL ~ AgeSex, data=adults))
summary(res); confint(res)

# Mean female size, LCL, UCL
res$coefficients[1]; confint(res)[1,1]; confint(res)[1,2]

# Mean male size, LCL, UCL
res$coefficients[1] + res$coefficients[2] # Male size
res$coefficients[1] + res$coefficients[2] + confint(res)[2,1] # Male LCL
res$coefficients[1] + res$coefficients[2] + confint(res)[2,2] # Male UCL
  ### Adult males and females do not differe in maximum carapace length
  ### this confirms suggestions in the literature (Savage 2002) 


##### Do males and females differ in length-mass ratios?

# Simple model of length-mass with males and females lumped
res1 <- lm(Mass ~ CL, data=adults)
summary(res1)

# Multivariate model where males differ in size than females (additive effect of sex)
res2 <- lm(Mass ~ CL + Male, data=adults)
summary(res2)
confint(res2)
	# carapace-mass relationship significant, females > males significant

# Saturated multivariate model; additive model now includes interaction term
res3 <- lm(Mass ~ CL*Male, data=adults)
summary(res3)
  ## since intereaction term was nonsignificant, we should remove it
  ## but first, let's test to see if res3 (complicated) is better than
  ## res2, a slightly more simple model

anova(res3,res2)
	## non significant, so more complicated model is not better

# thus, the best model explaining length mass is Mass ~ CL + Male
# discuss this model in the results
summary(res2)
confint(res2)


###### Does CL-PL relationship vary between sexes?
CLPLmod <- lm(PL ~ CL + Male + CL:Male, data=adults)
summary(CLPLmod)
# interaction non-significant; drop and rebuild model

CLPLmod2 <- lm(PL ~ CL + Male, data=adults)
summary(CLPLmod2)
confint(CLPLmod2)
# all variable significant; additive effect of sex on PL-CL relationship


##### Does tail length vary between sexes?
tailres <- lm(TL ~ CL*Male, data=adults)
summary(tailres)
  # Interaction not significant, so remove that term and rerun

tailres2 <- lm(TL ~ CL + Male, data=adults)
summary(tailres2)
confint(tailres2)
  # Males have longer tails than females, 
  # but no interaction between length*sex


###### Does plastron-vent distance vary between sexes?
pvmod <- lm(PV ~ CL*Male, data=adults)
summary(pvmod)
  # Interaction not significant, so remove that term and rerun

pvmod2 <- lm(PV ~ CL + Male, data=adults)
summary(pvmod2)
confint(pvmod2)
  # The vent is proportionally farther down the tail in males than females



###### Part II -- 
###### Plot morphometric results for graphs in manuscript
dev.off()

## A) Length to mass
summary(res2)

plot(Mass ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Mass (g)", axes=FALSE, xlim=c(20,36), ylim=c(1500,5000),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(Mass ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)

# line for females
lines(x=c(20,36), y=c(res2$coefficients[1]+20*res2$coefficients[2],res2$coefficients[1]+36*res2$coefficients[2]), type="l", lwd=2.5)
# line for males
lines(x=c(20,36), y=c(res2$coefficients[1]+20*res2$coefficients[2]+res2$coefficients[3],res2$coefficients[1]+36*res2$coefficients[2]+res2$coefficients[3]), type="l", lwd=2.5, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), pch=c(1,2), lty=c(1,2), cex=1.5, pt.lwd=2)


### B) Tail length
summary(tailres2)

# Custom plot
plot(TL ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Tail length (cm)", axes=FALSE, xlim=c(20,36), ylim=c(4,12),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(TL ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)

# line for females
lines(x=c(20,36), y=c(tailres2$coefficients[1]+20*tailres2$coefficients[2],tailres2$coefficients[1]+36*tailres2$coefficients[2]), type="l", lwd=2.5)
# line for males
lines(x=c(20,36), y=c(tailres2$coefficients[1]+20*tailres2$coefficients[2]+tailres2$coefficients[3],tailres2$coefficients[1]+36*tailres2$coefficients[2]+tailres2$coefficients[3]), type="l", lwd=2.5, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), pch=c(1,2), lty=c(1,2), cex=1.5, pt.lwd=2)


### C) Plastron to vent length
summary(pvmod2)

# Custom plot
plot(PV ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Plastron-vent length (cm)", axes=FALSE, xlim=c(20,36), ylim=c(2,8),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(PV ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)

# line for females
lines(x=c(20,36), y=c(pvmod2$coefficients[1]+20*pvmod2$coefficients[2],pvmod2$coefficients[1]+36*pvmod2$coefficients[2]), type="l", lwd=2.5)
# line for males
lines(x=c(20,36), y=c(pvmod2$coefficients[1]+20*pvmod2$coefficients[2]+pvmod2$coefficients[3],pvmod2$coefficients[1]+36*pvmod2$coefficients[2]+pvmod2$coefficients[3]), type="l", lwd=2.5, lty=2)

legend("topleft", inset=0.05, box.lwd=2, c("Female","Male"), 
	pch=c(1,2), lty=c(1,2), cex=1.5, pt.lwd=2)


### D) CL-PL relationship
summary(CLPLmod2)
dev.off()

plot(PL ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Plastron length (cm)", axes=FALSE, xlim=c(22,36), 
	ylim=c(20,30), cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(PL ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)

# line for females 
lines(x=c(22,36), y=c(CLPLmod2$coefficients[1]+22*CLPLmod2$coefficients[2],CLPLmod2$coefficients[1]+36*CLPLmod2$coefficients[2]), type="l", lwd=2.5)
# line for males
lines(x=c(22,36), y=c(CLPLmod2$coefficients[1]+22*CLPLmod2$coefficients[2]+CLPLmod2$coefficients[3],CLPLmod2$coefficients[1]+36*CLPLmod2$coefficients[2]+CLPLmod2$coefficients[3]), type="l", lwd=2.5, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), pch=c(1,2), lty=c(1,2), cex=1.5, pt.lwd=2)



#### Make two two-panel graph of the morphometric results

## (1) A two-panel graph of CL-Mass and CL-PL relationships
# Carapace length-mass relationsip
par(mfrow=c(2,1), oma=c(4,4,0,0)+0.5, mar=c(0,0,1,1)+0.5, cex.lab=1.7) 
plot(Mass ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
     ylab="Mass (g)", axes=FALSE, xlim=c(20,36), ylim=c(1500,5000),
     cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)
points(Mass ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)
lines(x=c(20,36), y=c(res2$coefficients[1]+20*res2$coefficients[2],res2$coefficients[1]+36*res2$coefficients[2]), type="l", lwd=2.5)
lines(x=c(20,36), y=c(res2$coefficients[1]+20*res2$coefficients[2]+res2$coefficients[3],res2$coefficients[1]+36*res2$coefficients[2]+res2$coefficients[3]), type="l", lwd=2.5, lty=2)
legend(20,4500, inset=0.05, box.lwd=2,
       c("Female","Male"), pch=c(1,2), lty=c(1,2), cex=1.5, pt.lwd=2)
text(21,4850, "(A)", cex=2)

# Carapace length-plastron length relationship
plot(PL ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
     ylab="Plastron length (cm)", axes=FALSE, xlim=c(22,36), 
     ylim=c(20,30), cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)
points(PL ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)
lines(x=c(22,36), y=c(CLPLmod2$coefficients[1]+22*CLPLmod2$coefficients[2],CLPLmod2$coefficients[1]+36*CLPLmod2$coefficients[2]), type="l", lwd=2.5)
lines(x=c(22,36), y=c(CLPLmod2$coefficients[1]+22*CLPLmod2$coefficients[2]+CLPLmod2$coefficients[3],CLPLmod2$coefficients[1]+36*CLPLmod2$coefficients[2]+CLPLmod2$coefficients[3]), type="l", lwd=2.5, lty=2)
text(22.8,29.3, "(B)", cex=2)

# Axis titles
title(xlab = "Carapace length (cm)", outer=TRUE, line=2.4, cex.sub=1.8, cex.lab=1.8)
title(ylab = "Plastron length (cm)                Mass (g)         ", outer=TRUE, line=2.4, cex.sub=1.6, cex.lab=1.6)


## (2)  A two-panel graph of CL-PV and CL-TL relationships
# Tail length by carapace length
par(mfrow=c(2,1), oma=c(4,4,0,0)+0.5, mar=c(0,0,1,1)+0.5, cex.lab=1.7) 
plot(TL ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
     ylab="Tail length (cm)", axes=FALSE, xlim=c(20,36), ylim=c(4,12),
     cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)
points(TL ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)
lines(x=c(20,36), y=c(tailres2$coefficients[1]+20*tailres2$coefficients[2],tailres2$coefficients[1]+36*tailres2$coefficients[2]), type="l", lwd=2.5)
lines(x=c(20,36), y=c(tailres2$coefficients[1]+20*tailres2$coefficients[2]+tailres2$coefficients[3],tailres2$coefficients[1]+36*tailres2$coefficients[2]+tailres2$coefficients[3]), type="l", lwd=2.5, lty=2)
legend("topleft", inset=0.05, box.lwd=2, c("Female","Male"), 
       pch=c(1,2), lty=c(1,2), cex=1.2, pt.lwd=2)
text(33,11.5, "(A)", cex=2)

# Plastron-vent length by carapace length
plot(PV ~ CL, data=adults, pch=as.numeric(Sex), xlab="Carapace length (cm)",
     ylab="Plastron-vent length (cm)", axes=FALSE, xlim=c(20,36), ylim=c(2,8),
     cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)
points(PV ~ CL, data=adults, pch=as.numeric(Sex), lwd=2, cex=1.5)
lines(x=c(20,36), y=c(pvmod2$coefficients[1]+20*pvmod2$coefficients[2],pvmod2$coefficients[1]+36*pvmod2$coefficients[2]), type="l", lwd=2.5)
lines(x=c(20,36), y=c(pvmod2$coefficients[1]+20*pvmod2$coefficients[2]+pvmod2$coefficients[3],pvmod2$coefficients[1]+36*pvmod2$coefficients[2]+pvmod2$coefficients[3]), type="l", lwd=2.5, lty=2)
text(33,7.6, "(B)", cex=2)

# Axis titles
title(xlab = "Carapace length (cm)", outer=TRUE, line=2.4, cex.sub=1.8, cex.lab=1.8)
title(ylab = "Plastron-vent length (cm)            Tail length (cm)          ", outer=TRUE, line=2.4, cex.sub=1.6, cex.lab=1.6)








