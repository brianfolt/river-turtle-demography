########################################################################
########  	Population demographics of the Black River Turtle   ########  	
########  	(Rhinoclemmys funerea) at La Selva, Costa Rica      ########  	
########  	              B Folt, April 2020                    ########  	
########################################################################

# Clear the workspace
rm(list=ls())

# Set the working directory
setwd("/Users/Brian/Dropbox/La Selva Turtles/Rhinoclemmys funerea project/GitHub")

# Load the table
datum <- read.csv("captures.csv", header = TRUE)
head(datum)

# Log transform length and mass measurements
logCL <- log(datum$CL)
logPL <- log(datum$PL)
logMass <- log(datum$Mass)
TailRatio <- datum$Plastron_to_vent/datum$Tail_length
datum <- cbind(datum,logCL,logPL,logMass,TailRatio)



#####################
### Morphometrics ###
#####################

morph <- subset(datum, datum$Status!="recap") #remove recaptures
head(morph)
morph

males <- subset(morph, morph$Sex=="M")
females <- subset(morph, morph$Sex=="F")

library(lawstat)

###
### Do males and females differ in body size?
###

#box-and-whiskers plots of adult males and females
plot(morph[-c(15),]$Sex, morph[-c(15),]$logCL, ylab="Carapace length (cm)", xlab="Sex")

#test for assumptions of parametric statistics
shapiro.test(males$CL)		#male carapace length is normal
shapiro.test(females$CL)	#female carapace length is normal
levene.test(morph[-c(15),]$CL, morph[-c(15),]$Sex)	#homoscedasticity is OK

#compare CL of adult males & females with a t-test
(results <- lm(CL ~ Sex, data=morph[-c(15),]))
summary(results)
confint(results)


###
### Do males and females differ in length-mass ratios?
###

# simple model of length-mass with males and females lumped
results2 <- lm(Mass ~ CL, data=morph[-c(15),])
summary(results2)

# multivariable model where males differ in size than females
results3 <- lm(Mass ~ CL + Male, data=morph[-c(15),])
summary(results3)
confint(results3)
	# carapace-mass relationship significant, females > males significant

# saturated multivariable model (with interaction term)
results4 <- lm(Mass ~ CL*Male, data=morph[-c(15),])
summary(results4)

# since intereaction term was nonsignificant, we should remove it
# but first, let's test to see if results4 (complicated) is better than
# results3, a slightly more simple model
anova(results4,results3)
	# non significant, so more complicated model is not better

# thus, the best model explaining length mass is Mass ~ CL + Male
# discuss this model in the results




###
### Plot length-mass relationships for adult males and females
####

plot(Mass ~ CL, data=morph[-c(15),], pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Mass (g)", axes=FALSE, xlim=c(22,36), ylim=c(1500,5000),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(Mass ~ CL, data=morph[-c(15),], pch=as.numeric(Sex),
	lwd=2, cex=2)

chicas <- lm(Mass ~ CL, data=females)
abline(chicas)
machos <- lm(Mass ~ CL, data=males)
abline(machos, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), 
	pch=c(1,3), cex=1.5, pt.lwd=2)


###
### Is sex externally distinguishable by tail length
### and position of the cloaca?
###

### Tail length

# Simple plot
plot(Tail_length ~ CL, pch=as.numeric(Sex), data=morph[-c(15),], 
	xlab="Carapace length (cm)", ylab="Tail length (cm)")

# Custom plot
plot(Tail_length ~ CL, data=morph[-c(15),], pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Tail length (cm)", axes=FALSE, xlim=c(22,36), ylim=c(6,12),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(Tail_length ~ CL, data=morph[-c(15),], pch=as.numeric(Sex),
	lwd=2, cex=2)

chicas <- lm(Tail_length ~ CL, data=females)
abline(chicas)
machos <- lm(Tail_length ~ CL, data=males)
abline(machos, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), 
	pch=c(1,3), cex=1.5, pt.lwd=2)
text(35,11.5,"A",cex=2)


# Let's do a multi-variable regression to see if males differ from females,
# and if there is an interaction between tail length and sex
model <- lm(Tail_length ~ CL*Male, data=datum)
summary(model)

# Interaction not significant, so remove that term and rerun
model2 <- lm(Tail_length ~ CL + Male, data=datum)
summary(model2)
confint(model2)
	# Males have longer tails than females, 
	# but no interaction between length*sex


### Plastron to vent 

# Simple plot
plot(Plastron_to_vent ~ CL, pch=as.numeric(Sex), data=morph[-c(15),])

# Custom plot
plot(Plastron_to_vent ~ CL, data=morph[-c(15),], pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Plastron-vent length (cm)", axes=FALSE, xlim=c(22,36), ylim=c(2,8),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(Plastron_to_vent ~ CL, data=morph[-c(15),], pch=as.numeric(Sex),
	lwd=2, cex=2)

chicas <- lm(Plastron_to_vent ~ CL, data=females)
abline(chicas)
machos <- lm(Plastron_to_vent ~ CL, data=males)
abline(machos, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), 
	pch=c(1,3), cex=1.5, pt.lwd=2)
text(35,7.4,"B",cex=2)

# Let's do a multi-variable regression to see if males differ from females,
# and if there is an interaction between tail length and sex
model <- lm(Plastron_to_vent ~ CL*Male, data=morph[-c(15),])
summary(model)

# Interaction not significant, so remove that term and rerun
model2 <- lm(Plastron_to_vent ~ CL + Male, data=morph[-c(15),])
summary(model2)
confint(model2)
	# The vent is proportionally farther down the tail
	# in males than females


### Ratio of plastron-vent/tail length
plot(TailRatio ~ CL, data=morph[-c(15),])

plot(TailRatio ~ CL, data=morph[-c(15),], pch=as.numeric(Sex), xlab="Carapace length (cm)",
	ylab="Plastron-vent/Tail length ratio", axes=FALSE, xlim=c(22,36), ylim=c(0.40,0.70),
	cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(TailRatio ~ CL, data=morph[-c(15),], pch=as.numeric(Sex),
	lwd=2, cex=2)

chicas <- lm(TailRatio ~ CL, data=females)
abline(chicas)
machos <- lm(TailRatio ~ CL, data=males)
abline(machos, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), 
	pch=c(1,3), cex=1.5, pt.lwd=2)


###
### Is sex externally distinguishable by CL-PL ratio?
###

plot(CL ~ PL, data=morph[-c(15),], pch=as.numeric(Sex), 
	xlab="Plastron length (cm)", ylab="Carapace length (cm)")
model <- lm(CL ~ PL + Male + PL:Male, data=morph[-c(15),])
summary(model)
confint(model)
	# all variable significant -- CL, Male, and interaction


plot(CL ~ PL, data=morph[-c(15),], pch=as.numeric(Sex), xlab="Plastron length (cm)",
	ylab="Carapace length (cm)", axes=FALSE, xlim=c(20,30), 
	ylim=c(22,36), cex.lab=1.3, type="n")
axis(1, cex.lab=1.3, cex.axis=1.3, lwd=3)
axis(2, cex.lab=1.3, cex.axis=1.3, lwd=3)

points(CL ~ PL, data=morph[-c(15),], pch=as.numeric(Sex),
	lwd=2, cex=2)

chicas <- lm(CL ~ PL, data=females)
abline(chicas)
machos <- lm(CL ~ PL, data=males)
abline(machos, lty=2)

legend("topleft", inset=0.05, box.lwd=2,
	c("Female","Male"), 
	pch=c(1,3), cex=1.5, pt.lwd=2)









