## Generate Tables of the paper: effect of early cases (and some other guys)
## By MW, GPLv3+, Oct. 2016

##
## ==== Load stuff
##
source("utilitaires.R")
load("./climat/161028_climep.RData") ## climep
load("./climat/161028_climep_lag.RData") ## climep.lag

## Let's take better AOT!!!
climep <- subset(climep, week < 16 | week > 40) ## Remain in the validity range of the AOT
climep.jfma <- subset(climep, week < 16)

## Only consider HC-year that report something
climep$season <- ifelse(climep$week > 30, as.numeric(climep$year)-1, as.numeric(climep$year)) ## Cut by seasons
ncases <- aggregate(climep$cas, by=list(climep$season, climep$nomfus), FUN=function(l){sum(l, na.rm=T)}) ## Get the number of cases per HC-year
sum(ncases==0) ## 1450 (/2848) -> 50% of the HC do not report cases
keep.hcy <- ncases[ncases$x>0,]
keep.hcy[,3] <- NULL
names(keep.hcy) <- c("season", "nomfus")
climep.reporting <- merge(keep.hcy, climep, by=c("season", "nomfus"), all.x=T, all.y=F)
####
####  ==== Table 1
####
## Table 1, line 1 (before resampling)
rr <-  glm(seu ~ aot, family=binomial(link=logit), data=climep)
summary(rr) ## Est, p-val: 1.920, <2e-16
exp(coef(rr)) ## OR: 6.823
exp(confint.default(rr))## CI: 4.906680636 9.487554055

## Same but focused on JFM
rr <-  glm(seu ~ aot, family=binomial(link=logit), data=climep.jfma)
summary(rr) ## Est, p-val: 0.98, 8.9e-7
exp(coef(rr)) ## OR: 2.7
exp(confint.default(rr))## CI: 1.8-3.9

## Table 1, line 2 (weeks matched)
rat <- 5
size <- sum(climep$seu, na.rm=T)
episamp <-  subset(climep, climep$seu==1)
nepisamp <- subset(climep, climep$seu==0)
distr <- melt(table(episamp$week))
distr <- distr[distr$value!=0,]
names(distr) <- c('week', 'count')

set.seed(41) ## Careful, this is not random anymore...
first <- TRUE
for (i in 1:nrow(distr)) {
    l <- distr[i,]
    su <- subset(nepisamp, as.numeric(week)==l$week)
    sa <- su[sample(1:nrow(su), l$count*rat),]
    if (first) {
        dres.dfw <- sa
        first <- FALSE
    } else {
        dres.dfw <- rbind(dres.dfw, sa)
    }
}

rr <- glm(seu ~ aot, family=binomial(link=logit), data=rbind(episamp, dres.dfw))
summary(rr)             ## -.08, p=0.73 NS
exp(coef(rr))           ## 0.92
exp(confint.default(rr))## CI: 0.5742524 1.4747406


## Table 1, line 3 (weeks + years matched)
rat <- 5
size <- sum(climep$seu, na.rm=T)
episamp <-  subset(climep, climep$seu==1)
nepisamp <- subset(climep, climep$seu==0)
distr <- melt(table(episamp$week, episamp$year))
distr <- distr[distr$value!=0,]
names(distr) <- c('week', 'year', 'count')

set.seed(42) ## Careful, this is not random anymore...
first <- TRUE
for (i in 1:nrow(distr)) {
    l <- distr[i,]
    su <- subset(nepisamp, year==l$year & as.numeric(week)==l$week)
    sa <- su[sample(1:nrow(su), l$count*rat),]
    if (first) {
        dres.dfwy <- sa
        first <- FALSE
    } else {
        dres.dfwy <- rbind(dres.dfwy, sa)
    }
}

rr <- glm(seu ~ aot, family=binomial(link=logit), data=rbind(episamp, dres.dfwy))
summary(rr)             ## 0.1107, p=0.653 NS
exp(coef(rr))           ## 1.12
exp(confint.default(rr))## CI: 0.6890787 1.8107303

####
####  ==== Table 2 -- model with daily AOT!!!!
####
climep.5 <- climep.lag[["5"]]

## Table 2, line 1 (before resampling)
rr <-  glm(seu ~ aot, family=binomial(link=logit), data=climep.5)
summary(rr)             ## Est, p-val: 1.948, <2e-16
exp(coef(rr))           ## OR: 7.017308821 
exp(confint.default(rr))## CI: 5.1529262522 9.556244486

## Table 2, line 2 (weeks matched)
rat <- 5
size <- sum(climep.5$seu, na.rm=T)
episamp <-  subset(climep.5, climep.5$seu==1)
nepisamp <- subset(climep.5, climep.5$seu==0)
distr <- melt(table(episamp$week))
distr <- distr[distr$value!=0,]
names(distr) <- c('week', 'count')

set.seed(40) ## Careful, this is not random anymore...
first <- TRUE
for (i in 1:nrow(distr)) {
    l <- distr[i,]
    su <- subset(nepisamp, as.numeric(week)==l$week)
    sa <- su[sample(1:nrow(su), l$count*rat),]
    if (first) {
        dres.5dfw <- sa
        first <- FALSE
    } else {
        dres.5dfw <- rbind(dres.5dfw, sa)
    }
}

rr <- glm(seu ~ aot, family=binomial(link=logit), data=rbind(episamp, dres.5dfw))
summary(rr)             ## -.184, p=0.407 NS
exp(coef(rr))           ## 0.8323376
exp(confint.default(rr))## CI: 0.5394317 1.284288


## Table 1, line 3 (weeks + years matched)
rat <- 5
size <- sum(climep.5$seu, na.rm=T)
episamp <-  subset(climep.5, climep.5$seu==1)
nepisamp <- subset(climep.5, climep.5$seu==0)
distr <- melt(table(episamp$week, episamp$year))
distr <- distr[distr$value!=0,]
names(distr) <- c('week', 'year', 'count')

set.seed(39) ## Careful, this is not random anymore...
first <- TRUE
for (i in 1:nrow(distr)) {
    l <- distr[i,]
    su <- subset(nepisamp, year==l$year & as.numeric(week)==l$week)
    sa <- su[sample(1:nrow(su), l$count*rat),]
    if (first) {
        dres.5dfwy <- sa
        first <- FALSE
    } else {
        dres.5dfwy <- rbind(dres.5dfwy, sa)
    }
}

rr <- glm(seu ~ aot, family=binomial(link=logit), data=rbind(episamp, dres.5dfwy))
summary(rr)             ## 0.153, p=0.509 NS
exp(coef(rr))           ## 1.16
exp(confint.default(rr))## CI: 0.7405298 1.8322634


####
####  ==== Table 3
####

## Table 3, line 1
rr <- glm(seu ~ precoces, family=binomial(link='logit'), subset(climep[!is.na(climep$seu),], as.Date(date)<as.Date("2014-08-01")))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))

## Table 3, line 2
rr <- glm(seu ~ precoces, family=binomial(link='logit'), subset(climep[!is.na(climep$seu),], !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))

## Table 3, line 3
rr <- glm(seu ~ aot.den + precoces, family=binomial(link='logit'), subset(climep[!is.na(climep$seu),], as.Date(date)<as.Date("2014-08-01")))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))

## Table 3, line 4
rr <- glm(seu ~ aot.den + precoces, family=binomial(link='logit'), subset(climep[!is.na(climep$seu),], !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))


###
### 3 ==== This models should assess for the correlation of AOT 
###        at the time of early cases

## with 2006
rr <- glm(seu ~ aot.precoces + precoces, family=binomial(link='logit'), subset(climep[!is.na(climep$seu),], as.Date(date)<as.Date("2014-08-01")))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))

## without 2006
rr <- glm(seu ~ aot.precoces + precoces, family=binomial(link='logit'), subset(climep[!is.na(climep$seu),], !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))

## Test the association between the number of early cases and early AOT (w/ 2006)
rr <- glm(precoces ~ aot.precoces + seu, family='poisson', subset(climep[!is.na(climep$seu),], as.Date(date)<as.Date("2014-08-01")))
summary(rr) ## neg. corr. between AOT and EC
exp(coef(rr))
exp(confint.default(rr))

## Test the association between the number of early cases and early AOT (no 2006)
rr <- glm(precoces ~ aot.precoces + seu, family='poisson', subset(climep[!is.na(climep$seu),], !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")))
summary(rr) ## neg. corr. between AOT and EC
exp(coef(rr))
exp(confint.default(rr))

###
### 3b ==== This models should assess for the correlation of AOT 
###        at the time of early cases
##  Same analysis as before but with climep.reporting instead of climep

## with 2006 (you can play to add/remove the aot.precoces ver)
rr <- glm(seu ~ aot.precoces + precoces, family=binomial(link='logit'), subset(climep.reporting[!is.na(climep.reporting$seu),], as.Date(date)<as.Date("2014-08-01")))
summary(rr) ## 0.13 (<0.0003) / 0.14 (<0.0002)
exp(coef(rr)) ## 1.14 / 1.15
exp(confint.default(rr)) ## 1.06-1.21 / 1.06-1.23

## without 2006 (you can play to add/remove the aot.precoces ver)
rr <- glm(seu ~ aot.precoces + precoces, family=binomial(link='logit'), subset(climep.reporting[!is.na(climep.reporting$seu),], !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")))
summary(rr) ## 0.18 / 0.20
exp(coef(rr)) ## 1.20 (1.94e-05) / 1.22 (7.37e-06)
exp(confint.default(rr)) ## 1.10-1.30 / 1.12-1.33

## Test the association between the number of early cases and early AOT (w/ 2006)
rr <- glm(precoces ~ aot.precoces + seu, family='poisson', subset(climep[!is.na(climep$seu),], as.Date(date)<as.Date("2014-08-01")))
summary(rr) ## neg. corr. between AOT and EC
exp(coef(rr))
exp(confint.default(rr))

## Test the association between the number of early cases and early AOT (no 2006)
rr <- glm(precoces ~ aot.precoces + seu, family='poisson', subset(climep[!is.na(climep$seu),], !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")))
summary(rr) ## neg. corr. between AOT and EC
exp(coef(rr))
exp(confint.default(rr))



## Plot the early AOT as a function of the number of early cases
pdf("./aot_early_cases/link_aot_ec.pdf")
par(mfrow=c(1,2))
boxplot(split(climep$aot.precoces, f=list(climep$precoces)),
        xlab="#early cases",
        ylab="Early AOT",
        notch=T, main="early AOT = f(#early cases)") ## No pool
boxplot(split(climep$aot.precoces, f=list(climep$precoces>0)),
        xlab="HC-week with early cases", ylab="Early AOT",
        notch=T, main="AOT = f(presence of early cases)") ## Pool
dev.off()

###
### ??? =====
###

## 2 (!?)
climep.save <- climep
climep <- climep.lag[["-5"]]

rat <- 5
size <- sum(climep$seu, na.rm=T)

episamp <-  subset(climep, climep$seu==1)
nepisamp <- subset(climep, climep$seu==0)
distr <- melt(table(episamp$week, episamp$year))
distr <- distr[distr$value!=0,]
names(distr) <- c('week', 'year', 'count')

## Resample! (5 days lag)
first <- TRUE
for (i in 1:nrow(distr)) {
    l <- distr[i,]
    su <- subset(nepisamp, year==l$year & as.numeric(week)==l$week)
    sa <- su[sample(1:nrow(su), l$count*rat),]
    if (first) {
        dres.df <- sa
        first <- FALSE
    } else {
        dres.df <- rbind(dres.df, sa)
    }
}

## Plot!
ep <- list(epi=episamp$aot.den, nonepi=dres.df$aot.den)
par(mfrow=c(3,1))
boxplot(ep)
hist(episamp$aot.den, breaks=50);hist(dres.df$aot.den, breaks=50)
climep <- climep.save


se <- 0.5

####
#### ====Table 2 (?)
####

## Table 2, line 1
rr <-  glm(seu ~ aot.den, family=binomial(link=logit), data=climep)
summary(rr)
exp(coef(rr))
exp(confint.default(rr))

## Je pense que ces valeurs viennent du script 17.160517_climat.R, l. 690 :-)

## Table 2, line 3 (no idea what I'm doing)
rr <- glm(seu ~ aot.den, family=binomial(link=logit), data=subset(rbind(episamp, dres.df), aot.den>se))
summary(rr)
exp(coef(rr))
exp(confint.default(rr))


summary(glm(seu ~ aot.den + precoces, family=binomial(link=logit), data=subset(rbind(episamp, dres.df), aot.den>se)))

summary(glm(seu ~ aot.den + precoces, family=binomial(link=logit), data=subset(rbind(episamp, dres.df), aot.den>se)))
climep <- climep.save
