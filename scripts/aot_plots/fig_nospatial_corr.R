## Generate Figure 9 of the paper: no spatial correlation btw AOT/meningitis
## By MW, GPLv3+, Oct. 2016
## The panels are the following:
## - i.   Colocalization matrix
## - ii.  Aligned reads with/without resampling + profiles

##
## ==== Load stuff
##
library(fields)
source("utilitaires.R")
load("./climat/161028_climep.RData") ## climep
load("./climat/161029_climat.RData") ## climat
load("./climat/161029_climday.RData") ## climday
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg
load("../data/foyers/160517_foyers.RData") # foyers, comp, compp

amh <- read.csv("../data/epi/150413_bdd_FRFE.csv")
amh$cas[amh$district=="titao" & amh$date=="2011-03-27"] <- 0 ## Correction d'un bug pesant
menin <- subset(amh, maladie=="menin")

fus <- read.csv("../data/geolocalisation/150413_fusion_fs_2500m.csv") ## Carte des fusions

menin.cas <- fusionner.fs(menin, fus, "cas", verbose=FALSE)
menin.pop <- fusionner.fs(menin, fus, "population", verbose=FALSE)
menin.inc <- menin.cas
idx <- 1:(ncol(menin.inc)-1)
menin.inc[,idx] <- menin.cas[,idx]/menin.pop[,idx]

##
## ===== Functions
##
getSeason <- function(da){
    ## Return the season of a date
    y <- as.numeric(substr(da, 1, 4))
    if (da < as.Date(paste(y, 07, 01, sep="-"))) {
        return(y)
    } else {
        return(y+1)
    }
}

plot.lag <- function(climep.epi, before, after,
                     plot=TRUE, plotA=TRUE, plotB=TRUE, ...) {
    time.span <- seq(-1*before-6, after, 1)
    climep.epi$year <- as.numeric(as.character(climep.epi$year))
    climep.epi$season <- unlist(lapply(climep.epi$date, getSeason))
    aligned <- matrix(NA, ncol=before+after+7, nrow=nrow(climep.epi))

    ##tss <- time.span<=0 & time.span>-10 ## this is cheating
    for (i in 1:nrow(climep.epi)) {
        ##print(i)
        l <- climep.epi[i,]
        t.nomfus <- as.character(l$nomfus)
        t.year <- as.character(l$season)
        t.date <- as.Date(as.character(l$date))
        t.span <- as.character(seq(t.date-before-6, t.date+after, by=1))
        c.names <- as.character(names(climday[[t.year]]))
        res <- as.numeric(climday[[t.year]][t.nomfus, c.names %in% t.span])
        if (length(res)==ncol(aligned)) {
            aligned[i,] <- res
        }
    }

    aot.means <- colMeans(aligned, na.rm=T)
    ##time.span[which.max(aot.means)] ## The time difference
    aot.smooth <- ksmooth(time.span, aot.means, bandwidth=7)
    idx <- order(apply(aligned[,], 1, mean, na.rm=TRUE))

    if (plot) {
        if (plotA) {
            image.plot(time.span, 1:nrow(aligned), t(aligned[idx,]),
                       col=rev(heat.colors(30)),
                       xlab='Position w.r.t epidemic onset',
                       ylab='localized epidemic id', ...)
            lines(c(0, 0), c(0,nrow(aligned)), col="red", lwd=2, lty=2)
        }
        if (plotB) {
            plot(time.span, aot.means, type='l', main="Mean of aligned 'reads'", ylab='mean AOT', xlab='Position w.r.t epidemic onset')
            lines(aot.smooth, col='orange')
            lines(c(0, 0), c(0,1), col="red", lwd=2, lty=2)
        }
    }
    
    return(list(time=time.span, aot=aot.means, aot.smooth=aot.smooth))
}


##
## ==== Panel 1
## /!\ Partial

## Build climate map
cli <- climat[[1]]
for (i in climat[2:length(climat)]) {
    cli <- merge(cli, i, by='nomfus', all.x=TRUE, all.y=TRUE)
    cli$annee.x <- NULL
    cli$annee.y <- NULL
}

## Match the two matrices: drop stuff on both sides
m.cli <- cli[,colnames(cli) %in% colnames(menin.cas)]
m.cas <- menin.cas[,colnames(menin.cas) %in% colnames(m.cli)]
m.inc <- menin.inc[,colnames(menin.inc) %in% colnames(m.cli)]
m.cli <- m.cli[rownames(m.cli) %in% rownames(m.cas),]
m.cas <- m.cas[rownames(m.cas) %in% rownames(m.cli),]
m.inc <- m.inc[rownames(m.inc) %in% rownames(m.cli),]
m.cli <- m.cli[,order(colnames(m.cli))]
m.cas <- m.cas[,order(colnames(m.cas))]
m.inc <- m.inc[,order(colnames(m.inc))]
m.cli <- m.cli[order(rownames(m.cli)),]
m.cas <- m.cas[order(rownames(m.cas)),]
m.inc <- m.inc[order(rownames(m.inc)),]

tst <- TRUE
if (tst) {
    all(colnames(m.cli)==colnames(m.cas))
    all(colnames(m.cli)==colnames(m.inc))
    all(rownames(m.cli)==rownames(m.cas))
    all(rownames(m.cli)==rownames(m.inc))
    dim(m.cas)
    dim(m.cli)
    dim(m.inc)
}

## Apply the same masks
m.cli[is.na(m.cas)] <- NA
m.cli[is.na(m.inc)] <- NA
m.cli[is.na(m.cli)] <- NA
m.cas[is.na(m.cli)] <- NA
m.cas[is.na(m.cas)] <- NA
m.inc[is.na(m.cli)] <- NA
m.inc[is.na(m.cas)] <- NA
m.inc[is.na(m.inc)] <- NA

z <- as.matrix(m.inc[,1:(ncol(m.cas)-1)]>40e-5)-2*as.matrix(m.cli[,1:(ncol(m.cas)-1)]>1)
y <- as.Date(colnames(z))
m <- as.numeric(substr(y, 6,7))
mm <- (m %in% c(1,2,3,4))
di <- as.factor(unlist(lapply(strsplit(m.inc$nomfus, "__"), function(el){el[[1]]})))
or <- order(as.numeric(di))
di <- sort(di)

## Compute scale bar (horizontal)
bar.r <- 20
bar <- matrix(rep(as.numeric(di)%%2+2, times=bar.r), ncol=bar.r)
bar <- cbind(bar, matrix(rep(2,times=nrow(bar)*4),nrow=nrow(bar)))
yy <- y[1]-rev((1:(ncol(bar))))*3

## Compute scale bar (vertical)
bar.r <- 20
bar.x <- matrix(rep(as.numeric(mm)%%2+1, times=bar.r), ncol=bar.r)
bar.x <- cbind(bar.x, matrix(rep(2,times=nrow(bar.x)*4),nrow=nrow(bar.x)))
bar.x <- t(rbind(matrix(NA, nrow=bar.r+4, ncol=bar.r+4), bar.x))
bar.x[bar.x==1] <- NA
bar.x[bar.x==2] <- -1
bar.x <- bar.x[1:2,]

## Compute scale bar (again)
v0 <- which(c(0,diff(as.numeric(di)),1)==1)
v1 <- as.integer(diff(c(0,v0))/2)+c(0,v0[1:(length(v0)-1)])

## Colors
cols <- c('blue', 'green', 'white', '#F5F5F5', "red","red", "black", "#F5F5F5")

pdf("../papers/woringer2/figures/aot_plots/high_dust_80.pdf", width=20, height=10)

image(y=c(yy,y),z=rbind(bar.x+2,cbind(bar,z[or,])), xaxt="n", ylab="Date",
      col=cols, main="Co-occurence of high AOT and high meningitis incidence")
axis(side=1, at=v1/nrow(z), labels=levels(di), tick=FALSE, line=-.75)
legend('topright', legend=c('aot', 'both', 'NA', 'None', 'meningitis'),
       pch=15, col=cols)

dev.off()

##
## ==== Supplementary fig 12 (various thresholds)
##
z.05 <- as.matrix(m.inc[,1:(ncol(m.cas)-1)]>05e-5)-2*as.matrix(m.cli[,1:(ncol(m.cas)-1)]>1)
z.20 <- as.matrix(m.inc[,1:(ncol(m.cas)-1)]>20e-5)-2*as.matrix(m.cli[,1:(ncol(m.cas)-1)]>1)
z.40 <- as.matrix(m.inc[,1:(ncol(m.cas)-1)]>40e-5)-2*as.matrix(m.cli[,1:(ncol(m.cas)-1)]>1)
z.80 <- as.matrix(m.inc[,1:(ncol(m.cas)-1)]>80e-5)-2*as.matrix(m.cli[,1:(ncol(m.cas)-1)]>1)

pdf("../papers/woringer2/figures/aot_plots/high_dust_cb.pdf", width=10, height=20)

par(mfrow=c(4,1))

image(y=c(yy,y),z=rbind(bar.x+2,cbind(bar,z.05[or,])), xaxt="n", ylab="Date",
      col=cols, main="Co-occurence of high AOT and high meningitis incidence (threshold 5/100000)")
axis(side=1, at=v1/nrow(z), labels=levels(di), tick=FALSE, line=-.75)
legend('topright', legend=c('aot', 'both', 'NA', 'None', 'meningitis'),
       pch=15, col=cols)

image(y=c(yy,y),z=rbind(bar.x+2,cbind(bar,z.20[or,])), xaxt="n", ylab="Date",
      col=cols, main="Co-occurence of high AOT and high meningitis incidence (threshold 20/100000)")
axis(side=1, at=v1/nrow(z), labels=levels(di), tick=FALSE, line=-.75)
legend('topright', legend=c('aot', 'both', 'NA', 'None', 'meningitis'),
       pch=15, col=cols)

image(y=c(yy,y),z=rbind(bar.x+2,cbind(bar,z.40[or,])), xaxt="n", ylab="Date",
      col=cols, main="Co-occurence of high AOT and high meningitis incidence (threshold 40/100000)")
axis(side=1, at=v1/nrow(z), labels=levels(di), tick=FALSE, line=-.75)
legend('topright', legend=c('aot', 'both', 'NA', 'None', 'meningitis'),
       pch=15, col=cols)

image(y=c(yy,y),z=rbind(bar.x+2,cbind(bar,z.80[or,])), xaxt="n", ylab="Date",
      col=cols, main="Co-occurence of high AOT and high meningitis incidence (threshold 80/100000)")
axis(side=1, at=v1/nrow(z), labels=levels(di), tick=FALSE, line=-.75)
legend('topright', legend=c('aot', 'both', 'NA', 'None', 'meningitis'),
       pch=15, col=cols)


dev.off()

##
## ==== Panel 2
##
## /!\ Partial

## The idea is to extract the daily AOT values from -20 days to + 10 days around the epidemic week. This allows to probe for the AOT environment, and also to detect potential patterns. This method is inspired from visualization often used in the RNA sequencing/genomics community.
before <- 30 # number of days to probe before the epidemic event
after <- 20 # number of days to probe after the epidemic event

climep.epi <- subset(climep, seu==TRUE) ## All epidemics
climep.epi.p <- plot.lag(climep.epi, before, after, plot=FALSE)


resamp.tc <- list() ## Resampled dates JFM (an old strategy)
for (i in 1:40) {
    d.y <- data.frame(year=climep.epi$year)
    d.y$year.n <- as.Date(paste(d.y$year, "01", "01", sep="-"))
    d.y$day.n <- sample(1*30+2*31, nrow(d.y), replace=TRUE)
    climep.epi.resample <- climep.epi
    climep.epi.resample$date <- d.y$year.n+d.y$day.n
    resamp.tc[[i]] <- plot.lag(climep.epi.resample, before, after, plot=FALSE)
}

resamp.sc <- list() ## Spatial scrambling
nomfus.l <- as.character(levels(as.factor(climep$nomfus)))
for (i in 1:40) {
    climep.epi.resample <- climep.epi
    ##climep.epi.resample$date <- sample(climep.epi.resample$date)
    climep.epi.resample$nomfus <- sample(nomfus.l, nrow(climep.epi.resample))
    resamp.sc[[i]] <- plot.lag(climep.epi.resample, before, after, plot=FALSE)
}

clust.nomfus <- data.frame( ## Clusters
    nomfus=names(comp$membership[!duplicated(comp$membership)]),
    component=as.numeric(comp$membership[!duplicated(comp$membership)]))
climep.epi.clust <- merge(compp, clust.nomfus, by="component")
climep.epi.clust$date <- climep.epi.clust$startdate
climep.epi.clust$nomfus <- unlist(lapply(strsplit(as.character(climep.epi.clust$nomfus), " "),
                                         function(l){l[[1]]}))
climep.epi.clust.p <- plot.lag(climep.epi.clust, before, after, plot=FALSE)


r.tc <- do.call("cbind", lapply(resamp.tc, function(l){l[["aot"]]}))
q.tc <- apply(r.tc, 1, quantile, c(.05,.95))
t.tc <- resamp.tc[[1]][["time"]]
r.sc <- do.call("cbind", lapply(resamp.sc, function(l){l[["aot"]]}))
q.sc <- apply(r.sc, 1, quantile, c(.05,.95))
t.sc <- resamp.sc[[1]][["time"]]
t.glo <- seq(-1*before-6, after, 1)
k.sc1 <- ksmooth(t.glo,q.sc[1,],bandwidth=7)
k.sc2 <- ksmooth(t.glo,q.sc[2,],bandwidth=7)
k.scm <- ksmooth(t.glo,rowMeans(r.sc),bandwidth=7)
k.tc1 <- ksmooth(t.glo,q.tc[1,],bandwidth=7)
k.tc2 <- ksmooth(t.glo,q.tc[2,],bandwidth=7)
k.tcm <- ksmooth(t.glo,rowMeans(r.tc),bandwidth=7)


## Actually make the plot
pdf("../papers/woringer2/figures/aot_plots/timelag.pdf", width=12, height=7)

plot(c(-1*before, after), c(0,1), pch="", xlab="Time lag (days)", ylab="AOT")
##lines(t.tc, rowMeans(r.tc), type="l", col="orange", lwd=2)
##lines(t.tc, q.tc[1,], col="orange")
##lines(t.tc, q.tc[2,], col="orange")
polygon(c(k.sc1$x, rev(k.sc2$x)), c(k.sc1$y, rev(k.sc2$y)),
        border="yellow", col="#EAEAEACC")
polygon(c(k.tc1$x, rev(k.tc2$x)), c(k.tc1$y, rev(k.tc2$y)),
        border="orange", col="#EAEAEACC")
lines(k.scm, col="yellow", lwd=2)
lines(k.tcm, col="orange", lwd=2)

lines(c(0, 0), c(0,1), col="black", lwd=2, lty=2)
lines(climep.epi.p$aot.smooth, col="red", lwd=2)
lines(climep.epi.clust.p$aot.smooth, col="blue", lwd=2)

legend("bottomleft", legend=c("AOT wrt. localized epidemics",
                              "AOT wrt. epidemic clusters",
                              "Temporal resampling",
                              "Spatial resampling"),
       col=c("red", "blue", "orange", "yellow"), lty=1, lwd=2)

dev.off()

pdf("../papers/woringer2/figures/aot_plots/timelag_matrix.pdf", width=8, height=8)

climep.epi <- subset(climep, seu==TRUE) ## All epidemics
climep.epi.p <- plot.lag(climep.epi, before, after, plot=TRUE)

dev.off()

####
#### ==== NO 2006 resampling
####

pdf("../papers/woringer2/figures/aot_plots/timelag_matrix_no2006.pdf", width=8, height=8)

climep.epi <- subset(climep, seu==TRUE & !(as.Date(date)>"2005-08-01" & as.Date(date)<"2006-06-01") & as.Date(date)<as.Date("2014-08-01")) ## All epidemics !!! REMoVE 2006
climep.epi.p <- plot.lag(climep.epi, before, after, plot=TRUE)

dev.off()

resamp.tc <- list() ## Resampled dates JFM (an old strategy)
for (i in 1:40) {
    d.y <- data.frame(year=climep.epi$year)
    d.y$year.n <- as.Date(paste(d.y$year, "01", "01", sep="-"))
    d.y$day.n <- sample(1*30+2*31, nrow(d.y), replace=TRUE)
    climep.epi.resample <- climep.epi
    climep.epi.resample$date <- d.y$year.n+d.y$day.n
    resamp.tc[[i]] <- plot.lag(climep.epi.resample, before, after, plot=FALSE)
}

resamp.sc <- list() ## Spatial scrambling
nomfus.l <- as.character(levels(as.factor(climep$nomfus)))
for (i in 1:40) {
    climep.epi.resample <- climep.epi
    ##climep.epi.resample$date <- sample(climep.epi.resample$date)
    climep.epi.resample$nomfus <- sample(nomfus.l, nrow(climep.epi.resample))
    resamp.sc[[i]] <- plot.lag(climep.epi.resample, before, after, plot=FALSE)
}

clust.nomfus <- data.frame( ## Clusters
    nomfus=names(comp$membership[!duplicated(comp$membership)]),
    component=as.numeric(comp$membership[!duplicated(comp$membership)]))
climep.epi.clust <- merge(compp, clust.nomfus, by="component")
climep.epi.clust$date <- climep.epi.clust$startdate
climep.epi.clust$nomfus <- unlist(lapply(strsplit(as.character(climep.epi.clust$nomfus), " "),
                                         function(l){l[[1]]}))
climep.epi.clust.p <- plot.lag(climep.epi.clust, before, after, plot=FALSE)


r.tc <- do.call("cbind", lapply(resamp.tc, function(l){l[["aot"]]}))
q.tc <- apply(r.tc, 1, quantile, c(.05,.95))
t.tc <- resamp.tc[[1]][["time"]]
r.sc <- do.call("cbind", lapply(resamp.sc, function(l){l[["aot"]]}))
q.sc <- apply(r.sc, 1, quantile, c(.05,.95))
t.sc <- resamp.sc[[1]][["time"]]
t.glo <- seq(-1*before-6, after, 1)
k.sc1 <- ksmooth(t.glo,q.sc[1,],bandwidth=7)
k.sc2 <- ksmooth(t.glo,q.sc[2,],bandwidth=7)
k.scm <- ksmooth(t.glo,rowMeans(r.sc),bandwidth=7)
k.tc1 <- ksmooth(t.glo,q.tc[1,],bandwidth=7)
k.tc2 <- ksmooth(t.glo,q.tc[2,],bandwidth=7)
k.tcm <- ksmooth(t.glo,rowMeans(r.tc),bandwidth=7)


## Actually make the plot
pdf("../papers/woringer2/figures/aot_plots/timelag_no2006.pdf", width=12, height=7)

plot(c(-1*before, after), c(0,1), pch="", xlab="Time lag (days)", ylab="AOT")
##lines(t.tc, rowMeans(r.tc), type="l", col="orange", lwd=2)
##lines(t.tc, q.tc[1,], col="orange")
##lines(t.tc, q.tc[2,], col="orange")
polygon(c(k.sc1$x, rev(k.sc2$x)), c(k.sc1$y, rev(k.sc2$y)),
        border="yellow", col="#EAEAEACC")
polygon(c(k.tc1$x, rev(k.tc2$x)), c(k.tc1$y, rev(k.tc2$y)),
        border="orange", col="#EAEAEACC")
lines(k.scm, col="yellow", lwd=2)
lines(k.tcm, col="orange", lwd=2)

lines(c(0, 0), c(0,1), col="black", lwd=2, lty=2)
lines(climep.epi.p$aot.smooth, col="red", lwd=2)
lines(climep.epi.clust.p$aot.smooth, col="blue", lwd=2)

legend("bottomleft", legend=c("AOT wrt. localized epidemics",
                              "AOT wrt. epidemic clusters",
                              "Temporal resampling",
                              "Spatial resampling"),
       col=c("red", "blue", "orange", "yellow"), lty=1, lwd=2)

dev.off()
