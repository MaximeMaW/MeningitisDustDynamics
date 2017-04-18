## Generate Figure 8 of the paper: link between AOT and meningitis
## By MW, GPLv3+, Oct. 2016
## The panels are the following:
## - i.   Sample maps of AOT vs. epi (incidence or clusters?)
## - ii.  AOT and epi pattern as a function of the calendar week
## - iii. Boxplot showing the AOT within epi and non-epi weeks
##        (shrinked to meningitis season)
## - iv.  Time course of AOT and incidence

##
## ==== Load stuff & prepare data
##
library(rgdal)
library(rgeos)
library(maptools)
library(viridis) # Colormap
source("utilitaires.R")
load("./climat/161028_climep_lag.RData") ## climep.lag (this is a heavy boy)
load("./climat/161028_climep.RData") ## climep
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg

basename <- "../data/carto"
bf.maps <- readOGR(dsn = paste(basename, "WEB", sep="/"), "BUF") # Main administrative info
bf.maps <- gBuffer(bf.maps, width=0, byid=TRUE) # The shapefile is somehow broken. Repair it.
bf.country <- unionSpatialPolygons(bf.maps, bf.maps$ADM0) # Subset the BF
bf.districts <- unionSpatialPolygons(bf.maps, bf.maps$ADM2)
bf.rivers <- readOGR(dsn = paste(basename, "BF", sep="/"), "RIVERS") # Main administrative info
bf.rivers <- gIntersection(bf.rivers, bf.country) # Make sure that we stay within the borders
bf.roads <- readOGR(dsn = paste(basename, "BF", sep="/"), "roads")
bf.roads <- subset(bf.roads, RDLNTYPE==1) # Keep only main roads
bf.roads <- gIntersection(bf.roads, bf.country) # Make sure that we stay within the borders
bf.cities <- readOGR(dsn = paste(basename, "BF", sep="/"), "airp") # Airports -> big cities

amh <- read.csv("../data/epi/150413_bdd_FRFE.csv")
amh$cas[amh$district=="titao" & amh$date=="2011-03-27"] <- 0 ## Correction d'un bug pesant
menin <- subset(amh, maladie=="menin")

fus <- read.csv("../data/geolocalisation/150413_fusion_fs_2500m.csv") ## Carte des fusions
load("../data/carto/fs_par_annee/merge_2500m.RData") # partition.merg
load("../data/foyers/160517_foyers.RData") # foyers, comp, compp

menin.cas <- fusionner.fs(menin, fus, "cas", verbose=FALSE)
menin.pop <- fusionner.fs(menin, fus, "population", verbose=FALSE)
menin.inc <- menin.cas
idx <- 1:(ncol(menin.inc)-1)
menin.inc[,idx] <- menin.cas[,idx]/menin.pop[,idx]

qq <- 40/1e5
seu <- is.epi.week.seuiladj(menin.cas[,1:(ncol(menin.cas)-1)],menin.pop[,1:(ncol(menin.cas)-1)],qq)
seu$annee <- menin.cas$annee
seu$nomfus <- menin.cas$nomfus

##
## ==== Functions
##
map.cols <- function(v, max.inc=NULL, colormap=viridis, reverse=FALSE,
                     overflow="#ff0000", missing="gray") {
    if (is.null(max.inc)) {
        max.inc <- as.integer(max(v, na.rm=TRUE))
    }
    v.int <- as.integer(v)
    cb <- colormap(max.inc+1)
    if (reverse) {
        cb <- rev(cb)
    }
    cb <- c(cb, overflow)
    v.int[v.int>max.inc] <- max.inc+1
    cs <- cb[v.int+1]
    cs[is.na(cs)] <- missing
    return(list(col=cs,
                cb=cb,
                min=0,
                max=max.inc))
    }
}

# Function to plot color bar (from: http://www.colbyimaging.com/wiki/statistics/color-bars)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', new=TRUE)
{
    scale = (length(lut)-1)/(max-min)
    if (new) {
        dev.new(width=1.75, height=5)
    }
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

## Format seasons data
seasons.cas <- to.seasons(menin.cas) ## split into seasons
seasons.cas <- lapply(seasons.cas, function(tt){tt <- tt[,sort(names(tt))]
                                                tt[,] <- tt[order(tt$nomfus),]
                                                rownames(tt) <- tt$nomfus
                                                return(tt)
                                            })
seasons.inc <- to.seasons(menin.inc) ## split into seasons
seasons.inc <- lapply(seasons.inc, function(tt){tt <- tt[,sort(names(tt))]
                                                tt[,] <- tt[order(tt$nomfus),]
                                                rownames(tt) <- tt$nomfus
                                                return(tt)
                                            })
seasons.seu <- to.seasons(seu) ## split into seasons
seasons.seu <- lapply(seasons.seu, function(tt){tt <- tt[,sort(names(tt))]
                                                tt[,] <- tt[order(tt$nomfus),]
                                                rownames(tt) <- tt$nomfus
                                                return(tt)
                                            })

##
## ==== Panel 1
## Choose 3 dates and plot maps of AOT and incidences
## 
plot.that.shit <- function(dat, max.inc=NULL, max.aot=NULL) {
    gray <- "#EAEAEA"
    yy <- substr(dat, 1, 4)
    inc.y <- seasons.inc[[yy]][,c("nomfus", dat)]
    names(inc.y) <- c("NAME", "inc")
    inc.y <- merge(partition.merg[[yy]], inc.y)

    aot.y <- subset(climep, as.Date(date)==(as.Date(dat)), c("nomfus", "aot"))
    names(aot.y) <- c("NAME", "aot")
    aot.y <- merge(partition.merg[[yy]], aot.y)

    colbar <- map.cols(as.integer(inc.y$inc*1e5), colormap=heat.colors, reverse=T,
                       max.inc=max.inc, missing=gray)
    colbar.aot <- map.cols(as.integer(aot.y$aot*1e3),
                           colormap=viridis, reverse=F, max.inc=max.aot*1e3, missing=gray)

    ## Plot the first guy
    plot(bf.districts, main=dat, col=gray, border=NA)
    plot(inc.y, col=colbar[["col"]], border=NA, add=T) # Plot the data
    plot(bf.districts, main=dat, col=NA, border="gray", add=T)
    ##plot(bf.rivers, add=T, col="#0093AF55")
    plot(bf.roads, add=T, col="#FF996655", lwd=2)
    plot(bf.country, lwd=1, add=T)
    plot(bf.cities, add=T, pch=20, cex=2, col="#ED2939")

    ## And the second one
    plot(bf.districts, col=gray, border=NA)
    plot(aot.y, col=colbar.aot[["col"]], border=NA, add=T) # Plot the data
    plot(bf.districts, main=dat, col=NA, border="gray", add=T)    
    ##plot(bf.rivers, add=T, col="#0093AF55")
    plot(bf.roads, add=T, col="#FF996655", lwd=2)
    plot(bf.country, lwd=1, add=T)
    plot(bf.cities, add=T, pch=20, cex=2, col="#ED2939")

    return(list(colbar, colbar.aot))
}

####
#### ==== A few experiments
####

## Extract the 3 highest AOT values (2010-04-25, 2007-04-08, 2012-03-25)
climep.jfma <- subset(climep, week <= 16)
head(climep.jfma[order(climep.jfma$aot, decreasing=T),], n=45) 

## Extract the dates with highest AOT quintile
climep.dat <- split(climep.jfma$aot, f=list(climep.jfma$date))
climep.quint <- do.call("rbind", lapply(climep.dat, FUN=quantile, 0.9, na.rm=T))
sort(climep.quint[!is.na(climep.quint),], decreasing=T)[1:10]

pdf("../papers/woringer2/figures/aot_plots/aot_map2.pdf", height=10, width=20)
par(mfcol=c(2,3))
off <- 0 ## days
cl1 <- plot.that.shit(as.character(as.Date("2010-04-25")+off), max.aot=2, max.inc=200)
cl2 <- plot.that.shit(as.character(as.Date("2007-04-08")+off), max.aot=2, max.inc=200)
cl3 <- plot.that.shit(as.character(as.Date("2012-03-25")+off), max.aot=2, max.inc=200)

par(mfcol=c(2,3))
off <- 7 ## days
cl1 <- plot.that.shit(as.character(as.Date("2010-04-25")+off), max.aot=2, max.inc=200)
cl2 <- plot.that.shit(as.character(as.Date("2007-04-08")+off), max.aot=2, max.inc=200)
cl3 <- plot.that.shit(as.character(as.Date("2012-03-25")+off), max.aot=2, max.inc=200)

par(mfcol=c(2,3))
off <- 14 ## days
cl1 <- plot.that.shit(as.character(as.Date("2010-04-25")+off), max.aot=2, max.inc=200)
cl2 <- plot.that.shit(as.character(as.Date("2007-04-08")+off), max.aot=2, max.inc=200)
cl3 <- plot.that.shit(as.character(as.Date("2012-03-25")+off), max.aot=2, max.inc=200)
dev.off()

pdf("../papers/woringer2/figures/aot_plots/aot_map2_cbinc.pdf", width=1.75, height=5)
colbar <- cl1[[1]]
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=F)      # And the colorbar
dev.off()

pdf("../papers/woringer2/figures/aot_plots/aot_map2_cbaot.pdf", width=1.75, height=5)
colbar <- cl1[[2]]
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=F)      # And the colorbar
dev.off()


####
#### ==== Current stuff
####
pdf("../papers/woringer2/figures/aot_plots/aot_map.pdf", height=10, width=20)
par(mfcol=c(2,3))
cl1 <- plot.that.shit("2006-02-19", max.aot=1, max.inc=200)
cl2 <- plot.that.shit("2008-03-09", max.aot=1, max.inc=200)
cl3 <- plot.that.shit("2012-02-12", max.aot=1, max.inc=200)
dev.off()

pdf("../papers/woringer2/figures/aot_plots/aot_map_cbinc.pdf", width=1.75, height=5)
colbar <- cl1[[1]]
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=F)      # And the colorbar
dev.off()

pdf("../papers/woringer2/figures/aot_plots/aot_map_cbaot.pdf", width=1.75, height=5)
colbar <- cl1[[2]]
color.bar(colbar[["cb"]], colbar[["min"]], colbar[["max"]], new=F)      # And the colorbar
dev.off()

##
## ==== Panel 2
## 
##
pdf("../papers/woringer2/figures/aot_plots/calendar_week.pdf", height=5, width=14)
boxplot(climep$aot.den ~ climep$week, main="AOT vs. calendar week", xlab='calendar week', ylab='AOT')
par(new=T)
plot(unlist(lapply(split(climep$inc, climep$week), mean)), col='red', axes=F, type='l', ylab=NA, xlab=NA)
mtext(side=4, line=3, 'Incidence')
axis(side=4)
dev.off()

##
## ==== Panel 3 
## Note that these figures are not very convincing in many respects
##

## Resample the data
rat <- 5
size <- sum(climep$seu, na.rm=T)
episamp <-  subset(climep, climep$seu==1)
nepisamp <- subset(climep, climep$seu==0)
distr <- melt(table(episamp$week, episamp$year))
distr <- distr[distr$value!=0,]
names(distr) <- c('week', 'year', 'count')

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

pdf("../papers/woringer2/figures/aot_plots/macroscopic_correlation1.pdf", height=6, width=4)

climep.jfma <- subset(climep, week < 16 | week > 40) ## Remain in the validity range of the AOT
probs <- c(0, 0.5, 0.75, 0.9, 0.95, 0.99, .999)
climep.jfma$aotq2 <- cut(climep.jfma$aot, quantile(climep.jfma$aot, probs=probs, na.rm=TRUE))

ss <- split(climep$aot, list(climep$seu)) ## Split by category
names(ss) <- c("Non-epi", "Epi")
ss[["Non-epi"]] <- subset(climep.jfma, seu=FALSE)$aot
ss[["Resampled"]] <-  dres.df$aot

boxplot(ss, notch=TRUE, xlab='', ylab='aot', main=paste("Threshold", qq))
dev.off()

pdf("../papers/woringer2/figures/aot_plots/macroscopic_correlation2.pdf", height=6, width=8)

plot(unlist(lapply(split(climep.jfma$inc, list(as.factor(climep.jfma$aotq2))), mean)),
     xlab='AOT class', ylab='mean incidence', col='red',
     xlim=c(0,7), ylim=c(0,7e-5))
text(1:length(levels(climep.jfma$aotq2)), unlist(lapply(split(climep.jfma$inc, list(as.factor(climep.jfma$aotq2))), mean))-3e-6, levels(climep.jfma$aotq2))
text(1:length(levels(climep.jfma$aotq2)), unlist(lapply(split(climep.jfma$inc, list(as.factor(climep.jfma$aotq2))), mean))-5e-6, paste(probs[1:(length(probs)-1)]*100, "-", probs[2:length(probs)]*100, "%", sep=""))

dev.off()

##
## ==== Panel 4
##
a1 <- unlist(lapply(split(climep$aot.den, list(as.factor(climep$date))), mean, na.rm=TRUE))
a2 <- data.frame(cas=colSums(menin.cas[,1:(ncol(menin.cas)-1)], na.rm=TRUE),
                 pop=colSums(menin.pop[,1:(ncol(menin.cas)-1)], na.rm=TRUE),
                 date=as.Date(names(menin.cas[,1:(ncol(menin.cas)-1)])))
a2$inc <- a2$cas/a2$pop

pdf("../papers/woringer2/figures/aot_plots/aot_timecourse.pdf", height=4, width=16)
plot(as.Date(names(a1)), a1, type="l")
par(new=T)
plot(inc ~ date, data=a2, type="l", col="red", lwd=2, xlab=NA, ylab=NA, axes=F)
mtext(side=4, line=3, 'Incidence')
axis(side=4)
dev.off()
