data <- read.csv("AllDerivs-my.csv", na.strings=c(9999, -9999))

# Lower bounds of the bins
LB <- data[,seq(4,55,3)]

# Upper bounds of the bins
UB <- data[,seq(5,56,3)]

# Probabilities of the bins
PR <- data[,seq(6,57,3)]

# Adjust lower and upper bounds
LB2 <- LB
LB2[,1] <- UB[,1] - (UB[,2] - UB[,1])
UB2 <- UB
for (i in 1:nrow(UB2)) { 
	for (j in 1:ncol(UB2)) { 
		if (is.na(UB2[i,j])) { 
			UB2[i,j] <- LB[i,j] + (LB[i,j] - LB[i,j-1])
		} 
	} 
}

# Central points of the bins
CP2 <- (UB2 + LB2) / 2

# Means of prediction markets
data$pm <- rowSums(CP2 * PR/100, na.rm=TRUE)

# Cummulative probabilities
PR.cum <- PR
for (i in 1:nrow(PR.cum)) { for (j in 2:ncol(PR.cum)) { PR.cum[i,j] <- PR.cum[i,j] + PR.cum[i,j-1] }}

# Quantiles sequence
taus <- 1:19*5/100

# Quantiles: artificially bounded sides, probabilities at the central points
Q <- matrix(NA, nrow(PR.cum), length(taus))
q <- 0
for (qq in taus) {
    q <- q + 1
    for (i in 1:nrow(PR.cum)) {
        j <- 1
        while (j <= ncol(PR.cum)) {
            if (PR.cum[i,j] >= qq*100) {
                if (j==1)
                    Q[i,q] <- LB2[i,j] + (CP2[i,j] - LB2[i,j]) * qq*100 / PR.cum[i,j]
                else
                    Q[i,q] <- CP2[i,j-1] + (CP2[i,j] - CP2[i,j-1]) * (qq*100 - PR.cum[i,j-1]) / (PR.cum[i,j] - PR.cum[i,j-1])
                break
            }
            j <- j + 1
        }
    }
}

# Remove unneeded columns
data[,4:57] <- list(NULL)

# Put quantiles into data set
q_names <- c()
for (i in 1:length(taus)) {
    q_names[i] <- paste("q", taus[i]*100, sep="")
}
Q <- data.frame(Q)
names(Q) <- q_names
data <- cbind(data, Q)

# Normalization of the data
ndata <- data
for (st in levels(data$stat)) {
    ndata[ndata$stat==st,2:ncol(ndata)] <- (data[data$stat==st,2:ncol(data)] - mean(data$rv[data$stat==st])) / sd(data$rv[data$stat==st])
}

# Load quantile regression library
library(quantreg)

# Quantile regression from means
# Expert survey
ex <- list()
for (st in levels(data$stat))
    ex[[st]] <- rq(rv ~ ex, data=data, tau=taus, subset=(stat==st))

# Prediction markets
pm <- list()
for (st in levels(data$stat))
    pm[[st]] <- rq(rv ~ pm, data=data, tau=taus, subset=(stat==st))

# Quantile regression from quantiles
qq <- list()
for (st in levels(data$stat)) {
    for (i in 1:length(taus)) { 
        qq[[st]][[paste('q', taus[i]*100, sep='')]] <- rq(rv ~ data[[4+i]], data=data, tau=taus[i], subset=(stat==st))
    }
    for (i in 1:length(taus)) { 
        qq[[st]][['rho']][i] <- qq[[st]][[paste('q', taus[i]*100, sep='')]]$rho
    }
}

# Loss function (rho) of interpolated quantiles from prediction markets
qp <- list()
for (st in levels(data$stat)) {
    qp[[st]] <- list()
    d <- data[data$stat==st,]
    for (i in 1:length(taus)) { 
        sum <- 0
        for (j in 1:nrow(d)) { 
            diff <- d$rv[j] - d[[4+i]][j]
            if (diff < 0)
                sum <- sum + diff*(i*5/100-1)
            else
                sum <- sum + diff*i*5/100
        }
        qp[[st]][['rho']][i] <- sum
    }
}

# Quantile regression on whole normalized data
ex[['norm']] <- rq(rv ~ ex, data=ndata, tau=taus)
pm[['norm']] <- rq(rv ~ pm, data=ndata, tau=taus)
for (i in 1:length(taus)) { 
    qq[['norm']][[paste('q', taus[i]*100, sep='')]] <- rq(rv ~ ndata[[4+i]], data=ndata, tau=taus[i])
}
for (i in 1:length(taus)) { 
    qq[['norm']][['rho']][i] <- qq[['norm']][[paste('q', taus[i]*100, sep='')]]$rho
}
qp[['norm']] <- list()
for (i in 1:length(taus)) { 
    sum <- 0
    for (j in 1:nrow(ndata)) { 
        diff <- ndata$rv[j] - ndata[[4+i]][j]
        if (diff < 0)
            sum <- sum + diff*(i*5/100-1)
        else
            sum <- sum + diff*i*5/100
    }
    qp[['norm']][['rho']][i] <- sum
}

# Restricted quantile regression (used for R1 calculations)
rr <- list()
for (st in levels(data$stat))
    rr[[st]] <- rq(rv ~ NULL, data=data, tau=taus, subset=(stat==st))
rr[['norm']] <- rq(rv ~ NULL, data=ndata, tau=taus)

# Calculations of R1
for (st in c(levels(data$stat), 'norm')) {
    ex[[st]]$R1 <- 1 - ex[[st]]$rho / rr[[st]]$rho
    pm[[st]]$R1 <- 1 - pm[[st]]$rho / rr[[st]]$rho
    qq[[st]]$R1 <- 1 - qq[[st]]$rho / rr[[st]]$rho
    qp[[st]]$R1 <- 1 - qp[[st]]$rho / rr[[st]]$rho
    qq.ar[[st]]$R1 <- 1 - qq.ar[[st]]$rho / rr[[st]]$rho
    qq.ex[[st]]$R1 <- 1 - qq.ex[[st]]$rho / rr[[st]]$rho
}

# Realised quantiles
# Calculate realised quantiles based on interpolated quantiles data
req <- vector()
for (i in 1:nrow(data)) {
    lb <- which(data[i,5:23] < data$rv[i])
    lbq <- lb[length(lb)] * 5 / 100
    lb <- data[i, 4+lb[length(lb)]]
    ub <- which(data[i,5:23] > data$rv[i])
    ubq <- ub[1] * 5 / 100
    ub <- data[i, 4+ub[1]]
    if (length(lbq) != 0 && lbq == 0.95)
        req[i] <- 1
    else if (length(ubq) != 0 && ubq == 0.05)
        req[i] <- 0
    else
        req[i] <- lbq + (ubq - lbq) * (data$rv[i] - lb) / (ub - lb)
#        req[i] <- lbq
}

# High and Low limits of uniformity test with 10 bins in histogram
ph <- 0.1 + 1.96 * sqrt(0.1*0.9/length(req))
pl <- 0.1 - 1.96 * sqrt(0.1*0.9/length(req))

# Histogram with lines
h <- hist(req, ylim=c(0,25))
abline(h=mean(h$counts)*pl/.1, lty=2)
abline(h=mean(h$counts)*ph/.1, lty=2)

# Calculate realised quantiles based on normalized fitted values of quantile regression on market quantiles
fq <- matrix(nrow=nrow(ndata))
for (i in 1:length(taus)) {
    fq <- cbind(fq, qq$norm[[i]]$fitted)
}
fq <- fq[,2:ncol(fq)]
req <- vector()
for (i in 1:nrow(ndata)) {
    lb <- which(fq[i,] < ndata$rv[i])
    lbq <- lb[length(lb)] * 5 / 100
    lb <- fq[i, lb[length(lb)]]
    ub <- which(fq[i,] >= ndata$rv[i])
    ubq <- ub[1] * 5 / 100
    ub <- fq[i, ub[1]]
    if (length(lbq) != 0 && lbq == 0.95)
        req[i] <- 1
    else if (length(ubq) != 0 && ubq == 0.05)
        req[i] <- 0
    else
        req[i] <- lbq + (ubq - lbq) * (ndata$rv[i] - lb) / (ub - lb)
}

# Calculate realised quantiles based on normalized fitted values of quantile regression on market means
fq <- pm$norm$fitted
req <- vector()
for (i in 1:nrow(ndata)) {
    lb <- which(fq[i,] < ndata$rv[i])
    lbq <- lb[length(lb)] * 5 / 100
    lb <- fq[i, lb[length(lb)]]
    ub <- which(fq[i,] >= ndata$rv[i])
    ubq <- ub[1] * 5 / 100
    ub <- fq[i, ub[1]]
    if (length(lbq) != 0 && lbq == 0.95)
        req[i] <- 1
    else if (length(ubq) != 0 && ubq == 0.05)
        req[i] <- 0
    else
        req[i] <- lbq + (ubq - lbq) * (ndata$rv[i] - lb) / (ub - lb)
}

# Calculate realised quantiles based on normalized fitted values of quantile regression on expert means
fq <- ex$norm$fitted
req <- vector()
for (i in 1:nrow(ndata)) {
    lb <- which(fq[i,] < ndata$rv[i])
    lbq <- lb[length(lb)] * 5 / 100
    lb <- fq[i, lb[length(lb)]]
    ub <- which(fq[i,] >= ndata$rv[i])
    ubq <- ub[1] * 5 / 100
    ub <- fq[i, ub[1]]
    if (length(lbq) != 0 && lbq == 0.95)
        req[i] <- 1
    else if (length(ubq) != 0 && ubq == 0.05)
        req[i] <- 0
    else
        req[i] <- lbq + (ubq - lbq) * (ndata$rv[i] - lb) / (ub - lb)
}

###################################################################################################################
# Graphs
###################################################################################################################

# Interpolated quantiles R1 compared to expert means quantiles R1
for (st in c(levels(data$stat), 'norm')) {
    dev.new()
    min.R1 <- min(ex[[st]]$R1, qp[[st]]$R1)
    max.R1 <- max(ex[[st]]$R1, qp[[st]]$R1)
    plot(taus, ex[[st]]$R1, ylim=c(min.R1, max.R1), 
        main=toupper(st), ylab='R1', xlab='Quantiles', lab=c(19,5,5), type='b')
    legend(x='bottom', 
        legend=c("Markets interpolated quantiles", "Expert mean quantile regression"), 
        col=c('black', 'black'), pch=c(20, 1), cex=0.8)
    points(taus, qp[[st]]$R1, pch=20, type='b')
}

# Interpolated quantiles R1 compared to market means quantiles R1
for (st in c(levels(data$stat), 'norm')) {
    dev.new()
    min.R1 <- min(pm[[st]]$R1, qp[[st]]$R1)
    max.R1 <- max(pm[[st]]$R1, qp[[st]]$R1)
    plot(taus, pm[[st]]$R1, ylim=c(min.R1, max.R1), 
        main=toupper(st), ylab='R1', xlab='Quantiles', lab=c(19,5,5), type='b')
    legend(x='bottom', 
        legend=c("Markets interpolated quantiles", "Markets mean quantile regression"), 
        col=c('black', 'black'), pch=c(20, 1), cex=0.8)
    points(taus, qp[[st]]$R1, pch=20, type='b')
}

# Market means quantiles R1 compared to expert means quantiles R1
for (st in c(levels(data$stat), 'norm')) {
    dev.new()
    min.R1 <- min(pm[[st]]$R1, ex[[st]]$R1)
    max.R1 <- max(pm[[st]]$R1, ex[[st]]$R1)
    plot(taus, ex[[st]]$R1, ylim=c(min.R1, max.R1), 
        main=toupper(st), ylab='R1', xlab='Quantiles', lab=c(19,5,5), type='b')
    legend(x='bottom', 
        legend=c("Markets mean quantile regression", "Experts mean quantile regression"), 
        col=c('black', 'black'), pch=c(20, 1), cex=0.8)
    points(taus, pm[[st]]$R1, pch=20, type='b')
}

# Interpolated quantiles R1 compared to market quantiles regression R1
for (st in c(levels(data$stat), 'norm')) {
    dev.new()
    min.R1 <- min(qq[[st]]$R1, qp[[st]]$R1)
    max.R1 <- max(qq[[st]]$R1, qp[[st]]$R1)
    plot(taus, qq[[st]]$R1, ylim=c(min.R1, max.R1), 
        main=toupper(st), ylab='R1', xlab='Quantiles', lab=c(19,5,5), type='b')
    legend(x='bottom', 
        legend=c("Markets interpolated quantiles", "Markets quantiles quantile regression"), 
        col=c('black', 'black'), pch=c(20, 1), cex=0.8)
    points(taus, qp[[st]]$R1, pch=20, type='b')
}

# Market means regression R1 compared to market quantiles regression R1
for (st in c(levels(data$stat), 'norm')) {
    dev.new()
    min.R1 <- min(qq[[st]]$R1, pm[[st]]$R1)
    max.R1 <- max(qq[[st]]$R1, pm[[st]]$R1)
    plot(taus, qq[[st]]$R1, ylim=c(min.R1, max.R1), 
        main=toupper(st), ylab='R1', xlab='Quantiles', lab=c(19,5,5), type='b')
    legend(x='top', 
        legend=c("Markets means quantile regression", "Markets quantiles quantile regression"), 
        col=c('black', 'black'), pch=c(20, 1), cex=0.8)
    points(taus, pm[[st]]$R1, pch=20, type='b')
}