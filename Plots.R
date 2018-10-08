## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## PLOTTING FOR LEAVE-ONE-OUT EVALUATION OF IN-SEASON BAYESIAN RUN ABUNDANCE UPDATING ##
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##

## ::::::::::::::::::::::: ##
## ----------------------- ##
## CODE BY: BEN STATON     ##
## ----------------------- ##
## LAST UPDATED: 7/31/2018 ##
## ----------------------- ##
## ::::::::::::::::::::::: ##

rm(list = ls(all = T))

# resolution of figures
ppi = 600

library(latex2exp)  # FOR GREEK LETTERS WITH SUBSCRIPTS ON PLOTS

# set the working directory to the location of this file
# setwd("C:/Users/bas0041/Dropbox/PhD Project/Manuscripts/Inseason Bayes Updating/Submission 3/)

# read in functions
source("1b_Function_Code.R")

# determine location of the output files
out_dir = paste(getwd(), "Output", sep = "/")

#### FIGURE 1: ABUNDANCE AND FORECAST ERROR TIME SERIES, DISTN OF FORECAST ERRORS ####
true.N = read.table("2e_Total Run_Data.txt", stringsAsFactors = F, header = T)
e = log(true.N[1:(nrow(true.N)-1),"N"]) - log(true.N[2:nrow(true.N),"N"])

# years corresponding to true run sizes
y = true.N$year

# years corresponding to forecast errors
ey = true.N$year[2:nrow(true.N)]

# years with long tick marks
long.x = seq(min(y), max(y), 8)

max(abs(e))
round(mean(e), 3)
round(sd(e), 3)

jpeg(paste(out_dir, "Figure1.jpg", sep = "/"), h = 5 * ppi, w = 3 * ppi, res = ppi)
par(mfrow = c(3,1), mar = c(1.5,4,0.5,1), oma = c(1,0,0,0), cex.lab = 1.2)

# (a) run size plot
plot(N ~ year, data = true.N, type = "o", pch = 16,
     yaxt = "n", xaxt = "n", ylim = c(90000, 450000), ylab = "Run Size (1,000s)")
axis(side = 2, at = seq(90000, 450000, 90000), labels = seq(90, 450, 90), las = 2)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
len = ifelse(y %in% long.x, ydiff * 0.05, ydiff * 0.025)
segments(y, usr[3], y, usr[3] - len, xpd = T)
text(long.x, usr[3] - len[len == max(len)] * 2, labels = long.x, xpd = T)
text(usr[1] + xdiff * 0.05, usr[4] - ydiff * 0.07, "(a)", cex = 1.3, font = 2)

# (b) forecast error time series plot
plot(e ~ ey, type = "n", ylab = TeX("$\\epsilon_{F,t}$"),
     xlim = range(y), las = 1, 
     ylim = max(abs(e)) * c(-1,1), xaxt = "n")
abline(h = 0, col = "grey")
lines(e ~ ey, type = "o", pch = 16)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
len = ifelse(y %in% long.x, ydiff * 0.05, ydiff * 0.025)
segments(y, usr[3], y, usr[3] - len, xpd = T)
text(long.x, usr[3] - len[len == max(len)] * 2, labels = long.x, xpd = T)
text(usr[1] + xdiff * 0.05, usr[4] - ydiff * 0.07, "(b)", cex = 1.3, font = 2)

# (c) histogram of errors
par(xaxs = "i", yaxs = "i")
hist(e, main = "", col = "grey80", xlim = c(-1, 1), las = 1, ylim = c(0,16), yaxt = "n")
axis(side = 2, at = seq(0, 16, 4), labels = seq(0, 16, 4), las = 2)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(usr[1] + xdiff * 0.05, usr[4] - ydiff * 0.07, "(c)", cex = 1.3, font = 2)
text(x = usr[2] - xdiff * 0.41, y = usr[4] - ydiff * 0.1, labels = TeX("$\\Mean(\\epsilon_{F,t}) = 0.009$"), pos = 4)
text(x = usr[2] - xdiff * 0.41, y = usr[4] - ydiff * 0.2, labels = TeX("$\\SD(\\epsilon_{F,t}) = 0.267$"), pos = 4)
box()
dev.off()

#### FIGURE 2: REGRESSION RELATIONSHIPS ####
N.eos.dat = read.table("2c_Regression_Data.txt", sep = " ", stringsAsFactors = F, header = T)

newx = seq(1, 1200, 1)
fit1 = lm(log(N.btf) ~ log(eos), data = subset(N.eos.dat, q == 1))
fit2 = lm(log(N.btf) ~ log(eos), data = subset(N.eos.dat, q == 2))

pred1 = exp(predict(fit1, newdata = data.frame(eos = newx)))
pred2 = exp(predict(fit2, newdata = data.frame(eos = newx)))

ppi = 600
jpeg(paste(out_dir, "Figure2.jpg", sep = "/"), h = 4 * ppi, w = 5 * ppi, res = ppi)
par(mar = c(4,4.5,0.5,1), xaxs = "i", yaxs = "i")
plot(N.btf ~ eos, type = "n", data = N.eos.dat, ylim = c(0, 400000), yaxt = "n",
     xlim = c(0, 1200), pch = N.eos.dat$q, xlab = expression("EOS"[t]), ylab = expression(paste("N"["vuln,t"], " (1,000s)")))
axis(side = 2, at = seq(0, 450000, 50000), labels = seq(0, 450, 50), las = 2)
text(x = N.eos.dat$eos, y = N.eos.dat$N.btf, labels = paste("'", substr(N.eos.dat$year, 3, 4), sep = ""),
     col = ifelse(N.eos.dat$q == 1, "grey", "black"))

lines(pred1 ~ newx, col = "grey", lwd =2)
lines(pred2 ~ newx, lwd = 2)
box()
dev.off()

#### FIGURE 3: ERROR SUMMARIES #####
load(paste(out_dir, "like.summ1", sep = "/"))  # summary stats of likelihood PDFs with timing forecast included
load(paste(out_dir, "like.summ2", sep = "/"))  # summary stats of likelihood PDFs without timing forecast
load(paste(out_dir, "post.summ1", sep = "/"))  # summary stats of posterior PDFs with timing forecast included
load(paste(out_dir, "post.summ2", sep = "/"))  # summary stats of posterior PDFs without timing forecast

# the dimensions of each of these arrays is [date,summary.stat,year]

# get the true abundance
true.N = read.table("2e_Total Run_Data.txt", header = T)

# extrat the evaluated dates and years
dates = dimnames(like.summ1)[[1]]; n.dates = length(dates)
years = as.numeric(dimnames(like.summ1)[[3]]); nyrs = length(years)

# get the summary stats for the priors each year
prior.mean = true.N$N[true.N$year %in% (years - 1)]
true.N = true.N$N[true.N$year %in% years]
prior.summ = matrix(NA, nyrs, 9)
for (y in 1:nyrs) {
  prior_samps = exp(rnorm(1e6, mean = log(prior.mean[y]) - 0.5 * 0.267^2, 0.267))
  prior.summ[y,] = summ(prior_samps)
}
dimnames(prior.summ) = list(years, names(summ(rnorm(10))))
prior.median = prior.summ[,"50%"]

# calculate different summaries of the errors
# x.method.pdf
  # x: type of error
    # e = error (median - truth)
    # ae = absolute error (abs(median - truth))
    # pe = percent error (e/truth)
    # ape = absolute percent error: (ae/truth)
  # method: type of estimator
    # f: run size forecast (this is the prior)
    # 1: with timing forecast
    # 2: without timing forecast
  # pdf: which pdf is used in error calculation
    # 1: likelihood
    # 2: posterior

# errors
ef = prior.median - true.N
e1.1 = apply(like.summ1[,"50%",], 1, function(x) x - true.N)  # method 1 likelihood
e1.2 = apply(post.summ1[,"50%",], 1, function(x) x - true.N)  # method 1 posterior
e2.1 = apply(like.summ2[,"50%",], 1, function(x) x - true.N)  # method 2 likelihood
e2.2 = apply(post.summ2[,"50%",], 1, function(x) x - true.N)  # method 2 posterior

# absolute errors
aef = abs(ef); ae1.1 = abs(e1.1); ae1.2 = abs(e1.2); ae2.1 = abs(e2.1); ae2.2 = abs(e2.2)

# percent errors
pef = ef/true.N
pe1.1 = apply(e1.1, 2, function(x) x/true.N)
pe1.2 = apply(e1.2, 2, function(x) x/true.N)
pe2.1 = apply(e2.1, 2, function(x) x/true.N)
pe2.2 = apply(e2.2, 2, function(x) x/true.N)

# absolute percent errors
apef = aef/true.N
ape1.1 = apply(ae1.1, 2, function(x) x/true.N)
ape1.2 = apply(ae1.2, 2, function(x) x/true.N)
ape2.1 = apply(ae2.1, 2, function(x) x/true.N)
ape2.2 = apply(ae2.2, 2, function(x) x/true.N)

# mean percent errors
mpef = mean(pef)
mpe1.1 = colMeans(pe1.1)
mpe1.2 = colMeans(pe1.2)
mpe2.1 = colMeans(pe2.1)
mpe2.2 = colMeans(pe2.2)

# mean absolute percent errors
mapef = mean(apef)
mape1.1 = colMeans(ape1.1)
mape1.2 = colMeans(ape1.2)
mape2.1 = colMeans(ape2.1)
mape2.2 = colMeans(ape2.2)

# variability of log-scale multiplicative errors
sigf = sd(log(prior.median/true.N))
sig1.1 = apply(apply(like.summ1[,"50%",], 1, function(x) log(x/true.N)), 2, sd)
sig2.1 = apply(apply(like.summ2[,"50%",], 1, function(x) log(x/true.N)), 2, sd)
sig1.2 = apply(apply(post.summ1[,"50%",], 1, function(x) log(x/true.N)), 2, sd)
sig2.2 = apply(apply(post.summ2[,"50%",], 1, function(x) log(x/true.N)), 2, sd)

# the coordinates of the x-axis
lab.days = seq(1, n.dates)

# plot
jpeg(paste(out_dir, "Figure3.jpg", sep = "/"), h = 4* ppi, w = 8 * ppi, res = ppi)

# (a): MAPE by day, method, and PDF type
par(mfrow = c(1,3), mar = c(3,3.5,2,0.25), cex.axis = 1.2, oma = c(0,0,0,0.5))
plot(mape1.1, type = "o", ylim = c(0.1, 0.6), lty = 2, lwd = 2, xaxt = "n", main = "MAPE", pch = 2, las = 1, cex = 1.5, cex.main = 1.5, ylab = "", xlab = "")
lines(mape1.2, type = "o", lty = 2, lwd = 2, pch = 1, cex = 1.5)
lines(mape2.1, type = "o", lty = 1, lwd = 2, pch = 17, cex = 1.5)
lines(mape2.2, type = "o", lty = 1, lwd = 2, pch = 16, cex = 1.5)
abline(h = mapef, col = "grey", lwd = 2)
axis(side = 1, at = (1:n.dates)[lab.days], labels = dates[lab.days], las = 1)
legend("topright", legend = c("Prior", "Likelihood (Hist. Timing)", "Likelihood (Fcst. Timing)", "Posterior (Hist. Timing)", "Posterior (Fcst. Timing)"),
       lwd = c(2,NA,NA,NA,NA), col = c("grey", "black", "black", "black", "black"),
       pch = c(NA, 17, 2, 16, 1), cex = 1.2, bty = "n",
       pt.cex = 1.5)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[1] + xdiff * 0.05, y = usr[4] - ydiff * 0.03, "(a)", font = 2, cex = 1.5)

# (b): MPE by day, method, and PDF type
plot(mpe1.1, type = "o", ylim = c(-0.25, 0.1), lty = 2, lwd = 2, xaxt = "n", main = "MPE", pch = 2, las = 1, cex = 1.5, cex.main = 1.5, ylab = "", xlab = "")
lines(mpe1.2, type = "o", lty = 2, lwd = 2, pch = 1, cex = 1.5)
lines(mpe2.1, type = "o", lty = 1, lwd = 2, pch = 17, cex = 1.5)
lines(mpe2.2, type = "o", lty = 1, lwd = 2, pch = 16, cex = 1.5)
abline(h = mpef, col = "grey", lwd = 2)
axis(side = 1, at = (1:n.dates)[lab.days], labels = dates[lab.days], las = 1)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[1] + xdiff * 0.05, y = usr[4] - ydiff * 0.03, "(b)", font = 2, cex = 1.5)

# (c) variability of errors by day, method, and PDF type
plot(sig1.1, type = "o", ylim = c(0.1, 0.7), lty = 2, lwd = 2, xaxt = "n", main = TeX("$\\sigma$"), pch = 2, las = 1, cex = 1.5, cex.main = 1.5, ylab = "", xlab = "")
lines(sig1.2, type = "o", lty = 2, lwd = 2, pch = 1, cex = 1.5)
lines(sig2.1, type = "o", lty = 1, lwd = 2, pch = 17, cex = 1.5)
lines(sig2.2, type = "o", lty = 1, lwd = 2, pch = 16, cex = 1.5)
abline(h = sigf, col = "grey", lwd = 2)
axis(side = 1, at = (1:n.dates)[lab.days], labels = dates[lab.days], las = 1)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[1] + xdiff * 0.05, y = usr[4] - ydiff * 0.03, "(c)", font = 2, cex = 1.5)
dev.off()

##### COVERAGE TABLE #####
lwrCL = rev(c("2.5%", "10%", "25%"))
uprCL = rev(c("97.5%", "90%", "75%"))

prior_coverage = numeric(3)
like1_coverage = matrix(NA, 3, n.dates)
like2_coverage = matrix(NA, 3, n.dates)
post1_coverage = matrix(NA, 3, n.dates)
post2_coverage = matrix(NA, 3, n.dates)
btwn = function(x, a, b) {
  ifelse (x >= a & x <= b, 1, 0)
}

for (l in 1:3) {
  prior_coverage[l] = round(mean(btwn(x = true.N, a = prior.summ[,lwrCL[l]], b = prior.summ[,uprCL[l]])),2) * 100
  for (d in 1:n.dates) {
    like1_coverage[l,d] = round(mean(btwn(x = true.N, a = like.summ1[dates[d],lwrCL[l],], b = like.summ1[dates[d],uprCL[l],])),2) * 100
    like2_coverage[l,d] = round(mean(btwn(x = true.N, a = like.summ2[dates[d],lwrCL[l],], b = like.summ2[dates[d],uprCL[l],])),2) * 100
    post1_coverage[l,d] = round(mean(btwn(x = true.N, a = post.summ1[dates[d],lwrCL[l],], b = post.summ1[dates[d],uprCL[l],])),2) * 100
    post2_coverage[l,d] = round(mean(btwn(x = true.N, a = post.summ2[dates[d],lwrCL[l],], b = post.summ2[dates[d],uprCL[l],])),2) * 100
  }
}

pretty_post1 = apply(post1_coverage[,dates %in% c("6/10", "6/24", "7/8")], 2, function(x) paste(x, collapse = ", "))
pretty_post2 = apply(post2_coverage[,dates %in% c("6/10", "6/24", "7/8")], 2, function(x) paste(x, collapse = ", "))
pretty_like1 = apply(like1_coverage[,dates %in% c("6/10", "6/24", "7/8")], 2, function(x) paste(x, collapse = ", "))
pretty_like2 = apply(like2_coverage[,dates %in% c("6/10", "6/24", "7/8")], 2, function(x) paste(x, collapse = ", "))
pretty_fcst = rep(paste(prior_coverage, collapse = ", "), 3)

coverage_out = rbind(pretty_fcst, pretty_like1, pretty_like2, pretty_post1, pretty_post2)
colnames(coverage_out) = c("6/10", "6/24", "7/8")

write.csv(coverage_out, paste(out_dir, "coverage_table.csv", sep = "/"))

##### Plot the %Error for Each Year Separately #####

# must have ran the error summary code first!!

# function to plot the within year time-series of errors for one year
within_year_error = function(y, pe1, pe2, pef) {
  ymax = 0.8
  
  if (y > 0) {
    # create plot
    plot(1,1, type = "n", xlim = c(0.5, 6.5), ylim = ymax * c(-1,1), axes = F, ann = F)
    abline(h = 0, col = "grey", lty = 1, lwd = 2)
    usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    
    # draw on errors
    lines(pe1[y,], type = "o", lty = 1, pch = 15, cex = 1.2, col = "blue")
    lines(pe2[y,], type = "o", lty = 1, pch = 15, cex = 1.2, col = "black")
    abline(h = pef[y], col = "red", lwd = 2)
    
    # add year label
    text(x = usr[1], y = usr[4] - ydiff * 0.1, pos = 4, labels = years[y], font = 2, cex = 1.2)
    box(lwd = 2)
    
    # draw an axis for the error if on outer margin
    if (y %in% outer) {
      axis(side = 2, at = round(seq(-ymax, ymax, 0.4), 1), labels = round(seq(-ymax, ymax, 0.4), 1), las = 2,
           cex.axis = 1.2, lwd = 2)
    }
    
    # draw an axis for the date if on bottom
    if (y %in% bottom) {
      axis(side = 1, at = 1:6, labels = dates, las = 2,
           cex.axis = 1.2, lwd = 2)
    }
    
  } else {  # if y is zero, draw a legend 
    plot.new()
    legend("bottom", legend = c("Prior", "Likelihood", "Posterior"),
           lty = c(1,1,1), pch = c(NA, 15, 15), bty = "n", col = c("red", "blue", "black"), cex = 1.2
           )
  }
}


outer = seq(1, nyrs, 4)  # these years are on the outer margin
bottom = (nyrs-3):nyrs   # these years are on the bottom margin

pdf("Output/Within Year Errors.pdf", h = 8, w = 6)

# with timing forecast
par(mfrow = c(6,4), mar = c(0.2,0,0.2,0), oma = c(5,5,2,1))
sapply(c(1:nyrs, 0), function(y) {
  within_year_error(y = y, pe1 = pe1.1, pe2 = pe1.2, pef = pef)
})
mtext(side = 1, outer = T, line = 4, "Date")
mtext(side = 2, outer = T, line = 3.5, "Proportional Error")
mtext(side = 3, outer = T, line = 0.5, "With Timing Forecast", font = 2)

# without timing forecast
par(mfrow = c(6,4), mar = c(0.2,0,0.2,0), oma = c(5,5,2,1))
sapply(c(1:nyrs, 0), function(y) {
  within_year_error(y = y, pe1 = pe2.1, pe2 = pe2.2, pef = pef)
})
mtext(side = 1, outer = T, line = 4, "Date")
mtext(side = 2, outer = T, line = 3.5, "Proportional Error")
mtext(side = 3, outer = T, line = 0.5, "Without Timing Forecast", font = 2)

dev.off()

##### PLOT THE TIME SERIES OF ERRORS BY DATE #####

# must have ran the error summary code first!!
d = 1
errors_ts = function(d) {
  ymax = 0.8
  
  ### with timing forecast
  # make empty plot with correct dimensions
  par(mar = c(2,2,1,0.1), mfrow = c(1,2), oma = c(2,2,1,2))
  plot(1,1, type = "n", axes = F, ann = F, ylim = ymax * c(-1,1), xlim = range(years))
  title("With Timing Forecast")
  abline(h = 0, lty = 1, lwd = 2, col = "grey")
  
  # draw on the errors for each year on this date
  points(pe1.1[,d] ~ years, col = "blue", pch = 15, type = "o")
  points(pe1.2[,d] ~ years, col = "black", pch = 15, type = "o")
  points(pef ~ years, col = "red", pch = 15, type = "o")
  
  # draw on the mean errors across years
  abline(h = mpe1.1[d], col = "blue", lty = 2)
  abline(h = mpe1.2[d], col = "black", lty = 2)
  abline(h = mpef, col = "red", lty = 2)
  
  axis(side = 1)
  axis(side = 2, at = seq(-ymax, ymax, 0.4), labels = seq(-ymax, ymax, 0.4), las = 1)
  box()

  ### without timing forecast
  # make empty plot with correct dimensions
  par(mar = c(2,0.1,1,2))
  plot(1,1, type = "n", axes = F, ann = F, ylim = ymax * c(-1,1), xlim = range(years))
  title("Without Timing Forecast")
  abline(h = 0, lty = 1, lwd = 2, col = "grey")
  
  # draw on the errors for each year on this date
  points(pe2.1[,d] ~ years, col = "blue", pch = 15, type = "o")
  points(pe2.2[,d] ~ years, col = "black", pch = 15, type = "o")
  points(pef ~ years, col = "red", pch = 15, type = "o")
  
  # draw on the mean errors across years
  abline(h = mpe2.1[d], col = "blue", lty = 2)
  abline(h = mpe2.2[d], col = "black", lty = 2)
  abline(h = mpef, col = "red", lty = 2)
  
  axis(side = 1)
  axis(side = 4, at = seq(-ymax, ymax, 0.4), labels = seq(-ymax, ymax, 0.4), las = 1)
  box()
  
  mtext(side = 1, outer = T, line = 1, "Year")
  mtext(side = 2, outer = T, line = 1, "Proportional Error")
  mtext(side = 3, outer = T, line = -0.5, dates[d], font = 2, cex = 1.5)
}

pdf("Output/Time Series Errors.pdf", h = 4, w = 8)
sapply(1:n.dates, errors_ts)
dev.off()
