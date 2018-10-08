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

# automatically sets wd to this location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read in functions
source("1b_Function_Code.R")

# determine location of the output files
out_dir = paste(getwd(), "Output", sep = "/")

#### FIGURE 1: ABUNDANCE AND FORECAST ERROR TIME SERIES, DISTN OF FORECAST ERRORS ####
true_N = read.table("2b_Total Run_Data.txt", stringsAsFactors = F, header = T)
e = log(true_N[1:(nrow(true_N)-1),"N"]) - log(true_N[2:nrow(true_N),"N"])

# years corresponding to true run sizes
y = true_N$year

# years corresponding to forecast errors
ey = true_N$year[2:nrow(true_N)]

# years with long tick marks
long.x = seq(min(y), max(y), 8)

max(abs(e))
round(mean(e), 3)
round(sd(e), 3)

jpeg(paste(out_dir, "Figure1.jpg", sep = "/"), h = 5 * ppi, w = 3 * ppi, res = ppi)
par(mfrow = c(3,1), mar = c(1.5,4,0.5,1), oma = c(1,0,0,0), cex.lab = 1.2)

# (a) run size plot
plot(N ~ year, data = true_N, type = "o", pch = 16,
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
btf = read.table("2a_BTF_Data.txt", stringsAsFactors = F, header = T)
true_N = read.table("2b_Total Run_Data.txt", stringsAsFactors = F, header = T)
rt_ests = read.table("2c_Run Timing_Data.txt", stringsAsFactors = F, header = T)

dates_eval = c("6/10", "6/24", "7/15")
lets = c("(a)", "(b)", "(c)")
d = 3

jpeg(paste(out_dir, "FigureX.jpg", sep = "/"), h = 5 * ppi, w = 3 * ppi, res = ppi)
par(mfrow = c(3,1), mar = c(0.75,2,1,1), oma = c(2.5,2,0,0), cex.lab = 1.2,
    tcl = -0.35, mgp = c(2,0.5,0), yaxs = "i")

for (d in 1:length(dates_eval)) {
  fit_data = prepare_fit_data(dt = dates_eval[d], yr = 1995, loo = F)
  fit = lm(log(N) ~ q * ccpue + d50, data = fit_data)
  
  min_d50 = min(fit_data$d50)
  mean_d50 = mean(fit_data$d50)
  max_d50 = max(fit_data$d50)
  
  min_ccpue = 1
  max_ccpue1 = max(fit_data$ccpue[fit_data$q == 1])
  max_ccpue2 = max(fit_data$ccpue[fit_data$q == 2])
  
  xccpue1 = seq(min_ccpue, max_ccpue1, length = 100)
  xccpue2 = seq(min_ccpue, max_ccpue2, length = 100)
  
  pred1 = expand.grid(ccpue = xccpue1,
                      d50 = c(min_d50, mean_d50, max_d50),
                      q = factor(1))
  pred2 = expand.grid(ccpue = xccpue2,
                      d50 = c(min_d50, mean_d50, max_d50),
                      q = factor(2))
  pred1$fit = exp(predict(fit, pred1))
  pred2$fit = exp(predict(fit, pred2))
  
  plot(N ~ ccpue, data = fit_data,
       pch = 16, xpd = T,
       xlim = c(0, max(btf$ccpue[btf$date == dates_eval[3]])),
       ylim = c(0, 500000),
       col = ifelse(q == 1, "grey", "black"),
       ylab = "", xlab = "", yaxt = "n")
  
  axis(side = 2, at = seq(0, 500000, 100000), labels = seq(0, 500, 100), las = 1)
  
  if (d == 1) {
    legend("topright", 
           legend = c("Early", "Average", "Late"),
           lty = c(2,1,3),
           title = "Run Timing Type", bty = "n", cex = 0.8)
    legend("right",
           legend = c("1984-2007", "2008-2017"),
           pch = 16, col = c("grey", "black"),
           lty = 1, title = "Catchability Period", bty = "n", cex = 0.8)
  }
  
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  
  text(x = usr[2], y = usr[3] + ydiff * 0.2, lets[d], font = 2, pos = 2, cex = 1.2)
  text(x = usr[2], y = usr[3] + ydiff * 0.08, dates_eval[d], font = 2, pos = 2, cex = 1.2)
  
  lines(fit ~ ccpue, data = subset(pred1, d50 == mean_d50), col = "grey", lty = 1)
  lines(fit ~ ccpue, data = subset(pred1, d50 == min_d50), col = "grey", lty = 2)
  lines(fit ~ ccpue, data = subset(pred1, d50 == max_d50), col = "grey", lty = 3)
  
  lines(fit ~ ccpue, data = subset(pred2, d50 == mean_d50), lty = 1)
  lines(fit ~ ccpue, data = subset(pred2, d50 == min_d50), lty = 2)
  lines(fit ~ ccpue, data = subset(pred2, d50 == max_d50), lty = 3)
  box()
}
mtext(side = 1, outer = T, "Cumulative CPE", line = 1, cex = 0.8)
mtext(side = 2, outer = T, "Run Size (1,000s)", line = 0.5, cex = 0.8)
dev.off()

#### RUN TIMING FIGURE ####

rt_ests = read.table("2c_Run Timing_Data.txt", header = T)

at.y = seq(156, 191, 5)
lab_dates = btf$date[btf$year == 2017 & btf$doy %in% at.y]

rt_ests = rt_ests[rt_ests$year >= 1995,]
rt_ests$fcst_lwr = rt_ests$fcst_d50 - 1.96 * rt_ests$fcst_se_d50
rt_ests$fcst_upr = rt_ests$fcst_d50 + 1.96 * rt_ests$fcst_se_d50

jpeg("Output/FigureY.jpg", h = 3 * ppi, w = 3.4 * ppi, res = ppi)
par(mar = c(2,3,0.5, 0.5))
plot(d50 ~ year, data = rt_ests, type = "o", pch = 16, 
     ylim = range(at.y),
     yaxt = "n", xaxt = "n")
lines(fcst_d50 ~ year, data = rt_ests, type = "o",
      lty = 1, pch = 16, col = "grey")
with(rt_ests, arrows(year, fcst_lwr, year, fcst_upr,
                     length = 0, col = "grey"))
axis(side = 2, at = at.y, labels = lab_dates, las = 1, tcl = -0.35, mgp = c(2,0.4,0), cex.axis = 0.8)
mtext(side = 2, expression("Date of D"[50]), cex = 0.8, line = 2)
axis(side = 1, at = seq(1995, 2015, 5), labels = seq(1995, 2015, 5), mgp = c(2,0.4,0), tcl = -0.35, cex.axis = 0.8)
axis(side = 1, at = seq(1995, 2017, 1), labels = rep("", nrow(rt_ests)),tcl = -0.35/2)
legend("bottomleft", lty = 1, pch = 16, col = c("black", "grey"),
       legend = c("Observed", "Forecast"), cex = 0.8, bty = "n", pt.cex = 1.05)
dev.off()


#### FIGURE 3: ERROR SUMMARIES #####
load(paste(out_dir, "like.summ1", sep = "/"))  # summary stats of likelihood PDFs with timing forecast included
load(paste(out_dir, "like.summ2", sep = "/"))  # summary stats of likelihood PDFs without timing forecast
load(paste(out_dir, "post.summ1", sep = "/"))  # summary stats of posterior PDFs with timing forecast included
load(paste(out_dir, "post.summ2", sep = "/"))  # summary stats of posterior PDFs without timing forecast

# the dimensions of each of these arrays is [date,summary.stat,year]

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

### with timing forecast
# maxe = max(abs(c(pe1.1, pe1.2, pef)))
maxe = 0.6
ppi = 600
png("Output/All Year PE_w_fcst.png", h = 8 * ppi, w = 6 * ppi, res = ppi)
par(mfrow = c(6,4), mar = c(0.2,0,0.2,0), oma = c(5,4,1,1))

outer = seq(1, nyrs, 4)
bottom = (nyrs-3):nyrs
sapply(1:nyrs, function(i) {
  
  plot(pe1.1[i,], type = "o", lty = 2, pch = 2, cex = 1.2, ylim = maxe * c(-1,1), yaxt = "n", xaxt = "n",
       xlim = c(0.5, 6.5))
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  
  lines(pe1.2[i,], type = "o", lty = 2, pch = 1, cex = 1.2)
  abline(h = 0, col = "blue", lty = 2)
  abline(h = pef[i], col = "grey", lwd = 2)
  text(x = usr[1], y = usr[4] - ydiff * 0.1, pos = 4, labels = years[i], font = 2, cex = 1.2)
  
  if (i %in% outer) {
    axis(side = 2, at = round(seq(-maxe, maxe, 0.2), 1), labels = round(seq(-maxe, maxe, 0.2), 1), las = 2,
         cex.axis = 1.2)
  }
  if (i %in% bottom) {
    axis(side = 1, at = 1:6, labels = dates, las = 2,
         cex.axis = 1.2)
  }
})

plot.new()

legend("bottom", legend = c("Prior", "Likelihood", "Posterior"),
       lty = c(1,2,2), pch = c(NA, 2, 1), bty = "n", col = c("grey", "black", "black"), cex = 1.2,
       title = "With Timing Fcst.")

mtext(side = 2, outer = T, line = 2.5, "% Error")
mtext(side = 1, outer = T, line = 4, "Date")
dev.off()

### without timing forecast
# maxe = max(abs(c(pe2.1, pe2.2, pef)))
maxe = 0.6
ppi = 600
png("Output/All Year PE_wo_fcst.png", h = 8 * ppi, w = 6 * ppi, res = ppi)
par(mfrow = c(6,4), mar = c(0.2,0,0.2,0), oma = c(5,4,1,1))

outer = seq(1, nyrs, 4)
bottom = (nyrs-3):nyrs
sapply(1:nyrs, function(i) {
  
  plot(pe2.1[i,], type = "o", lty = 1, pch = 17, cex = 1.2, ylim = maxe * c(-1,1), yaxt = "n", xaxt = "n",
       xlim = c(0.5, 6.5))
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  
  lines(pe2.2[i,], type = "o", lty = 1, pch = 16, cex = 1.2)
  abline(h = 0, col = "blue", lty = 2)
  abline(h = pef[i], col = "grey", lwd = 2)
  text(x = usr[1], y = usr[4] - ydiff * 0.1, pos = 4, labels = years[i], font = 2, cex = 1.2)
  
  if (i %in% outer) {
    axis(side = 2, at = round(seq(-maxe, maxe, 0.2), 1), labels = round(seq(-maxe, maxe, 0.2), 1), las = 2,
         cex.axis = 1.2)
  }
  if (i %in% bottom) {
    axis(side = 1, at = 1:6, labels = dates, las = 2,
         cex.axis = 1.2)
  }
})

plot.new()

legend("bottom", legend = c("Prior", "Likelihood", "Posterior"),
       lty = c(1,1,1), pch = c(NA, 17, 16), bty = "n", col = c("grey", "black", "black"), cex = 1.2,
       title = "Without Timing Fcst.")

mtext(side = 2, outer = T, line = 2.5, "% Error")
mtext(side = 1, outer = T, line = 4, "Date")
dev.off()

