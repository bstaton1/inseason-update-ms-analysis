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

jpeg(paste(out_dir, "Figure2.jpg", sep = "/"), h = 5 * ppi, w = 3 * ppi, res = ppi)
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

#### FIGURE 3: RUN TIMING FIGURE ####

rt_ests = read.table("2c_Run Timing_Data.txt", header = T)

at.y = seq(156, 191, 5)
lab_dates = btf$date[btf$year == 2017 & btf$doy %in% at.y]

rt_ests = rt_ests[rt_ests$year >= 1995,]
rt_ests$fcst_lwr = rt_ests$fcst_d50 - 1.96 * rt_ests$fcst_se_d50
rt_ests$fcst_upr = rt_ests$fcst_d50 + 1.96 * rt_ests$fcst_se_d50

jpeg("Output/Figure3.jpg", h = 3 * ppi, w = 3.4 * ppi, res = ppi)
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


#### FIGURE 4: ERROR SUMMARIES #####
load(file = "Output/prior_summ")
load(file = "Output/post_summ")
load(file = "Output/like_summ")

# the dimensions of each of these arrays is [date,summary.stat,year]

# dimensional stuff
dates = dimnames(post_summ)[[1]]; n_dates = length(dates)
years = as.numeric(dimnames(post_summ)[[3]]); n_years = length(years)
models = dimnames(post_summ)[[4]]; n_models = length(models)

true_N = read.table("2b_Total Run_Data.txt", header = T)
true_N = true_N[true_N$year %in% years,"N"]

d = 1
m = 1

prior_mpe = mean(calc_errors(prior_summ[,"50%"], true_N)$pe)
prior_mae = mean(calc_errors(prior_summ[,"50%"], true_N)$ape)
prior_sig = sd(log(calc_errors(prior_summ[,"50%"], true_N)$mult))
prior_cv = mean(prior_summ[,"sd"]/prior_summ[,"mean"])
like_mae = like_mpe = post_mae = post_mpe = matrix(NA, n_dates, n_models)
like_sig = post_sig = like_mae
like_cv = post_cv = like_mae

for (d in 1:n_dates) {
  for (m in 1:n_models) {
    like_mae[d,m] = mean(calc_errors(like_summ[d,"50%",,m], true_N)$ape)
    like_mpe[d,m] = mean(calc_errors(like_summ[d,"50%",,m], true_N)$pe)
    like_sig[d,m] = sd(log(calc_errors(like_summ[d,"50%",,m], true_N)$mult))
    post_mae[d,m] = mean(calc_errors(post_summ[d,"50%",,m], true_N)$ape)
    post_mpe[d,m] = mean(calc_errors(post_summ[d,"50%",,m], true_N)$pe)
    post_sig[d,m] = sd(log(calc_errors(post_summ[d,"50%",,m], true_N)$mult))
    
    post_cv[d,m] = mean(post_summ[d,"sd",,m]/post_summ[d,"mean",,m])
    like_cv[d,m] = mean(like_summ[d,"sd",,m]/like_summ[d,"mean",,m])
  }
}

pt.cex = 1.3
jpeg("Output/Figure4.jpg", h = 5 * ppi, w = 6 * ppi, res = ppi)
par(mfrow = c(2,2), mar = c(0.3,3.1,1,0.5), yaxs = "i", oma = c(4,0,0,0),
    tcl = -0.35, mgp = c(2,0.5,0))
plot(like_mae[,1], type = "o", pch = 17, ylim = c(0, 0.4), axes = F, ann = F, cex = pt.cex)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(a)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "MAE", font = 2, pos = 2, cex = 1.1)
lines(like_mae[,2], type = "o", pch = 2, lty = 2, cex = pt.cex)
lines(post_mae[,1], type = "o", pch = 16, cex = pt.cex)
lines(post_mae[,2], type = "o", pch = 1, lty = 2, cex = pt.cex)
abline(h = prior_mae, col = "grey")
axis(side = 1, at = 1:n_dates, labels = rep("", n_dates))
axis(side = 2, las = 2)
legend("bottomleft", title = "Likelihood",
       legend = c("w/ RTF", "w/o RTF"),
       pch = c(2, 17), lty = c(2,1), bty = "n", cex = 0.8)
legend("bottom", title = "Posterior",
       legend = c("w/ RTF", "w/o RTF"),
       pch = c(1, 16), lty = c(2,1), bty = "n", cex = 0.8)
legend(x = usr[2] - xdiff * 0.25, y = usr[3] + ydiff * 0.29, legend = "", title = "Prior", bty = "n", lty = 1, col = "grey", cex = 0.8)
box()

par(mar = c(0.3,0.5,1,3.1))
plot(like_mpe[,1], type = "o", pch = 17, ylim = c(-0.1, 0.1), axes = F, ann = F, cex = pt.cex)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(b)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "MPE", font = 2, pos = 2, cex = 1.1)
lines(like_mpe[,2], type = "o", pch = 2, lty = 2, cex = pt.cex)
lines(post_mpe[,1], type = "o", pch = 16, cex = pt.cex)
lines(post_mpe[,2], type = "o", pch = 1, lty = 2, cex = pt.cex)
abline(h = prior_mpe, col = "grey")
axis(side = 1, at = 1:n_dates, labels = rep("", n_dates))
axis(side = 4, las = 2)
box()

par(mar = c(0.3,3.1,1,0.5))
plot(like_sig[,1], type = "o", pch = 17, ylim = c(0, 0.4), axes = F, ann = F, cex = pt.cex)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(c)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = TeX("$\\sigma$"), font = 2, pos = 2, cex = 1.5)
lines(like_sig[,2], type = "o", pch = 2, lty = 2, cex = pt.cex)
lines(post_sig[,1], type = "o", pch = 16, cex = pt.cex)
lines(post_sig[,2], type = "o", pch = 1, lty = 2, cex = pt.cex)
abline(h = prior_sig, col = "grey")
axis(side = 1, at = 1:n_dates, labels = dates, las = 2)
axis(side = 2, las = 2)
box()

par(mar = c(0.5,0.5,1,3.1))
plot(like_cv[,1], type = "o", pch = 17, ylim = c(0, 0.4), axes = F, ann = F, cex = pt.cex)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(d)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "PDF CV", font = 2, pos = 2, cex = 1.1)
lines(like_cv[,2], type = "o", pch = 2, lty = 2, cex = pt.cex)
lines(post_cv[,1], type = "o", pch = 16, cex = pt.cex)
lines(post_cv[,2], type = "o", pch = 1, lty = 2, cex = pt.cex)
abline(h = prior_cv, col = "grey")
axis(side = 1, at = 1:n_dates, labels = dates, las = 2)
axis(side = 4, las = 2)
box()

mtext(side = 1, line = 2.5, "Date", outer = T)
dev.off()

#### FIGURE 5: EACH YEAR ERROR #####

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

