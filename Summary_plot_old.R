load(file = "C:/Users/bas0041/Desktop/Reviewer2 Method/prior_summ")
load(file = "C:/Users/bas0041/Dropbox/PhD Project/Manuscripts/Inseason Bayes Updating/Submission 3/Analysis/Original Method/Output/like.summ1")
load(file = "C:/Users/bas0041/Dropbox/PhD Project/Manuscripts/Inseason Bayes Updating/Submission 3/Analysis/Original Method/Output/like.summ2")
load(file = "C:/Users/bas0041/Dropbox/PhD Project/Manuscripts/Inseason Bayes Updating/Submission 3/Analysis/Original Method/Output/post.summ1")
load(file = "C:/Users/bas0041/Dropbox/PhD Project/Manuscripts/Inseason Bayes Updating/Submission 3/Analysis/Original Method/Output/post.summ2")

# dimensional stuff
dates = dimnames(like.summ1)[[1]]; n_dates = length(dates)
years = as.numeric(dimnames(like.summ1)[[3]]); n_years = length(years)
stats = dimnames(like.summ1)[[2]]

post_summ = array(NA, dim = c(dim(post.summ1), 2))
post_summ[,,,1] = post.summ2
post_summ[,,,2] = post.summ1
dimnames(post_summ) = list(dates, stats, years, c(1,2))
like_summ = array(NA, dim = c(dim(like.summ1), 2))
like_summ[,,,1] = like.summ2
like_summ[,,,2] = like.summ1
dimnames(like_summ) = list(dates, stats, years, c(1,2))
rm(like.summ1); rm(like.summ2); rm(post.summ1); rm(post.summ2)

n_models = 2

true_N = read.table("C:/Users/bas0041/Desktop/Reviewer2 Method/2e_Total Run_Data.txt", header = T)
true_N = true_N[true_N$year %in% years,"N"]

calc_errors = function(est, true) {
  e = est - true
  ae = abs(e)
  pe = e/true
  ape = ae/true
  mult = est/true
  list(
    e = e, ae = ae, pe = pe, ape = ape, mult = mult
  )
}

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

m1 = 1; m2 = 2
ppi = 600
jpeg("Summary_old.jpg", h = 5 * ppi, w = 6 * ppi, res = ppi)
par(mfrow = c(2,2), mar = c(1,3.1,1,0.5), yaxs = "i", oma = c(4,0,0,0))
plot(like_mae[,m1], type = "o", pch = 17, ylim = c(0, 0.55), axes = F, ann = F)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(a)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "MAE", font = 2, pos = 2, cex = 1.1)
lines(like_mae[,m2], type = "o", pch = 2, lty = 2)
lines(post_mae[,m1], type = "o", pch = 16)
lines(post_mae[,m2], type = "o", pch = 1, lty = 2)
abline(h = prior_mae, col = "grey")
axis(side = 1, at = 1:n_dates, labels = rep("", n_dates))
axis(side = 2, las = 2)
legend("bottomleft", title = "Likelihood",
       legend = c("w/o RTF", "RTF"),
       pch = c(2, 17), lty = c(2,1), bty = "n", cex = 0.8)
legend("bottom", title = "Posterior",
       legend = c("w/o RTF", "RTF"),
       pch = c(1, 16), lty = c(2,1), bty = "n", cex = 0.8)
legend("bottomright", legend = "", title = "Prior", bty = "n", lty = 1, col = "grey", cex = 0.8)
box()

par(mar = c(1,0.5,1,3.1))
plot(like_mpe[,m1], type = "o", pch = 17, ylim = c(-0.25, 0.06), axes = F, ann = F)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(b)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "MPE", font = 2, pos = 2, cex = 1.1)
lines(like_mpe[,m2], type = "o", pch = 2, lty = 2)
lines(post_mpe[,m1], type = "o", pch = 16)
lines(post_mpe[,m2], type = "o", pch = 1, lty = 2)
abline(h = prior_mpe, col = "grey")
axis(side = 1, at = 1:n_dates, labels = rep("", n_dates))
axis(side = 4, las = 2)
box()

par(mar = c(1,3.1,1,0.5))
plot(like_sig[,m1], type = "o", pch = 17, ylim = c(0, 0.5), axes = F, ann = F)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(c)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "sd(est/true)", font = 2, pos = 2, cex = 1.1)
lines(like_sig[,m2], type = "o", pch = 2, lty = 2)
lines(post_sig[,m1], type = "o", pch = 16)
lines(post_sig[,m2], type = "o", pch = 1, lty = 2)
abline(h = prior_sig, col = "grey")
axis(side = 1, at = 1:n_dates, labels = dates, las = 2)
axis(side = 2, las = 2)
box()

par(mar = c(1,0.5,1,3.1))
plot(like_cv[,m1], type = "o", pch = 17, ylim = c(0, 0.4), axes = F, ann = F)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.05, labels = "(d)", font = 2, pos = 2, cex = 1.1)
text(x = usr[2] + xdiff * 0.02, y = usr[4] - ydiff * 0.15, labels = "PDF CV", font = 2, pos = 2, cex = 1.1)
lines(like_cv[,m2], type = "o", pch = 2, lty = 2)
lines(post_cv[,m1], type = "o", pch = 16)
lines(post_cv[,m2], type = "o", pch = 1, lty = 2)
abline(h = prior_cv, col = "grey")
axis(side = 1, at = 1:n_dates, labels = dates, las = 2)
axis(side = 4, las = 2)
box()

mtext(side = 1, line = 2.5, "Date", outer = T)
dev.off()