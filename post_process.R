load(file = "C:/Users/bas0041/Desktop/Reviewer2 Method/prior_summ")
load(file = "C:/Users/bas0041/Desktop/Reviewer2 Method/post_summ")
load(file = "C:/Users/bas0041/Desktop/Reviewer2 Method/like_summ")

years = as.numeric(rownames(prior_summ))
true_N = read.table("C:/Users/bas0041/Desktop/Reviewer2 Method/2e_Total Run_Data.txt", header = T)
true_N = true_N[true_N$year %in% years,"N"]

Sobj = 92500
maxH = 110000

get_surplus = function(est) {
  surplus = est - Sobj
  surplus = ifelse(surplus > maxH, maxH, surplus)
  surplus = ifelse(surplus < 0, 0, surplus)
  surplus
}


true_surplus = get_surplus(true_N)
prior_surplus = get_surplus(prior_summ[,"50%"])

prior_meanH = sum(prior_surplus)/n_years
prior_p_meet = mean((true_N - prior_surplus) >= 65000)

post_meanH = matrix(NA, n_dates, n_models)
like_meanH = matrix(NA, n_dates, n_models)
post_p_meet = matrix(NA, n_dates, n_models)
like_p_meet = matrix(NA, n_dates, n_models)
for (d in 1:n_dates) {
  for (m in 1:n_models) {
    post_surplus = get_surplus(post_summ[d,"50%",,m])
    like_surplus = get_surplus(like_summ[d,"50%",,m])
    
    post_meanH[d,m] = sum(post_surplus)/n_years
    like_meanH[d,m] = sum(like_surplus)/n_years
    
    post_p_meet[d,m] = mean((true_N - post_surplus) >= 65000)
    like_p_meet[d,m] = mean((true_N - like_surplus) >= 65000)
  }
}

matplot(post_p_meet[,c(4,6)], type = "l", col = "black", ylim = c(0.7, 1))
matplot(like_p_meet[,c(4,6)], type = "l", col = "blue", add = T)
abline(h = prior_p_meet, col = "red")
matplot(post_meanH[,c(4,6)], type = "l", col = "black", ylim = c(80000, 90000))
matplot(like_meanH[,c(4,6)], type = "l", col = "blue", add = T)
abline(h = prior_meanH, col = "red")

prior_summ[1,"50%"]
post_summ[,"50%",1,1]
like_summ[,"50%",1,1]
true_N[1]

categorize = function(surplus) {
  x = factor(ifelse(surplus < 23000, 1, 
         ifelse(surplus < 60000, 2, 3)), levels = c(1,2,3))
}

weights = rbind(
  c( 1,    -0.5, -1),
  c(-0.5,   1,   -0.5),
  c(-1,   - 1,    1))
  
true_cats = true_cat
est_cats = post_cat
get_score = function(est_cats, true_cats, weights) {
  sum(table(est_cats, true_cats) * weights)/length(est_cats)
}

m = 6
date = "6/10"
dates = dimnames(like_summ)[[1]]; n_dates = length(dates)
years = as.numeric(dimnames(like_summ)[[3]]); n_years = length(years)
models = dimnames(like_summ)[[4]]; n_models = length(models)

post_scores = like_scores = matrix(NA, n_dates, n_models)

for (d in 1:n_dates) {
  for (m in 1:n_models) {
    post_scores[d,m] = get_score(categorize(post_summ[dates[d],"50%",,m] - Sobj), true_cat, weights)
    like_scores[d,m] = get_score(categorize(like_summ[dates[d],"50%",,m] - Sobj), true_cat, weights)
  }
}

prior_score = get_score(categorize(prior_summ[,"50%"] - Sobj), true_cat, weights)

matplot(post_scores[,c(4,6)], type = "l", col = "black", ylim = c(0.5, 1))
matplot(like_scores[,c(4,6)], type = "l", col = "blue", add = T)
abline(h = prior_score, col = "red")

prior_cat = categorize(prior_summ[,"50%"] - Sobj)
true_cat = categorize(true_N - Sobj)
post_cat = categorize()
like_cat = categorize(like_summ["6/10","50%",,4] - Sobj)





sum(table(like_cat, true_cat) * weights)/length(years)
sum(table(prior_cat, true_cat) * weights)/length(years)

calc_errors = function(est, true) {
  e = est - true
  ae = abs(e)
  pe = e/true
  ape = ae/true
  list(
    e = e, ae = ae, pe = pe, ape = ape
  )
}

length(like_summ[dates[d],"50%",,1])

like_mape = matrix(NA, n_models, n_dates)
like_mpe = matrix(NA, n_models, n_dates)
post_mape = matrix(NA, n_models, n_dates)
post_mpe = matrix(NA, n_models, n_dates)

for (m in 1:n_models) {
  like_mape[m,] = colMeans(sapply(dates, function(x) {
    calc_errors(like_summ[x,"50%",,names(fits)[m]], true = true_N)$ape
  }))
  
  like_mpe[m,] = colMeans(sapply(dates, function(x) {
    calc_errors(like_summ[x,"50%",,names(fits)[m]], true = true_N)$pe
  }))
  
  post_mape[m,] = colMeans(sapply(dates, function(x) {
    calc_errors(post_summ[x,"50%",,names(fits)[m]], true = true_N)$ape
  }))
  
  post_mpe[m,] = colMeans(sapply(dates, function(x) {
    calc_errors(post_summ[x,"50%",,names(fits)[m]], true = true_N)$pe
  }))
}

prior_mpe = mean(calc_errors(prior_summ[,"50%"], true_N)$pe)
prior_mape = mean(calc_errors(prior_summ[,"50%"], true_N)$ape)

matplot(post_mape, type = "l", ylim = c(0, 0.35), col = "blue")
matplot(like_mape, type = "l", add = T, col = "red")
abline(h = prior_mape)

matplot(post_mpe, type = "l", ylim = c(-0.1, 0.1), col = "blue")
matplot(like_mpe, type = "l", add = T, col = "red")
abline(h = prior_mpe)
