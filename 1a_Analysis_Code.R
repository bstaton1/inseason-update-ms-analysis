## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## LEAVE-ONE-OUT EVALUATION OF IN-SEASON BAYESIAN RUN ABUNDANCE UPDATING ##
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##

## ::::::::::::::::::::::: ##
## ----------------------- ##
## CODE BY: BEN STATON     ##
## ----------------------- ##
## LAST UPDATED: 7/31/2018 ##
## ----------------------- ##
## ::::::::::::::::::::::: ##

# this is the analysis presented in Staton and Catalano
 # "Bayesian information updating procedures for Pacific salmon run size indicators:
 # Evaluation in the presence and absence of auxiliary migration timing information"

# running this whole script should take approximately 40 minutes

# method codes
# 1 = with timing forecast
# 2 = without timing forecast (null model)

# automatically sets wd to this location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(mvtnorm)    # for drawing random regression coefficients from a bivariate normal dist
library(coda)       # for MCMC diagnostics

##### PART 0: SESSION SET UP #####
rm(list = ls(all = T))

starttime_all = Sys.time()

# directory for output
out_dir = paste(getwd(), "Output", sep = "/")
if(!dir.exists(out_dir)) dir.create(out_dir)

# load in other functions necessary for this analysis
source("1b_Function_Code.R")

##### PART 1: DATA AQUISITION #####
# get historical BTF data
  # day: day of sampling (1 = June 1 always, last is 85 = 8/24 always)
  # doy: day of year
  # date: M/DD formated date
  # year: year of sampling
  # cpue: daily cpue value
  # ccpue: cumulative cpue value
  # p.ccpue: the fraction of EOS ccpue that was caught by each date
btf = read.table("2a_BTF_Data.txt", stringsAsFactors = F, header = T)
head(btf)

# get observed total run
# mean estimated total run as presented by Smith and Liller (2018)
true_N = read.table("2b_Total Run_Data.txt", stringsAsFactors = F, header = T)
head(true_N)

# get historical and forecasted run timing values
  # year: the year of the timing estimates
  # d50: the doy of 50% of the run complete (observed)
  # fcst_d50: the mean run timing forecast that year
rt_ests = read.table("2c_Run Timing_Data.txt", stringsAsFactors = F, header = T)
head(rt_ests)

# add thd d50 to historical data
# N.eos.dat$hist_d50 = rt.ests[rt.ests$year %in% N.eos.dat$year,"hist_d50"]

##### PART 2: SET UP UPDATING CHARACTERISTICS #####
years = 1995:2017  # these are the years that will be evaluated for oos predictions
n_years = length(years)

# mcmc characteristics
ni = 1e5        # total samples per chain
nb = 1e4        # number of burnin samples
nt = 1          # thinning interval
prop_sig = 0.5  # sd of proposal dist.

((ni - nb)/nt) * 2

# maximum run size to calculate density for; runs beyond this will have prior and likelihood values of zero
dens_max = 5e6

# get run size forecasts
prior_mean = true_N[true_N$year %in% (years - 1),"N"]
e = log(true_N[1:(nrow(true_N)-1),"N"]) - log(true_N[2:nrow(true_N),"N"])
prior_sd = sd(e)

# when Monte Carlo draws from dists to make likelihood, how many to draw (not mcmc)
n_mc = 1e6

# dates to evaluate updates on
dates = c("6/10", "6/17", "6/24", "7/1", "7/8", "7/15")
n_dates = length(dates)

##### PART 3: LEAVE ONE OUT ANALYSIS #####
models = c("null", "fcst")
n_models = length(models)
n_stat = 9
like_summ = post_summ = array(NA, dim = c(n_dates, n_stat, n_years, n_models))
prior_summ = matrix(NA, n_years, n_stat)
accept = bgr = neff = array(NA, dim = c(n_dates, n_years, n_models))

y = 1
d = 1
m = 1

start = Sys.time()
pdf("Output/MCMC_diag.pdf", h = 8, w = 5)
for (y in 1:n_years) {
  cat("---", years[y], "---\n")
  
  # sample the prior
  prior_samps = exp(rnorm(n_mc, log(prior_mean[y]) - 0.5 * prior_sd^2, prior_sd))
  prior_fun = approxfun(density(prior_samps, from = 0, to = dens_max), yleft = 0, yright = 0)
  prior_summ[y,] = summ(prior_samps)
  
  for (d in 1:n_dates) {
    cat("  ", paste(dates[d], years[y], sep = "/"), "\n", sep = "")
    # data for fitting the model
    fit_data = prepare_fit_data(dt = dates[d], yr = years[y], loo = T)
    
    # fit the regression model to all historical data (except the one left out)
    fit_null = lm(log(N) ~ q * ccpue, data = fit_data)
    fit_fcst = lm(log(N) ~ q * ccpue + d50, data = fit_data)

    for (m in 1:n_models) {
      cat("    Model ", m, "\n", sep = "")
      
      # get the data for prediction from the regression model
      pred_data = prepare_predict_data(fit_data = fit_data, dt = dates[d], 
                                       yr = years[y], n_mc = n_mc,
                                       rt_type = models[m])
      
      # sample the likelihood
      cat("      Likelihood Sampling...")
      if (models[m] == "null") {
        like_samps = sample_likelihood_null(fit = fit_null, pred_data = pred_data)
      } else {
        like_samps = sample_likelihood_fcst(fit = fit_fcst, pred_data = pred_data)
      }
      like_fun = approxfun(density(like_samps, from = 0, to = dens_max), yleft = 0, yright = 0)
      cat("Done\n")
      
      # sample posterior
      cat("      MCMC Sampling: Chain 1...")
      chain1 = MH(
        init = sample(x = prior_samps, size = 1),
        ni = ni, nb = nb, nt = nt, prop.sig = prop_sig,
        like.fun = like_fun, prior.fun = prior_fun
      ); cat("Done\n")
      
      cat("      MCMC Sampling: Chain 2...")
      chain2 = MH(
        init = sample(x = prior_samps, size = 1),
        ni = ni, nb = nb, nt = nt, prop.sig = prop_sig,
        like.fun = like_fun, prior.fun = prior_fun
      ); cat("Done\n")
      
      # create diagnostic plots; they get dumped in a PDF file placed in /Output directory
      diag.plots(list(chain1$post.samp, chain2$post.samp),
                 paste(dates[d], years[y], models[m], sep = "/"))
      
      # summarize posterior and likelihood inference
      post_summ[d,,y,m] = summ(c(chain1$post.samp, chain2$post.samp))
      like_summ[d,,y,m] = summ(like_samps)
      
      # summarize diagnostics
      accept[d,y,m] = round(mean(c(chain1$accept.rate, chain2$accept.rate)),2)
      neff[d,y,m] = round(unname(effectiveSize(as.mcmc.list(lapply(data.frame(chain1 = chain1$post.samp, chain2 = chain2$post.samp), mcmc)))))
      bgr[d,y,m] = round(unname(gelman.diag(x = as.mcmc.list(lapply(data.frame(chain1 = chain1$post.samp, chain2 = chain2$post.samp), mcmc)), multivariate = F)[[1]][1,1]),2)
    }
  }
}
dev.off()
stop = Sys.time()

dimnames(like_summ) = dimnames(post_summ) = 
  list(dates, names(summ(rnorm(10))), years, models)
dimnames(prior_summ) = list(years, names(summ(rnorm(10))))

save(like_summ, file = "Output/like_summ")
save(post_summ, file = "Output/post_summ")
save(prior_summ, file = "Output/prior_summ")

min(neff)
range(accept)




