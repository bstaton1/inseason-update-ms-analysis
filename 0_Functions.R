## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## FUNCTIONS FOR LEAVE-ONE-OUT EVALUATION OF IN-SEASON BAYESIAN RUN ABUNDANCE UPDATING ##
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##

## :::::::::::::::::::::::: ##
## ------------------------ ##
## CODE BY: BEN STATON      ##
## ------------------------ ##
## LAST UPDATED: 10/20/2018 ##
## ------------------------ ##
## :::::::::::::::::::::::: ##

# these functions are split into three main categories:
  # 1.) those needed to prepare the x-y data to fit to and the x data to predict from
  # 2.) those needed to generate samples from the likelihood PDF
  # 3.) those needed to generate samples from the posterior PDF
  # 4.) those needed in post-processing (not plotting, just summaries)

##### 1.) FUNCTIONS TO PREPARE DATA #####

# PREPARE DATA FOR FITTING HISTORICAL REGRESSION MODEL
#  dt: the date (m/dd)
#  yr: the year (YYYY; could be any year between 1984-2017)
#  loo: do you wish to leave out the year you specified? if F, yr not used.

prepare_fit_data = function(dt, yr, loo = T) {
  # extract this day from btf for historical years
  btf_sub = subset(btf, date == dt)
  
  # merge in total run with btf_dat
  dat = merge(btf_sub, true_N, by = "year")
  
  # merge in run timing parameters
  dat = merge(rt_ests[,c("year", "d50")],
              dat, by = "year")
  
  # add a q period
  dat$q = factor(ifelse(dat$year < 2008, 1, 2))
  
  # if leaving-one-out, drop the year
  if (loo) {
    output = dat[dat$year != yr,]
  } else {
    output = dat
  }
  
  output[,c("year", "N", "q", "ccpue", "d50")]
}

# examples
# x = prepare_fit_data(dt = "6/10", yr = 2017, loo = T); nrow(x)
# y = prepare_fit_data(dt = "6/10", loo = F); nrow(y)

# PREPARE DATA FOR PREDICTING FROM A REGRESSION MODEL

# fit_data: a data.frame object created using prepare_fit_data
# dt: the date being estimated
# yr: the year being estimated (left out from fitting)
# n_mc: the number of monte carlo samples to draw
# rt_type: the run timing info type: "null" or "fcst"
# eiv: errors in variables - do you wish to include variability in run timing?

prepare_predict_data = function(fit_data, dt, yr, n_mc, rt_type, eiv = T) {
  # get the right q period
  q = factor(ifelse(yr < 2008, 1, 2))
  
  # error handler: rt_type
  if (!(rt_type %in% c("null", "fcst"))) {
    stop("'rt_type' must be one of 'null' or 'fcst'")
  }
  
  # obtain the correct d50 parameters
  if (rt_type == "null") {
    d50_mu = mean(fit_data$d50)
    d50_sig = sd(fit_data$d50)
  } else {
    d50_mu = subset(rt_ests, year == yr)$fcst_d50
    d50_sig = subset(rt_ests, year == yr)$fcst_se_d50
  }
  
  # obtain the correct ccpue to use
  ccpue = subset(btf, date == dt & year == yr)$ccpue
  
  # if doing errors in variables, generate mc samples of d50
  if (eiv) {
    output = data.frame(
      q = q,
      d50 = rnorm(n_mc, d50_mu, d50_sig),
      ccpue = ccpue
    )
  } else {
    output = data.frame(
      q = q,
      d50 = d50_mu,
      ccpue = ccpue
    )
  }
  
  # return the output data.frame
  output
}

##### 2.) FUNCTIONS TO SAMPLE FROM LIKELIHOOD #####

# DRAW SAMPLES FROM NULL REGRESSION MODEL PREDICTIVE DIST
sample_likelihood_null = function(fit, pred_data) {
  
  # extract the mean coefficients in the model
  coefs = coef(fit)
  
  # extract the covariance matrix on the coefficients
  Sigma = vcov(fit)
  
  # extract the residual standard error
  sigma = summary(fit)$sigma
  
  # the number of monte carlo draws
  n_mc = ifelse(nrow(pred_data) > 1, nrow(pred_data), n_mc)
  
  # turn q to be binary
  pred_data$q = ifelse(pred_data$q == 1, 0, 1) # turn to binary
  
  # create random coefficients
  b_rand = rmvnorm(n = n_mc, mean = coefs, sigma = Sigma)
  colnames(b_rand) = names(coefs)
  
  # create the log prediction
  pred_log_N = 
    b_rand[,"(Intercept)"] + 
    b_rand[,"q2"] * pred_data$q +
    b_rand[,"ccpue"] * pred_data$ccpue +
    b_rand[,"q2:ccpue"] * pred_data$q * pred_data$ccpue
  
  # add residual error and bring to natural scale
  exp(pred_log_N + rnorm(n_mc, -0.5 * sigma^2, sigma))
}

# DRAW SAMPLES FROM FCST REGRESSION MODEL PREDICTIVE DIST
sample_likelihood_fcst = function(fit, pred_data) {
  
  # extract the mean coefficients in the model
  coefs = coef(fit)
  
  # extract the covariance matrix on the coefficients
  Sigma = vcov(fit)
  
  # extract the residual standard error
  sigma = summary(fit)$sigma
  
  # the number of monte carlo draws
  n_mc = ifelse(nrow(pred_data) > 1, nrow(pred_data), n_mc)
  
  # turn q to be binary
  pred_data$q = ifelse(pred_data$q == 1, 0, 1) # turn to binary
  
  # create random coefficients
  b_rand = rmvnorm(n = n_mc, mean = coefs, sigma = Sigma)
  colnames(b_rand) = names(coefs)
  
  # create the log prediction
  pred_log_N = 
    b_rand[,"(Intercept)"] + 
    b_rand[,"q2"] * pred_data$q +
    b_rand[,"ccpue"] * pred_data$ccpue +
    b_rand[,"d50"] * pred_data$d50 +
    b_rand[,"q2:ccpue"] * pred_data$q * pred_data$ccpue
  
  # add residual error and bring to natural scale
  exp(pred_log_N + rnorm(n_mc, -0.5 * sigma^2, sigma))
}

##### 3.) FUNCTIONS TO SAMPLE FROM POSTERIOR #####

# GENERATE PROPOSAL FROM LOGNORMAL JUMPING DIST
  # it is centered on the current value
  # set the sig to control acceptance rate
generate.proposal = function(x, sig) {
  rlnorm(1, log(x) - 0.5 * sig^2, sig)
}

# CALCULATE MH CORRECTION FOR ASYMMETRIC JUMPING DIST
  # this is the ratio of P(c|p)/P(p|c)
  # where P() is the lognormal density function
  # c is the current state of MCMC
  # p is the proposed state of MCMC

calc.correction = function(current, proposal, sig) {
  # find the mu for dist centered on current
  mu.current = log(current) - 0.5 * sig^2
  
  # find the mu for dist centered on proposal
  mu.proposal = log(proposal) - 0.5 * sig^2
  
  # calculate density of drawing proposal from dist centered on current
  prop.given.curr = dlnorm(x = proposal, meanlog = mu.current, sdlog = sig)
  
  # calculate density of drawing current from dist centered on current
  curr.given.prop = dlnorm(x = current, meanlog = mu.proposal, sdlog = sig)
  
  # return the output
  return(curr.given.prop/prop.given.curr)
}

# FUNCTION TO PERFORM MH MCMC SAMPLING 
  # init: the initial value of the chain
  # ni: the number of MCMC iterations (total)
  # nb: the number of burnin interations
  # nt: the thinning interval
  # prop.sig: sd of the lognormal proposal dist
  # like.fun: the likelihood density function
  # prior.fun: the prior density function

MH = function(init, ni, nb, nt, prop.sig, like.fun, prior.fun) {
  # containers to store values
  current = numeric(ni)  # the current state in Markov chain
  accept = numeric(ni)   # an accept counter

  # step 1.) propose an initial value: is the current value for i = 1
  current[1] = init
  
  # sml = 1e-300
  sml = 0
  # repeat steps 2-4
  for (i in 1:(ni - 1)) {
    
    # step 2.) calculate the posterior prob of current value with bayes' theorem
    p.current = (like.fun(current[i]) + sml) * (prior.fun(current[i]) + sml)
    
    # step 3.) generate proposal and calculate posterior prob with bayes' theorem
    proposal = generate.proposal(current[i], sig = prop.sig)
    p.proposal = (like.fun(proposal) + sml) * (prior.fun(proposal) + sml)
    
    # step 4.) pick one to keep and one to discard
    # calculate MH correction factor
    correct = calc.correction(current = current[i], proposal = proposal, sig = prop.sig)
    # this is the metropolis jumping rule: accept proposal with probability equal to the ratio of post(p)/post(c)
    p.accept = min(c((p.proposal/p.current) * correct, 1))
    if(p.accept > runif(1, 0, 1)) {
      current[i+1] = proposal
      accept[i] = 1
    } else {
      current[i+1] = current[i]
      accept[i] = 0
    }
  }
  
  # calculate acceptance rate after burnin
  accept.rate = mean(accept[(nb+1):ni])
  
  # burnin and thin chain according to the settings
  keep.iter = numeric(ni)
  if(nt > 1) keep.iter[seq(nb + 1, ni, by = nt)] = 1 else keep.iter = rep(1, ni) 
  if(nb > 0) keep.iter[1:nb] = 0 else keep.iter[1] = 1
  post.samp = current[which(keep.iter == 1)]
  
  # summarize posterior
  summry = summ(post.samp)
  
  # create output
  output = list(
    accept.rate = accept.rate,  # fraction of all proposals accepted (before burnin/thin)
    summ = summry,              # posterior summary (mean, sd, and quantiles)
    post.samp = post.samp)      # all retained posterior samples
  
  # return the output
  output
}

# PLOT DIAGNOSTICS: TRACE PLOT AND DENSITY PLOTS
  # post.samp: a list with two elements; each storing the retained samples from an MH call
  # plot.main: the title of the plot

diag.plots = function(post.samp, plot.main = "") {
  
  if (length(post.samp) != 2) stop("This function can only plot two chains at a time!\n
                                   Supply each chain as an element from a list")
  
  # extract chains
  chain1 = post.samp[[1]]
  chain2 = post.samp[[2]]
  
  # set up device
  par(mfrow = c(2,1), mar = c(2,5,2,1))
  
  # create the density plot
  dens1 = density(chain1)
  dens2 = density(chain2)
  x.lim = range(c(dens1$x, dens2$x))
  y.lim = range(c(dens1$y, dens2$y))
  plot(dens1, type = "l", col = "blue", xlim = x.lim,
       ylim = y.lim, main = plot.main, lwd = 2, las = 1, ylab = "")
  lines(dens2, col = "red", lwd = 2)

  # create the trace plot
  ind = seq(1, length(chain1), by = 10)  # thin chains for traceplotting only - pdf loads too slowly otherwise
  y.lim = range(c(chain1, chain2)) # but have the y.lim show the full range in case bad samps were discarded
  plot(chain1[ind], type = "l", col = "blue", ylim = y.lim, las = 1, ylab = "")
  lines(chain2[ind], col = "red", ylim = y.lim)
  
}

##### 4.) POST-PROCESSING FUNCTIONS #####

# CALCULATE SUMMARY STATS
summ = function(x) {
  c(mean = mean(x), sd = sd(x),
    quantile(x, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
}

# CALCULATE DIFFERENT TYPES OF ERRORS
  # est and true are vectors of the estimated and true quantities
  # they must be in the same order and have the same length
calc_errors = function(est, true) {
  e = est - true     # 'absolute' error (in # fish)
  ae = abs(e)        # absolute 'absolute' error (in # fish; always positive)
  pe = e/true        # relative error (standardize size of error by size of abundance)
  ape = ae/true      # absolute relative error
  mult = est/true    # what multiple of the true value was the estimate?
  list(
    e = e, ae = ae, pe = pe, ape = ape, mult = mult
  )
}

