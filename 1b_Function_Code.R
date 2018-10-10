## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
## FUNCTIONS FOR LEAVE-ONE-OUT EVALUATION OF IN-SEASON BAYESIAN RUN ABUNDANCE UPDATING ##
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##

## ::::::::::::::::::::::: ##
## ----------------------- ##
## CODE BY: BEN STATON     ##
## ----------------------- ##
## LAST UPDATED: 7/31/2018 ##
## ----------------------- ##
## ::::::::::::::::::::::: ##

# these functions are split into two main categories:
  # those that are needed to generate samples from the likelihood PDF
  # those that are needed to generate samples from the posterior PDF

## ::::::::::::::::::::::::::::::::::: ##
## FUNCTIONS TO SAMPLE FROM LIKELIHOOD ##
## ::::::::::::::::::::::::::::::::::: ##

##### PREPARE DATA FOR FITTING HISTORICAL REGRESSION MODEL #####
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
# y = prepare_fit_data(dt = "6/10", yr = 2017, loo = F); nrow(y)

##### PREPARE DATA FOR PREDICTING FROM A REGRESSION MODEL #####

# fit_data: an object created using prepare_fit_data
# dt: the date
# yr: the year
# n_mc: the number of monte carlo samples to draw
# rt_type: the run timing info type: null or fcst
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
  
  output
}

##### SAMPLE THE LIKELIHOOD: DRAW SAMPLES FROM REGRESSION PREDICTIVE DIST #####
sample_likelihood = function(fit, pred_data) {
  
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

##### FUNCTION TO CALCULATE DIFFERENT TYPES OF ERRORS #####
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

## ::::::::::::::::::::::::::::::::::: ##
## FUNCTIONS TO SAMPLE FROM POSTERIOR  ##
## ::::::::::::::::::::::::::::::::::: ##

##### GENERATE PROPOSAL FROM LOGNORMAL JUMPING DIST #####
generate.proposal = function(x, sig) {
  rlnorm(1, log(x) - 0.5 * sig^2, sig)
}

##### CALCULATE METROPOLIS HASTINGS CORRECTION FOR ASYMMETRIC JUMPING DIST #####
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
##### FUNCTION TO PERFORM METROPOLIS-HASTINGS MCMC SAMPLING #####
## ARGUMENT LIST
  # init: the initial value of the chain
  # ni: the number of MCMC iterations (total)
  # nb: the number of burnin interations
  # nt: the thinning interval
  # prop.sig: sd of the lognormal proposal dist
  # like.fun: the likelihood density function
  # prior.fun: the prior density function

init = 5e7
ni = 10; nb = 1; nt = 1
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
    correct = calc.correction(current = current[i], proposal = proposal, sig = prop.sig)
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
  
  # burnin and thin chain
  keep.iter = numeric(ni)
  if(nt > 1) keep.iter[seq(nb + 1, ni, by = nt)] = 1 else keep.iter = rep(1, ni) 
  if(nb > 0) keep.iter[1:nb] = 0 else keep.iter[1] = 1
  post.samp = current[which(keep.iter == 1)]
  
  # summarize posterior
  summry = summ(post.samp)
  
  # create and return output
  output = list(accept.rate = accept.rate, summ = summry, post.samp = post.samp)
  return(output)
}

##### PLOT DIAGNOSTICS: TRACE PLOT AND DENSITY PLOTS #####
diag.plots = function(post.samp, plot.main = "") {
  
  if (length(post.samp) != 2) stop("This function can only plot two chains at a time!")
  
  chain1 = post.samp[[1]]
  chain2 = post.samp[[2]]
  
  par(mfrow = c(2,1), mar = c(2,5,2,1))
  dens1 = density(chain1)
  dens2 = density(chain2)
  x.lim = range(c(dens1$x, dens2$x))
  y.lim = range(c(dens1$y, dens2$y))
  plot(dens1, type = "l", col = "blue", xlim = x.lim, ylim = y.lim, main = plot.main, lwd = 2, las = 1, ylab = "")
  lines(dens2, col = "red", lwd = 2)

  ind = seq(1, length(chain1), by = 10)  # thin chains for traceplotting only - pdf loads too slowly otherwise
  y.lim = range(c(chain1, chain2))
  plot(chain1[ind], type = "l", col = "blue", ylim = y.lim, las = 1, ylab = "")
  lines(chain2[ind], col = "red", ylim = y.lim)
  
}


##### FUNCTION TO CALCULATE SUMMARY STATS #####
summ = function(x) {
  c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
}
