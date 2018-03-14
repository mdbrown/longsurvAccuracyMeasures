get.stats <- function(dc, risk_dat.si,
                       timevar.name,  status.name,meas.time.name,
                       risk.threshold,  pcf.threshold, pnf.threshold){

  #get d.pred
#browser()
  d.pred <- get.d.pred(dc = dc,
                       timevar = risk_dat.si[[timevar.name]],
                       eventvar =  risk_dat.si[[status.name]],
                       xi_var = risk_dat.si[[timevar.name]] - risk_dat.si[[meas.time.name]],
                       risk = risk_dat.si[, ncol(risk_dat.si)])

  #get roc
  tpfp <- get.tpfp(d.pred)
  roc <- tibble("risk.threshold" = c(1, d.pred$risk, 0),
                "risk.percentile" = c(1, rank(d.pred$risk), 0) /length(d.pred$risk),
                "TPF" = c(tpfp[,1, drop = FALSE]),
                "FPF" = c(tpfp[,2, drop = FALSE]),
                "PPV" = c(tpfp[,3, drop = FALSE]),
                "NPV" = c(tpfp[,4, drop = FALSE]),
                "landmark.time" = rep(dc$ti, nrow(tpfp)),
                "conditioning.time" =rep(dc$si, nrow(tpfp)),
                "prediction.time" = rep(dc$pred.time, nrow(tpfp)))


  #get stats

  stats <- tibble(
    measure = c("PE",
                "AUC",
                paste0("TPF(", risk.threshold, ")"),
                paste0("FPF(", risk.threshold, ")"),
                paste0("PPV(", risk.threshold, ")"),
                paste0("NPV(", risk.threshold, ")"),
                paste0("PCF(", pcf.threshold, ")"),
                paste0("PNF(", pnf.threshold, ")"),
                paste0("Proportion with risk > ", risk.threshold),
                "Prevalence"
    ),
    estimate =
      c(get.prediction.error(d.pred),
        auc <- get.auc(tpfp),
        get.tpfp(d.pred, risk.thr =risk.threshold)[2,1:4],
        get.pcf(d.pred, tpfp, prop = pcf.threshold),
        get.ppnf(d.pred, tpfp, prop = pnf.threshold),
        sum(d.pred$risk > risk.threshold)/length(d.pred$risk),
        sum(d.pred$hhat.inv.case)/sum(d.pred$hhat.inv)
      ))

     list(stats = stats, roc = roc)
  }




get.d.pred <- function(dc, timevar, xi_var, eventvar, risk){


  ix.case <- (xi_var <= dc$pred.time) & (eventvar == 1)
  ix.ctrl <- (xi_var > dc$pred.time)
  ix.cens <- (xi_var <= dc$pred.time) & (eventvar == 0)


  #get.censoring.weights is function below
  hhat.inv <- get.censoring.weights(dc = dc,
                                    timevar=timevar,
                                    eventvar=eventvar,
                                    xi_var=xi_var,
                                    ix.ctrl = ix.ctrl)
  hhat.inv <- matrix(hhat.inv, ncol=1)

  #browser()
  #print(paste("Cases:", sum(ix.case), "Controls:", sum(ix.ctrl), "Censored:", sum(ix.cens)))


  risk    <- risk[!ix.cens]
  #marker  <- data.s$marker[!ix.cens]
  ix.case <- ix.case[!ix.cens]
  ix.ctrl <- ix.ctrl[!ix.cens]

  if(!is.null(hhat.inv)){
    hhat.inv <- hhat.inv[!ix.cens]
  }else{
    hhat.inv <- rep(1, length(risk))
  }


  # check if risk is decreasing
  # if not, sort the data so that risk is decreasing (non-increasing)
  if(sum(diff(risk) > 0, na.rm = TRUE) != 0){
    ix.sort     <- sort(risk, decreasing = T, index.return = T)$ix
    risk        <- risk[ix.sort]
    #marker      <- marker[ix.sort]
    ix.case     <- ix.case[ix.sort]
    ix.ctrl     <- ix.ctrl[ix.sort]
    hhat.inv    <- hhat.inv[ix.sort]
  }
  #browser()
  d.pred <- list(risk         = risk,
                 risk.case    = risk[ix.case],
                 risk.ctrl    = risk[ix.ctrl],
                 #marker       = marker,
                 ix.case      = ix.case,
                 ix.ctrl      = ix.ctrl,
                 hhat.inv     = hhat.inv,
                 hhat.inv.case = hhat.inv[ix.case],
                 hhat.inv.ctrl = hhat.inv[ix.ctrl]
  )
  return(d.pred)
}


get.tpfp <- function(d.pred, risk.thr = NULL){

  #d.pred=dpred
  #risk.thr=NULL

  if((sum(d.pred$risk.case - sort(d.pred$risk.case, decreasing = T), na.rm =TRUE) != 0)
     | (sum(d.pred$risk.ctrl - sort(d.pred$risk.ctrl, decreasing = T), na.rm =TRUE) != 0)){
    print('Risk not sorted in decreasing order!')
    return(NULL)
  }

  if(is.null(risk.thr)){
    risk.thr <- d.pred$risk
  }

  tpfp <- matrix(NA, nrow = length(risk.thr), ncol = 4)

  for(i in 1:length(risk.thr)){
    tpfp[i, 1] <- sum(d.pred$hhat.inv.case[d.pred$risk.case > risk.thr[i]])
    tpfp[i, 2] <- sum(d.pred$risk.ctrl > risk.thr[i])

    #ppv
    tpfp[i, 3] = sum(d.pred$hhat.inv.case[d.pred$risk.case > risk.thr[i]]) /sum(d.pred$hhat.inv[d.pred$risk > risk.thr[i]])
    #npv
    tpfp[i, 4] = sum(d.pred$hhat.inv.ctrl[d.pred$risk.ctrl < risk.thr[i]]) /sum(d.pred$hhat.inv[d.pred$risk < risk.thr[i]])

  }
  tpfp[, 1] <- tpfp[, 1] / sum(d.pred$hhat.inv.case) #tpr
  tpfp[, 2] <- tpfp[, 2] / length(d.pred$risk.ctrl)  #fpr

  # sorting by increasing FPR & increasing TPR

  tpfp<-tpfp[order(tpfp[,2],tpfp[,1]),]
  tpfp<-rbind(c(0,0,NA,NA),tpfp,c(1,1,NA,NA)) #this is adding (0,0) and (1,1) for plotting and AUC
  return(tpfp)
}


# data = dataset on which survfit will be fit
# weights = weights for each individual in data
# new.times = times at which the censoring distribution is to be estimated
# ix.ctrl = true/false vector (true: survived beyond ti)
get.censoring.weights <- function(dc, timevar, xi_var, eventvar,weights = NULL,
                                  ix.ctrl){

  #note that this truncation step is not necessary (helps with computing time)
  ix <- xi_var > (dc$pred.time + 0.1) # the magic number is necessary, otherwise the KM estimate at dc$ti is off.
  xi_var[ix] <- (dc$pred.time + 0.1) #this truncates all times to (ti + .1) or less

  #this is modeling time to censoring, and S(t) is prob of not being censored past time t, given that you haven't censored up to time t
  cc  <- survfit(Surv(xi_var, eventvar == 0) ~ 1,
                 weights = weights, se.fit = F, type = 'kaplan-meier')
  new.times.case.ctrl<- xi_var
  new.times.case.ctrl[ix.ctrl] <- dc$pred.time #for controls, this truncates times to ti


  new.times.case.ctrl.sorted.incr <- sort(new.times.case.ctrl, decreasing = F, method = 'shell')
  recover.original.order <- rank(new.times.case.ctrl)

  cens.weights <- summary(cc, times = new.times.case.ctrl.sorted.incr)$surv[recover.original.order]
  #cens.weights are prob of not being censored at event times
  return(1/cens.weights) #this is inverse probability of not being censored
}





# ~~~~~~~~~~~~~~~
# ~~~ get.auc ~~~
# ~~~~~~~~~~~~~~~
get.auc <- function(tpfp){
  x <- c(0, tpfp[, 2])
  y <- c(tpfp[, 1], 1)
  n <- length(x)
  dx <- x[-1] - x[-n] #x[-1] is x with first value removed, x[-n] is x with last value removed
  #dx calculates the difference in x from t=1 to t=t+1
  mid.y <- (y[-n] + y[-1])/2 # y[-n] is last value removed, y[-1] is first value removed
  # mid.y calculates
  rm(x,y,n)
  return(sum(dx*mid.y))
}





# ~~~~~~~~~~~~~~~~
# ~~~ get.pcf ~~~
# ~~~~~~~~~~~~~~~~
# proportion of cases followed
get.pcf <- function(d.pred, tpfp, prop){
  weighted.prop <- prop * (sum(d.pred$hhat.inv))

  population.q <- 0
  weights.q.cases <- 0
  count <- 1

  while(population.q <= weighted.prop){
    population.q <- population.q + d.pred$hhat.inv[count]
    if(d.pred$ix.case[count] == T){
      weights.q.cases <- weights.q.cases + d.pred$hhat.inv[count]
    }
    count <- count + 1
  }
  return(weights.q.cases/sum(d.pred$hhat.inv.case))
}



# ~~~~~~~~~~~~~~~~
# ~~~ get.ppnf ~~~
# ~~~~~~~~~~~~~~~~
# proportion of the population needed to be followed (to capture prop of the cases)
#  we count the number of cases at highest risk to get to prop of cases
#  then we look at how many of all subjects are above that cutoff
#
# assumes that the subjects are sorted by decreasing risk
get.ppnf <- function(d.pred, tpfp, prop){

  reached.prop.cases <- prop * sum(d.pred$hhat.inv.case)

  # reals
  sum.q.case <- 0
  sum.q.ctrl <- 0

  # integer counters for cases, controls and the entire sample
  count.top.q.cases <- 0
  count.top.q.ctrls <- 0
  count.top.q.sample <- 0

  while(sum.q.case <= reached.prop.cases){
    count.top.q.sample <- count.top.q.sample + 1
    if(d.pred$ix.case[count.top.q.sample] == T){
      count.top.q.cases <- count.top.q.cases + 1
      sum.q.case <- sum.q.case + d.pred$hhat.inv.case[count.top.q.cases]
    }else{
      count.top.q.ctrls <- count.top.q.ctrls + 1
      sum.q.ctrl <- sum.q.ctrl + d.pred$hhat.inv.ctrl[count.top.q.ctrls]
    }
  }
  sum.num <- sum.q.case + sum.q.ctrl
  sum.den <- sum(d.pred$hhat.inv.case) + sum(d.pred$hhat.inv.ctrl)
  return(sum.num/sum.den)
}



get.prediction.error <- function(d.pred){

  mean.pe.case <- sum((1-d.pred$risk[d.pred$ix.case])^2 * d.pred$hhat.inv.case)
  mean.pe.ctrl <- sum((d.pred$risk[d.pred$ix.ctrl])^2 * d.pred$hhat.inv.ctrl)

  return((mean.pe.case + mean.pe.ctrl)/length(d.pred$risk))
}

