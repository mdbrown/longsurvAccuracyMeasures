#' Evaluate a partly conditional model
#'
#' Evaluate the predictive performance of a partly conditional (PC) model fit using the 'partlyconditional' R package.  Given validation data consisting of a time-to-event outcome and longitudinal marker information on a set of individuals up to a fixed time window (conditioning.time.window), this function estimates measures of prognostic accuracy for a specified future landmark prediction time ('prediction.time'). Bootstrap confidence intervals are provided for measures at the risk threshold provided.
#'
#' @param pc.object output from PC.Cox or PC.GLM functions in the 'partlyconditional' R package used to fit partly conditional Cox/GLM models.
#' @param newdata data.frame with new data for which to estimate summary measures. All variables used to fit the PC.Cox/PC.GLM model must be present. Observations with missing data will be removed.
#' @param prediction.time  numeric value of prediction time (from conditioning.time) to estimate future risk. Prediction time should be on the same scale as the measurement time and the survival times provided to fit the partly conditional model.
#' @param conditioning.time.window Time. All measurement times in newdata exceeding this conditioning time will be removed from analysis (and a message will be produced if silent = FALSE).
#' @param risk.threshold numeric threshold on the risk scale used to classify individuals as 'high-risk' for TPF, FPF, NPV, and PPV measures.
#' @param pnf.threshold  threshold q to estimate the proportion needed to follow PNF(q). Defaults to q = .5. PNF(q), is the proportion of the population at highest risk that one needs classify high risk in order that a proportion q of the cases will be identified.
#' @param pcf.threshold  threshold p to estimate the proportion of cases followed PCF(p) measure. Defaults to p = .25. PCF(p) is the proportion of cases included in the proportion p of individuals in the population at highest risk.
#' @param bootstraps  Number of bootstraps used for confidence intervals and standard error estimation. Default is 500. Set to 0 if no CI's are desired. See details below for further information on CI construction.
#' @param alpha  alpha level for confidence interval calculations. Default is 0.05 for 95\% confidence intervals.
#' @param silent set to TRUE to hide messages printed from function. Default is silent = FALSE.
#'
#' @details
#'
#' Summary measures are estimated using the set of observations with measurement time falling within the conditioning window set by 'conditioning.time.window'. If two observations from the same 'individual' fall within this time window, the observation with greater measurement time is used. These data can be accessed in the 'data.for.measures' element of the list returned from this function.
#'
#' Bootstrap standard error estimates and confidence intervals are constructed by resampling from this set of observations as well. The upper and lower CI bounds are calculated using the 1-alpha/2 and alpha/2 quantiles from the empirical bootstrap distribution for each measure. It is important to note that CI's are calculated conditional on the fitted model provided. An alternative bootstrapping routine should be run where the PC model is re-fit and evaluated in each bootstrap sample if the user wishes to account for the variation due to model fitting in standard error/CI estimates.
#'
#'
#' @return
#'
#' An object of class 'pc_evaluate' which is a list containing:
#'
#' \item{measures}{tibble consisting of estimates for prediction error, AUC, true positive fraction (TPF), false positive fraction (FPF), positive predictive value (PPV), negative predictive value (NPV), proportion of cases followed (PCF), proportion needed to follow-up (PNF), proportion high risk, and outcome prevalence.  }
#' \item{roc}{ tibble with components of the roc curve including, TPF, FPF, PPV, NPV, risk.threshold, and risk.percentile, at all risk thresholds observed.}
#' \item{data.for.measures}{data used to estimate summary measures. Consists of all observations with measurement time within the conditioning.time.window window. If two observations from the same 'individual' fall within this time window, the observation with greater measurement time is used.}
#' \item{call, bootstrap.info}{Inputs from function call. }
#'
#'
#'@examples
#'
#' library(partlyconditional)
#' library(longsurvAccuracyMeasures)
#' data(pc_data)
#'
#'pc.model.1 <-  PC.Cox(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  markers = c("marker", "marker_2"),
#'  data = pc_data,
#'  use.BLUP = c(FALSE, FALSE),
#'  knots.measurement.time = NA)
#'
#'pc.model.1
#'
#'result <- PC.evaluation(pc.cox.1,
#'                        newdata = pc_data,
#'                        conditioning.time.window = c(18, 24), #make predictions using data up to time 18-24
#'                        prediction.time = 12,  #one year predictions, conditional on historical data
#'                        risk.threshold = .4,
#'                        pnf.threshold = .5,
#'                        pcf.threshold = .18,
#'                        bootstraps = 25) #should use more bootstrap replicates in practice.
#'
#'result
#'
#' @import dplyr
#' @import survival
#' @import partlyconditional
#' @import tibble
#' @export

PC.evaluation <- function( pc.object,
                           newdata,
                           prediction.time = NULL,
                           conditioning.time.window = NULL,
                           risk.threshold,
                           pnf.threshold = .5,
                           pcf.threshold = .25,
                           bootstraps = 500,
                           alpha = .05,
                           silent = FALSE){

  call <- match.call()
  #write some checks
  stopifnot(is.element(class(pc.object), c("PC_cox", "PC_GLM")))
  if(class(pc.object) == "PC_cox"){
    stopifnot(is.numeric(prediction.time))
    stopifnot(length(prediction.time) ==1)
  }

  stopifnot(is.data.frame(newdata))
  #check if newdata has the right variables
  if(!all(is.element(pc.object$variable.names, names(newdata)) )){
    stop("variable(s): {", paste(pc.object$names[which(!is.element(pc.object$variable.names, names(newdata)))], collapse = ", "), "}, not found in 'newdata' ")
  }


  stopifnot(is.numeric(conditioning.time.window))
  stopifnot(length(conditioning.time.window) <= 2)
  if(length(conditioning.time.window) < 2) conditioning.time.window <- c(conditioning.time.window[1], conditioning.time.window[1])
  conditioning.time.window <- sort(conditioning.time.window)

  stopifnot(is.numeric(bootstraps))
  stopifnot(is.numeric(risk.threshold))
  stopifnot(is.logical(silent))
  stopifnot(risk.threshold < 1 & risk.threshold > 0 )
  stopifnot(pnf.threshold < 1 & pnf.threshold > 0)
  stopifnot(pcf.threshold < 1 & pcf.threshold > 0)

  #define some variables
  landmark.time <- conditioning.time.window[2] + prediction.time

  meas.time.name <- pc.object$variable.names[[4]]

  #filter out all measurement times greater than conditioning time. print a
  #message if silent is not TRUE
  newdata.si <- subset(newdata, newdata[[meas.time.name]] <= conditioning.time.window[2])

  if(!silent & nrow(newdata.si) < nrow(newdata) ){
   cat(paste0("... removing ",  nrow(newdata) - nrow(newdata.si) , " observations where ", meas.time.name, " is greater than conditioning.time.window = ", conditioning.time.window[2], ".\n"))
  }


  #get predicted risks on new data

  risk_dat  <- predict(pc.object, newdata = newdata.si, prediction.time = prediction.time)
  risk_dat.si <- subset(risk_dat,
                        risk_dat[[meas.time.name]] <= conditioning.time.window[2] &
                        risk_dat[[meas.time.name]] >= conditioning.time.window[1])
  if(!silent){
    cat(paste0("... calculating summary measures using ", nrow(risk_dat.si), " observations with measurement times within the conditioning time window  [", paste(conditioning.time.window, collapse = ", " ), "].\n"))
  }
  #
  dc <- list(ti = landmark.time,
             pred.time = prediction.time,
             si = conditioning.time.window[2])

  timevar.name <- pc.object$variable.names[[2]]
  status.name <- pc.object$variable.names[[3]]

  ss <- get.stats(dc = dc,
            risk_dat.si = risk_dat.si,
            timevar.name = timevar.name,
            status.name = status.name,
            meas.time.name = meas.time.name,
            risk.threshold = risk.threshold,
            pcf.threshold = pcf.threshold,
            pnf.threshold = pnf.threshold)
  stats = ss$stats
  ## bootstrap confidence intervals:

  if(bootstraps > 1){
    bootstrap.measures.dist <- matrix(ncol = bootstraps, nrow = nrow(stats) )
    for(b in 1:bootstraps){

        boot.ind <- sample.int(nrow(risk_dat.si), replace = TRUE)
        bootstrap.measures.dist[,b] <-
          c(unlist(get.stats(dc = dc,
                    risk_dat.si = risk_dat.si[boot.ind,],
                    timevar.name = timevar.name,
                    status.name = status.name,
                    meas.time.name = meas.time.name,
                    risk.threshold = risk.threshold,
                    pcf.threshold = pcf.threshold,
                    pnf.threshold = pnf.threshold)$stats[,2]))


    }
    stats$se <- apply(bootstrap.measures.dist, 1, sd)
    stats$lower <- apply(bootstrap.measures.dist, 1, quantile, type = 1, prob = alpha/2)
    stats$upper <- apply(bootstrap.measures.dist, 1, quantile, type = 1, prob = 1-alpha/2)
  }

  out <- list(measures = stats,
              roc = ss$roc,
              data.for.measures = as.tibble(risk_dat.si),
              call = call,
              bootstrap.info = list(alpha = alpha, bootstraps))
  class(out) = "pc_measures"
  return(out)

}


#' print function for PC.evaluation results
#' @export

print.pc_measures <- function(x, round = 3,  ...){

  cat("### Call:\n")
  print(x$call)
  cat("\n")
  cat("### Summary Measure Estimates:\n")
  ests <- x$measures
  ests[,-1] <- round(ests[,-1], round)
  print(ests)

  cat("\n")
  cat("### ROC curve estimates in pc.object$roc\n")
  cat("### showing top 5 rows...\n")
  print(head(x$roc, 5))
  cat("...\n")


  cat("### data used to estimate measures in pc.object$data.for.measures\n")
  }


