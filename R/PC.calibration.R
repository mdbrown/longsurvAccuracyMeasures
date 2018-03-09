#' Assess calibration of a partly conditional model
#'
#' Assess model calibration for a partly conditional (PC) model fit using the 'partlyconditional' R package. For 'newdata' consisting of a time-to-event outcome and longitudinal marker information on a set of individuals up to a fixed time window (conditioning.time.window), this function calculates the predicted risks from the model at a specified future landmark prediction time ('prediction.time'). It then bins observations using the percentiles of predicted risk (use 'n.groups' to set how many bins) and calculates observed risk using the Kaplan-Meier estimator. KM-based confidence intervals for observed risk are also provided.
#'
#' @param pc.object output from PC.Cox or PC.GLM functions in the 'partlyconditional' R package used to fit partly conditional Cox/GLM models.
#' @param newdata data.frame with new data for which to estimate summary measures. All variables used to fit the PC.Cox/PC.GLM model must be present. Observations with missing data will be removed.
#' @param prediction.time  numeric value of prediction time (from conditioning.time into the future) to estimate future risk. Prediction time should be on the same scale as the measurement time and the survival times provided to fit the partly conditional model.
#' @param conditioning.time.window Time. All measurement times in newdata exceeding this conditioning time will be removed from analysis (and a message will be produced if silent = FALSE).
#' @param n.groups number of groups to bin observations (using predicted risks) to calculate
#' @param alpha  alpha level for confidence interval calculations. Default is 0.05 for 95\% confidence intervals.
#' @param plot set to TRUE to generate a plot of
#' @param silent set to TRUE to hide messages printed from function. Default is silent = FALSE.
#'
#'
#' @return
#'
#' a list with two elements
#'
#' \item{observed}{tibble consisting of }
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
#'calibration.result <- PC.calibration( pc.model.1,
#'                                      newdata = pc_data,
#'                                      conditioning.time.window = c(18, 24),
#'                                      prediction.time = 12,
#'                                      n.groups = 8)
#' calibration.result$observed
#' calibration.result$predicted
#'
#' @import dplyr
#' @import survival
#' @import partlyconditional
#' @import tibble
#' @export
PC.calibration <- function( pc.object,
                            newdata,
                            prediction.time,
                            conditioning.time.window,
                            n.groups = 10,
                            alpha = 0.05,
                            plot = TRUE,
                            silent = FALSE){

  call <- match.call()
  #write some checks
  stopifnot(is.element(class(pc.object), c("PC_cox", "PC_GLM")))
  if(class(pc.object) == "PC_cox"){
    stopifnot(is.numeric(prediction.time))
    stopifnot(length(prediction.time) ==1)
  }
  stopifnot(is.logical(plot))
  stopifnot(is.numeric(n.groups))
  stopifnot(is.data.frame(newdata))
  #check if newdata has the right variables
  if(!all(is.element(pc.object$variable.names, names(newdata)) )){
    stop("variable(s): {", paste(pc.object$names[which(!is.element(pc.object$variable.names, names(newdata)))], collapse = ", "), "}, not found in 'newdata' ")
  }


  stopifnot(is.numeric(conditioning.time.window))
  stopifnot(length(conditioning.time.window) <= 2)
  if(length(conditioning.time.window) < 2) conditioning.time.window <- c(conditioning.time.window[1], conditioning.time.window[1])
  conditioning.time.window <- sort(conditioning.time.window)
  stopifnot(is.logical(silent))


  #define some variables
  landmark.time <- conditioning.time.window[2] + prediction.time

  meas.time.name <- pc.object$variable.names[[4]]

  #filter out all measurement times greater than conditioning time. print a
  #message if silent is not TRUE
  newdata.si <- subset(newdata, newdata[[meas.time.name]] <= conditioning.time.window)

  if(!silent & nrow(newdata.si) < nrow(newdata) ){
    cat(paste0("... removing ",  nrow(newdata) - nrow(newdata.si) , " observations where ", meas.time.name, " is greater than conditioning.time.window = ", conditioning.time.window, ".\n"))
  }


  #get predicted risks on new data

  risk_dat  <- predict(pc.object, newdata = newdata.si, prediction.time = prediction.time)
  risk_dat.si <- subset(risk_dat,
                        risk_dat[[meas.time.name]] <= conditioning.time.window[2] &
                          risk_dat[[meas.time.name]] >= conditioning.time.window[1])
  if(!silent){
    cat(paste0("... assessing calibration using ", nrow(risk_dat.si), " observations within conditioning time window  [", paste(conditioning.time.window, collapse = ", " ), "].\n"))
  }
  #
  dc <- list(ti = landmark.time,
             pred.time = prediction.time,
             si = conditioning.time.window[2])

  timevar.name <- pc.object$variable.names[[2]]
  status.name <- pc.object$variable.names[[3]]


  risk <- risk_dat.si[, ncol(risk_dat.si)]

  risk.tiles <- c(0,quantile(risk, probs=c(1:n.groups/n.groups)[-n.groups], na.rm = TRUE, type = 1), 1)

  if(length(unique(risk.tiles)) < length(risk.tiles)) stop("not enough data for n.groups to be formed. Please reduce n.groups.")
  risk.cut <- cut(risk, risk.tiles)

  #observed event rates from kaplan meier
  #note that Xi_star and event_gradevol are unique to PASS dataset labels
  tstar <- risk_dat.si[[timevar.name]] - risk_dat.si[[meas.time.name]]
  status <- risk_dat.si[[status.name]]
  obs.risk <-summary(survfit(Surv(tstar, status )~risk.cut, conf.int = 1-alpha), times = prediction.time, extend = TRUE)


  #average predicted risks by group
  obs.dat <- data.frame("observed.risk" = 1- obs.risk$surv,
                        "lower" = 1-obs.risk$upper,
                        "upper" = 1-obs.risk$lower,
                        "percentile" = c(c(1:n.groups/n.groups)[-n.groups], 1) - 1/(n.groups*2))


  risk.dat = data.frame("predicted.risk" = risk,
                        "percentile" = rank(risk)/length(risk))


    if(plot){
      order = order(risk.dat$risk)
    plot(risk.dat$percentile[order],risk.dat$risk[order],
         xlab="Percentile",
         ylab="Risk Probability",
         xlim=c(0,1),ylim=c(0,1),cex.lab=1.2, type = "l")
    points(obs.dat$percentile,obs.dat$observed,pch=15)
    segments(obs.dat$percentile,obs.dat$observed,obs.dat$percentile,obs.dat$lower, lty = 2)
    segments(obs.dat$percentile,obs.dat$observed,obs.dat$percentile,obs.dat$upper, lty = 2)
    }
    #plot(1, type="n", axes=F, xlab="", ylab="")
    #legend(.6,1,cex=1.3,c('Model Predicted','Observed'),pch=c(1,15),bty = "n")

    list(observed = as.tibble(obs.dat),
         predicted = as.tibble(risk.dat))

}
