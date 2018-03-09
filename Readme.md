### longsurvAccuracyMeasures

The R package `longsurvAccuracyMeasures` can be used to evaluate a partly conditional model fit using the R package  [`partlyconditional`](https://github.com/mdbrown/partlyconditional). Given validation data consisting of a time-to-event outcome and longitudinal marker information on a set of individuals up to a fixed time window, this package includes functions to assess model calibration (`PC.calibration`) and prognostic accuracy (`PC.evaluation`) at a specified future landmark prediction time. Point estimates, along with bootstrap confidence intervals, are provided for the following accuracy measures: AUC, ROC (TPF/FPF), PPV, NPV, PCF (proportion cases followed), and PNF (proportion needed to follow).

[A brief tutorial is available here](http://rpubs.com/mdbrown/longsurvAccuracyMeasures)
