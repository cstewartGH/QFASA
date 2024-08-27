---
output:
  html_document: default
  pdf_document: default
---
# News

## QFASA 1.2.1

* Made multiplicativeReplacement available to users of the package.
* Removed control specifications in the optimization function for p.MLE, 
p.MUFASA and likelihoodEstimates (for selection algorithms).  
Now using defaults.
* In p.SMUFASA, changed tolerance to "tol=1e-04".
* In p.SMUFASA, changed bounds of CCs to be 0 to 1 and then multiply estimated CCs by number of FAs. 
* In p.SMUFASA, changed starting values to average QFASA diet estimate.
* In p.SMUFASA, corrected description of prey.mat in Arguments in the help file.
* In conf.meth, added fat content to the function call in the example in the help file. 
* In testfordiff.ind.pval, added names to the list of outputs.


## QFASA 1.2.0

* Modified p.MLE to force prey species to be in alphabetical order. 
* Modified pseudo.pred (example and parameter description) and pseudo.pred.norm
(parameter description) to emphasize the proper order of the diet parameter.
* Fixed recent CRAN issues based on Kurt Hornik's emails.
* Added new functions for diet estimation based on forward selection and 
backward elimination algorithms.
* Changed split.prey, split.fatcont and mean.geometric to split_prey, 
split_fatcont and mean_geometric



## QFASA 1.1.2

* Added Tyler Rideout as a contributor.
* Fixed error with CS distance in conf.meth.
* Updated description of cal.mat in p.MUFASA.
* Removed ns1 as an argument in testfordiff.ind.pval (not needed).
* In MEANmeth, ensure that prey signatures sum to 1.  
* Added note to prey.mat description in p.QFASA.
* Updated p.MUFASA with a "See Also".
* Changed adonis to adonis2 in two.factor.rep and unbal.two.factor.rep.  Made corresponding necessary changes.
* Added function p.MLE
* Changed number of pseudo-predators generated to estimate error to 100 to help avoid switch error.


## QFASA 1.1.1

* In DESCRIPTION, changed order of names and added Jennifer as an author. Also modified wording in description.
* Added p.SMUFASA.
* Changed prey.mat description in p.MUFASA.
* Created Makefiles in folder src.
* Fixed title of QFASA_Workflow_Example vignette.
* Cleaned up MUFASA_Workflow_Example vignette.
* Added SMUFASA_Workflow_Example vignette.

## QFASA 1.1.0

* Made changes to DESCRIPTION including updating date, version, adding "cre" to Connie Stewart's role, removing "cre" from Justin Kamerman's role, adding Holly Steeves as an author, adding new imports and packages to LinkingTo. 
* Cleaned up p.QFASA: added reference, added to example, etc...
* Cleaned up pseudo.pred: changed description of parameters, added more explanation to "Value" and example.
* * Vignettes: 
  + Changed name of modelling workflow vignette to "QFASA_Workflow_Example" and updated it.  
  + Added a comment to the Parallel_Execution_for_Confidence_Intervals vignette.  
  + Added MUFASA_Workflow_Example vignette.
* beta.meths.CI is deprecated and replaced by conf.meth. Hid this function and bias.all from index. Made several corrections to conf.meth code.
* prey.cluster was missing code.  Added example and default values.
* Added functions to compute repeatability.
* Added functions to compute diet estimates and calibration coefficients simultaneously.
* Added functions to compute diet estimates via MUFASA (i.e. ML estimation)


Added the following functions:

## QFASA 1.0.3
* Added functions to test a difference between two independent samples of compositional data.
* Added functions to calculate confidence intervals on diet estimates.
* Added functions to perform hierarchical clustering on prey database.
* Added functions to calculate bias in diet estimates.
* Fixed: Add species names to diet estimates that are returned in p.QFASA
* Fixed: Set gamma=1 in p.QFASA
* Fixed: Add "Value" to p.QFASA function
* Fixed: Fix incorrect spelling of Shelley Lang and Chris Field
* Fixed: Package citation incorrect

## QFASA 1.0.2
* Fix bug whereby using data frames masks dimension incompatibilities of inputs.

## QFASA 1.0.1
* Fix documentation based on feedback from Connie Stewart. 
* Added separate tests for each distance measure.
* Fixed author list
* Added bibliographic references
* Add Shiny sample application


## QFASA 1.0.0
Initial release of code created in support of article

> Iverson, S. J., Field, C., Don Bowen, W. and Blanchard, W. (2004)
> QUANTITATIVE FATTY ACID SIGNATURE ANALYSIS: A NEW METHOD OF
> ESTIMATING PREDATOR DIETS Ecological Monographs 74: 211–235.
> doi:10.1890/02-4105
